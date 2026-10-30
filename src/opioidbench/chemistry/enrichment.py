"""Cached, rate-limited external chemical metadata enrichment."""

from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from typing import Any, Callable

from opioidbench.chemistry.identifiers import normalize_cas
from opioidbench.chemistry.structures import calculate_structure_features

PUBCHEM_PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class MetadataResolutionError(RuntimeError):
    """Raised when external metadata cannot be resolved conservatively."""


@dataclass(frozen=True, slots=True)
class EnrichmentRecord:
    query_type: str
    query: str
    source: str
    source_url: str
    retrieval_date: str
    pubchem_cid: str
    smiles: str
    canonical_smiles: str
    inchikey: str
    molecular_formula: str | None
    molecular_weight: float | None
    exact_mass: float | None = None

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class MetadataCache:
    """Filesystem JSON cache keyed by source/query, safe to persist outside Git."""

    def __init__(self, root: Path):
        self.root = root

    def _path(self, source: str, query_type: str, query: str) -> Path:
        digest = hashlib.sha256(f"{source}|{query_type}|{query}".encode("utf-8")).hexdigest()
        return self.root / source / f"{digest}.json"

    def get(self, source: str, query_type: str, query: str) -> dict[str, Any] | None:
        path = self._path(source, query_type, query)
        if not path.is_file():
            return None
        return json.loads(path.read_text(encoding="utf-8"))

    def put(self, source: str, query_type: str, query: str, payload: dict[str, Any]) -> Path:
        path = self._path(source, query_type, query)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return path


class PubChemResolver:
    """Conservative PubChem PUG REST resolver with local caching and backoff."""

    def __init__(
        self,
        cache: MetadataCache,
        *,
        allow_network: bool = False,
        min_interval_seconds: float = 0.25,
        max_retries: int = 3,
        opener: Callable[..., Any] | None = None,
        sleeper: Callable[[float], None] = time.sleep,
    ):
        self.cache = cache
        self.allow_network = allow_network
        self.min_interval_seconds = min_interval_seconds
        self.max_retries = max_retries
        self._opener = opener or urllib.request.urlopen
        self._sleeper = sleeper
        self._last_request = 0.0

    def resolve_cas(self, cas: str) -> EnrichmentRecord | None:
        normalized = normalize_cas(cas)
        if not normalized:
            raise ValueError(f"CAS is missing, placeholder, or invalid: {cas!r}")
        cached = self.cache.get("pubchem", "cas", normalized)
        if cached is not None:
            if cached.get("status") == "not_found":
                return None
            return EnrichmentRecord(**cached["record"])
        if not self.allow_network:
            return None
        record = self._fetch_registry_number(normalized)
        if record is None:
            self.cache.put(
                "pubchem",
                "cas",
                normalized,
                {
                    "status": "not_found",
                    "query": normalized,
                    "retrieval_date": date.today().isoformat(),
                },
            )
            return None
        self.cache.put(
            "pubchem",
            "cas",
            normalized,
            {"status": "resolved", "record": record.to_dict()},
        )
        return record

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_request
        remaining = self.min_interval_seconds - elapsed
        if remaining > 0:
            self._sleeper(remaining)

    def _request_json(self, url: str) -> dict[str, Any] | None:
        """Request one PubChem JSON endpoint with bounded retry/backoff."""
        payload: dict[str, Any] | None = None
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                with self._opener(url, timeout=30) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                self._last_request = time.monotonic()
                return payload
            except urllib.error.HTTPError as exc:
                self._last_request = time.monotonic()
                if exc.code == 404:
                    return None
                if exc.code not in {429, 500, 502, 503, 504} or attempt >= self.max_retries:
                    raise MetadataResolutionError(
                        f"PubChem request failed with HTTP {exc.code}"
                    ) from exc
            except urllib.error.URLError as exc:
                self._last_request = time.monotonic()
                if attempt >= self.max_retries:
                    raise MetadataResolutionError("PubChem request failed after retries") from exc
            self._sleeper(min(2**attempt, 8))
        return payload

    def _fetch_registry_number(self, query: str) -> EnrichmentRecord | None:
        encoded = urllib.parse.quote(query, safe="")
        cid_url = f"{PUBCHEM_PUG_REST}/compound/xref/RN/{encoded}/cids/JSON"
        cid_payload = self._request_json(cid_url)
        if cid_payload is None:
            return None
        cids = cid_payload.get("IdentifierList", {}).get("CID", [])
        if len(cids) != 1:
            if not cids:
                return None
            count = len(cids)
            raise MetadataResolutionError(
                f"PubChem registry-number query {query!r} returned {count} CIDs; refusing to guess"
            )
        cid = cids[0]
        properties = (
            "SMILES,ConnectivitySMILES,InChIKey,MolecularFormula,MolecularWeight,ExactMass"
        )
        url = f"{PUBCHEM_PUG_REST}/compound/cid/{cid}/property/{properties}/JSON"
        payload = self._request_json(url)
        properties_list = (payload or {}).get("PropertyTable", {}).get("Properties", [])
        if len(properties_list) != 1:
            if not properties_list:
                return None
            count = len(properties_list)
            raise MetadataResolutionError(
                f"PubChem query {query!r} returned {count} candidates; refusing to guess"
            )
        item = properties_list[0]
        smiles = (
            item.get("SMILES")
            or item.get("ConnectivitySMILES")
            or item.get("IsomericSMILES")
            or item.get("CanonicalSMILES")
        )
        inchikey = item.get("InChIKey")
        if not smiles or not inchikey:
            raise MetadataResolutionError(
                f"PubChem response for {query!r} lacks required structure fields"
            )

        calculated = calculate_structure_features(str(smiles))
        if calculated.inchikey != str(inchikey):
            raise MetadataResolutionError(
                f"PubChem/RDKit InChIKey mismatch for {query!r}; refusing enrichment"
            )
        molecular_weight = item.get("MolecularWeight")
        exact_mass = item.get("ExactMass")
        return EnrichmentRecord(
            query_type="cas",
            query=query,
            source="PubChem PUG REST",
            source_url=url,
            retrieval_date=date.today().isoformat(),
            pubchem_cid=str(cid),
            smiles=str(smiles),
            canonical_smiles=calculated.canonical_smiles,
            inchikey=calculated.inchikey,
            molecular_formula=(
                str(item["MolecularFormula"]) if item.get("MolecularFormula") else None
            ),
            molecular_weight=(float(molecular_weight) if molecular_weight is not None else None),
            exact_mass=(float(exact_mass) if exact_mass is not None else None),
        )


def enrich_compound_table(
    input_path: Path,
    *,
    output_dir: Path,
    cache_dir: Path,
    allow_network: bool = False,
    max_records: int | None = None,
) -> dict[str, int]:
    """Enrich a compound inventory by valid CAS, keeping sourced/calculated data separate."""
    import csv

    resolver = PubChemResolver(MetadataCache(cache_dir), allow_network=allow_network)
    output_dir.mkdir(parents=True, exist_ok=True)
    external_path = output_dir / "external_metadata.csv"
    descriptors_path = output_dir / "calculated_descriptors.csv"
    report_path = output_dir / "enrichment_report.json"

    with input_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if max_records is not None:
        rows = rows[:max_records]

    external_fields = [
        "compound_id",
        "query_type",
        "query",
        "status",
        "source",
        "source_url",
        "retrieval_date",
        "pubchem_cid",
        "smiles",
        "canonical_smiles",
        "inchikey",
        "molecular_formula",
        "molecular_weight",
        "exact_mass",
        "source_formula",
        "source_monoisotopic_mass",
        "source_exact_mass_delta",
        "error",
    ]
    descriptor_fields = [
        "compound_id",
        "canonical_smiles",
        "inchikey",
        "molecular_weight_rdkit",
        "clogp",
        "tpsa",
        "hbond_donors",
        "hbond_acceptors",
        "rotatable_bonds",
        "ring_count",
        "murcko_scaffold_smiles",
        "rdkit_source_smiles",
    ]
    counts = {
        "input_compounds": len(rows),
        "resolved": 0,
        "missing_or_invalid_cas": 0,
        "cache_miss_offline": 0,
        "not_found": 0,
        "source_formula_mismatch": 0,
        "source_exact_mass_mismatch": 0,
        "errors": 0,
    }
    external_rows: list[dict[str, Any]] = []
    descriptor_rows: list[dict[str, Any]] = []

    for row in rows:
        compound_id = row.get("compound_id", "")
        cas = normalize_cas(row.get("cas"))
        if not cas:
            counts["missing_or_invalid_cas"] += 1
            external_rows.append(
                {"compound_id": compound_id, "status": "missing_or_invalid_cas", "error": ""}
            )
            continue
        try:
            record = resolver.resolve_cas(cas)
        except (MetadataResolutionError, ValueError) as exc:
            counts["errors"] += 1
            external_rows.append(
                {
                    "compound_id": compound_id,
                    "query_type": "cas",
                    "query": cas,
                    "status": "error",
                    "error": str(exc),
                }
            )
            continue
        if record is None:
            status = "not_found" if allow_network else "cache_miss_offline"
            counts[status] += 1
            external_rows.append(
                {
                    "compound_id": compound_id,
                    "query_type": "cas",
                    "query": cas,
                    "status": status,
                    "error": "",
                }
            )
            continue

        source_formula = (row.get("molecular_formula") or "").strip() or None
        source_exact_mass = None
        raw_source_exact_mass = (row.get("source_monoisotopic_mass") or "").strip()
        if raw_source_exact_mass:
            try:
                source_exact_mass = float(raw_source_exact_mass)
            except ValueError:
                counts["errors"] += 1
                external_rows.append(
                    {
                        "compound_id": compound_id,
                        **record.to_dict(),
                        "status": "error",
                        "source_formula": source_formula or "",
                        "source_monoisotopic_mass": raw_source_exact_mass,
                        "error": "Invalid source_monoisotopic_mass in compound inventory",
                    }
                )
                continue

        if (
            source_formula
            and record.molecular_formula
            and source_formula != record.molecular_formula
        ):
            counts["source_formula_mismatch"] += 1
            external_rows.append(
                {
                    "compound_id": compound_id,
                    **record.to_dict(),
                    "status": "source_formula_mismatch",
                    "source_formula": source_formula,
                    "source_monoisotopic_mass": source_exact_mass or "",
                    "error": "PubChem molecular formula does not match SWGDRUG source formula",
                }
            )
            continue

        exact_mass_delta = None
        if source_exact_mass is not None and record.exact_mass is not None:
            exact_mass_delta = abs(source_exact_mass - record.exact_mass)
            if exact_mass_delta > 0.0001:
                counts["source_exact_mass_mismatch"] += 1
                external_rows.append(
                    {
                        "compound_id": compound_id,
                        **record.to_dict(),
                        "status": "source_exact_mass_mismatch",
                        "source_formula": source_formula or "",
                        "source_monoisotopic_mass": source_exact_mass,
                        "source_exact_mass_delta": exact_mass_delta,
                        "error": (
                            "PubChem exact mass differs from SWGDRUG source-derived exact mass "
                            "by more than 0.0001 Da; manual review required"
                        ),
                    }
                )
                continue

        counts["resolved"] += 1
        external_rows.append(
            {
                "compound_id": compound_id,
                **record.to_dict(),
                "status": "resolved",
                "source_formula": source_formula or "",
                "source_monoisotopic_mass": source_exact_mass or "",
                "source_exact_mass_delta": (
                    exact_mass_delta if exact_mass_delta is not None else ""
                ),
                "error": "",
            }
        )
        features = calculate_structure_features(record.smiles)
        descriptor_rows.append(
            {
                "compound_id": compound_id,
                "canonical_smiles": features.canonical_smiles,
                "inchikey": features.inchikey,
                "molecular_weight_rdkit": features.molecular_weight,
                "clogp": features.clogp,
                "tpsa": features.tpsa,
                "hbond_donors": features.hbond_donors,
                "hbond_acceptors": features.hbond_acceptors,
                "rotatable_bonds": features.rotatable_bonds,
                "ring_count": features.ring_count,
                "murcko_scaffold_smiles": features.murcko_scaffold_smiles,
                "rdkit_source_smiles": record.smiles,
            }
        )

    with external_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=external_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(external_rows)
    with descriptors_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=descriptor_fields)
        writer.writeheader()
        writer.writerows(descriptor_rows)
    report_path.write_text(json.dumps(counts, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return counts
