import json
from pathlib import Path

import pytest

pytest.importorskip("rdkit")

from opioidbench.chemistry.enrichment import MetadataCache, PubChemResolver


class _Response:
    def __init__(self, payload: dict):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self) -> bytes:
        return json.dumps(self.payload).encode("utf-8")


def test_pubchem_resolver_caches_validated_response(tmp_path: Path) -> None:
    cid_payload = {"IdentifierList": {"CID": [702]}}
    property_payload = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 702,
                    "SMILES": "CCO",
                    "ConnectivitySMILES": "CCO",
                    "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    "MolecularFormula": "C2H6O",
                    "MolecularWeight": "46.07",
                    "ExactMass": "46.041864812",
                }
            ]
        }
    }
    calls = []

    def opener(url, timeout):
        calls.append((url, timeout))
        return _Response(cid_payload if "/xref/RN/" in url else property_payload)

    cache = MetadataCache(tmp_path / "cache")
    resolver = PubChemResolver(
        cache,
        allow_network=True,
        min_interval_seconds=0,
        opener=opener,
        sleeper=lambda _: None,
    )
    first = resolver.resolve_cas("000064-17-5")
    second = resolver.resolve_cas("64-17-5")
    assert first is not None
    assert first.pubchem_cid == "702"
    assert first.canonical_smiles == "CCO"
    assert second == first
    assert len(calls) == 2
    assert "/xref/RN/64-17-5/cids/JSON" in calls[0][0]


def test_pubchem_resolver_is_offline_by_default(tmp_path: Path) -> None:
    resolver = PubChemResolver(MetadataCache(tmp_path / "cache"))
    assert resolver.resolve_cas("64-17-5") is None
