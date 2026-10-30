"""Optional RDKit-backed structure standardization and descriptor calculation."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass


class ChemistryDependencyError(RuntimeError):
    """Raised when an optional chemistry dependency is unavailable."""


@dataclass(frozen=True, slots=True)
class StructureFeatures:
    canonical_smiles: str
    inchikey: str
    molecular_weight: float
    clogp: float
    tpsa: float
    hbond_donors: int
    hbond_acceptors: int
    rotatable_bonds: int
    ring_count: int
    murcko_scaffold_smiles: str


def _rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise ChemistryDependencyError(
            "RDKit is required for structure calculations; install opioidbench[chem]."
        ) from exc
    return Chem, Crippen, Descriptors, Lipinski, rdMolDescriptors, MurckoScaffold


def calculate_structure_features(smiles: str) -> StructureFeatures:
    """Parse one SMILES and calculate publication-relevant structure descriptors."""
    Chem, Crippen, Descriptors, Lipinski, rdMolDescriptors, MurckoScaffold = _rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    inchikey = Chem.MolToInchiKey(mol)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)
    return StructureFeatures(
        canonical_smiles=canonical,
        inchikey=inchikey,
        molecular_weight=float(Descriptors.MolWt(mol)),
        clogp=float(Crippen.MolLogP(mol)),
        tpsa=float(rdMolDescriptors.CalcTPSA(mol)),
        hbond_donors=int(Lipinski.NumHDonors(mol)),
        hbond_acceptors=int(Lipinski.NumHAcceptors(mol)),
        rotatable_bonds=int(Lipinski.NumRotatableBonds(mol)),
        ring_count=int(Lipinski.RingCount(mol)),
        murcko_scaffold_smiles=scaffold_smiles,
    )


def morgan_fingerprint(
    smiles: str,
    *,
    radius: int = 2,
    n_bits: int = 2048,
    use_chirality: bool = True,
):
    """Return an explicit-configuration Morgan bit vector."""
    Chem, *_ = _rdkit()
    try:
        from rdkit.Chem import rdFingerprintGenerator
    except ImportError as exc:  # pragma: no cover
        raise ChemistryDependencyError("RDKit fingerprint generator is unavailable") from exc
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius,
        fpSize=n_bits,
        includeChirality=use_chirality,
    )
    return generator.GetFingerprint(mol)


def scaffold_id(murcko_scaffold_smiles: str) -> str:
    """Return a stable identifier for an exact Bemis-Murcko scaffold.

    Acyclic molecules share the explicit ``murcko:acyclic`` marker; downstream
    split code must handle that large group deliberately rather than silently
    assigning each acyclic molecule a different pseudo-scaffold.
    """
    if not murcko_scaffold_smiles:
        return "murcko:acyclic"
    digest = hashlib.sha256(murcko_scaffold_smiles.encode("utf-8")).hexdigest()[:20]
    return f"murcko:{digest}"
