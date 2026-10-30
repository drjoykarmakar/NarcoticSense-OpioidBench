# Methodology

## Current implemented scope

The implemented pipeline establishes traceable source import, EI-MS parsing, explicit raw/validated/normalized data layers, QC, conservative chemical identity normalization, structure-enrichment interfaces, RDKit calculations, taxonomy schema, and dataset manifests.

For SWGDRUG 3.14, `##TITLE` delimits concatenated source records and is preserved as the source record identifier. `##CAS NAME` is preserved as the raw compound name. Peak values are not normalized or binned at this stage.

Chemical identity resolution is deliberately conservative. A valid CAS Registry Number is preferred when available; otherwise an exact Unicode/whitespace-normalized source name is used to generate a stable compound key. Fuzzy name similarity is not used for merging. Placeholder or checksum-invalid CAS values normalize to missing and remain visible through QC flags.

Structure metadata are not inferred from formula or compound-name patterns. Optional PubChem enrichment operates only on valid CAS Registry Numbers in the current implementation, caches responses, and rejects ambiguous mappings. RDKit-derived descriptors are stored separately from externally sourced metadata.

Opioid taxonomy is evidence curated. No compound is labeled as an opioid solely because its name contains an opioid-related keyword, and absence from an opioid source never creates a negative label. `non-opioid` is an operational comparator label backed by an authoritative source-defined alternative class, not a claim of zero opioid-receptor activity. Structural similarity is used only to prioritize unresolved records for human review.

The curated subset can support controlled pilot experiments, but benchmark-freeze model evaluation remains deferred until the high-priority review queue and compound/scaffold leakage audits are complete.
