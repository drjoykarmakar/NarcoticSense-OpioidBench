# Data Sources

## SWGDRUG EI-MS

- Registry key: `swgdrug_eims`
- Public landing page: https://swgdrug.org/ms.htm
- Version used for local validation: 3.14, dated 2025-01-22.
- Supplied archive: `JCAMP_01222025.zip`
- Observed archive SHA-256: `da75d7c2e1383def05c691b94fbca37b83e244e9e333134b8a214c71905fd4a5`
- Archive member: `SWGDRUG 314.HPJ`
- Observed HPJ SHA-256: `c5703726b2a9a827bb0e4f737d1d3c995d324137d77c088948aa251f1d8e59fa`
- Local retrieval/import date for this development run: 2026-09-29.
- Explicit redistribution license: **not established**.

The public SWGDRUG page states that the library is available for download and provides a disclaimer, but this repository does not interpret that as a redistribution grant. Raw SWGDRUG spectra and processed spectrum payloads therefore remain outside Git. Automatic download is disabled until terms are documented.

The supplied HPJ is a concatenated JCAMP-DX stream containing 3,826 records. Records begin with `##TITLE=Library Entry ...` and do not use `##END=` separators. Chemical names are carried in `##CAS NAME`; `##TITLE` is treated as the source record identifier. The parser preserves peak order/values and verifies `NPOINTS` against the parsed peak count.

## ChEBI bulk chemical resolution

ChEBI is the preferred first external chemical resolver because its downloadable FULL SDF provides
curated chemical identifiers and structure fields in a bulk, reproducible form. The current
implementation streams a local `.sdf` or `.sdf.gz` once and retains only entries that touch an exact
SWGDRUG InChIKey, valid CAS Registry Number, or exact normalized source name/alias. Name-only
matches require molecular-formula support; structure conflicts, ambiguous candidates, and exact-mass
conflicts are quarantined rather than guessed.

The resolver does not automatically download ChEBI. Obtain a release from the official ChEBI
download page and pass it explicitly with `opioidbench metadata resolve --chebi-sdf ...`. The
2026-09-29 local audit used a user-supplied `chebi.sdf.gz` with SHA-256
`9e180253b1d6b57c197a06e5ea97c500488935c85a2df9ffcea6339e56ab3451`; its gzip header records
mtime `2026-09-09T20:35:50Z`. The generic filename does not itself encode a release identifier, so
the exact snapshot is anchored by the file hash rather than by an inferred release number. Official
downloads: https://www.ebi.ac.uk/chebi/downloads

That snapshot contains 199,606 parsable ChEBI entries. Conservative resolution against the 3,801
SWGDRUG compound identities yields 938 unique accepted ChEBI matches, seven ambiguous rows, 264
consistency-conflict rows, and 2,592 rows with no ChEBI candidate. Of the 938 accepted matches, 844
corroborate an already accepted NIST structure and 94 supply a previously missing structure. Final
unambiguous structure coverage is therefore 3,291/3,801 compounds (86.6%), with 510 compounds
still lacking a structure. The 94 newly structured compounds introduce 47 scaffold groups not
present in the NIST-only inventory, bringing the exact compound-level scaffold-group count to 792.
The resolver preserves source-linked NIST scaffold IDs rather than recalculating them during
external identity enrichment. ChEBI-derived fields remain attributable to ChEBI under CC BY 4.0.

## PubChem metadata enrichment

PubChem PUG REST is supported as an optional metadata source for compounds with a valid CAS Registry Number. Enrichment is cache-first and offline by default. Network access must be explicitly enabled; requests are rate-limited and retry with bounded backoff. Ambiguous registry-number results are rejected rather than guessed.

Externally sourced fields are written separately from RDKit-calculated descriptors. Every cached response records the query, source URL, retrieval date, and returned identifiers/structure fields.

The machine-readable source registry is the authority for acquisition policy and observed source checksums.

## SWGDRUG v3.14 companion compound list

The source-distributed workbook `Compound List 01222025.xlsx` was supplied alongside the
v3.14 spectra. Its observed SHA-256 is
`ab09287ff65ef993d01fffba9d393088b4f673c46cf13f02303d3b1e27bcdf61`.

The workbook contains one compound-name column with 3,826 data rows, exactly matching the
number and release order of records in `SWGDRUG 314.HPJ`. Red font is explicitly described by
the workbook header as indicating a new entry in v3.14; 229 rows carry that formatting.

Alignment against the HPJ yields:

- 3,805 exact name matches;
- 19 cases where the HPJ name is exactly 65 characters and is a prefix of the longer workbook
  name, consistent with a fixed-width HPJ truncation;
- 2 cases where the workbook name is a prefix of the longer HPJ name.

OpioidBench preserves both source strings. Only the observed 65-character HPJ-prefix pattern is
used as a high-confidence source-name repair for normalization. The two shorter-workbook cases
remain documented source disagreements and do not override the HPJ value. Any unrelated
row/name mismatch causes strict companion alignment to fail.

The workbook does not contain SMILES, InChI/InChIKey, PubChem CID, or a curated opioid taxonomy,
so it cannot by itself fill the remaining structure/taxonomy gaps.


## SWGDRUG v3.14 formula and exact-mass companion

The supplied `EXACTMASS_01222025.zip` has observed SHA-256
`cafebe04806cc90b73ecccf0822783f2ad4d87bbb21ac9f3496899c65c190bda` and contains:

- `SWGDrug 3.14 Excel list of names and formula.xlsx`
  (`17029ab729ede1c4f8aac6634b000d993b9292b0b3f42b6508198b1f814af911`);
- `SWGDrug 3.14.cdb`
  (`6347be2a9ce90c6452deefaa6cbe61d542bb28e893349cbd906d68c368a9b1b5`).

The workbook contains 3,826 rows after its header. Its compound names match the previously supplied
source companion list for all 3,826 rows, and its molecular formulas match the HPJ for all 3,826
rows. It supplies a monoisotopic exact mass for every row. The optional `M+proton` column contains
2,843 values and is blank for 983 rows. The sheet explicitly says that these masses were calculated
from the molecular formulas using Agilent PCDL software and are unverified. They must therefore be
treated as source-provided, formula-derived metadata rather than an independent measurement or
validation of molecular identity.

Two source records (mercuric bromide and mercuric chloride) have exact-mass values whose rounded
values do not equal the nominal `MW` in the HPJ. OpioidBench records this as a warning and preserves
both source values without choosing one as authoritative.

The `.cdb` member is an Agilent Personal Compound Database/Library file. Agilent documentation
indicates PCDLs can contain richer identifiers and structure fields, but the binary `.cdb` format is
not parsed by this open-source milestone. No SMILES, InChI/InChIKey, PubChem identifiers, or
structures are claimed to have been extracted from it. If a future legally permissible CSV/SDF
export is supplied, it should be imported as a separately provenance-tracked metadata source.

Relevant vendor documentation for the `.cdb` interpretation:

- https://www.agilent.com/Library/usermanuals/Public/G3336-90006_PCDL_QuickStart.pdf
- https://www.agilent.com/cs/library/usermanuals/public/G3336-90030_PCDLManager_QuickStart.pdf

## SWGDRUG v3.14 Agilent-format cross-check

The supplied `AGILENT_01222025.zip` has observed SHA-256
`5a41addddf725c78aa1ba61fc3bb2791e6d4a64e706703fabb9c44c1c3ede6be`.
It is another vendor-format representation of the same SWGDRUG v3.14 release, not an
independent data source. The embedded `SWGDRUG 3.14.L/HEADER` member has observed SHA-256
`d218828d7bd5cf13685a5e4a94c48abb7a0d8eb42c5902ceded915d4f72c9044`.

OpioidBench does not attempt to parse the proprietary Agilent spectral binary. It extracts only
printable metadata strings from the `HEADER` member and validates them against the already parsed
JCAMP release in source order. On the supplied release, all 3,826 records align. The HPJ source
comment is already identical for 1,779 records; for 2,047 records it is an exact prefix of a longer
Agilent-header comment. The longer value is retained as `SOURCE_COMMENT_FULL` while the raw HPJ
`SAMPLE_DESCRIPTION` remains unchanged.

This matters because SWGDRUG documents that contributors are cited in the comments field for each
entry. Cross-format comment recovery therefore improves spectrum-level provenance, but it does not
add statistically independent spectra.

## Opioid taxonomy evidence

`metadata/opioid_evidence.csv` is an open evidence seed, separate from SWGDRUG source data. Opioid
status/family labels come from traceable public pharmacology/classification sources including WHO
ATC/DDD, WHO ECDD, UNODC publications/EWA pages, and selected peer-reviewed or NLM evidence where
needed for exact high-priority cases. Exact CAS/name/alias matching is required; substring, fuzzy,
and scaffold-based label propagation are prohibited. Identity-only corroboration is recorded
separately from pharmacological evidence.

The current SWGDRUG 3.14 local build has 147 exact opioid-positive matches and 121 source-defined
non-opioid comparator matches. Another 222 records are excluded from the biological target: 217
isotope-labelled/TMS analytical variants and five exact target-policy exclusions. The other 3,311
compound identities remain unknown/pending. Absence from opioid evidence never creates a negative
label. The tracked `metadata/opioid_taxonomy.csv` contains curated and excluded rows; pending review
rows remain local/git-ignored. Exact target exclusions are tracked in
`metadata/taxonomy_exclusions.csv`.

A separate structure-similarity review command ranks unresolved compounds against structured
curated opioids using the documented Morgan configuration. The previous 51-compound queue was
adjudicated as 48 evidence-backed opioids and three N-benzyl target exclusions. Iterative
recomputation exposed eight additional near neighbors, which were also resolved with exact evidence.
At Tanimoto >= 0.75 the current priority queue is empty against 137 structured curated opioid
references. Similarity remains a review-order signal only and cannot assign opioid status.

The INCB Yellow List may corroborate official names/CAS identifiers, but it is not used by itself to
assign opioid status because it is a broader list of narcotic drugs under international control.

## SWGDRUG v3.14 NIST-format structure extraction

The supplied NIST-format archive is a same-release representation of SWGDRUG 3.14 and is therefore
not independent evidence. Its `USRSTRUC.DB` member contains 3,297 library-linked structure blocks.
The open pipeline parses 3,296 with RDKit; one linked molfile is quarantined after sanitization
fails. Mapping to the current JCAMP build is conservative: current library ID is preferred, followed
only by exact source-name + compatible-formula or valid-CAS + compatible-formula fallback.

On the supplied release, 3,207/3,826 spectra receive an accepted structure, covering 3,199/3,801
conservative compound identities. Formula relations are 3,088 exact and 119 `isotopologue_parent`
for deuterated source formulas whose NIST molfile encodes the corresponding unlabeled parent. One
spectrum has two distinct structure candidates and is excluded pending review. The accepted
structures produce 747 exact Bemis-Murcko scaffold IDs; acyclic structures use the explicit
`murcko:acyclic` marker so later split code must handle that group deliberately.

The build records both archive/member SHA-256 values and writes a local-only
`nist_structure_review.csv` for parse failures, unresolved formula/identity cases, competing
structures, and exact-mass QC warnings. No NIST/SWGDRUG structure payload is committed to Git.

The supplied Shimadzu-format archive adds no unique chemistry fields needed by the current build:
all 3,826 formulas agree with JCAMP and its 19 longer names correspond to the same HPJ truncations
already resolved by the official companion list. It is retained as provenance rather than counted
as an independent source. The supplied changelog is likewise retained as release-history provenance.

## ENFSI-DWG EI-MS 2025.11

- Registry key: `enfsi_dwg_eims`
- Public page: https://enfsi.eu/about-enfsi/structure/working-groups/drugs/
- Version: 2025.11.
- Publicly stated library size: 2,050 spectra.
- Supplied entry-list workbook: `2025ENFSI._DWG_MS_EntryList_2025_11.xlsx`
- Observed entry-list SHA-256: `55cafbf460edc96e88c5a499a95d04c71b1ba34c5d144cdc8af1acc26bee471b`
- Supplied JDX archive SHA-256: `90f5305c9bb431c2e562dd283a05e3fe001be57168d15a1a546f803f55de8946`
- JDX member SHA-256: `f7f96fc24dda46694e9383aea2b1d98025a6cd03f584479503a41213ea450146`
- Redistribution status: **restricted**.

ENFSI's public Data Usage Agreement states that its spectral-library data must not be shared,
published, or released publicly without prior approval; redistribution for profit is prohibited;
and redistribution to third parties is not allowed without prior approval. The open repository
therefore contains neither the ENFSI source files nor row-level source derivatives. Automatic
download is disabled. Local aggregate auditing is supported.

The entry list contains 2,050 rows, 2,035 unique names, 16 rows without a molecular formula, and 15
exact duplicate metadata rows beyond the first occurrence. A conservative metadata-only audit
against SWGDRUG v3.14 requires an explicit alias plus matching formula and rounded nominal molecular
weight. Under that rule, 938 ENFSI entry rows map uniquely to 855 SWGDRUG compounds; nine are
ambiguous and 1,103 remain unmatched. No fuzzy spelling correction is used for this aggregate count.

The supplied JDX contains 2,050 EI-MS records and all 2,050 pass current peak/schema validity checks.
It contains 18 exact duplicate-spectrum pairs (36 records), 931 source CAS values equal to `0` that
are treated as explicit missing-value placeholders, 613 other CAS strings that fail conservative CAS
validation, and two invalid/placeholder nominal molecular-weight values. The JDX aligns safely to the
entry list: all 2,050 `Library Entry N` identifiers and all 2,050 molecular formulas agree row-for-row;
2,048 nominal molecular weights agree exactly, with two placeholder-style values in the JDX. The
local JDX itself is not copied into repository artifacts.

For spectrum-level leakage auditing, OpioidBench represents each spectrum on integer nominal m/z bins
and compares raw intensity vectors after L2 normalization. This makes the comparison invariant to
the observed library-wide intensity scale difference (ENFSI base peak 999; SWGDRUG base peak 9,999).
The real-source audit found zero exact peak-list matches, 103 ENFSI records with identical m/z
support to at least one SWGDRUG spectrum, and 361 ENFSI records whose best SWGDRUG cosine similarity
is at least 0.9999. Threshold sensitivity is retained in the machine-readable audit: 95 records are
>=0.99999, 450 are >=0.999, 870 are >=0.99, and 1,350 are >=0.95. The 0.9999 flag is a deliberately
strict leakage-risk screen, not proof that two records are the same physical acquisition.

Of the 361 >=0.9999 matches, 271 also have a unique high-confidence metadata overlap, two have
ambiguous metadata overlap, and 88 were left unmatched by metadata. This demonstrates that a
metadata-only overlap audit misses materially important spectral overlap. After the >=0.9999 screen,
1,689 ENFSI records remain as *spectrum-disjoint candidates*. Of the 1,103 metadata-unmatched rows,
1,015 also fall below the 0.9999 spectral screen. Neither count establishes compound-level
independence: unresolved synonyms, salts/derivatives, or structures may still overlap.

Because SWGDRUG explicitly lists ENFSI/DWG as a library contributor, ENFSI 2025.11 must not be
described as an independent external test set without these checks. For unseen-compound evaluation,
structure/identity resolution and compound/scaffold leakage auditing are additionally required.

### Chemical/scaffold external-candidate audit

The external-candidate audit applies the metadata and spectral screens first, then resolves
structures and removes exact SWGDRUG chemical identity/scaffold overlap. Starting from 1,015
metadata-unmatched records below the 0.9999 spectral-overlap threshold, 66 receive conservative
ChEBI structures and three additional records receive public structures only after exact
name/alias-plus-formula agreement. Thus 69 base records have accepted structures.

After excluding exact InChIKey overlap, exact Murcko-scaffold overlap, and records belonging to an
ENFSI exact duplicate-spectrum group, 23 unique InChIKey/scaffold candidates remain. Exact public
taxonomy evidence labels two as opioids and 20 as non-opioid comparators; one remains pending and
there are zero label conflicts. Because only two opioid-positive candidates are available, the
project treats this as an external-test feasibility result and does not report classifier
performance on this cohort.

`metadata/external_structure_evidence.csv` contains only public supplemental structure identifiers
and structures; it does not contain ENFSI source rows. The resolver never selects those structures
by ENFSI CAS alone, which protects against source-identifier inconsistencies. ENFSI row-level
candidate decisions remain restricted/local.

