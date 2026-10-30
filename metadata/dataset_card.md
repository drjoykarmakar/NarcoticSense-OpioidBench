# Dataset Card — development build

## Status

A local, non-redistributed EI-MS development build has been validated against SWGDRUG Mass Spectral Library 3.14. A leakage-aware pilot benchmark is now frozen on the evidence-curated subset; this is not yet the final public benchmark release.

## Motivation

The dataset infrastructure is intended to support exact compound identification, opioid/non-opioid screening, family recognition, held-out analogue evaluation, open-set recognition, OOD detection, and uncertainty analysis as distinct tasks.

## Source and acquisition

The current source is SWGDRUG 3.14. Raw spectra are not committed because an explicit redistribution license has not been established. Local imports record source/version, retrieval date, archive/member names, and SHA-256 checksums.

## Local build audit

Current local build `opioidbench-eims-v0.1.0` contains 3,826 parsed EI-MS spectra resolving to 3,801 conservative compound identities. All 3,826 spectra pass current peak/schema validity checks. QC also records 1,502 placeholder CAS values, one non-placeholder CAS checksum failure, two exact duplicate-spectrum groups, 24 repeated conservative identity groups, and 19 identity groups with conflicting source names.

These counts are audit results only. Source metadata conflicts are not silently corrected.

## Chemical identity resolution

Valid CAS identifiers are preferred. When CAS is missing, placeholder, or invalid, an exact normalized source name is used as a conservative fallback identity. Fuzzy name matching is not used.

## Structures

The HPJ itself does not contain SMILES/InChI, but the same-release SWGDRUG NIST-format archive
contains library-linked molfiles. After conservative ID/name/CAS mapping plus molecular-formula
checks, the current local build assigns structures to 3,207/3,826 spectra, touching 3,199/3,801
conservative compound identities. At compound level, 3,197 identities have one unambiguous
NIST-derived canonical structure and 745 scaffold groups; two source-CAS identity collisions contain
conflicting structures and are excluded from compound/scaffold evaluation. Of the accepted spectrum
mappings, 3,088 have exact formula agreement and 119 are explicitly labeled as unlabeled parent
structures for deuterated isotopologues. One competing spectrum-level assignment is quarantined and
619 spectra remain unresolved. The accepted spectrum mappings contain 747 scaffold IDs.

A user-supplied ChEBI FULL SDF snapshot (SHA-256
`9e180253b1d6b57c197a06e5ea97c500488935c85a2df9ffcea6339e56ab3451`) was then streamed through
the conservative external resolver. It yields 938 accepted exact identity matches: 844 corroborate
existing NIST structures and 94 fill previously missing structures. Seven rows remain ambiguous and
264 are rejected for conflicts rather than forced. The combined unambiguous structure inventory is
therefore 3,291/3,801 compounds (86.6%), with 510 structureless identities and 792 exact
compound-level scaffold groups. The ChEBI additions introduce 47 scaffold groups beyond the
NIST-only inventory. PubChem CAS lookup remains a cache-first fallback and did not contribute new
matches in this offline audit.

## Opioid taxonomy

The current authoritative evidence seed produces 147 exact opioid-positive labels and 121
source-defined non-opioid comparator labels. Evidence is matched only by exact normalized name/alias
or validated CAS; no substring, fuzzy, spectral-similarity, scaffold, or fingerprint rule assigns a
class. Another 222 records are explicitly excluded from the biological target (217 analytical
variants plus five exact target-policy exclusions). The remaining 3,311 compound identities are
unknown/pending, not negative, and the current run has zero authoritative label conflicts.

The negative class is operational: an authoritative source explicitly classifies the substance in
a distinct comparator class or otherwise supports the non-opioid benchmark decision. It does not
assert zero activity at every opioid receptor.

A structure-only review-prioritization pass compares unresolved compounds with structured curated
opioids using radius-2, 2,048-bit chirality-aware Morgan fingerprints. The previous 51-compound
priority queue was closed with 48 evidence-backed opioid decisions and three N-benzyl target
exclusions; eight newly exposed near neighbors were then resolved with exact evidence. At Tanimoto
>= 0.75 the current queue is empty against 137 structured curated opioid references. Similarity
never assigns a label. Unresolved row-level review data remain local.

## Intended use

Research on analytical spectral identification, leakage-aware chemical ML benchmarking, open-set recognition, uncertainty, and related forensic-screening methodology.

## Inappropriate use

The software and dataset build are not independently validated for clinical, law-enforcement, evidentiary, diagnostic, or confirmatory analytical decision-making. They must not be represented as replacements for validated confirmatory analytical testing.

## Remaining limitations

Structure coverage remains incomplete (510/3,801 compound identities lack an accepted structure), and taxonomy coverage remains incomplete: 3,311 identities are still unknown/pending. Source metadata include placeholder identifiers and conflicting
names, including two source-CAS collisions with incompatible structures. Duplicate/replicate
interpretation still requires curation. Pilot model performance is reported only for the small
evidence-curated subset and must not be generalized to the full source inventory or external data.

## Formula/exact-mass companion metadata

SWGDRUG v3.14 also distributes a formula/exact-mass workbook. When explicitly supplied to the
build, OpioidBench requires row, name, and molecular-formula alignment before attaching its values.
The workbook states that monoisotopic masses were calculated from the library formulas using
Agilent PCDL software and are unverified. They are retained as source-provided derived metadata and
must not be interpreted as independent measurements, confirmatory identity evidence, or a substitute
for curated structures. The optional Agilent PCDL `.cdb` artifact is checksum-tracked but is not
parsed by the open-source pipeline.

### Additional source cross-checks

A same-release SWGDRUG Agilent-format archive may be used only to recover longer contributor/source
comments that are truncated in the JCAMP/HPJ representation. It is not counted as an independent
validation source.

ENFSI-DWG 2025.11 is tracked as a separate prospective EI-MS source. Its public Data Usage
Agreement restricts redistribution and public release without prior approval. Consequently, the
open benchmark repository stores only code and aggregate local cross-check behavior; ENFSI source
files and row-level derivatives must remain outside the public repository unless permission is
obtained.

A locally supplied ENFSI JDX library contains 2,050/2,050 valid spectra and aligns row-for-row with
the supplied entry list by library-entry number and molecular formula. It includes 18 exact
duplicate-spectrum pairs (36 records), which must not be treated as independent observations.
Against SWGDRUG 3.14, zero
peak lists are byte/exact identical because the two libraries use different intensity scales; 103
ENFSI records have identical m/z support and 361 have a best scale-invariant cosine similarity of
at least 0.9999. Eighty-eight of those 361 are not detected by the conservative metadata overlap
matcher. Therefore any ENFSI evaluation must use the spectrum-level audit, and the 1,689 records
below that screen are only spectrum-disjoint candidates, not proof of chemical independence.
Because ENFSI/DWG is also a documented SWGDRUG contributor, unseen-compound evaluation additionally
requires resolved chemical identities/structures plus compound/scaffold leakage checks.

## Pilot benchmark status

The current evidence seed curates 147 opioid positives and 121 source-defined non-opioid
comparators; 222 records are excluded from the biological target and 3,311 identities remain
unknown/pending. A frozen 5-fold pilot uses seed `20260929`.

The compound-held-out protocol contains 270 curated spectra. The exact Murcko-scaffold-held-out
protocol contains 257 structured curated spectra (139 opioid, 118 comparator) across 95 scaffold
groups. The stricter analogue-lineage protocol uses those same 257 spectra and additionally holds
out 17 broad curated opioid chemical lineages. All required leakage constraints are zero in every
fold for their respective protocols.

Random-forest mean balanced accuracy/AUROC are 0.938/0.987 for compound holdout, 0.850/0.878 for
scaffold holdout, and 0.616/0.651 for analogue-lineage holdout. The decline under stricter novelty
constraints is a pilot finding, not a final performance claim.

A separate ENFSI feasibility audit leaves 23 chemically de-leaked unique identity/scaffold
candidates, of which only two currently have source-backed opioid labels (20 comparator, one
pending). No external classifier performance metric is reported from this undersized positive set.
