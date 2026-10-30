# Results

Generated row-level outputs are ignored by Git. Real split assignments, predictions, taxonomy
review queues, and ENFSI candidate manifests remain local because they contain source-derived
record/compound identifiers. Aggregate-only pilot metrics and external-test feasibility counts are
documented under `benchmarks/`.

The current v0.4.0 pilot uses 147 curated opioids and 121 source-defined comparators. A further 222
records are excluded from the biological target rather than forced into either class. The previous
51-compound high-priority review queue is closed and the current Tanimoto >=0.75 queue is empty. The
pilot reports compound-, scaffold-, and analogue-lineage-held-out cross-validation, plus an
aggregate ENFSI chemical/scaffold de-leakage audit. The ENFSI candidate pool still has only two
source-backed opioid positives, so no external classifier metric is reported.

These are development results, not final OpioidBench performance claims.
