# Pilot benchmark results: interpretation and publication checklist

> **Status:** Preliminary, evidence-curated v0.4.0 pilot. Not a final benchmark, prospective forensic validation, or external classification result.
>
> **Data sharing:** This document contains aggregate research notes only. Row-level assignments, predictions, identifiers, spectra, and restricted-source records remain local, consistent with [results/README.md](../results/README.md).

## Evaluation design

The deterministic five-fold pilot uses seed `20260929`. The repository currently documents:

| Held-out protocol | Curated spectra | Opioid spectra | Comparator spectra | Evaluation constraint |
| --- | ---: | ---: | ---: | --- |
| Compound | 270 | 149 | 121 | No compound overlap |
| Exact Murcko scaffold | 257 | 139 | 118 | No compound or scaffold overlap |
| Analogue lineage | 257 | 139 | 118 | No compound, scaffold, or curated opioid-lineage overlap |

The analogue-lineage protocol groups exact scaffolds with 17 curated broad opioid chemical lineages. This is a generalization constraint, not a pharmacological relabeling. The random-held-out analysis discussed during exploratory review should not be treated as a frozen, independently verified public benchmark until its source artifacts and run provenance are confirmed.

## Frozen aggregate baseline reported in the repository README

The [main README](../README.md) reports the following **random forest** means for the fixed baseline using nominal integer m/z bins through 700, square-root intensity, and L2 normalization:

| Held-out protocol | Mean balanced accuracy | Mean AUROC |
| --- | ---: | ---: |
| Compound | 0.938 | 0.987 |
| Exact Murcko scaffold | 0.850 | 0.878 |
| Analogue lineage | 0.616 | 0.651 |

**Interpretation:** The marked reduction under analogue-lineage holdout suggests that conventional compound-held-out performance is not sufficient evidence of generalization to unseen opioid chemical families. It does **not** establish a causal effect of the splitting protocol: cohort composition and fold heterogeneity also matter.

### Exploratory screenshot-derived results (not yet reconciled with frozen aggregates)

Earlier local-terminal screenshots showed random-forest **ordinary accuracy**, F1, and sensitivity estimates of approximately:

| Protocol | Mean accuracy | Mean F1 | Mean opioid sensitivity |
| --- | ---: | ---: | ---: |
| Random | 0.937 | 0.941 | 0.926 |
| Compound | 0.941 | 0.946 | 0.933 |
| Scaffold | 0.864 | 0.795 | 0.764 |
| Analogue lineage | 0.621 | 0.321 | 0.256 |

These values are retained **only as exploratory notes**, not as independently verified publication tables. They must be reconciled against local `summary.json`, `fold_metrics.csv`, and `predictions.csv` artifacts before being cited. In particular, **ordinary accuracy and balanced accuracy are different metrics**, and screenshot-derived scaffold/analogue AUROC values differ from the frozen values in the main README. Do not silently combine metrics from distinct runs.

## Leakage and overlap

The repository README reports zero compound overlap, exact-spectrum overlap, and cross-fold raw-intensity cosine similarity pairs at or above 0.9999 for the compound-held-out pilot. Scaffold-held-out evaluation also reports zero exact-scaffold overlap; analogue-lineage-held-out evaluation additionally reports zero curated opioid-lineage overlap. These are protocol-specific audit findings, not a guarantee that all possible leakage pathways have been excluded.

## Five checks before manuscript submission

1. **Complete model table:** Extract all four models (spectral 1-NN, logistic regression, linear SVM, random forest) from one consistent, frozen run for every supported protocol. Report the metric definitions and fold mean/standard deviation.
2. **Per-fold errors:** Compute fold-specific TP, FP, TN, FN and pooled out-of-fold confusion matrices from local predictions; explicitly distinguish pooled metrics from fold means. Examine the highly variable analogue-lineage folds and test-set class counts.
3. **Leakage audit:** Preserve machine-readable, fold-specific audit summaries locally. Publish only aggregate, redistributable audit counts after verifying all folds and source restrictions.
4. **Comparable cohorts:** Identify exactly which 13 spectra are excluded from the structure-dependent protocols; where feasible, rerun easier splits on the same 257-spectrum subset to separate cohort effects from protocol effects.
5. **Limitations and provenance:** Record dataset snapshot, source licenses, taxonomy evidence version, code commit, random seed, feature representation, split assignments checksum, and execution environment. Distinguish pilot cross-validation from external validation.

## Scope and limitations

- The source inventory contains 3,801 conservative compound identities, but the curated classification pilot includes only 257–270 spectra depending on protocol. These are not equivalent sample sizes.
- The curated labels do not exhaust the source inventory; unknown/pending compounds must not be silently assigned to the negative class.
- Analogue-lineage folds can be small and chemically heterogeneous. Uncertainty across folds must be reported, not hidden behind a single mean.
- Spectra from curated reference libraries do not represent mixtures, instrument drift, real-world prevalence, or prospective forensic casework.
- The repository's ENFSI feasibility audit currently has only two source-backed opioid positives after its de-leakage filters. No external classifier metric should be claimed.

## Publication-safe artifact policy

Aggregate summary tables, methodology notes, source-attribution information, and figures based solely on approved aggregates can be reviewed for public inclusion. Do not commit source-derived row-level records, raw/restricted spectra, assignments, prediction CSVs, ENFSI candidate manifests, or identity review queues merely to reproduce a plot. Retain restricted artifacts locally and document how authorized users can regenerate results.

**Next requirement:** Retrieve and reconcile the original local experiment artifacts before producing the definitive 16-model/protocol comparison, exact confusion matrices, and final manuscript figures.
