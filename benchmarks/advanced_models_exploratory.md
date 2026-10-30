# Exploratory advanced-model experiment (not a frozen benchmark)

This branch adds **two predeclared nonlinear baselines**: an RBF-kernel support
vector machine and an ExtraTrees ensemble. Both use the existing EI-MS
preprocessing (nominal m/z bins, square-root intensity, L2 normalization) and
the *existing frozen fold assignments*. No hyperparameters are chosen on test
folds. Both models are evaluated with the same metrics as the classical pilot.

The additional nearest-training-spectrum cosine similarity is a descriptive
novelty proxy, **not** calibrated uncertainty, validated open-set recognition,
or a deployable rejection rule. A genuine OOD/abstention benchmark requires
separate unknown-class data, prespecified operating points, and train-only
calibration.

## Run locally with licensed source data

Use a local checkout of this branch with the user-managed normalized SWGDRUG
dataset and compounds/taxonomy tables. Do not upload or commit the source
spectra, row-level predictions, assignments, or restricted ENFSI data.

```python
from pathlib import Path
from opioidbench.experiments.advanced import evaluate_advanced
from opioidbench.experiments.splits import load_curated_samples, load_fold_assignments

dataset = Path("data/processed/normalized_records.jsonl")  # adjust to local file
compounds = Path("metadata/compounds.csv")
taxonomy = Path("metadata/opioid_taxonomy.csv")
for protocol in ("compound", "scaffold", "analogue"):
    samples, audit = load_curated_samples(
        dataset,
        compounds,
        taxonomy,
        require_scaffold=protocol in ("scaffold", "analogue"),
    )
    assignments_path = Path("results") / f"{protocol}_cv5_assignments.csv"
    assignments = load_fold_assignments(assignments_path)
    report = evaluate_advanced(
        samples,
        assignments,
        output_dir=Path("results/advanced") / protocol,
        mz_max=700,
        seed=20260929,
        input_hashes={
            "dataset": audit["dataset_sha256"],
            "compounds": audit["compounds_sha256"],
            "taxonomy": audit["taxonomy_sha256"],
        },
    )
    print(protocol, report["aggregate"])
```

The example file locations are **illustrative**; point them to the actual
local dataset and frozen assignment paths. Use the exact same cohort and
assignment SHA-256 as the reference run. Keep all generated row-level outputs
local and git-ignored.

## Scientific gate

1. Reconcile the unresolved frozen **m/z 700 scaffold** provenance mismatch
   before claiming a performance improvement over v0.4.0.
2. Confirm exact input hashes, sample counts, folds, class balance, and leakage
   audit for each run.
3. Compare mean balanced accuracy, AUROC, fold dispersion, sensitivity, and
   specificity on the **analogue-lineage** protocol first.
4. Report all three protocols and both new models, including unfavorable
   results; no cherry-picking or post-hoc test-set tuning.
5. Do not claim external validation or clinical/forensic deployment readiness.

This addition does **not** change the published frozen pilot table or taxonomy.
