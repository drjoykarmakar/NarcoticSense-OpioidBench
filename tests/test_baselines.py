from pathlib import Path

from opioidbench.experiments.baselines import evaluate_baselines
from opioidbench.experiments.splits import BenchmarkSample


def test_baselines_run_on_small_separable_fixture(tmp_path: Path):
    samples = []
    folds = {}
    for i in range(20):
        label = 1 if i % 2 else 0
        peak = 200.0 if label else 50.0
        sample = BenchmarkSample(
            spectrum_id=f"sp{i}",
            compound_id=f"c{i}",
            label=label,
            label_name="opioid" if label else "non-opioid",
            scaffold_id=f"s{i}",
            structure_conflict=False,
            peaks=((peak, 100.0),),
        )
        samples.append(sample)
        folds[sample.spectrum_id] = i % 5
    report = evaluate_baselines(samples, folds, output_dir=tmp_path / "bench", seed=3)
    assert report["n_samples"] == 20
    assert report["aggregate"]["spectral_1nn"]["balanced_accuracy_mean"] == 1.0
