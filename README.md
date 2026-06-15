# scPertEval publication figures

Reproduction of the figures from the **scPertEval** preprint:

- **DRF table** — distributional-rank-fidelity scores across the 7 benchmark datasets
  (`drf_table_figure.py` → `drf_table_mean.png`, `drf_table_median.png`).
- **Metric timing table** — per-metric wall-clock cost (`timing_table.py` →
  `metric_timing_table.png`).
- **DEG-Jaccard composite** — per-perturbation t-test vs MWU Spearman histograms, 12×12
  mean-Jaccard heatmaps, and median DEG set-size table
  (`deg_jaccard_figure.py` → `deg_jaccard_with_counts.png`).
- **DEG-concordance composite** — concordance@k (Squair AUCC) curves over the same heatmaps
  and count table (`deg_concordance_figure.py` → `deg_concordance.png`).
- **Overestim DEG variants** — the same two DEG composites, but comparing scanpy's
  conservative `t-test_overestim_var` against MWU instead of the standard Welch t-test
  (`deg_jaccard_overestim_figure.py`, `deg_concordance_overestim_figure.py` →
  `deg_jaccard_overestim_with_counts.png`, `deg_concordance_overestim.png`). See section (C).

The scPertEval code lives in a separate repo:
<https://github.com/Virtual-Cell-Research-Community/scPertEval>.

> **GCP is read-only here.** Every bucket access in this repo is download-only from public
> buckets; nothing is ever written back to GCS.

## (A) Reproduce everything — the notebook

Open and run `scperteval-publication-figures.ipynb` end-to-end. It:

1. installs scPertEval from GitHub
   (`pip install "scperteval @ git+https://github.com/Virtual-Cell-Research-Community/scPertEval.git"`),
2. downloads the 7 preprocessed datasets from the public bucket `gs://scperteval/processed`
   (needs the gcloud SDK / `gsutil` on `PATH`),
3. runs DRF for the 7 datasets and DE export for the 4 DEG datasets,
4. renders all four figures and displays them.

The dataset download and the scPertEval DRF run are **long-running** (Sinkhorn dominates).
Everything is deterministic (seed 42, 8192-cell subsample).

## (B) The two DEG figures, standalone

Both DEG figures run on their own with **no notebook and no gcloud auth**:

```bash
python deg_jaccard_figure.py
python deg_concordance_figure.py
```

Each script auto-downloads the per-gene DE HDF5s from the public bucket
`gs://scperteval/de_outputs/` over plain HTTPS (the GCS JSON API for listing, a direct object
URL for the download), caching them in `de_cache/`. No `gsutil`, no credentials. Outputs land
in `figures/`. Datasets not yet present in the bucket are skipped gracefully.

Requires Python with `h5py`, `numpy`, `pandas`, `scipy`, and `matplotlib` (all installed by the
scPertEval install in the notebook path).

## (C) DEG figures with the conservative t-test (`t-test_overestim_var` vs MWU)

scanpy's conservative-variance t-test (`t-test_overestim_var`) is a selectable scPertEval DE
backend, and the two DEG composites have **overestim variants** that compare it against MWU
instead of the standard Welch t-test:

```bash
python deg_jaccard_overestim_figure.py        # -> figures/deg_jaccard_overestim_with_counts.png
python deg_concordance_overestim_figure.py    # -> figures/deg_concordance_overestim.png
```

These auto-download the `t-test_overestim_var,MWU` DE HDF5s from the separate public folder
`gs://scperteval/de_outputs_overestim/` (cached in `de_cache_overestim/`), so they never
collide with the standard-t-test DE in `de_outputs/`. Criteria, heatmaps, and the count table
are identical to the standard figures; only the t-test backend differs (its rows are tagged
`t-test_ov`). The overestim t-test is deliberately more conservative — its `|t|` is smaller —
so the `|t|`-threshold rows select fewer genes than the standard variant.

To regenerate the underlying DE from scratch (then upload to your own bucket folder):

```bash
scperteval de <dataset>_processed_complete.h5ad --methods t-test_overestim_var,MWU \
    --subsample 8192 --seed 42
```

## Why is the wessels23 DEG concordance so low?

See [`why_concordance_is_low.md`](why_concordance_is_low.md) — the t-test vs MWU concordance
behaves sensibly on the strong datasets (arch1 0.66, replogle 0.34) but is near-zero on
wessels23. That writeup documents the per-dataset AUCC, the root cause (a weak combinatorial
Cas13 screen where the parametric t-test and rank-based MWU genuinely disagree), and the
artifact checks ruling out an implementation or precision bug.

## License

MIT — see `LICENSE`.
