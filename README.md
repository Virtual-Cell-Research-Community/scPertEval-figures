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

> **GCP is read-only here.** Every bucket access in this repo is download-only; nothing is ever
> written back to GCS.

## Downloading from the bucket

The data is in the [Requester Pays](https://cloud.google.com/storage/docs/requester-pays) bucket
`gs://scperteval`. Anyone can read it, but each download is billed to a Google Cloud project that
you provide. Anonymous requests and plain `https://storage.googleapis.com/...` links don't work.

You need a Google account, a Google Cloud project with billing enabled
([create a project](https://cloud.google.com/resource-manager/docs/creating-managing-projects),
[enable billing](https://cloud.google.com/billing/docs/how-to/modify-project)) and the
[gcloud CLI](https://cloud.google.com/sdk/docs/install). Then:

```bash
gcloud auth login
export SCPERTEVAL_BILLING_PROJECT=YOUR_PROJECT_ID
```

At [Google's rate](https://cloud.google.com/storage/pricing) of $0.12 per GiB to most destinations
(November 2026), the 7 datasets cost about $2 and each set of 4 DE files about $0.14.

## (A) Reproduce everything — the notebook

Open and run `scperteval-publication-figures.ipynb` end-to-end. It:

1. installs scPertEval from GitHub
   (`pip install "scperteval @ git+https://github.com/Virtual-Cell-Research-Community/scPertEval.git"`),
2. downloads the 7 preprocessed datasets from `gs://scperteval/processed`, billed to
   `SCPERTEVAL_BILLING_PROJECT` (see [Downloading from the bucket](#downloading-from-the-bucket)),
3. runs DRF for the 7 datasets and DE export for the 4 DEG datasets,
4. renders all four figures and displays them.

The dataset download and the scPertEval DRF run are **long-running** (Sinkhorn dominates).
Everything is deterministic (seed 42, 8192-cell subsample).

## (B) The two DEG figures, standalone

Both DEG figures run on their own with **no notebook**, once you've done the
[setup above](#downloading-from-the-bucket):

```bash
python deg_jaccard_figure.py
python deg_concordance_figure.py
```

Each script auto-downloads the per-gene DE HDF5s from `gs://scperteval/de_outputs/` over HTTPS
(the GCS JSON API for listing, a direct object URL for the download) using your gcloud credentials
and billing project, caching them in `de_cache/`. Once they're cached, re-runs need no credentials
and cost nothing. Outputs land in `figures/`. Datasets not yet present in the bucket are skipped
gracefully.

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

These auto-download the `t-test_overestim_var,MWU` DE HDF5s from the separate folder
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
