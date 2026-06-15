# Why the DEG concordance is low for some datasets (notably wessels23)

The DEG-concordance figure measures, per perturbation, how much two differential-expression
methods — a **Welch t-test** and the **Mann–Whitney U (MWU)** test — agree on *which genes
are most affected*. For each `k` it computes the overlap of the top-`k` genes ranked by each
method (genes ranked by `−log10(p_adj)`, tie-broken by `|statistic|`), and summarizes the
curve as the **AUCC** (Area Under the Concordance Curve, the Squair et al. 2021 construction).
AUCC = 1 means the two methods rank the same top genes; AUCC ≈ 0 means their top genes are
disjoint.

## What we observe

t-test ↔ MWU concordance (median AUCC across perturbations), GT-half vs the 8192-cell
all-perturbed pool:

| dataset | n perts | median AUCC | % perts AUCC > 0.3 |
|---|---:|---:|---:|
| arch1 | 150 | **0.66** | 83% |
| replogle22k562 | 1971 | **0.34** | 57% |
| nadig25hepg2 | 1818 | 0.10 | 10% |
| wessels23 | 157 | **0.002** | 5% |

So the metric behaves sensibly on strong datasets (arch1, replogle) and degrades with signal
strength. **wessels23 is the outlier**: for ~90% of its perturbations the two methods'
top-gene sets are essentially non-overlapping.

## Why wessels23 is low — the real reasons

1. **It is the weakest dataset.** wessels23 is a combinatorial Cas13 screen with small effect
   sizes and few cells per perturbation (~75–155 cells in the GT half). With little signal,
   each test is largely ranking noise, and two *different* tests rank noise differently.

2. **The two tests genuinely emphasize different genes.** The parametric Welch t-test ranks by
   `|t| = standardized mean shift`, which — against the large 8192-cell pool — is dominated by
   high-expression / low-variance genes that deviate from the pool baseline (largely
   *non-specific* to the perturbation). MWU (rank-based, Cliff's δ) surfaces genes with a
   consistent rank shift. On strong perturbations these coincide; on weak ones they don't.

3. **It is *below* random chance, which is diagnostic.** Two random top-500 sets out of
   ~16,775 genes would overlap ~15 genes (concordance ≈ 0.03). The typical wessels23
   perturbation sits at ~0 — the top sets are *actively disjoint*, not merely uncorrelated.
   That is the fingerprint of the t-test top list being **non-specific** (similar across
   perturbations, baseline-driven) while MWU recovers the actual perturbed genes.

## What it is NOT (ruled out)

- **Not a scPertEval implementation artifact.** scPertEval's Welch t-test is the *identical
  computation* as scanpy's standard `t-test` (verified: mutual top-gene concordance AUCC =
  **1.000**, i.e. the two rank genes identically). Running scanpy's standard t-test reproduces
  the same ~0.002 concordance for wessels23.

- **Not a p-value clipping artifact.** Significance is ranked with a `1e-300` floor only to
  avoid `−log10(0)`; it does not drop genes. The thousands of genes at `p_adj = 0` are tied
  *because the t-test p-values underflowed to exactly 0 at compute time* (`scipy.stats.t.sf`
  for large `|t|`), and they are correctly sub-ordered by `|t|`. We confirmed `scipy.stats.t.logsf`
  *also* underflows for the t-distribution at large `|t|`, so there is no higher-resolution
  log-space p-value to recover; the only robust alternative (a normal approximation) reduces to
  ranking by `|t|`, which is already what the tie-break does. So the saturated genes are ordered
  by their true significance, and the low concordance is real.

- **Not the same as the older drf figure (0.538 for wessels23).** That figure used scanpy's
  `t-test_overestim_var` — a deliberately conservative-variance variant that produces smaller
  `|t|`, saturates far less, and happens to rank more like MWU. That is a *different test*, not
  a difference in clipping or precision.

## Conclusion

The low wessels23 concordance is a faithful result: on this weak combinatorial screen the
parametric t-test and the rank-based MWU genuinely disagree about the top genes — the t-test's
ranking is dominated by non-specific high-expression genes relative to the all-perturbed pool,
while MWU recovers perturbation-specific genes. The metric and the implementation are correct;
wessels23 is simply the dataset where the two DE methods diverge most.
