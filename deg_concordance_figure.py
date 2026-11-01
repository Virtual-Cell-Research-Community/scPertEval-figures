"""DEG-method concordance figure from scPertEval-recomputed DE (self-contained).

Identical to ``deg_jaccard_figure.py`` EXCEPT the top band: instead of per-
perturbation Spearman histograms, this renders 4 concordance-at-k curves (one
per dataset) using the Squair AUCC construction -- concordance@k between the
t-test and MWU top-k gene rankings. The 12x12 mean-Jaccard heatmaps (middle) and
the median DEG count table (bottom) are reused verbatim by importing from
``deg_jaccard_figure``.

    python deg_concordance_figure.py    # standalone: auto-downloads DE from the bucket

``render()`` is parameterised by the t-test backend; ``deg_concordance_overestim_figure.py``
reuses it against the ``t-test_overestim_var`` DE in ``gs://scperteval/de_outputs_overestim/``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from deg_jaccard_figure import (
    build_criteria,
    compute_dataset,
    DESource,
    _ensure_de,
    _heatmap,
    _count_table,
    DATASETS,
    DE_CACHE,
    DE_DIR,
    FIG_DIR,
    GCS_DE_PREFIX,
    latest_de,
    TTEST_STAT_THR,
    MWU_CLIFF_THR,
    PADJ_THR,
)
import deg_jaccard_figure as _D

PADJ_FLOOR = 1e-300
CONCORDANCE_SUPTITLE = "Concordance@k: t-test vs MWU (Squair AUCC)"
HIST_COLOR = "#4477aa"


# --------------------------------------------------------------------------- #
# Concordance-at-k (Squair AUCC), numpy + scipy only
# --------------------------------------------------------------------------- #

def _rank_order(neglog, stat):
    """Descending rank order: primary -log10 padj, tie-break |statistic|."""
    return np.lexsort((-np.abs(stat), -neglog))


def _k_grid(n_genes, kmax=500):
    """Doubling grid (1,2,4,...) + fractional anchors (0.1..1.0), capped at min(kmax, n)."""
    limit = int(min(kmax, n_genes))
    if limit < 1:
        return np.array([], int)
    ks = set()
    k = 1
    while k <= limit:
        ks.add(k)
        k *= 2
    for frac in (0.1, 0.25, 0.5, 0.75, 1.0):
        ks.add(max(1, int(round(frac * limit))))
    ks.add(limit)
    return np.array(sorted(ks), int)


def _pert_concordance(arrays):
    """Return (ks, conc) for one perturbation, or None if <4 valid genes."""
    tt_padj, mw_padj = arrays["ttest_padj"], arrays["mwu_padj"]
    tt_stat, mw_stat = arrays["ttest_stat"], arrays["mwu_stat"]
    valid = (np.isfinite(tt_padj) & np.isfinite(mw_padj)
             & np.isfinite(tt_stat) & np.isfinite(mw_stat))
    n_genes = int(valid.sum())
    if n_genes < 4:
        return None
    tt_nl = -np.log10(np.clip(tt_padj[valid], PADJ_FLOOR, None))
    mw_nl = -np.log10(np.clip(mw_padj[valid], PADJ_FLOOR, None))
    tt_s = np.abs(tt_stat[valid])
    mw_s = np.abs(mw_stat[valid])
    tt_order = _rank_order(tt_nl, tt_s)
    mw_order = _rank_order(mw_nl, mw_s)
    tt_rankpos = np.empty(n_genes, dtype=np.int64)
    mw_rankpos = np.empty(n_genes, dtype=np.int64)
    tt_rankpos[tt_order] = np.arange(n_genes)
    mw_rankpos[mw_order] = np.arange(n_genes)
    ks = _k_grid(n_genes)
    conc = np.empty(ks.size, dtype=np.float64)
    for idx, k in enumerate(ks):
        inter = int(np.count_nonzero((tt_rankpos < k) & (mw_rankpos < k)))
        conc[idx] = inter / k
    return ks, conc


def compute_concordance(comp):
    """Per-dataset concordance curve: median + 25th/75th percentile across perts.

    Each per-pert curve is log-interpolated onto a shared grid = _k_grid(min over
    perts of their max k) before aggregation. AUCC per pert = mean(conc).
    """
    per_pert = []          # (ks, conc)
    aucc_vals = []
    for pert in comp.perturbations:
        arrays = comp.read_perturbation(pert)
        if arrays is None:
            continue
        res = _pert_concordance(arrays)
        if res is None:
            continue
        ks, conc = res
        per_pert.append((ks, conc))
        aucc_vals.append(float(np.mean(conc)))
    if not per_pert:
        return None
    shared_max = min(int(ks.max()) for ks, _ in per_pert)
    grid = _k_grid(shared_max)
    log_grid = np.log(grid)
    interp = np.empty((len(per_pert), grid.size), dtype=np.float64)
    for r, (ks, conc) in enumerate(per_pert):
        interp[r] = np.interp(log_grid, np.log(ks), conc, left=conc[0], right=conc[-1])
    return {
        "dataset": comp.dataset_name,
        "grid": grid,
        "median": np.median(interp, axis=0),
        "p25": np.percentile(interp, 25, axis=0),
        "p75": np.percentile(interp, 75, axis=0),
        "n_perturbations": len(per_pert),
        "median_aucc": float(np.median(aucc_vals)),
    }


# --------------------------------------------------------------------------- #
# Top-row plotting (replaces _hist)
# --------------------------------------------------------------------------- #

def _concordance(ax, curve, dataset, show_y):
    grid = curve["grid"]
    ax.fill_between(grid, curve["p25"], curve["p75"], color=HIST_COLOR, alpha=0.25, linewidth=0)
    ax.plot(grid, curve["median"], color=HIST_COLOR, lw=1.4)
    ax.set_xscale("log")
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(grid[0], grid[-1])
    ax.set_title(dataset, fontsize=11)
    ax.set_xlabel("k (top-ranked genes)", fontsize=10)
    if show_y:
        ax.set_ylabel("concordance@k", fontsize=10)
    ax.tick_params(labelsize=9)
    ax.annotate(f"median AUCC = {curve['median_aucc']:.2f}", xy=(0.04, 0.06),
                xycoords="axes fraction", fontsize=7.5, ha="left", va="bottom")


# --------------------------------------------------------------------------- #
# Figure (mirrors deg_jaccard_figure.plot_figure; only top row + top title differ)
# --------------------------------------------------------------------------- #

def plot_figure(results, curves, labels, out_path, suptitle=CONCORDANCE_SUPTITLE):
    """Three-band figure: concordance curves (top), 12x12 heatmaps (middle), set-size table (bottom)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    n = len(results)
    datasets = [res["dataset"] for res in results]
    size_matrix = np.column_stack([res["median_sizes"] for res in results])

    col_w = 3.6
    fig_w = col_w * n
    left, right = 0.10, 0.90
    heat_h = (right - left) * fig_w / n          # square heatmap content
    hist_h = 1.7
    fig_h = hist_h + heat_h + 9.5
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(2, n, figure=fig, height_ratios=[hist_h, heat_h],
                  left=left, right=right, top=0.92, bottom=0.50, hspace=0.95, wspace=0.30)

    first_hist_ax = None
    for col, res in enumerate(results):
        ax = fig.add_subplot(gs[0, col])
        _concordance(ax, curves[col], res["dataset"], show_y=(col == 0))
        first_hist_ax = first_hist_ax or ax

    last_im = last_heat_ax = None
    for col, res in enumerate(results):
        ax = fig.add_subplot(gs[1, col])
        last_im = _heatmap(ax, res["mean_jaccard"], labels, show_y=(col == 0))
        last_heat_ax = ax

    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    rend = fig.canvas.get_renderer()
    hist_row = gs[0, 0].get_position(fig)
    heat_row = gs[1, 0].get_position(fig)
    last_pos = last_heat_ax.get_position()
    xc = (left + right) / 2.0
    title_block = 0.05
    hist_bottom = first_hist_ax.get_tightbbox(rend).transformed(inv).y0
    label_bottom = last_heat_ax.get_tightbbox(rend).transformed(inv).y0
    gap = hist_bottom - (heat_row.y1 + title_block)
    table_top = label_bottom - gap - title_block

    ax_table = fig.add_axes([0.0, 0.03, 1.0, table_top - 0.03])
    table = _count_table(ax_table, labels, datasets, size_matrix)
    fig.canvas.draw()
    celld = table.get_celld()
    ncol = len(datasets) + 1
    data_w = 1.18 * max(celld[(0, j)].get_width() for j in range(1, ncol))
    table._autoColumns = [0]
    for j in range(1, ncol):
        for i in range(len(labels) + 1):
            celld[(i, j)].set_width(data_w)

    if last_im is not None:
        cb_h = last_pos.height * 0.6
        cax = fig.add_axes([right + 0.015, last_pos.y0 + (last_pos.height - cb_h) / 2.0, 0.012, cb_h])
        fig.colorbar(last_im, cax=cax, label="Mean Jaccard")

    title_kw = dict(ha="center", fontsize=14, fontweight="bold")
    sub_kw = dict(ha="center", va="top", fontsize=9)
    fig.text(xc, hist_row.y1 + 0.045, suptitle, va="bottom", **title_kw)
    fig.text(xc, heat_row.y1 + 0.022, "Mean Jaccard", va="bottom", **title_kw)
    fig.text(xc, heat_row.y1 + 0.010, "reference: all-perturbed 8192-cell pool", **sub_kw)
    fig.text(xc, table_top + 0.022, "Median number of DEGs selected per criterion",
             va="bottom", **title_kw)
    fig.text(xc, table_top + 0.010,
             "genes passing each threshold; rows match the heatmaps above", **sub_kw)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def render(ttest_key="t-test", de_prefix=GCS_DE_PREFIX, cache_dir=DE_CACHE,
           local_dirs=None, ttest_label="t-test", suptitle=CONCORDANCE_SUPTITLE,
           out_stem="deg_concordance") -> None:
    """Render the concordance composite for a chosen t-test backend.

    Defaults reproduce the standard (``t-test`` vs MWU) figure; the overestim variant
    overrides ``ttest_key``/``de_prefix``/``cache_dir`` (see
    ``deg_concordance_overestim_figure.py``).
    """
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    criteria = build_criteria(TTEST_STAT_THR, MWU_CLIFF_THR, PADJ_THR, ttest_label=ttest_label)
    labels = [c[0] for c in criteria]
    results, curves = [], []
    for ds in DATASETS:
        h5 = _ensure_de(ds, cache_dir=cache_dir, de_prefix=de_prefix, local_dirs=local_dirs)
        if h5 is None:
            print(f"[SKIP] {ds}: no DE available (local cache or bucket)")
            continue
        comp = DESource(h5, ds, ttest_key=ttest_key)
        res = compute_dataset(comp, criteria)
        if res["n_perturbations"] == 0:
            print(f"[SKIP] {ds}: no perturbations")
            continue
        curve = compute_concordance(comp)
        if curve is None:
            print(f"[SKIP] {ds}: no concordance (too few valid genes)")
            continue
        results.append(res)
        curves.append(curve)
        print(f"[OK]   {ds}: n_pert={res['n_perturbations']} median_aucc={curve['median_aucc']:.3f}")
    if not results:
        print("No datasets available; run `scperteval de` first.")
        return
    plot_figure(results, curves, labels, FIG_DIR / f"{out_stem}.png", suptitle=suptitle)
    rows = [(c["dataset"], int(k), float(m), float(lo), float(hi))
            for c in curves
            for k, m, lo, hi in zip(c["grid"], c["median"], c["p25"], c["p75"])]
    pd.DataFrame(rows, columns=["dataset", "k", "median_conc", "p25", "p75"]).to_csv(
        FIG_DIR / f"{out_stem}.csv", index=False)
    print(f"wrote {out_stem}.{{png,csv}} ({len(results)} datasets)")


def main() -> None:
    render()


if __name__ == "__main__":
    main()
