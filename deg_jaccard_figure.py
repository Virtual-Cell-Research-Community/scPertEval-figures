"""DEG-method Jaccard figure from scPertEval-recomputed DE (self-contained).

Consumes the per-gene DE HDF5 files written by ``scperteval de`` (t-test + MWU, GT
first-half vs the all-perturbed pool) and renders the composite figure: per-
perturbation Spearman histograms (top), 12x12 mean-Jaccard heatmaps (middle),
and a median DEG set-size table (bottom). All criteria / Jaccard / plotting code
is vendored here -- no external project dependency.

    python deg_jaccard_figure.py        # standalone: auto-downloads DE from the public bucket

scPertEval 0.1.0 writes **one method per HDF5** (``scperteval de --method t-test``, then
``--method MWU``), so a dataset's two channels normally come from two files; the exports
already in the public bucket predate that and carry both methods in one file. ``_ensure_de``
accepts either layout.

Standalone usage requires no auth: ``_ensure_de`` lists/downloads the per-dataset DE
HDF5s from the public bucket ``gs://scperteval/de_outputs/`` over plain HTTPS and caches
them in ``de_cache/`` next to this script.

``render()`` is parameterised by the t-test backend, so the companion
``deg_jaccard_overestim_figure.py`` reuses everything here against the
``t-test_overestim_var`` DE in ``gs://scperteval/de_outputs_overestim/``.
"""
from __future__ import annotations

import glob
import json
import urllib.parse
import urllib.request
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DE_DIR = ROOT / "de_outputs"
DE_CACHE = ROOT / "de_cache"
FIG_DIR = ROOT / "figures"

# Public GCS bucket holding the precomputed per-gene DE HDF5s (read-only, no auth).
GCS_BUCKET = "scperteval"
GCS_DE_PREFIX = "de_outputs"

DATASETS = ["wessels23", "replogle22k562", "nadig25hepg2", "arch1"]


# --------------------------------------------------------------------------- #
# DE loader: local ``scperteval de`` exports, else the public bucket (no gcloud / no auth)
# --------------------------------------------------------------------------- #

def _h5_methods(path) -> set:
    """The DE-method groups an export holds (the groups beside the genes/perturbations arrays)."""
    try:
        with h5py.File(str(path), "r") as f:
            return {k for k in f if isinstance(f[k], h5py.Group)}
    except OSError as exc:
        print(f"[WARN] unreadable DE export {path} ({exc})")
        return set()


def _file_method(name) -> str | None:
    """The method tagged in a ``scperteval de`` filename, or ``None`` for a two-method export.

    scPertEval 0.1.0 names each export ``<stem>__<method>__<stamp>__de.h5``; the files already in
    the public bucket predate that naming and hold both methods in ``<stem>__<stamp>__de.h5``.
    """
    parts = Path(name).name.split("__")
    return parts[1] if len(parts) == 4 else None


def _local_de(dataset: str, methods, dirs) -> dict:
    """``{method: path}`` for the methods covered by local ``<dataset>*__de.h5`` exports.

    Directories are searched in order (newest export first within each), so callers can put a
    dedicated cache ahead of the shared ones; each method takes the first export that holds it.
    """
    found: dict[str, Path] = {}
    for d in dirs:
        for path in sorted(glob.glob(str(Path(d) / f"{dataset}*__de.h5")), reverse=True):
            have = _h5_methods(path)
            for m in methods:
                if m not in found and m in have:
                    found[m] = Path(path)
            if len(found) == len(methods):
                return found
    return found


def _download_de(dataset: str, methods, cache_dir: Path, de_prefix: str) -> list:
    """Download the bucket exports covering ``methods`` for ``dataset``; return the local paths.

    Lists the public bucket under ``de_prefix`` via the GCS JSON API and picks, per method, the
    newest export tagged with it — falling back to the newest untagged (two-method) export.
    """
    list_url = (
        f"https://storage.googleapis.com/storage/v1/b/{GCS_BUCKET}/o"
        f"?prefix={urllib.parse.quote(f'{de_prefix}/{dataset}')}"
    )
    try:
        with urllib.request.urlopen(list_url, timeout=60) as resp:
            listing = json.load(resp)
    except Exception as exc:  # network / API failure
        print(f"[WARN] {dataset}: could not list public bucket ({exc})")
        return []
    names = [o["name"] for o in listing.get("items", []) if o["name"].endswith("__de.h5")]
    untagged = sorted(n for n in names if _file_method(n) is None)
    wanted = []
    for m in methods:
        tagged = sorted(n for n in names if _file_method(n) == m)
        pick = tagged[-1] if tagged else (untagged[-1] if untagged else None)
        if pick is not None and pick not in wanted:
            wanted.append(pick)
    if not wanted:
        print(f"[INFO] {dataset}: not present in gs://{GCS_BUCKET}/{de_prefix}/")
        return []

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out = []
    for obj_name in wanted:
        dest = cache_dir / Path(obj_name).name
        if not dest.exists():
            obj_url = f"https://storage.googleapis.com/{GCS_BUCKET}/{urllib.parse.quote(obj_name)}"
            print(f"[GET]  {dataset}: downloading gs://{GCS_BUCKET}/{obj_name} -> {dest}")
            tmp = dest.with_suffix(dest.suffix + ".part")
            try:
                urllib.request.urlretrieve(obj_url, tmp)
                tmp.replace(dest)
            except Exception as exc:
                if tmp.exists():
                    tmp.unlink()
                print(f"[WARN] {dataset}: download failed ({exc})")
                continue
        out.append(dest)
    return out


def _ensure_de(dataset: str, methods=("t-test", "MWU"), cache_dir: Path = DE_CACHE,
               de_prefix: str = GCS_DE_PREFIX, local_dirs=None):
    """Return ``{method: local path}`` covering every method, downloading what is missing.

    Local exports (in ``local_dirs``, default ``cache_dir`` + ``de_outputs/``) are reused first;
    whatever they do not cover is fetched from the public bucket into ``cache_dir``. Returns
    ``None`` if some method is available neither locally nor in the bucket.

    ``de_prefix`` selects the bucket folder (``de_outputs`` for the standard t-test,
    ``de_outputs_overestim`` for the overestim-variant DE); pair it with a dedicated
    ``cache_dir``/``local_dirs`` so the two never collide (same file-naming scheme).
    """
    search = list(local_dirs) if local_dirs is not None else [Path(cache_dir), DE_DIR]
    found = _local_de(dataset, methods, search)
    missing = [m for m in methods if m not in found]
    if missing:
        for path in _download_de(dataset, missing, cache_dir, de_prefix):
            for m in _h5_methods(path) & set(missing):
                found.setdefault(m, path)
        missing = [m for m in methods if m not in found]
    if missing:
        print(f"[INFO] {dataset}: no DE for {', '.join(missing)} "
              f"(local exports or gs://{GCS_BUCKET}/{de_prefix}/)")
        return None
    return found

# 12 DEG-calling criteria: t-test |t| / padj, then MWU |delta| / padj.
TTEST_STAT_THR = [2.0, 4.0, 8.0]
MWU_CLIFF_THR = [0.15, 0.33, 0.47]
PADJ_THR = [0.1, 0.05, 0.01]
PADJ_FLOOR = 1e-300
SUPTITLE = "corr(-log10 p_adj t-test, -log10 p_adj MWU)"


# --------------------------------------------------------------------------- #
# Criteria
# --------------------------------------------------------------------------- #

def build_criteria(ttest_stat_thr, mwu_cliff_thr, padj_thr, ttest_label="t-test"):
    """Ordered (label, method, kind, threshold); method in {ttest,mwu}, kind in {stat,padj}.

    ``ttest_label`` names the t-test backend in the row labels (e.g. ``t-test`` or
    ``t-test_ov`` for the overestim variant)."""
    criteria = []
    for thr in ttest_stat_thr:
        criteria.append((f"{ttest_label}, |t| > {thr:g}", "ttest", "stat", thr))
    for thr in padj_thr:
        criteria.append((f"{ttest_label}, padj < {thr:g}", "ttest", "padj", thr))
    for thr in mwu_cliff_thr:
        criteria.append((f"MWU, |δ| > {thr:g}", "mwu", "stat", thr))
    for thr in padj_thr:
        criteria.append((f"MWU, padj < {thr:g}", "mwu", "padj", thr))
    return criteria


def _criterion_mask(criterion, ttest_stat, ttest_padj, mwu_stat, mwu_padj):
    """Boolean per-gene mask for one criterion. NaNs -> False."""
    _, method, kind, thr = criterion
    stat, padj = (ttest_stat, ttest_padj) if method == "ttest" else (mwu_stat, mwu_padj)
    with np.errstate(invalid="ignore"):
        if kind == "stat":
            return np.where(np.isnan(stat), False, np.abs(stat) > thr)
        return np.where(np.isnan(padj), False, padj < thr)


def neg_log10_padj(padj):
    return -np.log10(np.clip(padj, PADJ_FLOOR, None))


# --------------------------------------------------------------------------- #
# Per-dataset accumulation: 12x12 mean Jaccard + per-pert Spearman + set sizes
# --------------------------------------------------------------------------- #

def _jaccard(mi, mj):
    union = int(np.count_nonzero(mi | mj))
    return int(np.count_nonzero(mi & mj)) / union if union else 0.0


def _spearman_corr(arrays):
    """Per-perturbation Spearman of -log10(p_adj) between t-test and MWU."""
    from scipy.stats import spearmanr

    tt_padj, mw_padj = arrays["ttest_padj"], arrays["mwu_padj"]
    valid = ~(np.isnan(tt_padj) | np.isnan(mw_padj))
    if int(valid.sum()) < 2:
        return np.nan
    tt_nl, mw_nl = neg_log10_padj(tt_padj[valid]), neg_log10_padj(mw_padj[valid])
    if np.ptp(tt_nl) <= 0 or np.ptp(mw_nl) <= 0:
        return np.nan
    with np.errstate(invalid="ignore"):
        s = spearmanr(tt_nl, mw_nl).correlation
    return float(s) if s is not None and np.isfinite(s) else np.nan


def compute_dataset(comp, criteria):
    """Single pass over a dataset: mean Jaccard matrix, per-pert Spearman, median set sizes."""
    n = len(criteria)
    sum_jacc = np.zeros((n, n), dtype=np.float64)
    count = 0
    corr_rows, size_rows = [], []
    for pert in comp.perturbations:
        arrays = comp.read_perturbation(pert)
        if arrays is None:
            continue
        count += 1
        corr_rows.append((pert, _spearman_corr(arrays)))
        masks = [_criterion_mask(c, arrays["ttest_stat"], arrays["ttest_padj"],
                                 arrays["mwu_stat"], arrays["mwu_padj"]) for c in criteria]
        size_rows.append([int(np.count_nonzero(m)) for m in masks])
        for i in range(n):
            for j in range(i, n):
                v = _jaccard(masks[i], masks[j])
                sum_jacc[i, j] += v
                if j != i:
                    sum_jacc[j, i] += v
    mean_jacc = sum_jacc / max(count, 1)
    median_sizes = (np.median(np.asarray(size_rows, dtype=np.float64), axis=0)
                    if size_rows else np.full(n, np.nan))
    return {
        "dataset": comp.dataset_name, "reference": comp.reference,
        "labels": [c[0] for c in criteria], "n_perturbations": count,
        "mean_jaccard": mean_jacc, "median_sizes": median_sizes,
        "corr_df": pd.DataFrame(corr_rows, columns=["perturbation", "spearman"]),
    }


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #

def _hist(ax, corr_df, dataset, show_y):
    vals = corr_df["spearman"].to_numpy(dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    if vals.size:
        ax.hist(vals, bins=20, range=(0.0, 1.0), color="#4477aa", edgecolor="white", linewidth=0.4)
        ax.axvline(float(np.median(vals)), color="crimson", ls="--", lw=1.0)
    ax.set_xlim(0.0, 1.0)
    ax.set_title(dataset, fontsize=11)
    ax.set_xlabel("spearman corr", fontsize=10)
    if show_y:
        ax.set_ylabel("count", fontsize=10)
    ax.tick_params(labelsize=9)


def _heatmap(ax, mat, labels, show_y):
    im = ax.imshow(mat, cmap="YlOrRd", vmin=0.0, vmax=1.0, aspect="equal")
    n = len(labels)
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, fontsize=8, rotation=90, ha="center", va="top")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels if show_y else [], fontsize=8)
    ax.tick_params(length=0)
    return im


def _count_table(ax, labels, datasets, size_matrix):
    ax.axis("off")
    # criterion label is a real first column; the whole block is sized to content and
    # centered in the (full-width) axes.
    cell_text = [[lab] + ["-" if not np.isfinite(size_matrix[r, c]) else f"{int(round(size_matrix[r, c]))}"
                          for c in range(len(datasets))]
                 for r, lab in enumerate(labels)]
    table = ax.table(cellText=cell_text, colLabels=[""] + list(datasets),
                     cellLoc="center", loc="upper center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.auto_set_column_width(list(range(len(datasets) + 1)))
    for (r, c), cell in table.get_celld().items():
        cell.set_facecolor("#ffffff")
        cell.set_edgecolor("#dddddd")
        cell.set_linewidth(0.4)
        cell.set_height(0.07)
        cell.PAD = 0.03
        if r == 0 or c == 0:
            cell.set_text_props(fontweight="bold")
        if c == 0:
            cell.get_text().set_ha("left")
    return table


def plot_figure(results, labels, out_path, suptitle=SUPTITLE):
    """Three-band figure: histograms (top), 12x12 heatmaps (middle), set-size table (bottom)."""
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
    # histograms + heatmaps share a 2-row grid in the upper portion; the table is placed
    # below by hand so the gap under the (tall) heatmap x-labels matches the hist->heat gap.
    gs = GridSpec(2, n, figure=fig, height_ratios=[hist_h, heat_h],
                  left=left, right=right, top=0.92, bottom=0.50, hspace=0.95, wspace=0.30)

    first_hist_ax = None
    for col, res in enumerate(results):
        ax = fig.add_subplot(gs[0, col])
        _hist(ax, res["corr_df"], res["dataset"], show_y=(col == 0))
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
    # even white gap between bands: measure each band's full bottom (incl. its x-labels), and
    # reserve a fixed block above each lower band for its section title.
    title_block = 0.05
    hist_bottom = first_hist_ax.get_tightbbox(rend).transformed(inv).y0   # below the "spearman corr" label
    label_bottom = last_heat_ax.get_tightbbox(rend).transformed(inv).y0   # below the vertical heatmap labels
    gap = hist_bottom - (heat_row.y1 + title_block)                       # hist -> jaccard-title gap
    table_top = label_bottom - gap - title_block                         # same gap above the table title

    # count table: auto-sized, centered, with its top at table_top
    ax_table = fig.add_axes([0.0, 0.03, 1.0, table_top - 0.03])
    table = _count_table(ax_table, labels, datasets, size_matrix)
    fig.canvas.draw()
    celld = table.get_celld()
    ncol = len(datasets) + 1
    data_w = 1.18 * max(celld[(0, j)].get_width() for j in range(1, ncol))  # uniform dataset cols
    table._autoColumns = [0]
    for j in range(1, ncol):
        for i in range(len(labels) + 1):
            celld[(i, j)].set_width(data_w)

    if last_im is not None:
        cb_h = last_pos.height * 0.6
        cax = fig.add_axes([right + 0.015, last_pos.y0 + (last_pos.height - cb_h) / 2.0, 0.012, cb_h])
        fig.colorbar(last_im, cax=cax, label="Mean Jaccard")

    # one section title per band, same font; each lower band's title sits in its reserved
    # title_block just above the plots, so the white gap above every title is equal.
    title_kw = dict(ha="center", fontsize=14, fontweight="bold")
    sub_kw = dict(ha="center", va="top", fontsize=9)
    fig.text(xc, hist_row.y1 + 0.045, suptitle, va="bottom", **title_kw)                  # histograms
    fig.text(xc, heat_row.y1 + 0.022, "Mean Jaccard", va="bottom", **title_kw)           # heatmaps
    fig.text(xc, heat_row.y1 + 0.010, "reference: all-perturbed 8192-cell pool", **sub_kw)
    fig.text(xc, table_top + 0.022, "Median number of DEGs selected per criterion",
             va="bottom", **title_kw)                                                    # count table
    fig.text(xc, table_top + 0.010,
             "genes passing each threshold; rows match the heatmaps above", **sub_kw)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# scPertEval DE source + render/main
# --------------------------------------------------------------------------- #

def _read_de(path, method):
    """``(genes, perturbations, statistic, pvalue_adj)`` for one method group of a DE export."""
    with h5py.File(str(path), "r") as f:
        return ([g.decode() for g in f["genes"][:]],
                [p.decode() for p in f["perturbations"][:]],
                f[f"{method}/statistic"][:], f[f"{method}/pvalue_adj"][:])


class DESource:
    """Reads ``scperteval de`` exports, yielding the {ttest,mwu}_{stat,padj} arrays per pert.

    ``sources`` is the ``{method: path}`` map from :func:`_ensure_de`: with scPertEval 0.1.0 the
    two methods come from two files, while a pre-0.1 export holds both (then both keys point at
    the same path). ``ttest_key`` selects which t-test group feeds the "ttest" channel --
    ``t-test`` (default) or ``t-test_overestim_var`` for the overestim variant; the other channel
    always reads ``MWU``. Two separate exports must agree on genes, and only the perturbations
    they share are kept (in the t-test export's order).
    """

    def __init__(self, sources: dict, dataset: str, ttest_key: str = "t-test"):
        self.dataset_name = dataset
        self.reference = "all_perturbed_8192"
        t_genes, perts, self._ts, self._tp = _read_de(sources[ttest_key], ttest_key)
        m_genes, m_perts, self._ms, self._mp = _read_de(sources["MWU"], "MWU")
        if t_genes != m_genes:
            raise ValueError(f"{dataset}: the {ttest_key} and MWU exports cover different genes "
                             f"({len(t_genes)} vs {len(m_genes)}) -- recompute both from one dataset")
        if perts != m_perts:
            m_idx = {p: i for i, p in enumerate(m_perts)}
            keep = [(i, m_idx[p]) for i, p in enumerate(perts) if p in m_idx]
            ti, mi = [i for i, _ in keep], [j for _, j in keep]
            print(f"[INFO] {dataset}: {ttest_key} and MWU exports differ in perturbations; "
                  f"keeping the {len(keep)} they share")
            perts = [perts[i] for i in ti]
            self._ts, self._tp = self._ts[ti], self._tp[ti]
            self._ms, self._mp = self._ms[mi], self._mp[mi]
        self.perturbations = perts
        self._idx = {p: i for i, p in enumerate(perts)}

    def read_perturbation(self, pert):
        i = self._idx[pert]
        return {"ttest_stat": self._ts[i], "ttest_padj": self._tp[i],
                "mwu_stat": self._ms[i], "mwu_padj": self._mp[i]}


def render(ttest_key="t-test", de_prefix=GCS_DE_PREFIX, cache_dir=DE_CACHE,
           local_dirs=None, ttest_label="t-test", suptitle=SUPTITLE,
           out_stem="deg_jaccard_with_counts") -> None:
    """Render the Jaccard composite for a chosen t-test backend.

    Defaults reproduce the standard (``t-test`` vs MWU) figure. The overestim variant
    passes ``ttest_key='t-test_overestim_var'``, ``de_prefix='de_outputs_overestim'`` and a
    dedicated ``cache_dir``/``local_dirs`` (see ``deg_jaccard_overestim_figure.py``).
    """
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    criteria = build_criteria(TTEST_STAT_THR, MWU_CLIFF_THR, PADJ_THR, ttest_label=ttest_label)
    labels = [c[0] for c in criteria]
    results = []
    for ds in DATASETS:
        sources = _ensure_de(ds, (ttest_key, "MWU"), cache_dir=cache_dir,
                             de_prefix=de_prefix, local_dirs=local_dirs)
        if sources is None:
            print(f"[SKIP] {ds}: no DE available (local cache or public bucket)")
            continue
        res = compute_dataset(DESource(sources, ds, ttest_key=ttest_key), criteria)
        if res["n_perturbations"] == 0:
            print(f"[SKIP] {ds}: no perturbations")
            continue
        results.append(res)
        print(f"[OK]   {ds}: n_pert={res['n_perturbations']}")
    if not results:
        print(f"No datasets available; run `scperteval de <dataset>.h5ad --method {ttest_key}` "
              "and `--method MWU` first.")
        return
    plot_figure(results, labels, FIG_DIR / f"{out_stem}.png", suptitle=suptitle)
    rows = [(r["dataset"], lab, int(round(sz)) if np.isfinite(sz) else "")
            for r in results for lab, sz in zip(labels, r["median_sizes"])]
    pd.DataFrame(rows, columns=["dataset", "criterion", "median_set_size"]).to_csv(
        FIG_DIR / f"{out_stem}.csv", index=False)
    print(f"wrote {out_stem}.{{png,csv}} ({len(results)} datasets)")


def main() -> None:
    render()


if __name__ == "__main__":
    main()
