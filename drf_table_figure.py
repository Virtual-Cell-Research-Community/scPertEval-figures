"""DRF table figure from scPertEval outputs, matching the DRF preprint figure's columns.

Reads the latest per-dataset DRF CSV in ``drf_outputs/`` (as written by
``scperteval calibrate --calibrator drf``), aggregates the
per-perturbation DRF to mean and median per protocol, and renders a
metric-row x dataset-column table (RdBu, signed 2-dp), one figure for each of mean
and median. The descriptor columns mirror the preprint exactly: Base Metric,
Reference Centering, Representation, Positive Control, Negative Control, Feature
Space, Weight Exponent.

    python drf_table_figure.py
"""
from __future__ import annotations

import glob
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path(__file__).resolve().parent
DRF_DIR = ROOT / "drf_outputs"
FIG_DIR = ROOT / "figures"

DATASETS = [
    ("wessels23", "Wessels23"),
    ("replogle22rpe1", "Replogle22\nRPE1"),
    ("replogle22k562", "Replogle22\nK562"),
    ("nadig25jurkat", "Nadig25\nJurkat"),
    ("nadig25hepg2", "Nadig25\nHepg2"),
    ("arch1", "ArcH1"),
    ("kaden25rpe1", "Kaden25\nRPE1"),
]

EM = "—"
APM, CTM, APC = "All-pert.\n(mean)", "Control\n(mean)", "All-pert.\n{8192 cells}"
# Columns mirror the DRF preprint figure exactly. Each row:
# (protocol, Base Metric, Reference Centering, Representation, Positive Control,
#  Negative Control, Feature Space, Weight Exponent). Protocol names are the current scPertEval
# names at their figure defaults (k=50 / padj=0.05).
ROWS = [
    ("pearson", "Pearson", EM, "Pseudobulk", "Interpolated TD", APM, EM, EM),
    ("pearson_ctrl", "Pearson", "Ctrl", "Pseudobulk", "Interpolated TD", APM, EM, EM),
    ("pearson_pert", "Pearson", "Pert", "Pseudobulk", "Interpolated TD", CTM, EM, EM),
    ("pearson_pert_degs_padj=0.05", "Pearson", "Pert", "Pseudobulk", "Interpolated TD", CTM, "Sig. DEGs", EM),
    ("pearson_pert_top_k=50", "Pearson", "Pert", "Pseudobulk", "Interpolated TD", CTM, "Top-50 DEGs", EM),
    ("mse", "MSE", EM, "Pseudobulk", "Interpolated TD", APM, EM, EM),
    ("mse_degs_padj=0.05", "MSE", EM, "Pseudobulk", "Interpolated TD", APM, "Sig. DEGs", EM),
    ("mse_top_k=50", "MSE", EM, "Pseudobulk", "Interpolated TD", APM, "Top-50 DEGs", EM),
    ("wmse_exp1", "WMSE", EM, "Pseudobulk", "Interpolated TD", APM, EM, "1.0"),
    ("wmse_exp2", "WMSE", EM, "Pseudobulk", "Interpolated TD", APM, EM, "2.0"),
    ("wmse_exp4", "WMSE", EM, "Pseudobulk", "Interpolated TD", APM, EM, "4.0"),
    ("transpose_rank", "TransposedRank", EM, "Pseudobulk", "Interpolated TD", APM, EM, EM),
    ("rank", "Rank", EM, "Pseudobulk", "Interpolated TD", APM, EM, EM),
    ("de_overlap_k=50", "DE Overlap @ 50", EM, "Single-cell", "Technical duplicate", APC, EM, EM),
    ("de_auprc", "DE AUPRC", EM, "Single-cell", "Technical duplicate", APC, EM, EM),
    ("de_auroc", "DE AUROC", EM, "Single-cell", "Technical duplicate", APC, EM, EM),
    ("unbiased_mmd_median_top_k=50", "MMD", EM, "Single-cell", "Technical duplicate", APC, "Top-50 DEGs", EM),
    ("unbiased_mmd_median_pca_k=50", "MMD", EM, "Single-cell", "Technical duplicate", APC, "PCA-50", EM),
    ("energy_distance_top_k=50", "ED", EM, "Single-cell", "Technical duplicate", APC, "Top-50 DEGs", EM),
    ("energy_distance_pca_k=50", "ED", EM, "Single-cell", "Technical duplicate", APC, "PCA-50", EM),
    ("sinkhorn_w2_top_k=50", "$W_2$ (Sinkhorn)", EM, "Single-cell", "Technical duplicate", APC, "Top-50 DEGs", EM),
    ("sinkhorn_w2_pca_k=50", "$W_2$ (Sinkhorn)", EM, "Single-cell", "Technical duplicate", APC, "PCA-50", EM),
]
DESC = ["Base\nMetric", "Reference\nCentering", "Representation", "Positive\nControl",
        "Negative\nControl", "Feature\nSpace", "Weight\nExponent"]


def latest_table(stem: str):
    """Newest ``drf_outputs/`` CSV matching ``stem`` (the filename timestamp orders the sort)."""
    cands = sorted(glob.glob(str(DRF_DIR / f"{stem}.csv")))
    return cands[-1] if cands else None


def latest_drf(ds: str):
    return latest_table(f"{ds}*__drf")


def aggregate() -> dict:
    """{dataset: {protocol: (mean_drf, median_drf)}}, mirroring scPertEval's drf calibrator."""
    out = {}
    for ds, _ in DATASETS:
        p = latest_drf(ds)
        if p is None:
            out[ds] = {}
            continue
        g = pd.read_csv(p).groupby("protocol")["drf"]
        out[ds] = {proto: (float(np.nanmean(v)), float(np.nanmedian(v))) for proto, v in g}
    return out


def render(agg: dict, stat_idx: int, path: Path) -> None:
    cmap, norm = plt.cm.RdBu, Normalize(vmin=-1.0, vmax=1.0)
    headers = DESC + [lbl for _, lbl in DATASETS]
    ncol = len(headers)
    cells, colors = [], []
    for proto, *desc in ROWS:
        vals, vcolors = [], []
        for ds, _ in DATASETS:
            mm = agg.get(ds, {}).get(proto)
            v = mm[stat_idx] if mm is not None else float("nan")
            vals.append("-" if not np.isfinite(v) else f"{v:+.2f}")
            vcolors.append("#f0f0f0" if not np.isfinite(v) else cmap(norm(v)))
        cells.append(list(desc) + vals)
        colors.append(["white"] * len(DESC) + vcolors)

    fig, ax = plt.subplots(figsize=(0.85 * ncol + 3.0, 0.55 * len(ROWS) + 1.6))
    ax.axis("off")
    tbl = ax.table(cellText=cells, colLabels=headers, cellColours=colors, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 2.0)   # tall rows: the two-line headers and Negative Control cells need two lines
    # Fit each column to its text, then give every dataset column a single uniform width
    # (= the widest dataset column) so all the DRF data cells are the same size.
    tbl.auto_set_column_width(list(range(ncol)))
    fig.canvas.draw()
    celld = tbl.get_celld()
    ds_w = max(celld[(0, j)].get_width() for j in range(len(DESC), ncol))
    tbl._autoColumns = list(range(len(DESC)))   # keep descriptor cols auto-fit; freeze dataset cols
    for j in range(len(DESC), ncol):
        for i in range(len(ROWS) + 1):
            celld[(i, j)].set_width(ds_w)
    for j in range(ncol):
        c = tbl[0, j]
        c.set_text_props(weight="bold")
        c.set_facecolor("#e8e8e8")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(cells, columns=[h.replace("\n", " ") for h in headers]).to_csv(
        str(path).replace(".png", ".csv"), index=False)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    agg = aggregate()
    render(agg, 0, FIG_DIR / "drf_table_mean.png")
    render(agg, 1, FIG_DIR / "drf_table_median.png")
    have = [ds for ds, _ in DATASETS if latest_drf(ds)]
    print(f"wrote drf_table_{{mean,median}}.{{png,csv}} ; datasets present: {have}")


if __name__ == "__main__":
    main()
