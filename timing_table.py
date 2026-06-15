"""Per-metric wall-clock timing table from EPPS ``--profile`` outputs.

Reads the latest per-dataset timing CSV in ``drf_outputs/`` and renders a table
with the same metric rows/order as the DRF table figure, one wall-clock column per
dataset plus a per-row Total and a per-dataset TOTAL row.

    python timing_table.py
"""
from __future__ import annotations

import glob
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from drf_table_figure import DATASETS, DESC, DRF_DIR, FIG_DIR, ROWS


def _is_dark(rgba, thresh=0.55):
    """White text reads better on a dark gray cell (relative luminance of the face)."""
    lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return lum < (1.0 - thresh)


def latest_timing(ds: str):
    cands = sorted(glob.glob(str(DRF_DIR / f"{ds}*__timing.csv")))
    return cands[-1] if cands else None


def collect() -> dict:
    """{dataset: {protocol: seconds}}"""
    out = {}
    for ds, _ in DATASETS:
        p = latest_timing(ds)
        if p is None:
            out[ds] = {}
            continue
        t = pd.read_csv(p)
        out[ds] = dict(zip(t.protocol, t.seconds))
    return out


def render(tim: dict, path: Path) -> None:
    headers = DESC + [lbl for _, lbl in DATASETS] + ["TOTAL (s)"]
    cells, cell_vals = [], []   # cell_vals: numeric seconds parallel to cells (NaN where "-")
    for proto, *desc in ROWS:
        secs = [tim.get(ds, {}).get(proto, float("nan")) for ds, _ in DATASETS]
        tot = np.nansum(secs) if np.any(np.isfinite(secs)) else float("nan")
        cells.append(list(desc)
                     + ["-" if not np.isfinite(s) else f"{s:.1f}" for s in secs]
                     + ["-" if not np.isfinite(tot) else f"{tot:.1f}"])
        cell_vals.append([float("nan")] * len(DESC) + list(secs) + [tot])
    # per-dataset TOTAL row
    totals = ["TOTAL (s)"] + [""] * (len(DESC) - 1)
    totals_vals = [float("nan")] * len(DESC)
    grand = 0.0
    for ds, _ in DATASETS:
        s = float(np.nansum(list(tim.get(ds, {}).values()))) if tim.get(ds) else float("nan")
        grand += 0.0 if not np.isfinite(s) else s
        totals.append("-" if not np.isfinite(s) else f"{s:.0f}")
        totals_vals.append(s)
    totals.append(f"{grand:.0f}")
    totals_vals.append(grand)
    cells.append(totals)
    cell_vals.append(totals_vals)
    cell_vals = np.array(cell_vals, dtype=float)

    ncol, nrow = len(headers), len(cells)
    d0, d1 = len(DESC), len(DESC) + len(DATASETS)   # per-dataset cols [d0, d1)
    total_col = ncol - 1

    # two separate grayscale (log) scales: per-cell times, and the much-larger TOTAL margins
    def _lognorm(mask):
        v = cell_vals[mask]
        v = v[np.isfinite(v) & (v > 0)]
        return LogNorm(vmin=v.min(), vmax=v.max())

    cell_mask = np.zeros_like(cell_vals, dtype=bool)
    cell_mask[: nrow - 1, d0:d1] = True             # metric rows × dataset cols
    total_mask = np.zeros_like(cell_vals, dtype=bool)
    total_mask[:, total_col] = True                 # right TOTAL column
    total_mask[nrow - 1, d0:total_col + 1] = True    # bottom TOTAL row
    cell_norm, total_norm = _lognorm(cell_mask), _lognorm(total_mask)
    cmap = plt.cm.Greys

    fig, ax = plt.subplots(figsize=(0.85 * ncol + 3.0, 0.55 * nrow + 1.6))
    ax.axis("off")
    tbl = ax.table(cellText=cells, colLabels=headers, cellLoc="center", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 2.0)   # tall rows for the two-line headers and Negative Control cells
    # Fit descriptor columns to text; give the numeric columns (datasets + Total) one uniform width.
    tbl.auto_set_column_width(list(range(ncol)))
    fig.canvas.draw()
    celld = tbl.get_celld()
    num_w = max(celld[(0, j)].get_width() for j in range(len(DESC), ncol))
    tbl._autoColumns = list(range(len(DESC)))
    for j in range(len(DESC), ncol):
        for i in range(nrow + 1):
            celld[(i, j)].set_width(num_w)
    for j in range(ncol):                            # bold gray header row
        c = tbl[0, j]
        c.set_text_props(weight="bold")
        c.set_facecolor("#e8e8e8")

    def shade(i, j, v, norm, bold=False):
        c = tbl[i, j]
        if bold:
            c.set_text_props(weight="bold")
        if np.isfinite(v) and v > 0:
            face = cmap(norm(v))
            c.set_facecolor(face)
            if _is_dark(face):
                c.set_text_props(weight="bold" if bold else "normal", color="white")
        elif bold:
            c.set_facecolor("#f3f3f3")

    # per-cell grayscale (cell_vals row r -> table row r+1)
    for r in range(nrow - 1):
        for j in range(d0, d1):
            shade(r + 1, j, cell_vals[r, j], cell_norm)
    # TOTAL margins on their own scale: right column then bottom row (bold)
    for r in range(nrow):
        shade(r + 1, total_col, cell_vals[r, total_col], total_norm, bold=True)
    for j in range(d0, total_col + 1):
        shade(nrow, j, cell_vals[nrow - 1, j], total_norm, bold=True)
    for j in range(d0):                              # bottom TOTAL row descriptor cells
        c = tbl[nrow, j]
        c.set_text_props(weight="bold")
        c.set_facecolor("#f3f3f3")

    fig.canvas.draw()   # title placed just above the table's actual top (the table overflows the axes)
    bb = tbl.get_window_extent(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text((bb.x0 + bb.x1) / 2, bb.y1 + 0.012, "Computation Time (Seconds)",
             ha="center", va="bottom", fontsize=15, fontweight="bold")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(cells, columns=[h.replace("\n", " ") for h in headers]).to_csv(
        str(path).replace(".png", ".csv"), index=False)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    render(collect(), FIG_DIR / "metric_timing_table.png")
    have = [ds for ds, _ in DATASETS if latest_timing(ds)]
    print(f"wrote metric_timing_table.{{png,csv}} ; datasets present: {have}")


if __name__ == "__main__":
    main()
