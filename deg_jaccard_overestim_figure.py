"""DEG-method Jaccard figure, overestim variant: scanpy ``t-test_overestim_var`` vs MWU.

Identical to ``deg_jaccard_figure.py`` but driven by the conservative-variance t-test
backend. It reuses ``deg_jaccard_figure.render()`` against the ``t-test_overestim_var`` DE in
the Requester Pays bucket ``gs://scperteval/de_outputs_overestim/`` (written by
``scperteval de --methods t-test_overestim_var,MWU``), cached separately in
``de_cache_overestim/`` so it never collides with the standard-t-test DE.

    python deg_jaccard_overestim_figure.py   # standalone: auto-downloads from the bucket
"""
from __future__ import annotations

from pathlib import Path

import deg_jaccard_figure as J

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "de_cache_overestim"


def main() -> None:
    J.render(
        ttest_key="t-test_overestim_var",
        de_prefix="de_outputs_overestim",
        cache_dir=CACHE,
        local_dirs=[CACHE],
        ttest_label="t-test_ov",
        suptitle="corr(-log10 p_adj t-test_overestim_var, -log10 p_adj MWU)",
        out_stem="deg_jaccard_overestim_with_counts",
    )


if __name__ == "__main__":
    main()
