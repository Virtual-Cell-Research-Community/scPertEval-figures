"""DEG-method concordance figure, overestim variant: scanpy ``t-test_overestim_var`` vs MWU.

Identical to ``deg_concordance_figure.py`` but driven by the conservative-variance t-test
backend. It reuses ``deg_concordance_figure.render()`` against the ``t-test_overestim_var`` DE
in the Requester Pays bucket ``gs://scperteval/de_outputs_overestim/`` (written by
``scperteval de --methods t-test_overestim_var,MWU``), cached separately in
``de_cache_overestim/`` so it never collides with the standard-t-test DE.

    python deg_concordance_overestim_figure.py   # standalone: auto-downloads from the bucket
"""
from __future__ import annotations

from pathlib import Path

import deg_concordance_figure as C

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "de_cache_overestim"


def main() -> None:
    C.render(
        ttest_key="t-test_overestim_var",
        de_prefix="de_outputs_overestim",
        cache_dir=CACHE,
        local_dirs=[CACHE],
        ttest_label="t-test_ov",
        suptitle="Concordance@k: t-test_overestim_var vs MWU (Squair AUCC)",
        out_stem="deg_concordance_overestim",
    )


if __name__ == "__main__":
    main()
