"""DEG-method concordance figure, overestim variant: scanpy ``t-test_overestim_var`` vs MWU.

Identical to ``deg_concordance_figure.py`` but driven by the conservative-variance t-test
backend. It reuses ``deg_concordance_figure.render()`` against the ``t-test_overestim_var`` DE
in the public bucket ``gs://scperteval/de_outputs_overestim/`` (written by
``scperteval de --method t-test_overestim_var``), cached separately in
``de_cache_overestim/`` so it never collides with the standard-t-test DE.

    python deg_concordance_overestim_figure.py   # standalone: auto-downloads from the public bucket
"""
from __future__ import annotations

import deg_concordance_figure as C
from deg_jaccard_overestim_figure import CACHE, LOCAL_DIRS   # same cache + search order


def main() -> None:
    C.render(
        ttest_key="t-test_overestim_var",
        de_prefix="de_outputs_overestim",
        cache_dir=CACHE,
        local_dirs=LOCAL_DIRS,
        ttest_label="t-test_ov",
        suptitle="Concordance@k: t-test_overestim_var vs MWU (Squair AUCC)",
        out_stem="deg_concordance_overestim",
    )


if __name__ == "__main__":
    main()
