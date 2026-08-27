"""reports/collateral_lethality_<indication>.pdf - the standalone G2c figure.

Generalises the POLR2A logic: a payload TARGET is a gene that is both recurrently
lost in the indication (TCGA) and more essential in the CN-loss lines (DepMap).
POLR2A is the mandatory positive control - it must recover before any novel hit
is trusted.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..config import load_gates  # noqa: E402

RRB_MAROON = "#6E1426"
_CAND = "#2C7FB8"
_NOT = "#9AA0A6"


def render_collateral_report(
    targets: pd.DataFrame,
    chronos: pd.DataFrame,
    copy_number: pd.DataFrame,
    *,
    indication: str,
    out_path: str | Path,
    positive_control: str = "POLR2A",
    config: dict | None = None,
) -> Path:
    cfg = config or load_gates()
    ccfg = cfg["g2"]["collateral"]
    thr = ccfg["cn_loss_threshold"]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle(f"WARHEAD - G2c Collateral Lethality ({indication})", color=RRB_MAROON,
                 fontsize=16, fontweight="bold", x=0.5, y=0.97)
    fig.text(0.5, 0.935,
             "Recurrent CN loss (TCGA) x stronger dependency on loss (DepMap)  ->  payload targets",
             ha="center", fontsize=10, color="#444")

    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.42, wspace=0.28,
                          left=0.12, right=0.96, top=0.90, bottom=0.10)

    # Panel A: target ranking by collateral score.
    axA = fig.add_subplot(gs[0, :])
    t = targets.sort_values("collateral_score").reset_index(drop=True)
    colors = [_CAND if c else _NOT for c in t["candidate"]]
    ypos = np.arange(len(t))
    axA.barh(ypos, t["collateral_score"], color=colors, edgecolor="#222", linewidth=0.4)
    axA.set_yticks(ypos)
    axA.set_yticklabels(t["gene"], fontsize=9)
    axA.set_xlabel("collateral score  =  -(Chronos shift on loss) x TCGA loss frequency")
    axA.set_title("A.  Candidate payload targets (blue = passes G2c)", fontsize=9, loc="left")

    # Panel B: POLR2A positive-control dependency distributions.
    axB = fig.add_subplot(gs[1, 0])
    g = chronos[chronos["gene"] == positive_control].merge(
        copy_number[copy_number["gene"] == positive_control][["ModelID", "cn_log2"]],
        on="ModelID", how="inner",
    )
    loss = g[g["cn_log2"] < thr]["chronos"].to_numpy()
    neutral = g[g["cn_log2"] >= thr]["chronos"].to_numpy()
    for i, (vals, lab, col) in enumerate([(neutral, "CN-neutral", _NOT), (loss, "CN-loss", RRB_MAROON)]):
        x = np.random.default_rng(0).normal(i, 0.06, size=vals.size)
        axB.scatter(x, vals, s=16, color=col, alpha=0.7, edgecolor="none")
        axB.hlines(np.median(vals), i - 0.2, i + 0.2, color="#222", lw=2)
    axB.set_xticks([0, 1])
    axB.set_xticklabels(["CN-neutral", "CN-loss"])
    axB.set_ylabel("Chronos dependency")
    axB.set_title(f"B.  Positive control {positive_control}: loss -> more dependent", fontsize=9, loc="left")
    axB.axhline(0, color="#bbb", lw=0.7, ls=":")

    # Panel C: volcano - dependency shift vs significance.
    axC = fig.add_subplot(gs[1, 1])
    q = targets["q"].clip(lower=1e-12)
    axC.scatter(targets["delta"], -np.log10(q), s=28,
                color=[_CAND if c else _NOT for c in targets["candidate"]],
                edgecolor="#222", linewidth=0.4, zorder=3)
    for _, r in targets.iterrows():
        if r["gene"] in ("POLR2A", "ME2"):
            axC.annotate(r["gene"], (r["delta"], -np.log10(max(r["q"], 1e-12))),
                         fontsize=8, color=RRB_MAROON, xytext=(4, 2), textcoords="offset points")
    axC.axhline(-np.log10(ccfg["fdr_alpha"]), color="#888", lw=0.8, ls="--")
    axC.set_xlabel("Chronos shift on loss (delta; <0 = more dependent)")
    axC.set_ylabel("-log10 q")
    axC.set_title("C.  Dependency shift vs significance", fontsize=9, loc="left")

    n_cand = int(targets["candidate"].sum())
    fig.text(0.12, 0.02,
             f"genes tested = {len(targets)}   |   candidate targets = {n_cand}"
             f"   |   positive control {positive_control} recovered = "
             f"{bool(targets.set_index('gene').loc[positive_control, 'candidate']) if positive_control in set(targets['gene']) else False}",
             fontsize=8, color="#555")

    fmt = out_path.suffix.lstrip(".").lower() or "pdf"
    fig.savefig(out_path, format=fmt, dpi=150)
    plt.close(fig)
    return out_path
