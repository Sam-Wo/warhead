"""Presentation figures for the exatecan (Top1i) combination-partner analysis."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

RRB = "#6E1426"
_CLASSES = [
    ("HSP90", ["HSP90"], "#6E1426"),
    ("MAPK / MEK / RAF", ["MAP2K", "MEK", "BRAF", "RAF1", "ERK", "MAPK"], "#C85A54"),
    ("PI3K / mTOR / ATR (DDR)", ["MTOR", "PIK3", "PI3K", "ATR", "AKT", "CHEK", "WEE1"], "#2E7D6B"),
    ("Proteasome", ["PSM", "PROTEAS"], "#7D5BA6"),
    ("HDAC / epigenetic", ["HDAC", "BRD", "EZH"], "#2C82C9"),
    ("Mitotic (tubulin/KIF/AURK/PLK)", ["TUBB", "TUBA", "KIF11", "AURK", "PLK", "BIRC5"], "#E08A1E"),
    ("Transcription / export", ["RELA", "XPO1", "POLR", "CDK9"], "#16A085"),
]
_OTHER = ("Other", "#AEB4BC")


def moa_class(target, moa=""):
    t = (str(target) + " " + str(moa)).upper()
    for name, keys, col in _CLASSES:
        if any(k in t for k in keys):
            return name, col
    return _OTHER


def render_combo_consensus(merged: pd.DataFrame, *, out_path, label_top=15) -> Path:
    """merged: one row per compound scored in BOTH screens, columns compound,
    target, combo_score_p (PRISM), combo_score_c (CTRP)."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    d = merged.copy()
    d["consensus"] = d["combo_score_p"] + d["combo_score_c"]
    d[["cls", "col"]] = d.apply(lambda r: pd.Series(moa_class(r["target"], r.get("moa", ""))), axis=1)

    fig, ax = plt.subplots(figsize=(11.5, 8.4))
    ax.axhline(0, color="#ddd", lw=1, zorder=0); ax.axvline(0, color="#ddd", lw=1, zorder=0)
    # quadrant wash: top-right = robust in both
    xmax, ymax = d["combo_score_p"].max() + .6, d["combo_score_c"].max() + .6
    xmin, ymin = d["combo_score_p"].min() - .4, d["combo_score_c"].min() - .4
    ax.axhspan(0, ymax, xmin=0.5, color="#f6eef0", zorder=0, alpha=.5)
    ax.scatter(d["combo_score_p"], d["combo_score_c"], s=26, c=d["col"], alpha=.55,
               edgecolor="white", linewidth=.4, zorder=2)

    top = d.sort_values("consensus", ascending=False).head(label_top)
    ax.scatter(top["combo_score_p"], top["combo_score_c"], s=95, c=top["col"],
               edgecolor="#222", linewidth=.7, zorder=4)
    for _, r in top.iterrows():
        ax.annotate(r["compound"], (r["combo_score_p"], r["combo_score_c"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=8.4, fontweight="bold",
                    color="#222", zorder=5)

    ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
    ax.set_xlabel("PRISM  combo score  (complementarity to Top1i)", fontsize=11)
    ax.set_ylabel("CTRP v2  combo score", fontsize=11)
    ax.set_title("Exatecan (Top1i) combination partners - reproducible across two screens",
                 color=RRB, fontsize=14, fontweight="bold", loc="left", pad=12)
    ax.text(0.985, 0.02, "combo score = potent on Top1i-RESISTANT lines  +  orthogonal response pattern.\n"
            "top-right = robust in both independent screens.  anchor = exatecan / SN-38 / irinotecan / "
            "topotecan / camptothecin (z-scored consensus).",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.8, color="#666")
    seen = list(dict.fromkeys(top["cls"]))
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=dict(
        [(n, c) for n, _, c in _CLASSES] + [_OTHER])[c], markersize=9, label=c) for c in seen]
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, frameon=False, title="mechanism",
              title_fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=160)
    plt.close(fig)
    return out_path


def render_two_paths(df: pd.DataFrame, *, out_path, wetlab=None) -> Path:
    """The narrative slide. df: compound, target, combo_score (orthogonality / coverage),
    ddr, arrest (mechanism scores). x = coverage (Path 1); y = mechanism axis (arrest
    below = antagonism risk, DDR above = potentiation/synergy). `wetlab`: optional dict
    {compound_norm: 'synergy'|'antagonism'} to overlay measured outcomes as ring markers."""
    import re
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    d = df.copy()
    d["y"] = d["ddr"] - d["arrest"]
    rng = np.random.RandomState(3)
    d["yj"] = d["y"] + rng.uniform(-.28, .28, len(d))

    def _c(r):
        if r["ddr"] > 0 and r["ddr"] >= r["arrest"]:
            return "DDR / checkpoint (synergy-plausible)", "#2E7D6B"
        if r["arrest"] > 0:
            return "cell-cycle arrest (antagonism risk)", "#C0392B"
        return "other", "#C9CDD3"
    d[["cls", "col"]] = d.apply(lambda r: pd.Series(_c(r)), axis=1)

    fig, ax = plt.subplots(figsize=(12, 8.2))
    ax.axhspan(0.5, d["y"].max() + 1, color="#e9f3ef", zorder=0)
    ax.axhspan(d["y"].min() - 1, -0.5, color="#fbecea", zorder=0)
    ax.axhline(0, color="#ccc", lw=1)
    for kind, sub in d.groupby("cls"):
        ax.scatter(sub["combo_score"], sub["yj"], s=np.where(sub["cls"] == "other", 14, 46),
                   c=sub["col"], alpha=.6 if kind == "other" else .85,
                   edgecolor="white", linewidth=.3, zorder=2)
    # label the exemplars: top DDR (by ddr then coverage) and high-coverage arrest
    lab = pd.concat([d[d["ddr"] > 0].sort_values(["ddr", "combo_score"], ascending=False).head(9),
                     d[(d["arrest"] > 0) & (d["combo_score"] > 2)].sort_values("combo_score", ascending=False).head(8)])
    for _, r in lab.drop_duplicates("compound").iterrows():
        ax.annotate(r["compound"], (r["combo_score"], r["yj"]), xytext=(4, 3), textcoords="offset points",
                    fontsize=8.2, fontweight="bold", color="#222", zorder=6)
    if wetlab:
        d["k"] = d["compound"].map(lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower()))
        for _, r in d[d["k"].isin(wetlab)].iterrows():
            oc = wetlab[r["k"]]
            ax.scatter([r["combo_score"]], [r["yj"]], s=230, facecolor="none",
                       edgecolor="#1b7a3d" if oc == "synergy" else "#B01818", linewidth=2.4, zorder=7)

    ax.text(0.015, 0.97, "Path 2 — mechanistic SYNERGY\nDDR / checkpoint inhibition potentiates the damage",
            transform=ax.transAxes, fontsize=10, color="#2E7D6B", va="top", fontweight="bold")
    ax.text(0.985, 0.03, "Path 1 — orthogonality TRAP\ncoverage looks great, but cell-cycle arrest pulls cells\n"
            "out of S-phase and antagonises the S-phase-specific Top1i",
            transform=ax.transAxes, fontsize=10, color="#C0392B", va="bottom", ha="right", fontweight="bold")
    ax.set_xlabel("orthogonality  /  complementary coverage   (Path 1 score)   →", fontsize=11)
    ax.set_ylabel("←  cell-cycle arrest  (antagonism)        mechanism        DDR / checkpoint  (synergy)  →", fontsize=10.5)
    ax.set_title("Two paths to an exatecan dual-payload partner", color=RRB, fontsize=15, fontweight="bold",
                 loc="left", pad=12)
    ax.text(0.5, -0.115, "x = potent on Top1i-resistant lines + orthogonal response (PRISM).  y = target engagement of the "
            "DDR/replication-stress checkpoint (up) vs the cell-cycle-arrest machinery (down).\nKey point: the orthogonality "
            "axis alone favours the arrest agents (mean coverage +0.20) and buries the DDR synergisers (mean -0.20) — you "
            "need the mechanism axis to surface Path 2.", transform=ax.transAxes, ha="center", va="top", fontsize=7.9, color="#666")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=10, label=n)
               for n, c in [("DDR / checkpoint (synergy-plausible)", "#2E7D6B"),
                            ("cell-cycle arrest (antagonism risk)", "#C0392B"), ("other", "#C9CDD3")]]
    if wetlab:
        handles += [Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#1b7a3d",
                           markeredgewidth=2, markersize=12, label="wetlab: synergy"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#B01818",
                           markeredgewidth=2, markersize=12, label="wetlab: antagonism")]
    ax.legend(handles=handles, loc="lower left", fontsize=8.4, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout(rect=[0, 0.05, 1, 1])
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=160)
    plt.close(fig)
    return out_path


def render_moa_orthogonality(matched: pd.DataFrame, *, anchor_dist, out_path, label_top=16) -> Path:
    """matched: compound, target, combo (consensus complementarity score), moa_distance
    (Tahoe distance from Top1i). Shows top partners engage a transcriptional program
    distinct from Top1i (right of the anchor) - non-redundant combinations."""
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)
    d = matched.copy()
    d[["cls", "col"]] = d.apply(lambda r: pd.Series(moa_class(r["target"], r.get("moa", ""))), axis=1)

    fig, ax = plt.subplots(figsize=(11, 7.6))
    ax.axvspan(0.35, anchor_dist + .05, color="#f6eef0", zorder=0)
    ax.axvline(anchor_dist, color=RRB, lw=1.3, ls="--", zorder=1)
    ax.text(anchor_dist - .004, d["combo"].max(), "Top1i\n(topotecan / irinotecan)", color=RRB,
            fontsize=8.5, ha="right", va="top", fontweight="bold")
    ax.scatter(d["moa_distance"], d["combo"], s=70, c=d["col"], alpha=.85,
               edgecolor="#222", linewidth=.5, zorder=3)
    top = d.sort_values("combo", ascending=False).head(label_top)
    for _, r in top.iterrows():
        ax.annotate(r["compound"], (r["moa_distance"], r["combo"]), xytext=(5, 3),
                    textcoords="offset points", fontsize=8.4, fontweight="bold", color="#222", zorder=5)
    ax.set_xlabel("MOA distance from Top1i   (Tahoe-100M signatures)   →  orthogonal", fontsize=11)
    ax.set_ylabel("complementarity  combo score  (PRISM + CTRP)", fontsize=11)
    ax.set_title("MOA orthogonality: complementary partners also engage a different program than Top1i",
                 color=RRB, fontsize=13, fontweight="bold", loc="left", pad=10)
    ax.text(0.995, 0.02, "x = 1 - directional signature concordance with the Top1i anchor (Topotecan + Irinotecan) "
            "in Tahoe-100M.\nEvery top complementary partner is transcriptionally orthogonal to Top1i "
            "(well right of the anchor) = non-redundant.", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.6, color="#666")
    seen = list(dict.fromkeys(top["cls"]))
    cmap = dict([(n, c) for n, _, c in _CLASSES] + [_OTHER])
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=cmap[c], markersize=9, label=c)
                       for c in seen], loc="lower left", fontsize=8.5, frameon=False, title="mechanism",
              title_fontsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=160)
    plt.close(fig)
    return out_path


def render_slfn11_orthogonality(dep: pd.DataFrame, *, anchor_norms, partner_names, out_path) -> Path:
    """Lollipop of SLFN11 dependence: Top1i anchors (SLFN11-dependent) vs the top combo
    partners (~0 = orthogonal to the Top1i resistance axis). dep = slfn11_dependence()."""
    import re
    out_path = Path(out_path); out_path.parent.mkdir(parents=True, exist_ok=True)

    def _n(s):
        return re.sub(r"[^a-z0-9]", "", str(s).lower())
    dep = dep.copy(); dep["k"] = dep["compound"].map(_n)
    anc = dep[dep["k"].isin(set(anchor_norms))].sort_values("slfn11_corr")
    pnorm = [_n(p) for p in partner_names]
    par = dep[dep["k"].isin(set(pnorm))].copy()
    par["ord"] = par["k"].map({k: i for i, k in enumerate(pnorm)})
    par = par.sort_values("ord")

    rows = []
    for _, r in anc.iterrows():
        rows.append((r["compound"], r["slfn11_corr"], RRB, "Top1i anchor", 82, "D"))
    for _, r in par.iterrows():
        _, col = moa_class(r["target"], r.get("moa", ""))
        rows.append((r["compound"], r["slfn11_corr"], col, "partner", 70, "o"))

    fig, ax = plt.subplots(figsize=(10.5, 0.42 * len(rows) + 2.2))
    # context: all compounds' SLFN11 corr as a faint strip
    ax.scatter(dep["slfn11_corr"], np.full(len(dep), -1.1) + np.random.RandomState(0).uniform(-.28, .28, len(dep)),
               s=6, color="#D9DCE1", alpha=.5, zorder=0)
    ax.text(dep["slfn11_corr"].median(), -1.1, "all compounds", fontsize=7.5, color="#9aa0a8", va="center", ha="left")

    ax.axvspan(-0.06, 0.06, color="#EAF2EF", zorder=0)
    ax.axvline(0, color="#bbb", lw=1)
    for i, (name, corr, col, kind, sz, mk) in enumerate(rows):
        ax.plot([0, corr], [i, i], color=col, lw=1.6, alpha=.5, zorder=1)
        ax.scatter([corr], [i], s=sz, c=col, marker=mk, edgecolor="#222", linewidth=.6, zorder=3)
        ax.text(corr + (0.006 if corr >= 0 else -0.006), i, name, va="center",
                ha="left" if corr >= 0 else "right", fontsize=8.6,
                fontweight="bold" if kind == "Top1i anchor" else "normal", color="#222")
    ax.set_yticks([]); ax.set_ylim(-1.9, len(rows) + 1.1)
    ax.set_xlabel("SLFN11 dependence   =   Spearman( potency , SLFN11 expression )", fontsize=11)
    ax.set_title("Resistance orthogonality: Top1i needs SLFN11; good partners do not",
                 color=RRB, fontsize=13.5, fontweight="bold", loc="left", pad=10)
    ax.text(0.0, len(rows) + .5, "orthogonal\n(SLFN11-independent)", fontsize=8.5, color="#2E7D6B", ha="center", va="center")
    ax.text(0.22, len(rows) + .5, "Top1i-dependent\n(fails on SLFN11-low cells) →", fontsize=8.5, color=RRB, ha="center", va="center")
    ax.text(0.995, 0.01, "Top1i drugs sit in the positive tail (need SLFN11 for their damage to kill); a partner near 0 "
            "covers the\nSLFN11-low, Top1i-resistant cells - orthogonal resistance. DepMap SLFN11 (24Q2) × PRISM IC50.",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.6, color="#666")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, format=out_path.suffix.lstrip(".").lower() or "pdf", dpi=160)
    plt.close(fig)
    return out_path
