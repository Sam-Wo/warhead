"""G2 - Does the mechanism survive the delivery bottleneck?

G2a  efflux dependence         - sensitivity vs ABCB1/ABCG2 expression
G2b  proliferation independence - sensitivity vs DepMap doubling time (HCC lever)
G2c  collateral-lethality scan  - CN-loss shifts the dependency distribution

G2b is the quantitative backbone of the MASH-HCC ADC argument: HCC has a lower
proliferative index than CRC, so a payload whose potency is independent of
doubling time is worth far more there than one that only kills fast cyclers.
Expected controls: auristatins / maytansinoids show a strong positive slope,
Top1i intermediate, transcription / translation / spliceosome / degrader agents
approximately flat. We KEEP the flat ones.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sstats

from ..config import load_gates
from ..stats import balance_weights, benjamini_hochberg, weighted_linregress
from .base import GateResult


# ---------------------------------------------------------------------------
# Sensitivity metric
# ---------------------------------------------------------------------------
def sensitivity_from_fits(
    fits: pd.DataFrame,
    *,
    metric: str = "log10_ic50",
    compound_col: str = "compound_id",
    model_col: str = "ModelID",
    use_qc_pass: bool = True,
) -> pd.DataFrame:
    """Build the per-(compound x line) axis G2b regresses on doubling time.

      * ``log10_ic50``     - log10(IC50 [M]); a potency-LOSS axis. Antimitotics
                             rise with doubling time (positive slope), matching
                             the spec's expected sign.
      * ``neg_log10_ic50`` - -log10(IC50); flips the sign (higher = more potent).
      * ``emax``           - completeness of kill (1 - Emax; higher = more kill).

    The gate keys off |slope|, so pass/fail is sign-invariant; the metric only
    sets which way "mitotic-like" points.
    """
    df = fits
    if use_qc_pass and "qc_pass" in df.columns:
        df = df[df["qc_pass"]]
    df = df.copy()
    if metric == "log10_ic50":
        df["sensitivity"] = df["log10_ic50_M"]
    elif metric == "neg_log10_ic50":
        df["sensitivity"] = -df["log10_ic50_M"]
    elif metric == "emax":
        df["sensitivity"] = 1.0 - df["emax"]
    else:
        raise ValueError(f"unknown sensitivity metric: {metric}")
    return df[[compound_col, model_col, "sensitivity"]]


# ---------------------------------------------------------------------------
# G2b - proliferation independence
# ---------------------------------------------------------------------------
def _classify(std_slope: float, q: float, alpha: float) -> str:
    """Label consistent with the gate: a compound is proliferation-independent
    iff its slope is not significant. Significant positives split by effect size
    into mitotic-like (steep, auristatin/maytansinoid regime) vs intermediate
    (Top1i regime)."""
    if not np.isfinite(std_slope) or not np.isfinite(q):
        return "undetermined"
    if q >= alpha:
        return "proliferation-independent"
    if std_slope >= 0.55:
        return "mitotic-like"
    if std_slope > 0:
        return "intermediate"
    return "significant-negative"


def proliferation_stats(
    sensitivity: pd.DataFrame,
    model_meta: pd.DataFrame,
    *,
    compound_col: str = "compound_id",
    model_col: str = "ModelID",
    config: dict | None = None,
) -> pd.DataFrame:
    """Per-compound regression of sensitivity on doubling time.

    Weights correct pooled-PRISM under-representation of slow-growing lines by
    balancing doubling-time quantile bins before the fit.
    """
    cfg = (config or load_gates())["g2"]["proliferation"]
    dt_col = cfg["doubling_time_col"]
    meta = model_meta[[model_col, dt_col]].dropna(subset=[dt_col])
    merged = sensitivity.merge(meta, on=model_col, how="inner")

    weighting = cfg.get("weighting", "quantile_balance")
    rows = []
    for comp, g in merged.groupby(compound_col, sort=False):
        x = g[dt_col].to_numpy(float)
        y = g["sensitivity"].to_numpy(float)
        if weighting == "balance":
            w = balance_weights(x, n_bins=cfg.get("weight_bins", 5))
        else:
            w = np.ones_like(x)
        fit = weighted_linregress(x, y, w)
        rows.append(
            {
                compound_col: comp,
                "slope": fit.slope,
                "std_slope": fit.std_slope,
                "slope_se": fit.slope_se,
                "p": fit.p,
                "n_lines": fit.n,
                "r": fit.r,
            }
        )
    stats = pd.DataFrame(rows)
    if len(stats):
        stats["q"] = benjamini_hochberg(stats["p"].to_numpy())
        alpha = cfg["fdr_alpha"]
        stats["prolif_class"] = [
            _classify(s, q, alpha) for s, q in zip(stats["std_slope"], stats["q"])
        ]
    return stats


def gate_g2b(
    sensitivity: pd.DataFrame,
    model_meta: pd.DataFrame,
    *,
    compound_col: str = "compound_id",
    model_col: str = "ModelID",
    config: dict | None = None,
) -> GateResult:
    cfg = (config or load_gates())["g2"]["proliferation"]
    stats = proliferation_stats(
        sensitivity, model_meta, compound_col=compound_col, model_col=model_col, config=config
    )
    if not len(stats):
        empty = stats.assign(g2b_reason=pd.Series(dtype=str))
        return GateResult("G2b", empty, empty, "g2b_reason", cfg, {"n_compounds": 0})

    # Primary criterion (WARHEAD.md G2b): keep compounds whose slope is NOT
    # significantly different from zero. std_slope (the standardised effect size,
    # == the weighted correlation for a simple regression) is reported and used
    # in classification, but the pass/fail is the significance test so we do not
    # silently drop a genuinely flat-but-slightly-noisy compound.
    alpha = cfg["fdr_alpha"]
    min_lines = cfg["min_lines"]

    reasons, passed_mask = [], []
    for _, row in stats.iterrows():
        if row["n_lines"] < min_lines:
            passed_mask.append(False)
            reasons.append(f"only {int(row['n_lines'])} lines (need {min_lines})")
        elif not (row["q"] > alpha):
            passed_mask.append(False)
            reasons.append(
                f"proliferation-dependent slope (q={row['q']:.3g} <= {alpha}, "
                f"std slope {row['std_slope']:+.2f})"
            )
        else:
            passed_mask.append(True)
            reasons.append("")
    stats = stats.assign(passed=passed_mask, g2b_reason=reasons)

    passed = stats[stats["passed"]].drop(columns="passed").reset_index(drop=True)
    failed = stats[~stats["passed"]].drop(columns="passed").reset_index(drop=True)
    return GateResult(
        gate="G2b",
        passed=passed,
        failed=failed,
        reason_col="g2b_reason",
        config=cfg,
        summary={
            "n_compounds": len(stats),
            "n_pass": len(passed),
            "class_counts": stats["prolif_class"].value_counts().to_dict(),
        },
    )


# ---------------------------------------------------------------------------
# G2a - efflux dependence
# ---------------------------------------------------------------------------
def gate_g2a(
    sensitivity: pd.DataFrame,
    expression: pd.DataFrame,
    *,
    compound_col: str = "compound_id",
    model_col: str = "ModelID",
    gene_col: str = "gene",
    expr_col: str = "expression",
    config: dict | None = None,
) -> GateResult:
    """Regress the resistance axis (log10_ic50) on ABCB1/ABCG2 expression. A
    strong positive slope on either transporter means IC50 rises with transporter
    level = efflux substrate = fails the way MMAE / DM1 fail. Pass the same
    ``sensitivity_from_fits(metric='log10_ic50')`` frame used elsewhere."""
    cfg = (config or load_gates())["g2"]["efflux"]
    genes = cfg["genes"]
    expr = expression[expression[gene_col].isin(genes)]

    rows = []
    for comp, g in sensitivity.groupby(compound_col, sort=False):
        rec = {compound_col: comp}
        best_slope, best_p = -np.inf, np.nan
        for gene in genes:
            e = expr[expr[gene_col] == gene][[model_col, expr_col]]
            m = g.merge(e, on=model_col, how="inner")
            fit = weighted_linregress(m[expr_col].to_numpy(float), m["sensitivity"].to_numpy(float))
            rec[f"std_slope_{gene}"] = fit.std_slope
            rec[f"p_{gene}"] = fit.p
            # Track the transporter with the strongest POSITIVE dependence.
            if np.isfinite(fit.std_slope) and fit.std_slope > best_slope:
                best_slope, best_p = fit.std_slope, fit.p
        rec["max_efflux_std_slope"] = best_slope if np.isfinite(best_slope) else np.nan
        rec["max_efflux_p"] = best_p
        rows.append(rec)
    stats = pd.DataFrame(rows)

    smax = cfg["std_slope_max"]
    alpha = cfg.get("fdr_alpha", 0.05)
    # Significance guard: a large slope must also be significant, or a null
    # compound trips the gate ~5% of the time by chance.
    stats["max_efflux_q"] = benjamini_hochberg(stats["max_efflux_p"].to_numpy())
    is_substrate = (stats["max_efflux_std_slope"] >= smax) & (stats["max_efflux_q"] < alpha)
    stats["passed"] = ~is_substrate
    stats["g2a_reason"] = np.where(
        stats["passed"],
        "",
        "efflux substrate: std slope "
        + stats["max_efflux_std_slope"].round(2).astype(str)
        + " (q=" + stats["max_efflux_q"].map(lambda v: f"{v:.2g}") + ")",
    )
    passed = stats[stats["passed"]].drop(columns="passed").reset_index(drop=True)
    failed = stats[~stats["passed"]].drop(columns="passed").reset_index(drop=True)
    return GateResult("G2a", passed, failed, "g2a_reason", cfg,
                      {"n_compounds": len(stats), "n_pass": len(passed)})


# ---------------------------------------------------------------------------
# G2c - collateral-lethality scan (generalises the POLR2A logic)
# ---------------------------------------------------------------------------
def collateral_lethality_scan(
    chronos: pd.DataFrame,
    copy_number: pd.DataFrame,
    *,
    gene_col: str = "gene",
    model_col: str = "ModelID",
    chronos_col: str = "chronos",
    cn_col: str = "cn_log2",
    indication_models: list[str] | None = None,
    config: dict | None = None,
) -> pd.DataFrame:
    """For each gene, test whether hemizygous CN loss shifts the Chronos
    dependency distribution leftward (more essential) using a one-sided
    Mann-Whitney U. BH-controlled across genes.

    POLR2A on 17p (co-deleted with TP53) is the positive control and must recover
    before any novel hit is trusted.
    """
    cfg = (config or load_gates())["g2"]["collateral"]
    thr = cfg["cn_loss_threshold"]
    cn = copy_number
    ch = chronos
    if indication_models is not None:
        cn = cn[cn[model_col].isin(indication_models)]
        ch = ch[ch[model_col].isin(indication_models)]

    merged = ch.merge(cn, on=[gene_col, model_col], how="inner")
    rows = []
    for gene, g in merged.groupby(gene_col, sort=False):
        loss = g[g[cn_col] < thr][chronos_col].to_numpy(float)
        neutral = g[g[cn_col] >= thr][chronos_col].to_numpy(float)
        if loss.size < 3 or neutral.size < 3:
            continue
        # One-sided: loss lines MORE dependent (lower Chronos) than neutral.
        u, p = sstats.mannwhitneyu(loss, neutral, alternative="less")
        rows.append(
            {
                gene_col: gene,
                "n_loss": int(loss.size),
                "n_neutral": int(neutral.size),
                "median_chronos_loss": float(np.median(loss)),
                "median_chronos_neutral": float(np.median(neutral)),
                "delta": float(np.median(loss) - np.median(neutral)),
                "p": float(p),
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out["q"] = benjamini_hochberg(out["p"].to_numpy())
        out["significant"] = out["q"] < cfg["fdr_alpha"]
        out = out.sort_values("q").reset_index(drop=True)
    return out


def rank_collateral_targets(
    scan: pd.DataFrame,
    tcga_recurrence: pd.DataFrame,
    *,
    indication: str,
    gene_col: str = "gene",
    config: dict | None = None,
) -> pd.DataFrame:
    """Combine the DepMap dependency-shift scan with TCGA loss recurrence to
    produce the payload-TARGET list.

    A collateral-lethality target must (WARHEAD.md G2c):
      * show a significant leftward DepMap dependency shift on CN loss,
      * be recurrently hemizygously lost in the indication (TCGA frequency), and
      * not be a common essential (pan-dependent regardless of CN).

    Ranked by ``collateral_score = -delta * loss_frequency`` (bigger differential
    dependency in a more frequently deleted gene). Steps 2-3 of the spec (ChEMBL
    chemical matter + a substitutable linker position) are the chemistry follow-up
    and are left to the G4/G5 gates.
    """
    cfg = (config or load_gates())["g2"]["collateral"]
    rec = tcga_recurrence[tcga_recurrence["indication"] == indication][
        [gene_col, "loss_frequency", "co_deleted"]
    ]
    out = scan.merge(rec, on=gene_col, how="left")
    out["loss_frequency"] = out["loss_frequency"].fillna(0.0)

    out["common_essential"] = out["median_chronos_neutral"] < cfg["common_essential_chronos"]
    out["recurrent"] = out["loss_frequency"] >= cfg["recurrence_min_frequency"]
    out["candidate"] = (
        out["significant"] & out["recurrent"] & ~out["common_essential"] & (out["delta"] < 0)
    )
    out["collateral_score"] = (-out["delta"]).clip(lower=0) * out["loss_frequency"]
    return out.sort_values(["candidate", "collateral_score"], ascending=[False, False]).reset_index(drop=True)
