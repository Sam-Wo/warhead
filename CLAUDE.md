# CLAUDE.md — WARHEAD working context

Trimmed operating context for work in this repo. The full spec is `WARHEAD.md`
(read it before touching a gate); this file is the day-to-day rules.

## What this is

An **explicit conjunctive filter cascade** for prioritising novel ADC payload
classes (CRC/HCC). NOT a supervised model — ~20 positives, so any classifier just
memorises auristatin/maytansinoid/camptothecin. Every gate is a stated threshold
in `config/gates.yaml` that must be independently defensible to a reviewer.

Two program questions the cascade must answer directly:
1. Partner warhead for exatecan (resistance orthogonal to Top1i; SLFN11 axis) — G3b.
2. Non-mitotic payload for HCC (potency independent of doubling time) — **G2b**.

## Hard rules (do not violate)

- **No composite score.** Gates are conjunctive; a candidate failing any gate is
  dead regardless of other dimensions.
- Every gate returns a `GateResult(passed, failed, reason_col, ...)` — nothing is
  silently dropped (`gates/base.py`).
- Compound identity → **InChIKey** at ingest; cell line → DepMap **ModelID**
  (ACH-######). Cross-screen name matching loses 20–30% otherwise.
- All thresholds live in `config/gates.yaml`. Never hard-code one in a module.
  Threshold sensitivity analysis is a required deliverable.
- Free-drug IC50 ≠ ADC IC50 (linker/DAR/internalisation move it 1–2 logs). Every
  output is a hypothesis **for conjugation**, not a ranked list of ADCs.
- G1: refit curves; keep IC50, Emax, Hill **separate**; never use AUC. Handle the
  assay floor/ceiling as **interval-censored** (Tobit) — do not impute the lowest
  tested dose (it flattens the ranking).
- G2b: 2D monolayer screening over-rates antimitotics; G2b exists to correct that.
  Do not reintroduce raw potency rank at the end.
- Tahoe is 24 h → reads MOA/early stress, **not** cell killing. Not a viability proxy.

## Current state (build order §5)

Wired + tested: identity, G1 refit (+censoring) & potency gate, G2a efflux,
**G2b proliferation independence (+ PDF report)**, G2c collateral lethality.
Deferred (contracts only in `gates/g3..g6.py`): G3 MOA/orthogonal-resistance,
G4 bystander, G5 conjugation, G6 window.

Loaders: DepMap + PRISM wired to documented schemas; others stubbed in `io/`.
The synthetic fixture (`fixtures/synth.py`) matches the real loader schemas, so
downstream code is source-agnostic.

## Conventions

- Sign: G2b/G2a regress the **`log10_ic50`** (resistance) axis, so antimitotic /
  efflux dependence reads as a **positive** slope (matches the spec's wording).
  `std_slope` == the weighted correlation; the gate keys off **significance**
  (BH-adjusted `q > alpha`), with `std_slope_max` advisory (report band + class).
- Pooled PRISM under-represents slow lines → `stats.balance_weights` up-weights
  sparse doubling-time (slow) bins before the G2b regression.

## Commands

```
PYTHONPATH=src python -m warhead demo     # G1->G2b on synthetic data -> reports/
PYTHONPATH=src python -m warhead info     # print resolved thresholds
PYTHONPATH=src python -m pytest -q        # tests (refit ~slow; ~45s total)
```

Python is available here as the `py` launcher (3.13). Use `py -m pytest`, etc.

## Deliverables (§7)

`reports/warhead_shortlist.xlsx`, `reports/proliferation_independence.pdf` (G2b —
done), `reports/exatecan_partner.pdf` (G3b), `reports/collateral_lethality_crc_hcc.xlsx`
(G2c), + threshold sensitivity across all gates.
