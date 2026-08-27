# WARHEAD

**Public-data cascade for prioritising novel ADC payload classes** (colorectal &
hepatocellular carcinoma). Owner: PFR / RedRidge Bio.

This is an **explicit conjunctive filter cascade**, not a supervised model. There
are ~15–25 clinically validated ADC payloads — too few to learn from — so every
gate is a stated, arguable threshold that must stand on its own to a reviewer. See
[`WARHEAD.md`](WARHEAD.md) for the full framing, data sources, and gate design,
and [`CLAUDE.md`](CLAUDE.md) for the trimmed working context.

## Status

Build order per `WARHEAD.md` §5. What is wired end-to-end today:

| Piece | State |
|---|---|
| Identity resolution (InChIKey, DepMap ModelID) | ✅ implemented |
| **G1** — 4PL refit with interval-censored (Tobit) likelihood + QC | ✅ implemented, tested |
| **G1** — potency gate (sub-nM in ≥20% lines, median Emax < 0.15) | ✅ implemented, tested |
| **G2a** — efflux dependence (ABCB1/ABCG2), FDR-guarded | ✅ implemented, tested |
| **G2b** — proliferation independence (the HCC lever) | ✅ implemented, tested, **report** |
| **G2c** — collateral-lethality scan (POLR2A positive control) | ✅ implemented, tested |
| G3–G6 | 🚧 interfaces + contracts only (deferred; see `src/warhead/gates/`) |
| Real loaders (DepMap, PRISM, CTRP, GDSC, …) | 🟡 DepMap/PRISM wired to documented schemas; rest stubbed |

Nothing in the cascade requires the real datasets to *run*: the synthetic fixture
(`warhead.fixtures`) emits frames in the **same tidy schema** the real loaders
produce, so the identical downstream code runs on synthetic or real data.

## Quickstart

```bash
# Python 3.10+; scientific stack (numpy/scipy/pandas/matplotlib/pyyaml/openpyxl).
pip install -e .            # or: pip install -e '.[chem,dev]'  for RDKit + pytest

# Run the G1 -> G2b slice on the synthetic fixture and write the deliverables:
warhead demo               # -> reports/proliferation_independence.{pdf,xlsx}

# Or without installing:
PYTHONPATH=src python -m warhead demo

# Inspect the resolved gate thresholds:
warhead info

# Once real data is in data/raw/{depmap,prism}/ :
warhead g2b --raw-dir data/raw
```

### What `warhead demo` demonstrates

The fixture plants a known per-compound proliferation dependence in each cell
line's IC50. The pipeline refits every curve, then G2b regresses recovered
`log10(IC50)` on DepMap doubling time and re-detects the planted slope:

- recovered vs planted slope correlation ≈ **0.99**
- antimitotic controls (auristatin/maytansinoid) → **steep positive** slope → fail G2b
- Top1i controls (exatecan/SN38) → **intermediate** → fail G2b
- non-mitotic classes (RNAPII/translation/spliceosome/degrader) → **flat** → pass G2b

exactly the ordering `WARHEAD.md` §G2b predicts. `reports/proliferation_independence.pdf`
is the standalone figure for the MASH-HCC ADC argument.

## Layout

```
config/gates.yaml     every arguable threshold (single source of truth)
src/warhead/
  identity.py         InChIKey + ModelID resolution
  curves/refit.py     4PL + interval-censored (Tobit) refit   <- G1a
  stats.py            weighted OLS, BH-FDR, slow-line balance weights
  gates/g1_potency.py G1
  gates/g2_delivery.py G2a efflux / G2b proliferation / G2c collateral
  gates/g3..g6        contracts for the deferred gates
  io/                 one loader per source (tidy frames)
  fixtures/synth.py   schema-faithful synthetic data with planted truth
  reporting/          deliverable figures
  cascade.py          orchestration + per-gate provenance
tests/                pytest (refit recovery, censoring, G2b/G2a/G2c behaviour)
```

## Design rules (non-negotiable)

- **No composite score.** Gates are conjunctive; a fatal gate is fatal.
- **Every gate emits `(passed, failed, reason)`** — nothing is silently dropped.
- **Free-drug IC50 ≠ ADC IC50.** Every output is a hypothesis *for conjugation*.
- All thresholds live in `config/gates.yaml`; threshold sensitivity analysis is a
  required deliverable, not optional.

## License

Proprietary — RedRidge Bio. All rights reserved. Not for distribution.
