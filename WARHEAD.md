# WARHEAD

**Public-data cascade for prioritising novel ADC payload classes.**
Indications in scope: colorectal carcinoma (CRC) and hepatocellular carcinoma (HCC).
Owner: PFR / RedRidge Bio.

---

## 0. Framing — read this before writing any code

There are roughly 15–25 clinically validated ADC payloads. That is the entire positive
set. **This is not a supervised learning problem.** Any model trained on that n will
memorise the auristatin/maytansinoid/camptothecin scaffolds and tell you nothing.

WARHEAD is therefore an **explicit filter cascade** with stated, arguable thresholds at
each gate. Every gate must be independently defensible to a reviewer. A candidate that
fails G2 is dead regardless of how it scores at G5 — no weighted composite score that
lets a strong dimension mask a fatal one.

Output is a ranked shortlist of **payload targets and chemical classes**, with a
provenance trail per candidate, not a single scalar.

Two program-specific questions the cascade must answer directly:

1. **Partner warhead for exatecan.** What payload class has a resistance profile
   orthogonal to Top1i, given SLFN11 as the sensitising axis? (Current lead partner is
   ATRi — the cascade should either corroborate that or produce alternatives.)
2. **Non-mitotic payload for HCC.** HCC has a lower proliferative index than CRC.
   Which payload classes retain potency independent of doubling time? This is the
   quantitative backbone of the MASH-HCC ADC argument.

---

## 1. The constraint vector

A payload must clear all six gates simultaneously:

| Gate | Constraint | Primary evidence |
|---|---|---|
| G1 | Sub-nM potency **with complete kill** (Emax → 0) | NCI-60, CTRP, GDSC, PRISM |
| G2 | Mechanism survives the delivery bottleneck (no efflux dependence, no mitosis requirement) | DepMap, Tahoe-100M, CCLE |
| G3 | Novel or under-exploited MOA space | LINCS, JUMP-CP, Tahoe-100M |
| G4 | Bystander-competent physicochemistry | RDKit over COCONUT/NPAtlas/ChEMBL |
| G5 | Conjugatable handle at an SAR-tolerant position | SMARTS + ChEMBL SAR |
| G6 | Payload target not enriched in recurring DLT organs | GTEx, HPA, FAERS |

---

## 2. Data sources

Priority: **P0** = build the first pass on these. **P1** = second wave. **P2** = nice to have.

### Compound sensitivity

| Source | Content | Access | Pri | Notes |
|---|---|---|---|---|
| **NCI-60 / DTP** | GI50, TGI, LC50 for >50,000 compounds, 60 lines | `dtp.cancer.gov`, CellMiner at `discover.nci.nih.gov/cellminer` | **P0** | The only public source that routinely reaches the sub-nM regime at scale, and heavily enriched for natural products. Every founding payload class came out of this screening lineage. Vintage QC is heterogeneous — hypothesis generation only. |
| **CTRP v2** | 481 compounds × ~860 lines, dense dose-response, MOA-annotated | CTD² portal; mirrored on DepMap | **P0** | Best MOA annotation density. Raw curves available — use them. |
| **GDSC1 / GDSC2** | ~500 compounds × ~1000 lines | `cancerrxgene.org` | **P0** | Raw dose-response point data downloadable. GDSC2 supersedes GDSC1 for overlapping compounds; do not pool naively. |
| **PRISM Repurposing** | Barcoded pooled viability, thousands of compounds × ~900 lines | DepMap portal (`depmap.org/portal/download`) | **P0** | Breadth over depth. Primary screen is single-dose; the secondary screen has 8-pt dose-response — **only the secondary is usable for G1**. Pooled format means slow-growing lines are underrepresented; correct for this before G2. |
| **CellMinerCDB** | Integration layer across NCI-60 / GDSC / CTRP / CCLE | `discover.nci.nih.gov` | P1 | Use for cross-screen reconciliation, not as a primary source. |

### Perturbation response / MOA

| Source | Content | Access | Pri | Notes |
|---|---|---|---|---|
| **Tahoe-100M** | >100M single-cell transcriptomes, 50 cancer lines × ~1,100 small molecules, ~60k conditions | HuggingFace `tahoebio/Tahoe-100M`; Arc Virtual Cell Atlas | **P0** | The single most valuable new addition. Single-cell resolution means you can detect **response heterogeneity within a line** — which is the persister-fraction signal that bulk Emax only approximates. Caveat: 24 h timepoint, so it reads MOA and early stress response, **not** killing kinetics. Use pseudobulk for the first pass; go single-cell only for G1b and G2b. Large — plan storage before pulling. |
| **LINCS L1000** | ~1.3M level-5 signatures, ~30k perturbagens | `clue.io` | **P0** | Broad compound coverage, coarse readout (978 landmark genes). Complements Tahoe's depth-over-breadth. |
| **JUMP Cell Painting** | Morphological profiles, ~116k compounds | `cellpainting-gallery` on AWS Open Data | P1 | Orthogonal modality to transcriptomics. **Agreement between JUMP and LINCS/Tahoe is a much stronger novelty signal than either alone** — that concordance filter is the point of including it. Confounded by dose; normalise. |
| **scBaseCount** | ~500M cells, 75 tissues, AI-curated public scRNA-seq | Arc Institute | P2 | Already in use on the target-discovery side. Here it is only useful for G6 normal-tissue expression of payload targets. |

### Genetics and dependency

| Source | Content | Access | Pri | Notes |
|---|---|---|---|---|
| **DepMap 26Q1** | CRISPR (Chronos), expression, CN, doubling time, common essentials | `depmap.org/portal/download` | **P0** | Already local — reuse the loader from `depmap_adc_payload_scoring.py`. Doubling-time metadata is the input to G2b. |
| **TCGA / GDC (COAD, READ, LIHC)** | GISTIC2 copy number, expression | `portal.gdc.cancer.gov`, cBioPortal | **P0** | Feeds the collateral-lethality scan (G2c). |
| **Sanger Project Score** | Independent CRISPR dependency | `score.depmap.sanger.ac.uk` | P1 | Orthogonal replication of dependency calls. |

### Chemistry

| Source | Content | Access | Pri | Notes |
|---|---|---|---|---|
| **ADCdb** | 6,572 ADCs; 359 approved or clinical, 501 preclinical; payload structures + drug-like properties | `adcdb.idrblab.net` | **P0** | The labelled positive set. Everything the cascade produces gets benchmarked against it. |
| **ChEMBL** | Bioactivity, SAR, structures | `ebi.ac.uk/chembl` | **P0** | Also the source for "does chemical matter exist for this target at sub-nM cellular potency". |
| **COCONUT** | ~400k natural products | `coconut.naturalproducts.net` | **P0** | Payload chemical space is overwhelmingly natural-product derived. |
| **NPAtlas** | Microbial natural products | `npatlas.org` | P1 | Higher curation quality than COCONUT, smaller. |

### Tolerability

| Source | Content | Access | Pri | Notes |
|---|---|---|---|---|
| **GTEx** | Bulk normal-tissue expression | `gtexportal.org` | **P0** | G6 input. |
| **Human Protein Atlas** | Protein-level tissue and single-cell expression | `proteinatlas.org` | **P0** | Better than GTEx for the specific DLT compartments (cornea, GI crypt, alveolar type II). |
| **FAERS** | Adverse event reports | FDA quarterly data extract files | P1 | Disproportionality analysis restricted to ADC regimens, stratified by payload class. |
| **ClinicalTrials.gov / Drugs@FDA** | Trial design, labels, DLTs | public APIs | P1 | Ground truth for the payload-class → DLT-organ map. |

---

## 3. The cascade

### G1 — Potency, done properly

**G1a. Refit every curve. Abandon AUC.**

AUC and IC50 conflate potency with completeness of kill, and for payloads those are not
interchangeable. A systemic small molecule with an Emax plateau at 0.35 can be dosed
harder; an ADC cannot — delivery caps intracellular concentration, so there is no
headroom to out-dose a persister fraction.

- Pull raw dose-response point data from CTRP v2, GDSC1/2, PRISM secondary.
- Refit per compound × line with a 4-parameter logistic. Extract **IC50, Emax, Hill
  slope** as separate features. Do not collapse them.
- **Handle left-censoring properly.** Most screens bottom out at 1–10 nM, which is
  exactly where payload-relevant activity begins. Treat "no response at lowest dose"
  and "full response at lowest dose" as interval-censored observations and fit with a
  Tobit-style likelihood. Imputing the lowest tested dose is the single most common
  error here and it systematically flattens the ranking.
- NCI-60 GI50 already extends lower — bring it in as an independent potency axis,
  cross-check against CTRP/GDSC on overlapping compounds, and report the disagreement
  rather than averaging it away.

**Gate:** IC50 < 1 nM in ≥ 20% of lines **and** median Emax < 0.15.

**G1b. Persister check (Tahoe).**
For compounds present in Tahoe-100M, test whether the single-cell response distribution
at 24 h is unimodal or shows a non-responding subpopulation. A bimodal response in a
sensitive line is a red flag that bulk Emax missed.

### G2 — Does the mechanism survive delivery?

**G2a. Efflux dependence.**
Regress per-compound sensitivity on `ABCB1` and `ABCG2` expression across CCLE lines.
Strong positive dependence = the compound is an efflux substrate = it will fail the way
MMAE and DM1 fail. This is the shared failure mode across all three incumbent classes,
so escaping it is a differentiator in itself.

**G2b. Proliferation independence — the HCC lever.**
Regress compound sensitivity against DepMap doubling time across lines. Expected
behaviour: auristatins and maytansinoids show a strong positive slope, Top1 inhibitors
intermediate, transcription / translation / spliceosome / protein-degradation agents
approximately flat.

Keep compounds where the slope is not significantly different from zero. This produces a
quantitative argument for a non-mitotic payload class in a low-proliferative-index
indication, which is a much stronger claim than the usual mechanistic hand-waving.

Note: pooled PRISM underrepresents slow-growing lines by construction. Weight or
restrict accordingly before running this regression, or the effect will be attenuated.

**G2c. Collateral-lethality scan (generalises the POLR2A logic).**
For every gene:
1. Does hemizygous loss in TCGA-COAD/READ or TCGA-LIHC shift the DepMap Chronos
   dependency distribution leftward? (Mann-Whitney on CN-loss vs CN-neutral lines,
   FDR-controlled.)
2. Does chemical matter exist in ChEMBL with sub-nM cellular activity against the
   product?
3. Is there a substitutable position for linker attachment?

Known anchor: POLR2A on 17p, co-deleted with TP53. Candidates worth explicit
interrogation: **ME2 on 18q** (co-deleted with SMAD4, high frequency in CRC), 1p36 and
8p losses in HCC. Output is a payload *target* list — the right granularity to hand to
chemistry.

### G3 — MOA novelty and orthogonality

**G3a. Embed and locate the holes.**
Build a joint embedding of LINCS L1000 signatures, Tahoe pseudobulk profiles, and JUMP
morphological profiles. Annotate the regions occupied by known payload classes using
ADCdb. Score candidates on distance to the nearest known-payload centroid **conditional
on having passed G1** — novelty without potency is worthless.

**G3b. Orthogonal-resistance search — the exatecan partner question.**
Take the PRISM/CTRP sensitivity matrix. Regress out the Top1i component (or use SLFN11
expression as the covariate) and rank compounds by potency in the **residual** space,
i.e. on lines that exatecan does not handle. Stratify by ABCB1 status so you are not
just rediscovering efflux.

This generates partner-warhead candidates empirically rather than by mechanistic
argument, and gives you a way to check whether ATRi is genuinely the best occupant of
that slot or merely the most obvious one.

### G4 — Bystander competence

Do **not** rebuild the bystander predictor from scratch. Prior art exists: Guo et al.,
*Adv. Sci.* 2024 (10.1002/advs.202306309) trained a graph attention model on membrane
permeability, validated it against 80+ clinical and development payloads from ADCdb, and
set a B-score threshold of 1.5 for bystander vs non-bystander. They also report a
correlation between IC50 and calculated cLogD across the ADCdb payload set.

Layer on top of that:
- **Charge state at lysosomal pH ~4.8 vs cytosol 7.2.** The MMAE/MMAF split is
  essentially this — MMAF's charged C-terminal phenylalanine is why it does not travel.
- cLogD(7.4), TPSA, MW.

Run RDKit over the G1-passing set intersected with COCONUT / NPAtlas / ChEMBL.

Note that bystander is a *design choice*, not a universal good — for a heterogeneous CRC
tumour it is essential; for a narrow-window target it is a liability. Tag rather than
filter, and let the antigen program decide.

### G5 — Conjugatability

SMARTS matching for a primary/secondary amine, hydroxyl, thiol, or carboxylic acid, at a
position that published SAR says tolerates substitution. Cross-reference ChEMBL SAR
series for the scaffold — a handle in the pharmacophore is not a handle.

### G6 — Therapeutic window

Payload class determines DLT organ far more than antigen does: ocular for MMAF, ILD for
DXd, thrombocytopenia for maytansinoids, peripheral neuropathy for MMAE.

1. Run disproportionality analysis (ROR / PRR, with the usual shrinkage) on FAERS
   restricted to ADC-containing regimens, stratified by payload class. This builds an
   empirical class → toxicity map rather than relying on the review-article version.
2. For each candidate payload target, score GTEx/HPA expression across the five
   recurring DLT compartments: **HSC/bone marrow, GI crypt, cornea, alveolar type II,
   peripheral nerve.** Low expression across all five is the window signal.

---

## 4. Repo layout

```
warhead/
  CLAUDE.md                    # this file, trimmed to working context
  data/
    raw/                       # never edited, checksummed on download
    interim/
    processed/
  src/warhead/
    io/                        # one loader per source, all returning tidy frames
      depmap.py                # reuse from depmap_adc_payload_scoring.py
      nci60.py
      ctrp.py  gdsc.py  prism.py
      tahoe.py                 # pseudobulk first; streaming for single-cell
      lincs.py  jump.py
      adcdb.py  chembl.py  coconut.py
      gtex.py  hpa.py  faers.py
      tcga.py
    curves/
      refit.py                 # 4PL + interval-censored likelihood
      qc.py
    gates/
      g1_potency.py
      g2_delivery.py           # efflux, proliferation, collateral lethality
      g3_moa.py                # embedding + orthogonal-resistance search
      g4_bystander.py
      g5_conjugation.py
      g6_window.py
    cascade.py                 # orchestration, provenance, per-gate audit trail
  notebooks/
  reports/                     # RRB maroon template #6E1426
  tests/
```

Conventions:
- Every gate emits `(passed_df, failed_df, reason_column)`. Nothing is silently dropped.
- Compound identity resolved to InChIKey at ingest. Cross-screen name matching will
  otherwise quietly lose 20–30% of overlaps.
- Cell line identity resolved to DepMap `ModelID` (ACH-######) at ingest.
- All thresholds live in one `config/gates.yaml`. Sensitivity analysis over thresholds
  is a required deliverable, not an optional one.

---

## 5. Build order

1. **Ingest + identity resolution** (InChIKey, ModelID). Nothing works until this does.
2. **G1 curve refitting** across CTRP + GDSC2 + PRISM secondary. Validate the refit by
   confirming that known payloads in ADCdb land where they should.
3. **G2b proliferation regression.** Fastest route to a defensible, program-relevant
   result. This is the one that feeds the HCC deck.
4. **G2c collateral-lethality scan.** Confirm POLR2A recovers as a positive control
   before trusting any novel hit.
5. **G3b orthogonal-resistance search.** Answers the exatecan-partner question.
6. **G3a embedding**, once Tahoe and JUMP are local.
7. **G4/G5 chemical gates.**
8. **G6 window scoring.**

Steps 3 and 5 are independently publishable/presentable before the full cascade exists.
Build toward them first.

---

## 6. Non-negotiables and known failure modes

- **Free-drug IC50 does not predict ADC IC50.** Conjugation site, DAR, linker chemistry
  and internalisation rate move it by one to two logs, and no public dataset carries
  linker context. Every output of this cascade is a hypothesis for conjugation, not a
  ranked list of ADCs.
- **2D monolayer screening systematically overrates antimitotics.** G2b exists to
  correct for exactly this bias — do not undo it by reintroducing raw potency rank at
  the end.
- **NCI-60 spans decades of variable QC.** Never let an NCI-60-only hit through to a
  deck without independent confirmation.
- **Cell Painting MOA clusters confound with potency.** Normalise by dose before
  computing distances, or "novel MOA" will just mean "tested at a different
  concentration".
- **Tahoe is 24 h.** It reads mechanism and early stress response. It does not read
  cell killing. Do not use it as a viability proxy.
- **No composite score.** Gates are conjunctive. If a reviewer can only argue with a
  weighted sum, the analysis is not defensible.
- With ~20 positives, nested cross-validation will not rescue a supervised model. If a
  model appears anywhere in this repo, it is a filter with a stated prior, and its
  threshold is in `gates.yaml` where someone can argue with it.

---

## 7. Deliverables

- `reports/warhead_shortlist.xlsx` — ranked candidates with per-gate pass/fail and
  provenance.
- `reports/proliferation_independence.pdf` — G2b, standalone, for the HCC/MASH argument.
- `reports/exatecan_partner.pdf` — G3b, standalone, for the dual-payload program.
- `reports/collateral_lethality_crc_hcc.xlsx` — G2c scan output with POLR2A as positive
  control.
- Threshold sensitivity analysis across all gates.
