# data/raw — never edited, checksummed on download

Raw source files are **not** committed (see `.gitignore`). Drop each source into
its own subdirectory; loaders in `src/warhead/io/` read from here and raise an
actionable error naming the expected file if it is absent.

## Expected files (G1 -> G2b slice)

### `depmap/`  (depmap.org/portal/download, release 26Q1)
- `Model.csv` — model metadata; must carry a doubling-time column
  (`doubling_time_hours` / `DoublingTime` / …). If doubling time ships separately,
  add `DoublingTime.csv` with `ModelID` + a doubling-time column.
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` — expression (G2a: ABCB1/ABCG2).
- `CRISPRGeneEffect.csv` — Chronos dependency (G2c).

### `prism/`  (DepMap portal, PRISM Repurposing **secondary** screen)
- `prism_secondary_long.parquet` (or `.csv`) — tidy long dose-response with
  columns `[compound_id, ModelID, dose_M, viability]`. Only the secondary screen
  (8-pt dose-response) is usable for G1; the primary screen is single-dose.

### `gdsc/`  (cancerrxgene.org bulk download, release 8.5)
- `GDSC2_fitted_dose_response.xlsx` — fitted LN_IC50 + AUC per drug x line
  (from `cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/`). Drives `warhead gdsc`
  (EC90 + selectivity) and `warhead g2b-real`.
- `Cell_Lines_Details.xlsx`, `screened_compounds.csv` — annotation (optional).

### `depmap/`  (figshare mirror — the DepMap PORTAL is behind Cloudflare)
The portal download API is gated by a Cloudflare Turnstile challenge, so pull the
release from its **figshare** article instead (public API, no gate). E.g. DepMap
24Q2 Public = figshare article `25880521`; list files via
`https://api.figshare.com/v2/articles/25880521/files` and download by
`https://ndownloader.figshare.com/files/<id>`. For `warhead g2b-real` you need
only `Model.csv` (0.6 MB) and `CRISPRInferredModelGrowthRate.csv` (0.04 MB); for a
real G2c you also need `CRISPRGeneEffect.csv` (~419 MB) and `OmicsCNGene.csv`
(~817 MB), plus TCGA GISTIC2 for recurrence.

### `prism/`  (PRISM Repurposing secondary, figshare)
- `secondary-screen-dose-response-curve-parameters.csv` — figshare article 9393293
  (`ndownloader.figshare.com/files/20237739`, ~252 MB). Fitted ic50/ec50/slope +
  a FREE lower_limit (Emax) + name/moa/target/phase inline. Powers `warhead prism`.

### `pdxe/`  (Novartis PDX Encyclopedia, Gao et al. 2015)
- `pdxe_MOESM10.xlsx` — Nat Med supplement
  (`static-content.springer.com/esm/art%3A10.1038%2Fnm.3954/MediaObjects/41591_2015_BFnm3954_MOESM10_ESM.xlsx`,
  ~122 MB). Sheet `PCT curve metrics` = per model x treatment BestAvgResponse /
  ResponseCategory / target; `PCT raw data` = Model -> Tumor Type. In-vivo tumour
  response (no IC50/EC90; no HCC arm). Powers `warhead pdxe`.

### CTRP v2  (BLOCKED without R)
The classic `v20.*` files are gone from NCI; the only stable copy is a Zenodo
PharmacoGx PharmacoSet **.rds** (record 3905470, ~43 MB) that needs R/Bioconductor
to read (S4 object; not readable by pyreadr). The R-free Argonne IMPROVE CSA
`response.txt` mirror carries CTRP's ic50/ec50/hs/einf but only for the ~43-line /
78-drug cross-dataset intersection - too small for tissue selectivity. Add CTRP by
exporting the PharmacoSet slots to the canonical schema from an R session.

### NCI-60  (BLOCKED here; no liver lines)
CellMiner's `nci60_Drug_act.zip` endpoint returned an empty body (needs the portal
download form). NCI-60 also has NO liver cell line and only ~7 colon lines, so it
cannot support HCC selectivity regardless.

## Other sources (deferred gates)
See `WARHEAD.md` §2 and the docstring of each loader in `src/warhead/io/` for the
source URL, access route, and target tidy schema (CTRP, GDSC, NCI-60, Tahoe,
LINCS, JUMP, ADCdb, ChEMBL, COCONUT, NPAtlas, GTEx, HPA, FAERS, TCGA).

## Integrity
`warhead.io.download_file(url, dest, sha256=...)` verifies a checksum on download
and refuses a rotated/corrupt file. Record checksums when you pull a release.
