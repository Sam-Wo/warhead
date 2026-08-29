# Read PharmacoSet slots as raw attributes (S4 slots are stored as attributes),
# so no PharmacoGx class definition is required.
pset <- readRDS("data/raw/ctrp/CTRPv2.rds")
at <- attributes(pset)
cat("slot names:\n"); print(setdiff(names(at), "class"))
sens <- at$sensitivity
cat("\n== sensitivity names ==\n"); print(names(sens))
cat("\n-- profiles cols --\n"); print(colnames(sens$profiles)); print(head(sens$profiles, 3))
cat("\n-- info cols --\n"); print(colnames(sens$info)); print(head(sens$info, 2))
cat("\n-- raw dims/dimnames --\n"); print(dim(sens$raw)); print(lapply(dimnames(sens$raw), head, 3))
cell <- if (!is.null(at$cell)) at$cell else at$sample
cat("\n== cell cols ==\n"); print(colnames(cell))
tissue_cols <- grep("tissue|site|hist|lineage|primary", colnames(cell), ignore.case=TRUE, value=TRUE)
cat("tissue-ish cols:", tissue_cols, "\n"); print(head(cell[, tissue_cols, drop=FALSE], 3))
drug <- if (!is.null(at$drug)) at$drug else at$treatment
cat("\n== drug cols ==\n"); print(colnames(drug))
tgt_cols <- grep("target|gene|activity|moa", colnames(drug), ignore.case=TRUE, value=TRUE)
cat("target-ish cols:", tgt_cols, "\n"); print(head(drug[, c(colnames(drug)[1], tgt_cols), drop=FALSE], 3))
