# Export CTRPv2 PharmacoSet -> a tidy CSV for the WARHEAD Python loader.
# Base R only (reads S4 slots as attributes); no PharmacoGx required.
pset <- readRDS("data/raw/ctrp/CTRPv2.rds")
at <- attributes(pset); sens <- at$sensitivity

prof <- as.data.frame(sens$profiles); prof$exp <- rownames(prof)
info <- as.data.frame(sens$info);     info$exp <- rownames(info)
m <- merge(prof, info[, c("exp","cellid","drugid","chosen.min.range","chosen.max.range",
                          "Number of Doses Tested")], by="exp")

cell <- as.data.frame(at$cell)
cell_sub <- unique(cell[, c("cellid","ccl_name","ccle_primary_site","tissueid")])
drug <- as.data.frame(at$drug)
drug_sub <- unique(drug[, c("drugid","cpd_name","gene_symbol_of_protein_target",
                            "target_or_activity_of_compound","FDA")])

m <- merge(m, cell_sub, by="cellid")
m <- merge(m, drug_sub, by="drugid")

out <- data.frame(
  cellid = m$cellid, ccl_name = m$ccl_name,
  drug = m$cpd_name, target = m$gene_symbol_of_protein_target,
  moa = m$target_or_activity_of_compound, fda = m$FDA,
  primary_site = m$ccle_primary_site, tissueid = m$tissueid,
  ic50_uM = m$ic50_recomputed, ec50_uM = m$EC50, HS = m$HS, E_inf = m$E_inf,
  min_conc_uM = m$chosen.min.range, max_conc_uM = m$chosen.max.range,
  stringsAsFactors = FALSE)
dir.create("data/interim", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "data/interim/ctrp_export.csv", row.names = FALSE)
cat("wrote", nrow(out), "rows;",
    length(unique(out$drug)), "drugs,", length(unique(out$cellid)), "lines\n")
cat("liver lines:", length(unique(out$cellid[out$primary_site=="liver"])),
    "| large_intestine lines:", length(unique(out$cellid[out$primary_site=="large_intestine"])), "\n")
