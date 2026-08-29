# Export raw measured dose-response for the top CTRP compounds (base R).
pset <- readRDS("data/raw/ctrp/CTRPv2.rds")
at <- attributes(pset); sens <- at$sensitivity
raw <- sens$raw
info <- as.data.frame(sens$info); info$exp <- rownames(info)
drug <- as.data.frame(at$drug)
top <- read.csv("data/interim/ctrp_top_drugs.csv", stringsAsFactors=FALSE)$compound
name2id <- setNames(drug$drugid, drug$cpd_name)
top_ids <- unname(name2id[top]); top_ids <- top_ids[!is.na(top_ids)]
sel <- info[info$drugid %in% top_ids, c("exp","drugid","cellid")]
cat("experiments for top drugs:", nrow(sel), "\n")
dn3 <- dimnames(raw)[[3]]
di <- match("Dose", dn3); vi <- match("Viability", dn3)
rows <- vector("list", nrow(sel))
for (k in seq_len(nrow(sel))) {
  e <- sel$exp[k]
  d <- raw[e, , di]; v <- raw[e, , vi]
  ok <- is.finite(d) & is.finite(v)
  if (!any(ok)) next
  rows[[k]] <- data.frame(drugid=sel$drugid[k], cellid=sel$cellid[k],
                          dose_uM=as.numeric(d[ok]), viability=as.numeric(v[ok]))
}
long <- do.call(rbind, rows)
id2name <- setNames(drug$cpd_name, drug$drugid)
long$drug <- id2name[as.character(long$drugid)]
write.csv(long, "data/interim/ctrp_curves.csv", row.names=FALSE)
cat("wrote", nrow(long), "points;", length(unique(long$drug)), "drugs; viability range",
    round(range(long$viability),2), "\n")
