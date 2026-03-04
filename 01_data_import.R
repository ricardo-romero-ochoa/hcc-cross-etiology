source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr"),
  bioc = c("GEOquery","Biobase")
)

library(GEOquery); library(Biobase)
library(readr); library(dplyr)

gse_ids <- c("GSE121248","GSE41804","GSE83148","GSE38941")

get_eset_cached <- function(gse_id) {
  out_rds <- file.path(ROOT, "data", "raw", paste0(gse_id, "_eset.rds"))
  if (file.exists(out_rds)) return(readRDS(out_rds))
  message("Downloading ", gse_id, " (SeriesMatrix)...")
  gset_list <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
  eset <- gset_list[[1]]
  saveRDS(eset, out_rds)
  eset
}

for (gse in gse_ids) {
  eset <- get_eset_cached(gse)
  ph <- Biobase::pData(eset) |> as.data.frame()
  ph$gsm <- rownames(ph)
  write_csv(ph |> dplyr::select(gsm, dplyr::everything()),
            file.path(ROOT, "data", "raw", paste0(gse, "_pheno_raw.csv")))
}

message("Done. Cached ExpressionSets to data/raw/.")
message("If using a pre-assembled GSE14520 matrix, place it at: data/external/GSE14520_expression_matrix_reordered.csv")