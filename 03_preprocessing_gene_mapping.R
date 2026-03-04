source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("dplyr","tibble","matrixStats"),
  bioc = c("Biobase","AnnotationDbi","hgu133plus2.db")
)

library(Biobase)
library(AnnotationDbi)
library(hgu133plus2.db)
library(matrixStats)
library(dplyr); library(tibble)

probe_to_symbol <- function(probes) {
  suppressMessages(
    AnnotationDbi::select(hgu133plus2.db, keys = unique(probes),
                          columns = c("SYMBOL"), keytype = "PROBEID")
  ) |>
    as_tibble() |>
    filter(!is.na(SYMBOL), SYMBOL != "") |>
    distinct(PROBEID, SYMBOL)
}

to_symbol_matrix_gpl570 <- function(eset) {
  expr <- Biobase::exprs(eset)
  probes <- rownames(expr)
  map <- probe_to_symbol(probes)
  symbols <- map$SYMBOL[match(probes, map$PROBEID)]
  collapse_to_genes_by_variance(expr, symbols)
}

for (gse in c("GSE121248","GSE41804","GSE83148","GSE38941")) {
  eset <- readRDS(file.path(ROOT, "data", "raw", paste0(gse, "_eset.rds")))
  message("Mapping probes→genes: ", gse)
  expr_g <- to_symbol_matrix_gpl570(eset)
  saveRDS(expr_g, file.path(ROOT, "results", "matrices", paste0(gse, "_expr_geneSymbol.rds")))
}

message("Done. Gene matrices saved in results/matrices/.")