source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","tibble"),
  bioc = c("limma","GSVA","msigdbr","SummarizedExperiment")
)

library(readr); library(dplyr); library(tibble)
library(limma); library(GSVA); library(msigdbr)

expr121 <- readRDS(file.path(ROOT,"results/matrices/GSE121248_expr_geneSymbol.rds"))
expr418 <- readRDS(file.path(ROOT,"results/matrices/GSE41804_expr_geneSymbol.rds"))
meta121 <- read_csv(file.path(ROOT,"data/metadata/GSE121248_metadata.csv"), show_col_types=FALSE)
meta418 <- read_csv(file.path(ROOT,"data/metadata/GSE41804_metadata.csv"), show_col_types=FALSE)

a121 <- align_expr_meta(expr121, meta121); expr121 <- a121$expr; meta121 <- a121$meta
a418 <- align_expr_meta(expr418, meta418); expr418 <- a418$expr; meta418 <- a418$meta

hallmark_df <- msigdbr(species="Homo sapiens", category="H") |>
  dplyr::select(gs_name, gene_symbol) |> distinct()
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

gsva121 <- run_gsva_safe(expr121, hallmark_sets)
gsva418 <- run_gsva_safe(expr418, hallmark_sets)

saveRDS(gsva121, file.path(ROOT,"results/gsva/GSE121248_Hallmark_GSVA.rds"))
saveRDS(gsva418, file.path(ROOT,"results/gsva/GSE41804_Hallmark_GSVA.rds"))

write_csv(as.data.frame(gsva121) |> tibble::rownames_to_column("pathway"),
          file.path(ROOT,"results/gsva/GSE121248_Hallmark_GSVA.csv"))
write_csv(as.data.frame(gsva418) |> tibble::rownames_to_column("pathway"),
          file.path(ROOT,"results/gsva/GSE41804_Hallmark_GSVA.csv"))

# limma contrasts
meta121$tissue_type <- factor(meta121$tissue_type, levels=c("adjacent_or_nontumor","tumor"))
k121 <- which(!is.na(meta121$tissue_type))
fit121 <- eBayes(lmFit(gsva121[,k121,drop=FALSE], model.matrix(~ tissue_type, data=meta121[k121,])))
tt121 <- topTable(fit121, coef="tissue_typetumor", number=Inf, sort.by="P") |> rownames_to_column("pathway")
write_csv(tt121, file.path(ROOT,"results/gsva/GSE121248_GSVA_tumor_vs_adjacent.csv"))

meta418$tissue_type <- factor(meta418$tissue_type, levels=c("adjacent_or_nontumor","tumor"))
meta418$il28b_genotype <- factor(meta418$il28b_genotype, levels=c("TG_GG","TT"))
k418 <- which(!is.na(meta418$tissue_type) & !is.na(meta418$il28b_genotype))
design418 <- model.matrix(~ tissue_type * il28b_genotype, data=meta418[k418,])
fit418 <- eBayes(lmFit(gsva418[,k418,drop=FALSE], design418))

tt418_main <- topTable(fit418, coef="tissue_typetumor", number=Inf, sort.by="P") |> rownames_to_column("pathway")
tt418_int  <- topTable(fit418, coef="tissue_typetumor:il28b_genotypeTT", number=Inf, sort.by="P") |> rownames_to_column("pathway")

write_csv(tt418_main, file.path(ROOT,"results/gsva/GSE41804_GSVA_tissue_main.csv"))
write_csv(tt418_int,  file.path(ROOT,"results/gsva/GSE41804_GSVA_interaction.csv"))

# Paper Table 1A (top 10 each)
topN <- 10
tab1A <- bind_rows(
  tt121 |> arrange(adj.P.Val) |> slice(1:topN) |> mutate(dataset="GSE121248 (HBV-HCC)"),
  tt418_main |> arrange(adj.P.Val) |> slice(1:topN) |> mutate(dataset="GSE41804 (HCV-HCC)")
) |> select(dataset, pathway, logFC, P.Value, adj.P.Val)

write_csv(tab1A, file.path(ROOT,"paper_package/tables/Table_1A_TopHallmarks.csv"))

message("Done. Hallmark scoring + contrasts written to results/gsva/ and paper_package/tables/.")