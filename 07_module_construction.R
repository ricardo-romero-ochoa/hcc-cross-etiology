source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","tibble","ggplot2","tidyr","matrixStats"),
  bioc = c()
)

library(readr); library(dplyr); library(tibble); library(ggplot2); library(tidyr); library(matrixStats)

meta_out <- read_csv(file.path(ROOT,"meta/tables/Meta_DE_full.csv"), show_col_types=FALSE)
expr121 <- readRDS(file.path(ROOT,"results/matrices/GSE121248_expr_geneSymbol.rds"))
expr418 <- readRDS(file.path(ROOT,"results/matrices/GSE41804_expr_geneSymbol.rds"))
meta121 <- read_csv(file.path(ROOT,"data/metadata/GSE121248_metadata.csv"), show_col_types=FALSE)
meta418 <- read_csv(file.path(ROOT,"data/metadata/GSE41804_metadata.csv"), show_col_types=FALSE)

a121 <- align_expr_meta(expr121, meta121); expr121 <- a121$expr; meta121 <- a121$meta
a418 <- align_expr_meta(expr418, meta418); expr418 <- a418$expr; meta418 <- a418$meta

cons <- meta_out |> filter(conserved_FDR05) |> arrange(FDR_meta)
if (nrow(cons) < 200) stop("Too few conserved genes for stable module construction.")

max_genes <- 1500
cons <- cons |> slice(1:min(max_genes, nrow(cons)))
genes <- intersect(cons$gene, intersect(rownames(expr121), rownames(expr418)))

hub_connectivity <- function(expr, power=6) {
  X <- t(scale(t(expr)))
  C <- cor(t(X), use="pairwise.complete.obs")
  A <- abs(C)^power
  diag(A) <- 0
  rowSums(A)
}

k121 <- hub_connectivity(expr121[genes,,drop=FALSE], power=6)
k418 <- hub_connectivity(expr418[genes,,drop=FALSE], power=6)

hub <- tibble(
  gene=genes,
  k_121 = as.numeric(k121),
  k_418 = as.numeric(k418)
) |>
  mutate(rank_121 = rank(-k_121),
         rank_418 = rank(-k_418),
         rank_sum = rank_121 + rank_418) |>
  left_join(meta_out |> select(gene, beta, FDR_meta, I2), by="gene") |>
  arrange(rank_sum, FDR_meta)

write_csv(hub, file.path(ROOT,"meta/tables/HubCandidates_consensus_coexpression.csv"))

panelA <- hub |> filter(beta > 0, FDR_meta < 1e-6) |> arrange(rank_sum) |> slice(1:20) |> pull(gene)
panelB <- meta_out |> filter(conserved_FDR05, beta < 0) |> arrange(FDR_meta) |> slice(1:20) |> pull(gene)

writeLines(panelA, file.path(ROOT,"results/modules/PanelA_prolif_hubs_TOP20.txt"))
writeLines(panelB, file.path(ROOT,"results/modules/PanelB_hepatocyte_loss_TOP20.txt"))

# score discovery cohorts
make_score_df <- function(expr, meta, label) {
  meta$tissue_type <- factor(meta$tissue_type, levels=c("adjacent_or_nontumor","tumor"))
  tibble(
    gsm=colnames(expr),
    tissue_type=meta$tissue_type[match(colnames(expr), meta$gsm)],
    ProlifHubScore=score_module(expr, panelA),
    HepLossScore=score_module(expr, panelB)
  ) |>
    mutate(HCCStateScore = ProlifHubScore - HepLossScore, dataset = label)
}

df121 <- make_score_df(expr121, meta121, "GSE121248")
df418 <- make_score_df(expr418, meta418, "GSE41804")

write_csv(df121, file.path(ROOT,"meta/plots/GSE121248_module_scores.csv"))
write_csv(df418, file.path(ROOT,"meta/plots/GSE41804_module_scores.csv"))

plot_box <- function(df, out) {
  long <- df |> pivot_longer(cols=c(ProlifHubScore,HepLossScore,HCCStateScore), names_to="score", values_to="value")
  p <- ggplot(long, aes(x=tissue_type, y=value)) +
    geom_boxplot(outlier.size=0.4) +
    geom_jitter(width=0.15, alpha=0.35, size=0.6) +
    facet_wrap(~score, scales="free_y", ncol=3) +
    theme_classic(base_size=12) +
    labs(title=paste0(df$dataset[1], ": module scores"), x=NULL, y="mean z-score")
  ggsave(out, p, width=11, height=4.2, device=cairo_pdf)
}

plot_box(df121, file.path(ROOT,"meta/plots/GSE121248_module_scores_boxplots.pdf"))
plot_box(df418, file.path(ROOT,"meta/plots/GSE41804_module_scores_boxplots.pdf"))

message("Done. Panels + discovery module score plots generated.")