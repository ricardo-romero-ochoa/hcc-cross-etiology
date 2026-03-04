source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","tibble","ggplot2","metafor"),
  bioc = c("limma")
)

library(readr); library(dplyr); library(tibble); library(ggplot2); library(metafor)
library(limma)

expr121 <- readRDS(file.path(ROOT,"results/matrices/GSE121248_expr_geneSymbol.rds"))
expr418 <- readRDS(file.path(ROOT,"results/matrices/GSE41804_expr_geneSymbol.rds"))
meta121 <- read_csv(file.path(ROOT,"data/metadata/GSE121248_metadata.csv"), show_col_types=FALSE)
meta418 <- read_csv(file.path(ROOT,"data/metadata/GSE41804_metadata.csv"), show_col_types=FALSE)

a121 <- align_expr_meta(expr121, meta121); expr121 <- a121$expr; meta121 <- a121$meta
a418 <- align_expr_meta(expr418, meta418); expr418 <- a418$expr; meta418 <- a418$meta

run_limma_gene <- function(expr, meta, name) {
  meta$tissue_type <- factor(meta$tissue_type, levels=c("adjacent_or_nontumor","tumor"))
  keep <- which(!is.na(meta$tissue_type))
  m <- meta[keep, , drop=FALSE]
  x <- expr[, m$gsm, drop=FALSE]
  fit <- eBayes(lmFit(x, model.matrix(~ tissue_type, data=m)))
  topTable(fit, coef="tissue_typetumor", number=Inf, sort.by="none") |>
    rownames_to_column("gene") |>
    transmute(dataset=name, gene, logFC, t, P.Value, adj.P.Val, AveExpr,
              SE = abs(logFC / t))
}

tt121 <- run_limma_gene(expr121, meta121, "GSE121248_HBV")
tt418 <- run_limma_gene(expr418, meta418, "GSE41804_HCV")

write_csv(tt121, file.path(ROOT,"meta/tables/GSE121248_gene_limma.csv"))
write_csv(tt418, file.path(ROOT,"meta/tables/GSE41804_gene_limma.csv"))

comb <- tt121 |> select(gene, logFC_121=logFC, SE_121=SE) |>
  inner_join(tt418 |> select(gene, logFC_418=logFC, SE_418=SE), by="gene") |>
  filter(is.finite(SE_121), is.finite(SE_418), SE_121>0, SE_418>0)

meta_one <- function(b1,se1,b2,se2) {
  yi <- c(b1,b2); vi <- c(se1^2,se2^2)
  fit <- tryCatch(metafor::rma(yi=yi, vi=vi, method="REML"),
                  error=function(e) metafor::rma(yi=yi, vi=vi, method="FE"))
  beta <- as.numeric(fit$b); se <- as.numeric(fit$se)
  z <- beta/se; p <- 2*pnorm(-abs(z))
  c(beta=beta, se=se, z=z, p=p, tau2=as.numeric(fit$tau2), I2=as.numeric(fit$I2))
}

meta_mat <- mapply(meta_one, comb$logFC_121, comb$SE_121, comb$logFC_418, comb$SE_418)
meta_df <- as.data.frame(t(meta_mat)) |> mutate(gene=comb$gene) |> relocate(gene)
meta_df$FDR_meta <- p.adjust(meta_df$p, method="BH")

meta_out <- comb |>
  left_join(meta_df, by="gene") |>
  mutate(direction_concordant = sign(logFC_121) == sign(logFC_418),
         conserved_FDR05 = direction_concordant & (FDR_meta < 0.05)) |>
  arrange(FDR_meta)

write_csv(meta_out, file.path(ROOT,"meta/tables/Meta_DE_full.csv"))
write_csv(meta_out |> filter(conserved_FDR05), file.path(ROOT,"meta/tables/Meta_DE_conserved_FDR05.csv"))

p <- meta_out |> mutate(neglog10=-log10(p)) |>
  ggplot(aes(beta, neglog10)) +
  geom_point(alpha=0.30, size=1.1) +
  geom_point(data=\(d) d |> filter(conserved_FDR05), alpha=0.9, size=1.1) +
  theme_classic(base_size=12) +
  labs(title="Meta-analysis: tumor vs adjacent/non-tumor (HBV+HCV)", x="meta logFC", y=expression(-log[10](p)))

ggsave(file.path(ROOT,"meta/plots/Meta_volcano.pdf"), p, width=7, height=5, device=cairo_pdf)

message("Done. Meta-analysis outputs in meta/.")