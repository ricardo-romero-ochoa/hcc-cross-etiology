source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","tibble","stringr","ggplot2"),
  bioc = c("Biobase","limma","AnnotationDbi","hgu133plus2.db")
)

library(readr); library(dplyr); library(tibble); library(stringr); library(ggplot2)
library(Biobase); library(limma); library(AnnotationDbi); library(hgu133plus2.db)

set.seed(1)

eset831 <- readRDS(file.path(ROOT,"data/raw/GSE83148_eset.rds"))
meta831 <- read_csv(file.path(ROOT,"data/metadata/GSE83148_metadata.csv"), show_col_types=FALSE)

expr121 <- readRDS(file.path(ROOT,"results/matrices/GSE121248_expr_geneSymbol.rds"))
meta121 <- read_csv(file.path(ROOT,"data/metadata/GSE121248_metadata.csv"), show_col_types=FALSE)
gsva121 <- readRDS(file.path(ROOT,"results/gsva/GSE121248_Hallmark_GSVA.rds"))

a121 <- align_expr_meta(expr121, meta121); expr121 <- a121$expr; meta121 <- a121$meta

probe_to_symbol <- function(probes) {
  suppressMessages(
    AnnotationDbi::select(hgu133plus2.db, keys=unique(probes),
                          columns=c("SYMBOL"), keytype="PROBEID")
  ) |> as.data.frame()
}

tt_trait <- function(expr_probe, md, flag_col, trait_name) {
  idx <- which(md$cohort=="hbv_hepatitis" & !is.na(md[[flag_col]]))
  x <- expr_probe[, idx, drop=FALSE]
  m <- md[idx, , drop=FALSE]
  m[[flag_col]] <- factor(m[[flag_col]], levels=c(0,1))
  design <- model.matrix(~ m[[flag_col]])
  colnames(design) <- c("Intercept", paste0(trait_name,"_high"))
  tt <- topTable(eBayes(lmFit(x, design)), coef=paste0(trait_name,"_high"), number=Inf, sort.by="P") |>
    rownames_to_column("probe_id")
  map <- probe_to_symbol(tt$probe_id)
  tt |> left_join(map, by=c("probe_id"="PROBEID")) |> filter(!is.na(SYMBOL), SYMBOL!="")
}

pick_up_down <- function(tt, fdr=0.05, min_abs=0.3, max_genes=200) {
  sig <- tt |> filter(adj.P.Val < fdr, abs(logFC) >= min_abs)
  list(
    up = sig |> filter(logFC>0) |> arrange(adj.P.Val) |> pull(SYMBOL) |> unique() |> head(max_genes),
    down = sig |> filter(logFC<0) |> arrange(adj.P.Val) |> pull(SYMBOL) |> unique() |> head(max_genes)
  )
}

expr831 <- Biobase::exprs(eset831)
tt_alt <- tt_trait(expr831, meta831, "alt_high", "ALT")
tt_ast <- tt_trait(expr831, meta831, "ast_high", "AST")

sig_alt <- pick_up_down(tt_alt)
sig_ast <- pick_up_down(tt_ast)

inj_up <- intersect(sig_alt$up, sig_ast$up)
inj_dn <- intersect(sig_alt$down, sig_ast$down)
if (length(inj_up) < 20) inj_up <- unique(c(sig_alt$up, sig_ast$up))
if (length(inj_dn) < 20) inj_dn <- unique(c(sig_alt$down, sig_ast$down))

writeLines(inj_up, file.path(ROOT,"results/programs/HBV_INJURY_UP_genes.txt"))
writeLines(inj_dn, file.path(ROOT,"results/programs/HBV_INJURY_DN_genes.txt"))

hbv_injury <- signed_zscore(expr121, inj_up, inj_dn, min_genes=10)

df <- tibble(gsm=colnames(expr121), HBV_INJURY=hbv_injury) |>
  left_join(meta121 |> select(gsm,tissue_type), by="gsm") |>
  mutate(tissue_type=factor(tissue_type, levels=c("adjacent_or_nontumor","tumor")))

# adjustment with E2F/G2M Hallmarks
gsva_mat <- gsva121
df$E2F <- as.numeric(gsva_mat["HALLMARK_E2F_TARGETS", match(df$gsm, colnames(gsva_mat))])
df$G2M <- as.numeric(gsva_mat["HALLMARK_G2M_CHECKPOINT", match(df$gsm, colnames(gsva_mat))])
df2 <- df |> filter(!is.na(tissue_type), is.finite(HBV_INJURY), is.finite(E2F), is.finite(G2M))

m0 <- lm(HBV_INJURY ~ tissue_type, data=df2)
m1 <- lm(HBV_INJURY ~ tissue_type + E2F + G2M, data=df2)

tab1C <- tibble(
  model=c("unadjusted","adjusted_E2F_G2M"),
  estimate=c(coef(m0)["tissue_typetumor"], coef(m1)["tissue_typetumor"]),
  p_value=c(summary(m0)$coefficients["tissue_typetumor","Pr(>|t|)"],
            summary(m1)$coefficients["tissue_typetumor","Pr(>|t|)"])
)
write_csv(tab1C, file.path(ROOT,"paper_package/tables/Table_1C_Injury_Adjustment_E2F_G2M.csv"))

# residualized for plotting
m_adj <- lm(HBV_INJURY ~ E2F + G2M, data=df2)
df2$HBV_INJURY_resid <- resid(m_adj)
write_csv(df2, file.path(ROOT,"results/programs/GSE121248_HBV_INJURY_residualized.csv"))

p_raw <- ggplot(df2, aes(x=tissue_type, y=HBV_INJURY)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.15, alpha=0.45, size=0.7) +
  theme_classic(base_size=12) + labs(title="GSE121248: HBV_INJURY (raw)", x=NULL, y="score")
p_res <- ggplot(df2, aes(x=tissue_type, y=HBV_INJURY_resid)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.3) +
  geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.15, alpha=0.45, size=0.7) +
  theme_classic(base_size=12) + labs(title="GSE121248: HBV_INJURY residualized by E2F+G2M", x=NULL, y="residual")

ggsave(file.path(ROOT,"results/figures/Fig_GSE121248_injury_raw.pdf"), p_raw, width=5.8, height=4.4, device=cairo_pdf)
ggsave(file.path(ROOT,"results/figures/Fig_GSE121248_injury_residualized.pdf"), p_res, width=5.8, height=4.4, device=cairo_pdf)

# Table 1B (score effects in HBV-HCC: injury only here; module score effects are added later)
score_mat <- matrix(df2$HBV_INJURY, nrow=1); rownames(score_mat) <- "HBV_INJURY"; colnames(score_mat) <- df2$gsm
fit <- eBayes(lmFit(score_mat, model.matrix(~ tissue_type, data=df2)))
tab1B <- topTable(fit, coef="tissue_typetumor", number=Inf, sort.by="P") |> rownames_to_column("score")
write_csv(tab1B, file.path(ROOT,"paper_package/tables/Table_1B_TopScoreEffects_GSE121248.csv"))

message("Done. Injury axis outputs written to results/ and paper_package/tables/.")