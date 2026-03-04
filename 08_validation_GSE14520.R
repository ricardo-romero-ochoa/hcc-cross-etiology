source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","tibble","ggplot2","tidyr","pROC"),
  bioc = c()
)

library(readr); library(dplyr); library(tibble); library(ggplot2); library(tidyr); library(pROC)

panelA <- readLines(file.path(ROOT,"results/modules/PanelA_prolif_hubs_TOP20.txt")) |> unique()
panelB <- readLines(file.path(ROOT,"results/modules/PanelB_hepatocyte_loss_TOP20.txt")) |> unique()

mat_path <- file.path(ROOT,"data/external/GSE14520_expression_matrix_reordered.csv")
stopifnot(file.exists(mat_path))

dat <- read.csv(mat_path, check.names=FALSE, stringsAsFactors=FALSE)
cls <- dat[dat$`#NAME`=="#CLASS", -1, drop=FALSE] |> unlist() |> as.character()
dat_gene <- dat[dat$`#NAME`!="#CLASS", , drop=FALSE]

expr <- as.matrix(dat_gene[, -1, drop=FALSE])
expr <- apply(expr, 2, as.numeric)
rownames(expr) <- dat_gene$`#NAME`
colnames(expr) <- colnames(dat_gene)[-1]

group <- cls; names(group) <- colnames(expr)
meta <- tibble(sample=colnames(expr), group=group[match(colnames(expr), names(group))]) |>
  mutate(group = factor(group, levels=c("non-tumor","tumor")))
stopifnot(!any(is.na(meta$group)))

scores <- meta |>
  mutate(ProlifHubScore = score_module(expr, panelA),
         HepLossScore = score_module(expr, panelB),
         HCCStateScore = ProlifHubScore - HepLossScore)

write_csv(scores, file.path(ROOT,"meta/validation/GSE14520_module_scores.csv"))

stats_one <- function(score_name) {
  x_t <- scores |> filter(group=="tumor") |> pull(.data[[score_name]])
  x_n <- scores |> filter(group=="non-tumor") |> pull(.data[[score_name]])
  tibble(score=score_name,
         mean_non=mean(x_n, na.rm=TRUE),
         mean_tumor=mean(x_t, na.rm=TRUE),
         delta=mean(x_t, na.rm=TRUE)-mean(x_n, na.rm=TRUE),
         p_ttest=t.test(x_t, x_n)$p.value,
         cohens_d=cohens_d(x_t, x_n))
}

stats <- bind_rows(stats_one("ProlifHubScore"), stats_one("HepLossScore"), stats_one("HCCStateScore"))
roc_obj <- pROC::roc(as.integer(scores$group=="tumor"), scores$HCCStateScore, quiet=TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))
stats <- stats |> mutate(AUC_HCCState = if_else(score=="HCCStateScore", auc_val, NA_real_))
write_csv(stats, file.path(ROOT,"meta/validation/GSE14520_module_score_stats.csv"))

long <- scores |> pivot_longer(cols=c(ProlifHubScore,HepLossScore,HCCStateScore), names_to="score", values_to="value")
p <- ggplot(long, aes(x=group, y=value)) +
  geom_boxplot(outlier.size=0.4) +
  geom_jitter(width=0.15, alpha=0.35, size=0.6) +
  facet_wrap(~score, scales="free_y", ncol=3) +
  theme_classic(base_size=12) +
  labs(title="GSE14520 validation: module scores", x=NULL, y="mean z-score")

ggsave(file.path(ROOT,"meta/validation/plots/GSE14520_module_scores_boxplots.pdf"), p,
       width=11, height=4.2, device=cairo_pdf)

pdf(file.path(ROOT,"meta/validation/plots/GSE14520_HCCState_ROC.pdf"), width=5, height=5)
plot(roc_obj, main=paste0("GSE14520 HCCStateScore ROC (AUC=", round(auc_val,3), ")"))
dev.off()

message("Done. GSE14520 validation outputs generated.")