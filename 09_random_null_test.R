source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","ggplot2","pROC"),
  bioc = c()
)

library(readr); library(dplyr); library(ggplot2); library(pROC)

set.seed(1)

scores <- read_csv(file.path(ROOT,"meta/validation/GSE14520_module_scores.csv"), show_col_types=FALSE)

dat <- read.csv(file.path(ROOT,"data/external/GSE14520_expression_matrix_reordered.csv"),
                check.names=FALSE, stringsAsFactors=FALSE)
dat_gene <- dat[dat$`#NAME`!="#CLASS", , drop=FALSE]
expr <- as.matrix(dat_gene[, -1, drop=FALSE])
expr <- apply(expr, 2, as.numeric)
rownames(expr) <- dat_gene$`#NAME`
colnames(expr) <- colnames(dat_gene)[-1]

panelA <- readLines(file.path(ROOT,"results/modules/PanelA_prolif_hubs_TOP20.txt")) |> unique()
panelB <- readLines(file.path(ROOT,"results/modules/PanelB_hepatocyte_loss_TOP20.txt")) |> unique()
nA <- length(panelA); nB <- length(panelB)

obs_auc <- as.numeric(pROC::auc(pROC::roc(as.integer(scores$group=="tumor"), scores$HCCStateScore, quiet=TRUE)))

all_genes <- rownames(expr)
B <- 500

null_auc <- replicate(B, {
  gA <- sample(all_genes, nA)
  gB <- sample(all_genes, nB)
  s <- score_module(expr, gA) - score_module(expr, gB)
  as.numeric(pROC::auc(pROC::roc(as.integer(scores$group=="tumor"), s, quiet=TRUE)))
})

emp_p <- mean(null_auc >= obs_auc)
write_csv(data.frame(null_auc=null_auc), file.path(ROOT,"meta/validation/GSE14520_random_null_auc.csv"))

p <- ggplot(data.frame(null_auc=null_auc), aes(x=null_auc)) +
  geom_histogram(bins=40, color="black", fill="grey80") +
  geom_vline(xintercept=obs_auc, linetype="dashed", linewidth=0.8) +
  annotate("text", x=obs_auc, y=Inf,
           label=paste0("Observed AUC = ", round(obs_auc,3), "\nEmpirical p = ", format(emp_p, scientific=TRUE)),
           vjust=1.2, hjust=1.05, size=3.6) +
  theme_classic(base_size=12) +
  labs(title="GSE14520 matched-size random gene-set null", x="AUC (random matched gene sets)", y="Count")

ggsave(file.path(ROOT,"meta/validation/plots/GSE14520_random_null_AUC.pdf"), p,
       width=6.5, height=4.8, device=cairo_pdf)

message("Done. Empirical p = ", emp_p)