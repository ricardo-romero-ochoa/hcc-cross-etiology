source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("readr","dplyr","tibble","ggplot2","tidyr","qpdf"),
  bioc = c()
)

library(readr); library(dplyr); library(tibble); library(ggplot2); library(tidyr); library(qpdf)

# Figure 2: Hallmark cross-etiology dot plot
tt121 <- read_csv(file.path(ROOT,"results/gsva/GSE121248_GSVA_tumor_vs_adjacent.csv"), show_col_types=FALSE) |>
  mutate(dataset="HBV-HCC (GSE121248)")
tt418 <- read_csv(file.path(ROOT,"results/gsva/GSE41804_GSVA_tissue_main.csv"), show_col_types=FALSE) |>
  mutate(dataset="HCV-HCC (GSE41804)")

topN <- 15
pw <- unique(c(tt121 |> arrange(adj.P.Val) |> slice(1:topN) |> pull(pathway),
               tt418 |> arrange(adj.P.Val) |> slice(1:topN) |> pull(pathway)))

df <- bind_rows(tt121, tt418) |>
  filter(pathway %in% pw) |>
  mutate(pathway = factor(pathway, levels=rev(pw)))

p2 <- ggplot(df, aes(x=dataset, y=pathway)) +
  geom_point(aes(size=-log10(adj.P.Val), color=logFC), alpha=0.9) +
  theme_classic(base_size=12) +
  labs(title="Cross-etiology Hallmark program shifts (tumor vs adjacent/non-tumor)",
       x=NULL, y=NULL, color="logFC", size=expression(-log[10](FDR)))

ggsave(file.path(ROOT,"paper_package/figures/Figure_2_Hallmark_cross_etiology.pdf"),
       p2, width=8.5, height=6, device=cairo_pdf)

# Figure 3: injury axis (combine raw + residualized)
inj_raw <- file.path(ROOT,"results/figures/Fig_GSE121248_injury_raw.pdf")
inj_res <- file.path(ROOT,"results/figures/Fig_GSE121248_injury_residualized.pdf")
if (file.exists(inj_raw) && file.exists(inj_res)) {
  qpdf::pdf_nup(c(inj_raw, inj_res),
                output=file.path(ROOT,"paper_package/figures/Figure_3_HBV_injury_axis.pdf"),
                nup="2x1")
}

# Figure 4: discovery module score boxplots (combine two PDFs)
b121 <- file.path(ROOT,"meta/plots/GSE121248_module_scores_boxplots.pdf")
b418 <- file.path(ROOT,"meta/plots/GSE41804_module_scores_boxplots.pdf")
if (file.exists(b121) && file.exists(b418)) {
  qpdf::pdf_combine(c(b121,b418),
                    output=file.path(ROOT,"paper_package/figures/Figure_4_Module_scores_discovery.pdf"))
}

# Figure 5: GSE14520 validation (A,B,C sequential pages)
vA <- file.path(ROOT,"meta/validation/plots/GSE14520_module_scores_boxplots.pdf")
vB <- file.path(ROOT,"meta/validation/plots/GSE14520_HCCState_ROC.pdf")
vC <- file.path(ROOT,"meta/validation/plots/GSE14520_random_null_AUC.pdf")
keep <- c(vA,vB,vC)[file.exists(c(vA,vB,vC))]
if (length(keep) > 0) {
  qpdf::pdf_combine(keep,
                    output=file.path(ROOT,"paper_package/figures/Figure_5_GSE14520_validation.pdf"))
}

# Copy/ensure key tables exist in paper_package/tables
# (Table_1A created in 04; Table_1B + 1C in 05)
message("Done. Figures assembled in paper_package/figures/.")
writeLines(capture.output(sessionInfo()), file.path(ROOT,"paper_package/text/sessionInfo.txt"))