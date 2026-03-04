source("R/_shared.R")
make_dirs()

bootstrap_packages(
  cran = c("dplyr","readr","stringr","tibble","tidyr","purrr"),
  bioc = c("Biobase")
)

library(Biobase)
library(dplyr); library(readr); library(stringr); library(tibble); library(tidyr); library(purrr)

load_eset <- function(gse_id) readRDS(file.path(ROOT, "data", "raw", paste0(gse_id, "_eset.rds")))

# ---- GSE83148 ----
standardize_83148 <- function(eset) {
  ph <- pheno_parsed_from_eset(eset, "GSE83148")
  ph |>
    mutate(
      study="GSE83148", platform="GPL570", etiology="HBV",
      cohort = case_when(str_detect(str_to_lower(source_name_ch1), "healthy") ~ "healthy_control",
                         TRUE ~ "hbv_hepatitis"),
      alt = na_if(alt, "NON"),
      ast = na_if(ast, "NON"),
      hbv_dna = na_if(hbv_dna, "NON"),
      alt_group = case_when(alt == ">40" ~ "high", alt == "<=40" ~ "normal", TRUE ~ NA_character_),
      ast_group = case_when(ast == ">35" ~ "high", ast == "<=35" ~ "normal", TRUE ~ NA_character_),
      hbv_dna_group = case_when(hbv_dna == ">10E6" ~ "high", hbv_dna == "<=10E6" ~ "low", TRUE ~ NA_character_),
      alt_high = case_when(alt_group=="high" ~ 1L, alt_group=="normal" ~ 0L, TRUE ~ NA_integer_),
      ast_high = case_when(ast_group=="high" ~ 1L, ast_group=="normal" ~ 0L, TRUE ~ NA_integer_),
      hbv_dna_high = case_when(hbv_dna_group=="high" ~ 1L, hbv_dna_group=="low" ~ 0L, TRUE ~ NA_integer_)
    )
}

# ---- GSE38941 ----
standardize_38941 <- function(eset) {
  ph <- pheno_parsed_from_eset(eset, "GSE38941")
  ph |>
    mutate(
      study="GSE38941", platform="GPL570", etiology="HBV",
      condition = case_when(
        disease_state == "normal" ~ "donor_normal",
        str_detect(str_to_lower(disease_state), "acute liver failure") ~ "hbv_alf",
        TRUE ~ "other"
      ),
      patient_id = readr::parse_number(subject),
      piece_id = suppressWarnings(as.integer(str_match(title, "patient\\s*\\d+-(\\d+)$")[,2])),
      block_id = if_else(condition=="hbv_alf", paste0("ALF_P", patient_id), paste0("DONOR_", gsm)),
      necrosis_group = case_when(
        condition != "hbv_alf" ~ NA_character_,
        patient_id %in% c(241,31) ~ "MHN",
        patient_id %in% c(219,32) ~ "SHN",
        TRUE ~ NA_character_
      )
    )
}

# ---- GSE121248 ----
standardize_121248 <- function(eset) {
  ph <- pheno_parsed_from_eset(eset, "GSE121248")
  blob <- row_blob(ph)
  ph |>
    mutate(
      study="GSE121248", platform="GPL570", etiology="HBV",
      tissue_type = vapply(blob, classify_tissue, character(1))
    )
}

# ---- GSE41804 ----
standardize_41804 <- function(eset) {
  ph <- pheno_parsed_from_eset(eset, "GSE41804")
  blob <- row_blob(ph)
  ph |>
    mutate(
      study="GSE41804", platform="GPL570", etiology="HCV",
      tissue_type = vapply(blob, classify_tissue, character(1)),
      il28b_genotype = vapply(blob, classify_il28b, character(1))
    )
}

meta831 <- standardize_83148(load_eset("GSE83148"))
meta389 <- standardize_38941(load_eset("GSE38941"))
meta121 <- standardize_121248(load_eset("GSE121248"))
meta418 <- standardize_41804(load_eset("GSE41804"))

write_csv(meta831, file.path(ROOT, "data", "metadata", "GSE83148_metadata.csv"))
write_csv(meta389, file.path(ROOT, "data", "metadata", "GSE38941_metadata.csv"))
write_csv(meta121, file.path(ROOT, "data", "metadata", "GSE121248_metadata.csv"))
write_csv(meta418, file.path(ROOT, "data", "metadata", "GSE41804_metadata.csv"))

message("Done. Standardized metadata written to data/metadata/.")