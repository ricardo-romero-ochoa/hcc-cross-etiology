# ============================================================
# Shared utilities for HCC cross-etiology pipeline
# ============================================================

options(stringsAsFactors = FALSE)

bootstrap_packages <- function(cran = character(), bioc = character()) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  for (p in bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

ROOT <- normalizePath(".", winslash = "/", mustWork = FALSE)

make_dirs <- function() {
  dirs <- c(
    "data/raw", "data/metadata", "data/external",
    "results/matrices", "results/gsva", "results/programs", "results/modules", "results/figures",
    "meta/tables", "meta/plots", "meta/validation/plots",
    "paper_package/figures", "paper_package/tables", "paper_package/text"
  )
  for (d in dirs) dir.create(file.path(ROOT, d), showWarnings = FALSE, recursive = TRUE)
}

# -------- phenotype parsing helpers --------
parse_kv <- function(x) {
  x <- stringr::str_squish(as.character(x))
  if (is.na(x) || x == "") return(tibble::tibble(key = NA_character_, value = NA_character_))
  if (stringr::str_detect(x, ":")) {
    parts <- stringr::str_split_fixed(x, ":", 2)
    key <- stringr::str_squish(parts[,1])
    val <- stringr::str_squish(parts[,2])
  } else {
    key <- "characteristics_misc"
    val <- x
  }
  key <- key |>
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("^_|_$", "")
  tibble::tibble(key = key, value = val)
}

pheno_parsed_from_eset <- function(eset, gse_id) {
  ph <- Biobase::pData(eset) |> tibble::as_tibble(rownames = "gsm") |> dplyr::mutate(gse = gse_id, .before = 1)
  char_cols <- grep("^characteristics_ch1", colnames(ph), value = TRUE)
  if (length(char_cols) == 0) return(ph)

  char_long <- ph |>
    dplyr::select(gse, gsm, dplyr::all_of(char_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(char_cols), names_to = "char_field", values_to = "char_value") |>
    dplyr::filter(!is.na(char_value), stringr::str_squish(char_value) != "") |>
    dplyr::mutate(tmp = purrr::map(char_value, parse_kv)) |>
    tidyr::unnest(tmp) |>
    dplyr::filter(!is.na(key), key != "")

  char_wide <- char_long |>
    dplyr::group_by(gse, gsm, key) |>
    dplyr::summarise(value = paste(unique(value), collapse = " | "), .groups = "drop") |>
    tidyr::pivot_wider(names_from = key, values_from = value)

  core_cols <- intersect(c("title","source_name_ch1","organism_ch1"), colnames(ph))

  ph |>
    dplyr::select(gse, gsm, dplyr::any_of(core_cols), dplyr::everything()) |>
    dplyr::select(-dplyr::all_of(char_cols)) |>
    dplyr::left_join(char_wide, by = c("gse","gsm"))
}

row_blob <- function(df) {
  apply(df, 1, function(r) {
    r <- r[!is.na(r)]
    paste(r, collapse = " | ") |> stringr::str_to_lower()
  })
}

classify_tissue <- function(blob) {
  b <- stringr::str_to_lower(blob)
  if (stringr::str_detect(b,
    "non[- ]?tumou?r|nontumou?r|non[- ]?tumoral|non[- ]?cancer|noncancerous|adjacent|surrounding|paired\\s*normal|normal\\s*(liver|tissue)|non[- ]?neoplastic"
  )) return("adjacent_or_nontumor")
  if (stringr::str_detect(b,
    "tumou?r\\s*tissue|cancerous\\s*tissue|carcinoma\\s*tissue|\\bhcc\\b|\\btumou?r\\b|\\bcancer\\b"
  )) return("tumor")
  NA_character_
}

classify_il28b <- function(blob) {
  b <- stringr::str_to_lower(blob)
  if (stringr::str_detect(b, "\\btt\\b")) return("TT")
  if (stringr::str_detect(b, "tg/?gg|\\btg\\b|\\bgg\\b")) return("TG_GG")
  NA_character_
}

align_expr_meta <- function(expr, meta, id_col = "gsm") {
  common <- intersect(colnames(expr), meta[[id_col]])
  expr2 <- expr[, common, drop = FALSE]
  meta2 <- meta |> dplyr::filter(.data[[id_col]] %in% common) |> dplyr::arrange(match(.data[[id_col]], common))
  stopifnot(all(meta2[[id_col]] == colnames(expr2)))
  list(expr = expr2, meta = meta2)
}

# -------- expression mapping --------
collapse_to_genes_by_variance <- function(expr, symbols) {
  ok <- !is.na(symbols) & symbols != ""
  expr <- expr[ok, , drop = FALSE]
  symbols <- symbols[ok]

  v <- matrixStats::rowVars(expr)
  df <- tibble::tibble(probe = rownames(expr), gene = symbols, var = v) |>
    dplyr::group_by(gene) |>
    dplyr::slice_max(order_by = var, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  expr_g <- expr[df$probe, , drop = FALSE]
  rownames(expr_g) <- df$gene
  expr_g
}

# -------- scoring --------
signed_zscore <- function(expr_sym, up_genes, down_genes, min_genes = 10) {
  up <- intersect(up_genes, rownames(expr_sym))
  dn <- intersect(down_genes, rownames(expr_sym))
  if (length(up) < min_genes && length(dn) < min_genes) return(rep(NA_real_, ncol(expr_sym)))
  z <- t(scale(t(expr_sym)))
  up_score <- if (length(up) >= min_genes) colMeans(z[up, , drop = FALSE]) else 0
  dn_score <- if (length(dn) >= min_genes) colMeans(z[dn, , drop = FALSE]) else 0
  as.numeric(up_score - dn_score)
}

score_module <- function(expr_gene, genes, min_genes = 5) {
  g <- intersect(genes, rownames(expr_gene))
  if (length(g) < min_genes) return(rep(NA_real_, ncol(expr_gene)))
  Z <- t(scale(t(expr_gene[g, , drop = FALSE])))
  colMeans(Z, na.rm = TRUE)
}

cohens_d <- function(x_tumor, x_non) {
  x_tumor <- x_tumor[is.finite(x_tumor)]
  x_non <- x_non[is.finite(x_non)]
  m1 <- mean(x_tumor); m0 <- mean(x_non)
  s1 <- sd(x_tumor); s0 <- sd(x_non)
  sp <- sqrt(((length(x_tumor)-1)*s1^2 + (length(x_non)-1)*s0^2) / (length(x_tumor)+length(x_non)-2))
  (m1 - m0) / sp
}

# -------- GSVA compatibility --------
gsva_compat <- function(expr, gsets, method = c("gsva","ssgsea")) {
  method <- match.arg(method)
  expr <- as.matrix(expr)

  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    if (method == "gsva") {
      param <- GSVA::gsvaParam(exprData = expr, geneSets = gsets)
    } else {
      param <- GSVA::ssgseaParam(exprData = expr, geneSets = gsets)
    }
    res <- tryCatch(GSVA::gsva(param, verbose = FALSE), error = function(e) GSVA::gsva(param))
    if (is.matrix(res)) return(res)
    m <- tryCatch(SummarizedExperiment::assay(res), error = function(e) NULL)
    if (!is.null(m)) return(as.matrix(m))
    return(as.matrix(res))
  }

  out <- GSVA::gsva(expr, gsets, method = method, kcdf = "Gaussian", verbose = FALSE)
  as.matrix(out)
}

run_gsva_safe <- function(expr, gsets) {
  tryCatch(gsva_compat(expr, gsets, method="gsva"),
           error = function(e) gsva_compat(expr, gsets, method="ssgsea"))
}
