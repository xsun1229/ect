#' Run PLINK-based SNP association analyses in parallel
#'
#' @description
#' Performs SNP-level association analyses for a given list of SNPs by invoking
#' PLINK in parallel. Each SNP is tested using its own covariate file, and
#' association results can optionally be grouped by phenotype factors.
#'
#' @param snp_list Character vector of SNP IDs to analyze.
#' @param folder_covar Path to the folder containing SNP-specific covariate files, in PLINK format
#'   (e.g., `"snp123.cov.txt"`).
#' @param folder_assoc Path to the output folder for PLINK association results.
#' @param plink_path Path to the PLINK executable.
#' @param file_bfile PLINK binary file prefix (`.bed/.bim/.fam`).
#' @param file_phe Path to the phenotype file to be analyzed by PLINK.
#'   This file should contain the phenotype columns (i.e., the factors)
#'   that are later referenced by `group_by_factor` to group and summarize
#'   association results.
#' @param single_factor Optional; a single phenotype name (column name from
#'   `file_phe`) to test individually. If `NULL`, all phenotypes in
#'   `file_phe` are analyzed (`--all-pheno` mode).
#' @param add_MAF Logical; if `TRUE`, the function merges association results
#'   with a MAF table. When `add_MAF = TRUE`, the argument `df_maf` **must**
#'   be provided. If `FALSE`, MAF information is ignored and `df_maf` can be
#'   left as `NULL`.
#' @param df_maf Data frame containing columns `SNP`, `A2`, and `MAF`.
#'   Required only when `add_MAF = TRUE`; otherwise leave it as `NULL`.
#' @param ncores Number of CPU cores to use for parallel processing.
#' @param cleanup Logical; if `TRUE`, removes intermediate `.assoc` files
#'   after grouping.
#' @param group_by_factor Logical; if `TRUE`, the function combines all SNP-level
#'   results and groups them by phenotype factors defined in `file_phe`.
#'
#' @return
#' Writes PLINK association results into `folder_assoc/` and (optionally)
#' grouped per-factor summary files (e.g., `"factor_pvalues.txt"`). Returns no
#' explicit object.
#'
#' @details
#' Each SNP is processed in its own subfolder under `folder_assoc/` to avoid
#' overwriting during parallel execution. Only additive test results
#' (`TEST == "ADD"`) are retained.
#'
#' @examples
#' \dontrun{
#' run_snps_assoc(
#'   snp_list = c("rs123", "rs456"),
#'   folder_covar = "covariates/",
#'   folder_assoc = "assoc_results/",
#'   plink_path = "/usr/bin/plink",
#'   file_bfile = "data/genotypes",
#'   file_phe = "data/pheno.txt",
#'   add_MAF = TRUE,
#'   df_maf = maf_table,
#'   ncores = 4
#' )
#' }
#'
#' @importFrom data.table fread fwrite rbindlist
#' @importFrom parallel mclapply
#' @export
run_snps_assoc <- function(
    snp_list,
    folder_covar,
    folder_assoc,
    plink_path,
    file_bfile,
    file_phe,
    single_factor = NULL,
    add_MAF = TRUE,
    df_maf = NULL,
    ncores = 4,
    cleanup = FALSE,
    group_by_factor = TRUE) {

  run_single <- function(snp) {
    cat(paste0("Processing ", snp, "\n"))

    file_covar <- file.path(folder_covar, paste0(snp, ".cov.txt"))
    folder_assoc_snp <- file.path(folder_assoc, snp)
    dir.create(folder_assoc_snp, recursive = TRUE, showWarnings = FALSE)
    file_out_snp <- file.path(folder_assoc_snp, snp)

    # Construct PLINK command
    if (is.null(single_factor)) {
      cmd <- sprintf(
        '%s --bfile %s --snp %s --covar %s --linear --allow-no-sex --pheno %s --all-pheno --out %s',
        plink_path, file_bfile, snp, file_covar, file_phe, file_out_snp
      )
    } else {
      cmd <- sprintf(
        '%s --bfile %s --snp %s --covar %s --linear --allow-no-sex --pheno %s --pheno-name %s --out %s',
        plink_path, file_bfile, snp, file_covar, file_phe, single_factor, file_out_snp
      )
    }

    try(system(cmd), silent = TRUE)

    files <- list.files(folder_assoc_snp, pattern = "\\.assoc\\.linear$", full.names = TRUE)
    if (length(files) == 0) {
      message(sprintf("No PLINK output for %s, skipping.", snp))
      return(NULL)
    }

    all_assoc <- tryCatch({
      data.table::rbindlist(lapply(files, function(f) {
        dt <- data.table::fread(f)
        dt[, factor := sub("\\.assoc\\.linear$", "", basename(f))]
        dt
      }))
    }, error = function(e) {
      message(sprintf("Error reading PLINK output for %s: %s", snp, e$message))
      return(NULL)
    })

    all_assoc <- all_assoc[TEST == "ADD"]

    if (!is.null(single_factor)) {
      all_assoc$factor <- single_factor
    } else {
      split_factor <- strsplit(all_assoc$factor, "\\.")
      all_assoc$factor <- sapply(split_factor, `[`, 2)
    }

    out_file <- file.path(folder_assoc, paste0(snp, ".assoc"))
    data.table::fwrite(all_assoc, file = out_file, sep = "\t", quote = FALSE)
    unlink(folder_assoc_snp, recursive = TRUE)
  }

  # Parallel execution
  system.time({
    parallel::mclapply(snp_list, run_single, mc.cores = ncores)
  })

  message("SNP association analyses complete.")

  # Group by factor
  if (group_by_factor) {
    assoc_files <- list.files(folder_assoc, pattern = "\\.assoc$", full.names = TRUE)
    all_data <- data.table::rbindlist(lapply(assoc_files, data.table::fread))
    factors <- unique(all_data$factor)

    for (f in factors) {
      dt_factor <- all_data[factor == f, .(SNP, CHR, BP, A1, BETA, STAT, P, NMISS)]
      dt_factor[, SE := ifelse(STAT == 0, NA, abs(BETA / STAT))]

      if (add_MAF) {
        if (is.null(df_maf))
          stop("df_maf must be provided when add_MAF = TRUE.")
        dt_factor <- merge(dt_factor, df_maf[, .(SNP, A2, MAF)],
                           by = "SNP", all.x = TRUE)
      }

      out_file <- file.path(folder_assoc, paste0(f, "_pvalues.txt"))
      data.table::fwrite(dt_factor, out_file, sep = "\t", quote = FALSE)
    }

    message("Grouped results by factor.")
  }

  # Cleanup
  if (cleanup) {
    assoc_files <- list.files(folder_assoc, pattern = "\\.assoc$", full.names = TRUE)
    suppressWarnings(file.remove(assoc_files))
    message("Clean-up done.")
  }
}
