#' Prune SNPs based on clumping and LD thresholds using PLINK
#'
#' @description
#' Performs SNP pruning in two stages: (1) PLINK-based clumping to select
#' approximately independent lead SNPs, and (2) additional distance- and
#' LD-based pruning to further reduce correlated variants.
#'
#' The function generates temporary files under a subdirectory of `out.dir`,
#' invokes PLINK commands, and returns a pruned GWAS summary statistics table.
#'
#' @param file_sumstats Path to the GWAS summary statistics file containing
#'   columns `SNP`, `CHR`, `BP`, and `P`.
#' @param plink_path Path to the PLINK executable.
#' @param bfile PLINK binary file prefix (`.bed/.bim/.fam`) used as the LD reference.
#' @param out.dir Output directory where intermediate and result files will be written.
#' @param out.pref Prefix for temporary PLINK output and final file naming.
#' @param clean_up Logical; if `TRUE`, temporary files and subdirectories are deleted
#'   after successful completion or if errors occur.
#'
#' @return
#' A `data.table` of pruned GWAS summary statistics containing the same columns as
#' the input file but restricted to SNPs that remain after distance- and LD-based
#' pruning. Returns `NULL` if PLINK clumping fails to produce output.
#'
#' @details
#' The pruning proceeds in the following steps:
#' 1. **PLINK clumping:** using thresholds `--clump-p1 5e-8`, `--clump-p2 1e-3`,
#'    `--clump-r2 0.1`, and `--clump-kb 1000`.
#' 2. **Distance-based pruning:** removes SNPs within ±1 Mb of a more significant SNP.
#' 3. **LD-based pruning:** calculates pairwise LD (`r²`) using PLINK and removes
#'    SNPs with `r² > 0.1` relative to previously retained variants.
#' 4. Returns GWAS summary statistics restricted to the remaining independent SNPs.
#'
#' @examples
#' \dontrun{
#' pruned_gwas <- prune_snps(
#'   file_sumstats = "assoc/gwas.txt",
#'   plink_path = "/usr/bin/plink",
#'   bfile = "ref/1000G_EUR",
#'   out.dir = "pruning_results/",
#'   out.pref = "GWAS1",
#'   clean_up = TRUE
#' )
#' }
#'
#' @importFrom data.table fread
#' @export
prune_snps <- function(file_sumstats,
                       plink_path,
                       bfile,
                       out.dir,
                       out.pref,
                       clean_up = FALSE) {

  #---------------------------------------------------------------
  # Setup
  #---------------------------------------------------------------
  tmp.dir <- file.path(out.dir, out.pref)
  if (!dir.exists(tmp.dir)) dir.create(tmp.dir, recursive = TRUE)
  tmp.out.pref <- file.path(tmp.dir, out.pref)

  #---------------------------------------------------------------
  # Step 1: PLINK clumping
  #---------------------------------------------------------------
  cmd <- sprintf(
    '%s --bfile %s --clump %s --clump-snp-field SNP --clump-field P --clump-p1 5e-8 --clump-p2 1e-3 --clump-r2 0.1 --clump-kb 1000 --out %s',
    plink_path, bfile, file_sumstats, tmp.out.pref
  )
  system(cmd)

  clumped_file <- paste0(tmp.out.pref, ".clumped")

  if (!file.exists(clumped_file)) {
    warning("PLINK clumping did not produce output. Check log file: ", tmp.out.pref, ".log")

    if (clean_up) {
      unlink(list.files(tmp.dir, full.names = TRUE))
      unlink(tmp.dir, recursive = TRUE)
    }

    return(NULL)
  }

  #---------------------------------------------------------------
  # Step 2: Distance-based pruning (1 Mb window)
  #---------------------------------------------------------------
  clumped <- read.table(clumped_file, header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  clumped$P  <- as.numeric(clumped$P)
  clumped$BP <- as.numeric(clumped$BP)
  clumped <- clumped[order(clumped$CHR, clumped$P), ]

  keep <- c()
  for (chr in unique(clumped$CHR)) {
    sub <- clumped[clumped$CHR == chr, ]
    selected <- c()
    while (nrow(sub) > 0) {
      top <- sub[1, ]
      selected <- rbind(selected, top)
      sub <- sub[abs(sub$BP - top$BP) > 1e6, ]
    }
    keep <- rbind(keep, selected)
  }
  final_snps <- keep

  #---------------------------------------------------------------
  # Step 3: Compute LD matrix using PLINK
  #---------------------------------------------------------------
  snp_list_file <- file.path(tmp.dir, "snp_list_distance_pruned.txt")
  write.table(final_snps$SNP, file = snp_list_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  ld_tmp.out.pref <- file.path(tmp.dir, "snp_ld_distance_pruned.txt")
  cmd <- sprintf(
    '%s --bfile %s --extract %s --r2 square yes-really --out %s',
    plink_path, bfile, snp_list_file, ld_tmp.out.pref
  )
  system(cmd)

  ld_file <- paste0(ld_tmp.out.pref, ".ld")
  ld <- as.matrix(read.table(ld_file, header = FALSE))
  rownames(ld) <- colnames(ld) <- final_snps$SNP

  #---------------------------------------------------------------
  # Step 4: LD-based pruning (R² threshold)
  #---------------------------------------------------------------
  pvals <- setNames(final_snps$P, final_snps$SNP)
  keep <- character(0)
  snps <- names(sort(pvals))  # ascending order by P-value

  while (length(snps) > 0) {
    s <- snps[1]
    keep <- c(keep, s)
    correlated <- names(which(ld[s, ] > 0.1))
    snps <- setdiff(snps, correlated)
  }

  final_pruned <- final_snps[final_snps$SNP %in% keep, ]

  #---------------------------------------------------------------
  # Step 5: Return pruned GWAS summary statistics
  #---------------------------------------------------------------
  gwas <- data.table::fread(file_sumstats)
  gwas <- gwas[SNP %in% final_pruned$SNP]
  gwas <- gwas[order(as.numeric(P))]

  #---------------------------------------------------------------
  # Cleanup
  #---------------------------------------------------------------
  if (clean_up) {
    unlink(list.files(tmp.dir, full.names = TRUE))
    unlink(tmp.dir, recursive = TRUE)
  }

  return(gwas)
}
