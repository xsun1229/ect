#' Effect Consistency Test (ECT)
#'
#' @description
#' Performs an adaptive resampling-based Effect Consistency Test (ECT) to evaluate
#' whether SNP effect sizes on an outcome (e.g.,GWAS trait) are consistent with those
#' on an exposure (e.g., pathway-level factor). The test compares the observed R² from
#' the exposure–outcome regression against a null distribution obtained by
#' permutation of outcome effects.
#'
#' @param assoc_pair Data frame containing association summary statistics for a
#'   given factor. Must include columns specified by `col_beta_exposure`,
#'   `col_beta_outcome`, and `col_pval_exposure`.
#' @param factor_name Character string indicating the name of the factor or pathway.
#' @param num_snp_cutoff Minimum number of supporting SNPs required to perform the
#'   test (default: `3`).
#' @param snp_fdr_cutoff FDR threshold used to select supporting SNPs based on
#'   `col_pval_exposure` (default: `0.2`).
#' @param seeds Base random seed for reproducibility (default: `1000`).
#' @param resample_stages Integer vector specifying the number of resampling
#'   iterations for each adaptive stage (default: `c(1e3, 1e4, 1e6)`).
#' @param p_thresholds Numeric vector of stopping thresholds for adaptive
#'   resampling (default: `c(0.1, 0.001)`).
#' @param col_beta_exposure Column name for exposure effect sizes (default: `"beta.exposure"`).
#' @param col_beta_outcome Column name for outcome effect sizes (default: `"beta.outcome"`).
#' @param col_pval_exposure Column name for exposure p-values (default: `"pval.exposure"`).
#'
#' @return
#' A data frame containing the following columns:
#' - `factor_name`: factor/pathway tested
#' - `p_std`: p-value from the observed linear regression
#' - `r_std`: R² from the observed regression
#' - `p_ECT`: empirical p-value from adaptive resampling
#' - `num_supp_SNP`: number of SNPs used in the test
#'
#' If the number of supporting SNPs is below `num_snp_cutoff`, the function
#' returns a placeholder row indicating insufficient SNPs.
#'
#' @details
#' The function fits a no-intercept linear model:
#' \deqn{beta_outcome ~ 0 + beta_exposure}
#' and compares its R² to those obtained after permuting outcome effect sizes.
#' Resampling proceeds through progressively larger sample sizes until the
#' empirical p-value exceeds the corresponding threshold in `p_thresholds` or
#' all stages are completed.
#'
#' @examples
#' \dontrun{
#' ECT(
#'   assoc_pair = mydata,
#'   factor_name = "pathway_PC1",
#'   snp_fdr_cutoff = 0.1,
#'   resample_stages = c(1e3, 1e4),
#'   p_thresholds = c(0.05, 0.001)
#' )
#' }
#'
#' @importFrom parallel mclapply detectCores
#' @importFrom gtools combinations
#' @export
ECT <- function(assoc_pair,
                factor_name,
                num_snp_cutoff = 3,
                snp_fdr_cutoff = 0.2,
                seeds = 1000,
                resample_stages = c(1e3, 1e4, 1e6),
                p_thresholds = c(0.1, 0.001),
                col_beta_exposure = "beta.exposure",
                col_beta_outcome  = "beta.outcome",
                col_pval_exposure = "pval.exposure") {

  # Check for required columns
  required_cols <- c(col_beta_exposure, col_beta_outcome, col_pval_exposure)
  if (!all(required_cols %in% colnames(assoc_pair))) {
    stop("Some required columns are missing from assoc_pair: ",
         paste(setdiff(required_cols, colnames(assoc_pair)), collapse = ", "))
  }

  # Select SNPs with significant exposure association
  assoc_pair$fdr.exposure <- p.adjust(assoc_pair[[col_pval_exposure]], method = "fdr")
  snp_select <- assoc_pair[assoc_pair$fdr.exposure < snp_fdr_cutoff, ]
  num_resamp <- nrow(snp_select)

  # Return early if too few SNPs
  if (num_resamp < num_snp_cutoff) {
    return(data.frame(
      factor_name = factor_name,
      p_std = "Not enough supporting SNPs",
      r_std = NA,
      p_ECT = NA,
      num_supp_SNP = num_resamp
    ))
  }

  # Fit main regression
  fit <- lm(snp_select[[col_beta_outcome]] ~ 0 + snp_select[[col_beta_exposure]], data = snp_select)
  p_std <- summary(fit)$coefficients[1, 4]
  r_std <- summary(fit)$r.squared

  # Adaptive resampling
  p_rsp_pair_final <- NA

  for (i in seq_along(resample_stages)) {
    times_resam <- resample_stages[i]
    set.seed(seeds + i)

    r_rsp_pair <- unlist(parallel::mclapply(
      1:times_resam, function(s) {
        snp_select_s <- snp_select
        snp_select_s[[col_beta_outcome]] <- sample(
          assoc_pair[[col_beta_outcome]],
          size = num_resamp,
          replace = FALSE
        )
        summary(lm(snp_select_s[[col_beta_outcome]] ~ 0 + snp_select_s[[col_beta_exposure]]))$r.squared
      },
      mc.cores = max(1, parallel::detectCores() - 1)
    ))

    p_rsp_pair <- round((sum(r_rsp_pair > r_std) + 1) / (length(r_rsp_pair) + 1), 7)

    # Stop early if threshold satisfied
    if (i == length(resample_stages) || p_rsp_pair >= p_thresholds[i]) {
      p_rsp_pair_final <- p_rsp_pair
      break
    }
  }

  data.frame(
    factor_name = factor_name,
    p_std = p_std,
    r_std = r_std,
    p_ECT = p_rsp_pair_final,
    num_supp_SNP = num_resamp
  )
}
