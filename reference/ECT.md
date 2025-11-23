# Effect Consistency Test (ECT)

Performs an adaptive resampling-based Effect Consistency Test (ECT) to
evaluate whether SNP effect sizes on an outcome (e.g.,GWAS trait) are
consistent with those on an exposure (e.g., pathway-level factor). The
test compares the observed R² from the exposure–outcome regression
against a null distribution obtained by permutation of outcome effects.

## Usage

``` r
ECT(
  assoc_pair,
  factor_name,
  num_snp_cutoff = 3,
  snp_fdr_cutoff = 0.2,
  seeds = 1000,
  resample_stages = c(1000, 10000, 1e+06),
  p_thresholds = c(0.1, 0.001),
  col_beta_exposure = "beta.exposure",
  col_beta_outcome = "beta.outcome",
  col_pval_exposure = "pval.exposure"
)
```

## Arguments

- assoc_pair:

  Data frame containing association summary statistics for a given
  factor. Must include columns specified by `col_beta_exposure`,
  `col_beta_outcome`, and `col_pval_exposure`.

- factor_name:

  Character string indicating the name of the factor or pathway.

- num_snp_cutoff:

  Minimum number of supporting SNPs required to perform the test
  (default: `3`).

- snp_fdr_cutoff:

  FDR threshold used to select supporting SNPs based on
  `col_pval_exposure` (default: `0.2`).

- seeds:

  Base random seed for reproducibility (default: `1000`).

- resample_stages:

  Integer vector specifying the number of resampling iterations for each
  adaptive stage (default: `c(1e3, 1e4, 1e6)`).

- p_thresholds:

  Numeric vector of stopping thresholds for adaptive resampling
  (default: `c(0.1, 0.001)`).

- col_beta_exposure:

  Column name for exposure effect sizes (default: `"beta.exposure"`).

- col_beta_outcome:

  Column name for outcome effect sizes (default: `"beta.outcome"`).

- col_pval_exposure:

  Column name for exposure p-values (default: `"pval.exposure"`).

## Value

A data frame containing the following columns:

- `factor_name`: factor/pathway tested

- `p_std`: p-value from the observed linear regression

- `r_std`: R² from the observed regression

- `p_ECT`: empirical p-value from adaptive resampling

- `num_supp_SNP`: number of SNPs used in the test

If the number of supporting SNPs is below `num_snp_cutoff`, the function
returns a placeholder row indicating insufficient SNPs.

## Details

The function fits a no-intercept linear model: \$\$beta_outcome ~ 0 +
beta_exposure\$\$ and compares its R² to those obtained after permuting
outcome effect sizes. Resampling proceeds through progressively larger
sample sizes until the empirical p-value exceeds the corresponding
threshold in `p_thresholds` or all stages are completed.

## Examples

``` r
if (FALSE) { # \dontrun{
ECT(
  assoc_pair = mydata,
  factor_name = "pathway_PC1",
  snp_fdr_cutoff = 0.1,
  resample_stages = c(1e3, 1e4),
  p_thresholds = c(0.05, 0.001)
)
} # }
```
