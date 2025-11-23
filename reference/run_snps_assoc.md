# Run PLINK-based SNP association analyses in parallel

Performs SNP-level association analyses for a given list of SNPs by
invoking PLINK in parallel. Each SNP is tested using its own covariate
file, and association results can optionally be grouped by phenotype
factors.

## Usage

``` r
run_snps_assoc(
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
  group_by_factor = TRUE
)
```

## Arguments

- snp_list:

  Character vector of SNP IDs to analyze.

- folder_covar:

  Path to the folder containing SNP-specific covariate files, in PLINK
  format (e.g., `"snp123.cov.txt"`).

- folder_assoc:

  Path to the output folder for PLINK association results.

- plink_path:

  Path to the PLINK executable.

- file_bfile:

  PLINK binary file prefix (`.bed/.bim/.fam`).

- file_phe:

  Path to the phenotype file to be analyzed by PLINK. This file should
  contain the phenotype columns (i.e., the factors) that are later
  referenced by `group_by_factor` to group and summarize association
  results.

- single_factor:

  Optional; a single phenotype name (column name from `file_phe`) to
  test individually. If `NULL`, all phenotypes in `file_phe` are
  analyzed (`--all-pheno` mode).

- add_MAF:

  Logical; if `TRUE`, the function merges association results with a MAF
  table. When `add_MAF = TRUE`, the argument `df_maf` **must** be
  provided. If `FALSE`, MAF information is ignored and `df_maf` can be
  left as `NULL`.

- df_maf:

  Data frame containing columns `SNP`, `A2`, and `MAF`. Required only
  when `add_MAF = TRUE`; otherwise leave it as `NULL`.

- ncores:

  Number of CPU cores to use for parallel processing.

- cleanup:

  Logical; if `TRUE`, removes intermediate `.assoc` files after
  grouping.

- group_by_factor:

  Logical; if `TRUE`, the function combines all SNP-level results and
  groups them by phenotype factors defined in `file_phe`.

## Value

Writes PLINK association results into `folder_assoc/` and (optionally)
grouped per-factor summary files (e.g., `"factor_pvalues.txt"`). Returns
no explicit object.

## Details

Each SNP is processed in its own subfolder under `folder_assoc/` to
avoid overwriting during parallel execution. Only additive test results
(`TEST == "ADD"`) are retained.

## Examples

``` r
if (FALSE) { # \dontrun{
run_snps_assoc(
  snp_list = c("rs123", "rs456"),
  folder_covar = "covariates/",
  folder_assoc = "assoc_results/",
  plink_path = "/usr/bin/plink",
  file_bfile = "data/genotypes",
  file_phe = "data/pheno.txt",
  add_MAF = TRUE,
  df_maf = maf_table,
  ncores = 4
)
} # }
```
