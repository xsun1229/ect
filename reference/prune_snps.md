# Prune SNPs based on clumping and LD thresholds using PLINK

Performs SNP pruning in two stages: (1) PLINK-based clumping to select
approximately independent lead SNPs, and (2) additional distance- and
LD-based pruning to further reduce correlated variants.

The function generates temporary files under a subdirectory of
`out.dir`, invokes PLINK commands, and returns a pruned GWAS summary
statistics table.

## Usage

``` r
prune_snps(
  file_sumstats,
  plink_path,
  bfile,
  out.dir,
  out.pref,
  clean_up = FALSE
)
```

## Arguments

- file_sumstats:

  Path to the GWAS summary statistics file containing columns `SNP`,
  `CHR`, `BP`, and `P`.

- plink_path:

  Path to the PLINK executable.

- bfile:

  PLINK binary file prefix (`.bed/.bim/.fam`) used as the LD reference.

- out.dir:

  Output directory where intermediate and result files will be written.

- out.pref:

  Prefix for temporary PLINK output and final file naming.

- clean_up:

  Logical; if `TRUE`, temporary files and subdirectories are deleted
  after successful completion or if errors occur.

## Value

A `data.table` of pruned GWAS summary statistics containing the same
columns as the input file but restricted to SNPs that remain after
distance- and LD-based pruning. Returns `NULL` if PLINK clumping fails
to produce output.

## Details

The pruning proceeds in the following steps:

1.  **PLINK clumping:** using thresholds `--clump-p1 5e-8`,
    `--clump-p2 1e-3`, `--clump-r2 0.1`, and `--clump-kb 1000`.

2.  **Distance-based pruning:** removes SNPs within ±1 Mb of a more
    significant SNP.

3.  **LD-based pruning:** calculates pairwise LD (`r²`) using PLINK and
    removes SNPs with `r² > 0.1` relative to previously retained
    variants.

4.  Returns GWAS summary statistics restricted to the remaining
    independent SNPs.

## Examples

``` r
if (FALSE) { # \dontrun{
pruned_gwas <- prune_snps(
  file_sumstats = "assoc/gwas.txt",
  plink_path = "/usr/bin/plink",
  bfile = "ref/1000G_EUR",
  out.dir = "pruning_results/",
  out.pref = "GWAS1",
  clean_up = TRUE
)
} # }
```
