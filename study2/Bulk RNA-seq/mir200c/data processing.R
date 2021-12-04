library(dplyr)
library(ggplot2)
library(DESeq2)
library(SummarizedBenchmark)
library(BiocParallel)
library(recount)


## load helper functions
for (f in list.files("../R", "\\.(r|R)$", full.names = TRUE)) {
  source(f)
}

## project data/results folders
datdir <- "data"
resdir <- "results"
sbdir <- "../../results/RNAseq"
dir.create(datdir, showWarnings = FALSE)
dir.create(resdir, showWarnings = FALSE)
dir.create(sbdir, showWarnings = FALSE)

## intermediary files we create below
count_file <- file.path(datdir, "rse_gene.Rdata")
result_file <- file.path(resdir, "mir200c-results.rds")
bench_file <- file.path(sbdir, "mir200c-benchmark.rds")
bench_file_uninf <- file.path(sbdir, "mir200c-uninf-benchmark.rds")

## set up parallel backend
cores <- as.numeric(Sys.getenv("SLURM_NTASKS"))
multicoreParam <- SerialParam()

if (!file.exists(count_file)) {
  download_study('SRP030475', outdir = datdir)
}
load(count_file)
dsd <- scale_counts(rse_gene)

dsd <- dsd[, grepl("WT|200c", colData(dsd)$title)]
colData(dsd)$mir200c <- factor(ifelse(grepl("WT", colData(dsd)$title), "WT", "KO"))
dsd <- as(dsd, "DESeqDataSet")
storage.mode(assays(dsd)[["counts"]]) <- "integer"

save(dsd, file = "E://github_local//MixTwice//study2//Bulk RNA-seq//mir200c//mir200c_raw data.RData")

if (file.exists(result_file)) {
  res <- readRDS(result_file)
} else {
  design(dsd) <- ~ mir200c
  dds <- DESeq(dsd)
  res <- results(dds, independentFiltering = FALSE) %>% 
    as.data.frame() %>%
    na.omit() %>% 
    dplyr::select(pvalue, baseMean, log2FoldChange, lfcSE, stat) %>%
    dplyr::rename(pval = pvalue,
                  ind_covariate = baseMean, 
                  effect_size = log2FoldChange,
                  SE = lfcSE, 
                  test_statistic = stat)
  saveRDS(res, file = result_file)
}

res$rand_covar <- rnorm(nrow(res))

save(res, file = "E://github_local//MixTwice//study2//Bulk RNA-seq//mir200c//mir200c_processed data.RData")
