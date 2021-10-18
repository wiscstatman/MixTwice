### workspace setup

library(data.table)
library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(BiocParallel)
## load helper functions
for (f in list.files("../R", "\\.(r|R)$", full.names = TRUE)) {
  source(f)
}
# data and results directories
datdir <- "data"
resdir <- "results"
sbdir <- "../../results/GWAS"
dir.create(datdir, showWarnings = FALSE)
dir.create(resdir, showWarnings = FALSE)
dir.create(sbdir, showWarnings = FALSE)
# results files
resfile_N <- file.path(sbdir, paste0("bmi-samplesize-benchmark.rds"))
resfile_AF <- file.path(sbdir, paste0("bmi-maf-benchmark.rds"))
resfile_uninf <- file.path(sbdir, paste0("bmi-uninf-benchmark.rds"))
# set up parallel backend
cores <- 8
multicoreParam <- MulticoreParam(workers = cores)


if (!file.exists(file.path(datdir, "BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt"))) {
  download.file(url = "http://portals.broadinstitute.org/collaboration/giant/images/3/3a/BMI.SNPadjSMK.zip", 
                destfile = file.path(datdir, "BMI.SNPadjSMK.zip")) 
  unzip(file.path(datdir, "BMI.SNPadjSMK.zip"), exdir = datdir)
  file.remove(file.path(datdir,"BMI.SNPadjSMK.zip"))
  
  dfiles <- list.files(path = datdir, pattern = "BMI.SNPadjSMK.*.txt", 
                       full.names = TRUE)
  dfiles <- dfiles[!grepl("BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt", dfiles)]
  file.remove(dfiles)
}

reffile <- file.path(datdir, "1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01")
if (!file.exists(paste0(reffile, ".fam"))) {
  download.file("http://neurogenetics.qimrberghofer.edu.au/iSECA/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01.zip", 
                destfile = paste0(reffile, ".zip"))
  unzip(paste0(reffile, ".zip"), exdir = datdir)
  file.remove(paste0(reffile, ".zip"))
}

bmi <- fread(file.path(datdir, "BMI.SNPadjSMK.CombinedSexes.EuropeanOnly.txt"),
             header = TRUE)

save(bmi, file = "E://github_local//MixTwice//study2//GWAS//bmi.RData")

bmi <- bmi %>% 
  dplyr::rename(pval = p_value,
                SE = stderr,
                effect_size = effect,
                ind_covar_N = N,
                ind_covar_AF = Freq_Allele1_HapMapCEU) %>%
  mutate(test_statistic = effect_size / SE) %>%
  select(pval, SE, effect_size, ind_covar_N, ind_covar_AF, test_statistic)

set.seed(39580)
bmi <- bmi %>% mutate(ind_covar_uninf = rnorm(nrow(bmi)))

ggplot(data=bmi, aes(effect_size)) +
  geom_histogram(bins=30)