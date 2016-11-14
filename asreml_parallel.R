## Run arseml in parallel over n number of jobs for a given phenotype and region.

#### NB!!! ####################################
# Set working dir to analysis folder.         #
# Script will create subdirectories           #
###############################################

## phenofile format.
# animal    nobs    dyd_*           dyd_*
# 1810      3865    0.4201330846    593.5648206
# 1893      8320    0.1953620621    593.9099203

## genotype file format (no header)
# animalID_1   0     1     2
# animalID_2   0.2   1.1   1.8

# Set parameters ----------------------------------------------------------
n_jobs = 4
phenotype = "C10"
region <- "Chr13:65000000-65500000"

phenofile <- "~tikn/Projects/Fatty_acids_bovine/GWAS/asreml/Data/AM_dyd_20_obs.txt"
pedigree <- "/mnt/users/tikn/Projects/R-packages/asremlParallel/data/testdata/pedigree/fa_20_daughters_Pedigree_asreml.txt.SRT"
# No change required below ------------------------------------------------


# setup packages, data and functions ---------------------------------------------
library(asremlParallel)
RLinuxModules::moduleInit()
RLinuxModules::module("load slurm")
# suppressPackageStartupMessages(require(plyr))
# suppressPackageStartupMessages(require(data.table))
# suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(dplyr))
library(data.table)

# Divide region string
chrom_r <- unlist(strsplit(region, ":|-"))[1]
start_r <- as.integer(unlist(strsplit(region, ":|-"))[2])
end_r <- as.integer(unlist(strsplit(region, ":|-"))[3])

# Set directories
old_workdir <- getwd()
analysis_dir <-  paste(phenotype, chrom_r, start_r, end_r, sep = "_")
dir.create(paste(phenotype, chrom_r, start_r, end_r, sep = "_"))
setwd(sprintf("%s/%s/", old_workdir, analysis_dir))

# # READ PHENOTYPIC INFO
pheno <- data.table::fread(phenofile)
names(pheno)[1]<-"animal" ## Suboptimal

# Read mapfile for genotypes ----------------------------------------------
mapfile <-
  sprintf(
    "~tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/%s_map_info.txt",
    chrom_r
  )
variant_map <-
  data.table::fread(mapfile, verbose = FALSE, data.table = FALSE)
variant_map <-
  mutate(variant_map, variant_id = paste(V1, V2,  V3, V4, sep = "_"))

#variant_map <- dplyr::mutate(variant_map, variant_id = paste(V1, V2,  V3, V4, sep = "_"))

# Set columns to read from geno-file from region.

region_final <- dplyr::filter(variant_map, V2 >= start_r & V2 <= end_r & V1 == chrom_r)
region_ids <- region_final$variant_id
region_final <- c(1, which(variant_map$variant_id %in% region_final$variant_id) + 1)

# READ GENOTYPIC INFO from dosage genoype format file.
genofile <-
  sprintf("~tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/%s_finalseqimputemergedNRF_dosage.transposed.sampleID_merge.txt", chrom_r)
geno <- read_genotypes(genofile, markers = region_final)
names(geno) <- c("animal", region_ids)

# format genotypes data frame
geno <- subset_common(geno) ## Change to just use common animals geno/pheno, and report Diff like gcta

# setup jobs
runs <- job_setup(snplist = region_ids, n_jobs = n_jobs)

# Run jobs ---------------------------------------------------------------------
plyr::l_ply(1:n_jobs, .fun = split_n_run, runs, phenotype, phenofile, pedigree)

# return to old dir -------------------------------------------------------
setwd(old_workdir)

