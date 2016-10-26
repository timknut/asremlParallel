## Run arseml in parallel for n number of jobs for a given jobname.

#### NB!!! ####################################
# Set working dir to analysis folder          #
###############################################

setwd("~/Projects/R-packages/asremlParallel/tests/")

## phenofile format.
# animal  pheno1       pheno2        n_obs_1   n_obs_2
# 1810    0.4201330846    593.5648206     3865    3865
# 1893    0.1953620621    593.9099203     8320    8319

# Set parameters ----------------------------------------------------------
n_jobs = 4
jobname = "C14"
phenofile <- "~tikn/Projects/R-packages/asremlParallel/data/testdata/AM_FA_dyd_asreml_20_daughters.txt"
#genofile <- "~tikn/Projects/R-packages/asremlParallel/data/testdata/fa_subset_final_march2016_updateids_updateparents_me_removed_tenperc_chr1_bp1_bp200000.raw" ## plink --recode A format
genofile <- "~tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/Chr13_finalseqimputemergedNRF_dosage.transposed.sampleID_merge.txt" ## plink --recode A format
pedigree <- "/mnt/users/tikn/Projects/R-packages/asremlParallel/data/testdata/pedigree/fa_20_daughters_Pedigree_asreml.txt.SRT"
mapfile <- "~tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/Chr13_map_info.txt"
# delim = " "

# No change required below ------------------------------------------------

# setup packages, data and functions ---------------------------------------------
suppressPackageStartupMessages(require(plyr))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(dplyr))
library(asremlParallel)
library(RLinuxModules) # move these three lines main script.
moduleInit()
module("load slurm")

dir.create(jobname)
dircheck(jobname)

# # READ PHENOTYPIC INFO
pheno <- data.table::fread(phenofile)
names(pheno)[1]<-"animal" ## Suboptimal


# Read mapfile for genotypes ----------------------------------------------
variant_map <- data.table::fread(mapfile, data.table = FALSE, verbose = FALSE)
# Subset region
region <- "Chr13:1-5000"
chrom_r <- unlist(stringr::str_split(region, ":|-"))[1]
start_r <- as.numeric(unlist(stringr::str_split(region, ":|-"))[2])
end_r <- as.numeric(unlist(stringr::str_split(region, ":|-"))[3])
region_final <- filter(variant_map, V2 > start_r & V2 < end_r & V1 == chrom_r)
region_final <- 1:nrow(region_final)
# End ---------------------------------------------------------------------


# READ GENOTYPIC INFO from plink --recode raw A file
#geno <- read_genotypes(genofile, markers = region_final)
system.time(geno <- read_genotypes(genofile, markers = region_final))
names(geno) <- stringr::str_replace(names(geno), "_.$", "") ## remove _C or _2 etc from snp name.

# format genotypes data frame
geno <- subset_common(geno) ## Change to just use common animals geno/pheno, and report Diff like gcta

# Get SNP names
snp_names <- names(geno)

# setup runs
runs <- job_setup(snplist = snp_names, n_jobs = n_jobs)


#  Below does not WORK!! ---------------------------------------------------------
# Run_1
# Error in eval(expr, envir, enclos) :
#   'IID' column not found in rhs, cannot join
# Calls: %>% ... semi_join -> semi_join.tbl_df -> semi_join_impl -> .Call
# Execution halted


# Run jobs ---------------------------------------------------------------------
plyr::l_ply(1:n_jobs, .fun = split_n_run, runs, jobname, phenofile, pedigree)


