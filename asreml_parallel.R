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
genofile <- "~tikn/Projects/R-packages/asremlParallel/data/testdata/fa_subset_final_march2016_updateids_updateparents_me_removed_tenperc_chr1_bp1_bp200000.raw" ## plink --recode A format
pedigree <- "/mnt/users/tikn/Projects/R-packages/asremlParallel/data/testdata/pedigree/fa_20_daughters_Pedigree_asreml.txt.SRT"
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

# READ GENOTYPIC INFO from plink --recode raw A file
geno <- data.table::fread(genofile, data.table = F, verbose = FALSE)
names(geno) <- stringr::str_replace(names(geno), "_.$", "") ## remove _C or _2 etc from snp name.

# format genotypes data frame
geno <- subset_common(geno) ## Change to just use common animals geno/pheno, and report Diff like gcta

# Get SNP names
snp_names <- names(geno)

# setup runs
runs <- job_setup(snplist = snp_names, n_jobs = n_jobs)

# Run jobs ---------------------------------------------------------------------
plyr::l_ply(1:n_jobs, .fun = split_n_run, runs, jobname, phenofile, pedigree)
