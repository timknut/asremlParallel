## Run arseml in parallel for n number of jobs for a given jobname.
## Consider using plink2 to subset the datafiles if slow or memory hungry. Not very portable.But efficient.

#### NB!!! ####################################
# Set wd to analysis folder with scripts      #
###############################################
setwd("~/Projects/R-packages/asremlParallel/tests/")

# Set parameters ----------------------------------------------------------
n_jobs = 2
jobname = "milk"
phenofile <- "~tikn/Projects/Mastitis/GWAS/Data/fenos_melkeegsk_GC.txt"
genofile <- "~tikn/Projects/Mastitis/GWAS/Data/genos_777_dmu_chr_6.raw" ## plink --recode A format
pedigree <- "/mnt/users/tikn/Projects/Mastitis/GWAS/Data/hgo_ped3_sort.ped"
delim = "\t"
# No change required below ------------------------------------------------

# setup packages, data and functions ---------------------------------------------

dir.create(jobname)
dirs <- list.dirs(recursive = F)
dircheck <- unlist(strsplit(dirs, "/"))
dircheck <- dircheck[match(jobname, dircheck)]


if(!grepl(jobname, dircheck)){
	stop("Did you set working directory?") } else {
		# Load different funcitons.
		source("~tikn/Projects/R-packages/asremlParallel/asreml_utils/asreml_utils_2.R")
		suppressPackageStartupMessages(require(plyr))
		suppressPackageStartupMessages(require(data.table))
		suppressPackageStartupMessages(require(stringr))
		suppressPackageStartupMessages(require(dplyr))
	}

# setwd(jobname) # Should not setwd!

# # READ PHENOTYPIC INFO
pheno <- readr::read_delim(phenofile, delim = delim, col_names = TRUE)
names(pheno)[1]<-"animal" ## Suboptimal

# READ GENOTYPIC INFO from plink --recode raw A file
geno <- data.table::fread(genofile, data.table = F, verbose = FALSE)
names(geno) <- str_replace(names(geno), "_.$", "") ## remove _C or _2 etc from snp name.

# format genotypes data frame
geno <- subset_common(geno) ## Change to just use common animals geno/pheno, and report Diff like gcta

# Get SNP names
snp_names <- names(geno)

# setup runs
runs <- job_setup(snplist = snp_names, n_jobs = n_jobs)

# Run jobs ---------------------------------------------------------------------
plyr::l_ply(1:n_jobs, .fun = split_n_run, runs, jobname, phenofile, pedigree)
