## Run arseml in parallel for n number of jobs for a given trait.
## Consider using plink2 to subset the datafiles if slow or memory hungry. Not very portable.But efficient.

#### NB!!! ####################################
# Set wd to analysis folder with scripts      #
###############################################

# Set parameters ----------------------------------------------------------
n_jobs = 10
jobname = "milk"
phenofile <- "~tikn/Projects/Mastitis/GWAS/Data/fenos_melkeegsk_GC.txt"
genofile <- "~tikn/Projects/Mastitis/GWAS/Data/genos_777_dmu_chr_6.raw" ## plink --recode A format
# No change required below ------------------------------------------------

# setup packages, data and functions ---------------------------------------------

dir.create(trait)
dirs <- list.dirs(recursive = F)
dircheck <- unlist(strsplit(dirs, "/"))
dircheck <- dircheck[match(trait, dircheck)]


if(!grepl(trait, dircheck)){
	stop("Did you set working directory?") } else {
		# Load different funcitons.
		source("~tikn/Projects/R-packages/asreml_utils/asreml_utils_2.R")
		suppressPackageStartupMessages(require(plyr))
		suppressPackageStartupMessages(require(data.table))
		suppressPackageStartupMessages(require(stringr))
		suppressPackageStartupMessages(require(dplyr))
	}

setwd(trait)
# READ PEDIGREE
# ped  <- read.table("Data/pedigree/fa_20_daugters_Pedigree.txt", header=FALSE)[1:3]

# # READ PHENOTYPIC INFO
pheno <- readr::read_delim(phenofile, delim = "\t", col_names = TRUE)
names(pheno)[1]<-"animal" ## Not generic. Should not matter what col 1 is names. 

# READ GENOTYPIC INFO from plink --recode raw A file
geno <- data.table::fread(genofile, data.table = F, verbose = FALSE)
names(geno) <- str_replace(names(geno), "_.$", "") ## remove _C or _2 etc from snp name.

# format genotypes data frame
geno <- genosubset(geno, multicore = F) ## Change to just use common animals geno/pheno, and report Diff like gcta

# Get SNP names
names_snp=names(geno)

# setup runs
runs <- data_frame(marker = names_snp[2:length(names_snp)])
runs <- mutate(runs, run = paste("run", ntile(seq_along(marker), n_jobs), sep = "_"))
runs <- split(runs$marker, runs$run)
dir.create("runfolder")

# Run jobs ---------------------------------------------------------------------
l_ply(1:n_jobs, .fun = split_n_run)
