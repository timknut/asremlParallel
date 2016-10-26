# setup packages, data and functions ---------------------------------------------
lib_loc <- "/mnt/users/tikn/R/x86_64-pc-linux-gnu-library/3.2"
suppressPackageStartupMessages(require(data.table, lib.loc = lib_loc))
suppressPackageStartupMessages(require(dplyr, lib.loc = lib_loc))
suppressPackageStartupMessages(require(stringr, lib.loc = lib_loc))
library(asremlParallel)
Sys.unsetenv("DISPLAY")  # To prevent error when running interacively or via shell.
work_dir <- getwd()

# functions for processing the genotypes and ANOVA (Don't require changes)
# source("~tikn/Projects/R-packages/asremlParallel/asreml_utils/asreml_utils_2.R") # not needed. Remove after tests.
# source("~tikn/Projects/R-packages/asremlParallel/asreml_utils/parse_results_Tim.R") # same.

# Parse command line args.
args <- commandArgs(TRUE)
run <- args[1]
phenotype <- args[2]
phenofile <- args[3]
pedigree <- args[4]

# Test --------------------------------------------------------------------
# run <- 1
# phenotype <- "milk"
# phenofile <- "~tikn/Projects/Mastitis/GWAS/Data/fenos_melkeegsk_GC.txt"
# pedigree <- "/mnt/users/tikn/Projects/Mastitis/GWAS/Data/hgo_ped3_sort.ped"
# end ---------------------------------------------------------------------


cat(paste0("Run_", run, "\n"))
genofile <- sprintf("runfolder/run_%s/run_%s_genos.txt", run, run)
geno <- data.table::fread(genofile, verbose = FALSE, data.table = FALSE)

# READ PHENOTYPIC INFO
pheno <- read.table(phenofile, header=TRUE)
names(pheno)[1]<-"animal"

## select singel phenotype keeping only common geno/pheno animals
pheno <- semi_join(pheno, select(geno, 1), by = c("animal" = "V1")) %>%
	select(1,ends_with(phenotype), nobs)


if(ncol(pheno) > 3){
  stop("more than 3 columns in phenofile")
}

# Get SNP names
names_snp=names(geno)

# Singel SNP loop ---------------------------------------------------------
setwd(sprintf("runfolder/run_%s", run))
temp_folder <- sprintf("temp_%s", run)
if (!file.exists(temp_folder)) dir.create(temp_folder) # create runfolder

for (i in 2:ncol(geno)) {
	SNP <- names_snp[i] # snp[i]
	marker_matrix <- select(geno, animal = V1, snp = i) # c(1,i)])
	data_loop <- dplyr::inner_join(pheno, marker_matrix, by = "animal")
	readr::write_delim(data_loop, path = sprintf('%s/data_loop.dat',temp_folder),
					col_names=TRUE, delim = " ")
	pedline <- pedigree
	dataline <- sprintf("%s/runfolder/run_%s/%s/data_loop.dat", work_dir, run, temp_folder)
	stopifnot(file.exists(dataline))
	# .as file
	# need to be adjusted for your own model and parameters included in the as file
	as.file <- paste(SNP,".as",sep="")
	as.file <- sprintf("%s/%s", temp_folder, as.file)
	# fixed: sprintf("%s !SKIP1 !MAXIT 20 !MVINCLUDE !AISING !FCON !DDF", dataline),
	# random: !AISING !MAXIT 20 !EXTRA 10 !DDF
	cat("!NODISP !WORKSPACE 2000",
		 "GWAS",
		 "animal !P",
		 sprintf("dyd_%s", phenotype),
		 #sprintf("n_%s", phenotype), # the weightings column name.
		 sprintf("nobs"),             # weigtings column name
		 "SNP !D-1",  # SHOULD BE CODED AS 0, 1, and 2. Missing (-1) will be deleted
		 sprintf("%s !ALPHA !MAKE", pedline),
		 sprintf("%s !SKIP1 !AISING !MAXIT 20 !EXTRA 5 !FCON !DDF", dataline),
		 sprintf("dyd_%s !WT nobs ~ mu SNP !r animal",phenotype)
		 ,file = as.file, sep="\n")


# Uncomment to keep log. --------------------------------------------------
		log_asreml <-
	  system(paste(
	    "/local/genome/packages/asreml/3.0.22.2-vb/bin/asreml ",
	    as.file
	  ),
	  intern = T)
#  system(paste("/local/genome/packages/asreml/3.0.22.2-vb/bin/asreml ", as.file))
# End ---------------------------------------------------------------------

	# results  <- parse_results(data_loop, multicore = TRUE)
	# asr_file <- sprintf("%s/%s", temp_folder, SNP) # use when snp as random regression.
	results <- parse_results_Tim(data_loop, multicore = TRUE)

	readr::write_csv(results, path = sprintf('summary_results_%s.csv', run), append=TRUE, col_names=FALSE)
	system(sprintf("rm %s/%s.*", temp_folder, SNP)) # Uncomment to keep temp log-files
}

# If you restart some anlysis, make sure that you delete the previous
# summary_results file because it will keep writing the results in the same file
