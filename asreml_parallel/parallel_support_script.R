# setup packages, data and functions ---------------------------------------------
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(stringr))
Sys.unsetenv("DISPLAY")  # To prevent error when running interacively or via shell.
work_dir <- getwd()
setwd(getwd())
#setwd("~hanngo/ASREML-melk-mastitt-HD_multicore/")
# functions for processing the genotypes and ANOVA (Don't require changes)
source("~tikn/Projects/R-packages/asreml_utils/asreml_utils_2.R")
source("~tikn/Projects/R-packages/asreml_utils/parse_results_Tim.R")

# Parse command line args.
args <- commandArgs(TRUE)
run <- args[1]
phenotype <- args[2]
phenofile <- args[3

cat(run, "\n")
genofile <- sprintf("runfolder/%s/%s_genos.txt", run, run)
geno <- data.table::fread(genofile)

# READ PHENOTYPIC INFO
pheno <- read.table(phenofile, header=TRUE)
names(pheno)[1]<-"animal"

## select singel phenotype keeping only common geno/pheno animals
pheno <- semi_join(pheno, select(geno, 1), by = c("animal" = "IID")) %>%
	select(1,ends_with(phenotype))

# Get SNP names
names_snp=names(geno)

# Singel SNP loop ---------------------------------------------------------
setwd(sprintf("runfolder/%s", run))
temp_folder <- sprintf("temp_%s", run)
if (!file.exists(temp_folder)) dir.create(temp_folder) # create runfolder

for (i in 2:ncol(geno)) {
	SNP=names_snp[i] # snp[i]
	marker <- select(geno, animal = IID, snp = i) # c(1,i)])
	data_loop <- merge(pheno, marker, by='animal')
	write.table(data_loop, sprintf('%s/data_loop.dat',temp_folder),
					col.names=TRUE, row.names=FALSE, quote=FALSE)
	pedline <- "/mnt/users/tikn/Projects/Mastitis/GWAS/Data/hgo_ped3_sort.ped"
	dataline <- sprintf("%s/runfolder/%s/%s/data_loop.dat", work_dir, run, temp_folder)
	# .as file
	# need to be adjusted for your own model and parameters included in the as file
	as.file <- paste(SNP,".as",sep="")
	as.file <- sprintf("%s/%s", temp_folder, as.file)

	# fixed: sprintf("%s !SKIP1 !MAXIT 20 !MVINCLUDE !AISING !FCON !DDF", dataline),
	# radom: !AISING !MAXIT 20 !EXTRA 10 !DDF
	cat("!NODISP !WORKSPACE 2000",
		 "GWAS",
		 "animal !P",
		 sprintf("dyd_%s", phenotype),
		 sprintf("n_%s", phenotype),
		 "snp",  # SHOULD BE CODED AS 0, 1, and 2
		 sprintf("%s !ALPHA !MAKE", pedline),
		 sprintf("%s !SKIP1 !AISING !MAXIT 20 !EXTRA 5 !FCON !DDF", dataline),
		 sprintf("dyd_%s !WT n_%s ~ mu snp !r animal",phenotype, phenotype)
		 ,file = as.file, sep="\n")

	system(paste("/local/genome/packages/asreml/3.0.22.2-vb/bin/asreml ", as.file))

	# results  <- parse_results(data_loop, multicore = TRUE)
	# asr_file <- sprintf("%s/%s", temp_folder, SNP) # use when snp as random regression.
	results <- parse_results_Tim(data_loop, multicore = TRUE)
	system(sprintf("rm %s/%s.*", temp_folder, SNP))

	write.table(results, sprintf('summary_results_%s.csv', run), row.names=FALSE,
					sep=",", quote=FALSE, append=TRUE, col.names=FALSE)
}

# If you restart some anlysis, make sure that you delete the previous
# summary_results file because it will keep writing the results in the same file
