# setup packages, data and functions ---------------------------------------------
lib_loc <- "/mnt/users/tikn/R/x86_64-pc-linux-gnu-library/3.3"
suppressPackageStartupMessages(require(data.table, lib.loc = lib_loc))
suppressPackageStartupMessages(require(dplyr, lib.loc = lib_loc))
suppressPackageStartupMessages(require(stringr, lib.loc = lib_loc))
library(asremlParallel)
Sys.unsetenv("DISPLAY")  # To prevent error when running interacively or via shell.
work_dir <- getwd()
cat("Working dir is: ", work_dir, "\n")
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


cat(run, "\n")
genofile <- sprintf("%s_genos.txt", run)
geno <- data.table::fread(genofile, verbose = FALSE, data.table = FALSE)

# READ PHENOTYPIC INFO
pheno <- read.table(phenofile, header=TRUE)
names(pheno)[1]<-"animal"

## select singel phenotype keeping only common geno/pheno animals
pheno <- semi_join(pheno, select(geno, 1), by = "animal" ) %>%
	select(1,ends_with(phenotype), nobs)

if(ncol(pheno) > 3){
  stop("more than 3 columns in phenofile")
}

# Get SNP names
names_snp=names(geno)

# Singel SNP loop ---------------------------------------------------------
# Make tempfolder at /work
#temp_folder <- tempdir()
user <- system("whoami", intern = T)
rundir <- basename(tempfile(pattern = "rundir"))
temp_folder <- paste0("/work/users/", user, "/", rundir)
system(paste0("mkdir -p ", temp_folder))
stopifnot(file.exists(temp_folder))
cat(sprintf("Using tempfolder: %s", temp_folder), "\n")
#setwd(sprintf("runfolder/%s", run))
#setwd(run_dir)
#temp_folder <- sprintf("temp_%s", run)
#if (!file.exists(temp_folder)) dir.create(temp_folder) # create runfolder


for (i in 2:ncol(geno)) {
  if (i == 2)
    timer_start <- proc.time() # Simple time left estimator.
  SNP <- names_snp[i] # snp[i]
  marker_matrix <- select(geno, animal, SNP = i) # c(1,i)])
  data_loop <- dplyr::inner_join(pheno, marker_matrix, by = "animal")
  data_loop_path <- sprintf('%s/data_loop.dat', temp_folder)
  readr::write_delim(
    data_loop,
    path = data_loop_path,
    col_names = TRUE, delim = " ", na = '-9'
  )
  pedline <- pedigree
  dataline <-  data_loop_path
  stopifnot(file.exists(dataline))
  # .as file need to be adjusted for your own model
  # and parameters included in the as file
  as.file <- paste0(temp_folder,"/",SNP, ".as")
  # fixed: sprintf("%s !SKIP1 !MAXIT 20 !MVINCLUDE !AISING !FCON !DDF", dataline),
  # random: !AISING !MAXIT 20 !EXTRA 10 !DDF
  cat(
    "!NODISP !WORKSPACE 2000",
    "GWAS",
    "animal !P",
    sprintf("dyd_%s", phenotype),
    #sprintf("n_%s", phenotype), # the weightings column name.
    sprintf("nobs"),
    # weigtings column name
    "SNP !D-9",
    # SHOULD BE CODED AS 0, 1, and 2. Missing (-9) will be deleted
    sprintf("%s !ALPHA !MAKE", pedline),
    sprintf("%s !SKIP1 !AISING !MAXIT 20 !EXTRA 5 !FCON !DDF", dataline),
    sprintf("dyd_%s !WT nobs ~ mu SNP !r animal", phenotype)
    ,
    file = as.file,
    sep = "\n"
  )
  stopifnot(file.exists(as.file))

  # Uncomment to keep log. --------------------------------------------------
  log_asreml <-
    system(
      paste(
        "/local/genome/packages/asreml/3.0.22.2-vb/bin/asreml ",
        as.file
      ),
      intern = T
    )
  # system(paste("/local/genome/packages/asreml/3.0.22.2-vb/bin/asreml ", as.file))
   cat(log_asreml, sep = "\n")
  # End ---------------------------------------------------------------------

  # results  <- parse_results(data_loop, multicore = TRUE)
  # asr_file <- sprintf("%s/%s", temp_folder, SNP) # use when snp as random regression.
  results <- parse_results_Tim(x = data_loop, multicore = TRUE,
                               tempfolder = temp_folder)

  readr::write_csv(
    results,
    path = sprintf('summary_results_%s.csv', run),
    append = TRUE, col_names = FALSE
  )
  # cat(log_asreml, file = sprintf('%s.log', SNP), sep = '\n')  # uncomment to write separate logs for each variant.
  system(sprintf("rm -rf %s/%s.*", temp_folder, SNP))
  if (i == 2) {
    timer_end <- proc.time() - timer_start
    cat(paste0((timer_end[3] * length(names_snp))/60, " estimated minutes  left."), sep = "\n")
  }
}
cat(names(results), sep = "\t", file = "summary_results_headers.txt")
# If you restart some anlysis, make sure that you delete the previous
# summary_results file because it will keep writing the results in the same file
