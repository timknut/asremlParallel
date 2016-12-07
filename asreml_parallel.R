## Run arseml in parallel over n number of jobs for a given phenotype and region.

#### NB!!! ####################################
# Set working dir to analysis folder.         #
# Script will create subdirectories           #
# and uses tempdir() or ASreml analysis       #
# CAUTION!: geno and mapfiles are FA-project  #
# spesific for the moment                     #
###############################################

## phenofile format.
# animal    nobs    dyd_*           dyd_*
# 1810      3865    0.4201330846    593.5648206
# 1893      8320    0.1953620621    593.9099203

## genotype file format (no header)
# animalID_1   0     1     2
# animalID_2   0.2   1.1   1.8

# For testing -------------------------------------------------------------
#setwd("~/Projects/R-packages/asremlParallel/tests")
setwd("~/fettsyrepaper_2/fullsekvens_GWAS/")

# Set parameters ----------------------------------------------------------
n_jobs = 30
phenotype = "C18"
region <- "Chr1:120000000-150000000"
phenofile <-
  "/mnt/users/tikn/Projects/Fatty_acids_bovine/GWAS/asreml/Data/AM_dyd_20_obs.txt"
pedigree <-
  "/mnt/users/tikn/Projects/R-packages/asremlParallel/data/testdata/pedigree/fa_20_daughters_Pedigree_asreml.txt.SRT"
# No change required below ------------------------------------------------

# setup packages, data and functions ---------------------------------------------
library(asremlParallel)
RLinuxModules::moduleInit()
RLinuxModules::module("load slurm")
suppressPackageStartupMessages(require(dplyr))


run_asremlParallel <-
  function(n_jobs,
           phenotype,
           region,
           phenofile,
           pedigree) {
    # # Set region variables. -------------------------------------------------
    if (length(unlist(strsplit(region, ":|-"))) != 3) {
      stop(sprintf(
        "Region string: %s should be in format: 'Chr:from_bp-to_bp'",
        region
      ))
    }
    chrom_r <- unlist(strsplit(region, ":|-"))[1]
    start_r <- as.integer(unlist(strsplit(region, ":|-"))[2])
    end_r <- as.integer(unlist(strsplit(region, ":|-"))[3])

    # Check filenames ---------------------------------------------------------
    if (!file.exists(phenofile))
      stop(sprintf("phenofile: %s does not exist", phenofile))
    if (!file.exists(pedigree))
      stop(sprintf("pedigreefile: %s does not exist", pedigree))

    # Set directories ---------------------------------------------------------
    old_workdir <- getwd()
    analysis_dir <-
      paste(phenotype, chrom_r, start_r, end_r, sep = "_")
    dirlist <- list("old" = old_workdir, "new" = analysis_dir)
    dir.create(paste(phenotype, chrom_r, start_r, end_r, sep = "_"))
    setwd(sprintf("%s/%s/", old_workdir, analysis_dir))

    # # READ PHENOTYPIC ANIMAL COLUMN
    pheno <- data.table::fread(phenofile, showProgress = FALSE)
    names(pheno)[1] <- "animal" ## Suboptimal
    if(length(grep(
      pattern = paste0(phenotype,"$"),
      x = names(pheno),
      ignore.case = TRUE
    )) != 1)
      stop("Did not find unique phenotype in phenofile"
      )
    analysed_pheno <- names(pheno)[grep(pattern = paste0(phenotype,"$"), names(pheno), ignore.case = T)]
    message(sprintf('Will analyse phenotype "%s" from phenofile.', analysed_pheno))
    ## Set genofile path
    genofile <-
      sprintf(
        "~tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/%s_finalseqimputemergedNRF_dosage.transposed.sampleID_merge.txt",
        chrom_r
      )

    # Read mapfile for genotypes ----------------------------------------------
    mapfile <-
      sprintf(
        "~tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/%s_map_info.txt",
        chrom_r
      )

    # Check file paths --------------------------------------------------------
    if (!file.exists(genofile))
      stop(sprintf("genofile: %s does not exist", genofile))
    if (!file.exists(mapfile))
      stop(sprintf("mapfile: %s does not exist", pedigree))

    variant_map <-
      data.table::fread(
        mapfile,
        data.table = FALSE,
        showProgress = FALSE
      )
    variant_map <-
      mutate(variant_map, variant_id = paste(V1, V2,  V3, V4, sep = "_"))


      genolist <- asremlParallel::read_region(
        variant_map = variant_map,
        start_r = start_r,
        end_r = end_r,
        chrom_r = chrom_r,
        genofile = genofile
      )


    # format genotypes data frame
    geno <-
      subset_common(genolist$genos, pheno = pheno) ## Change to just use common animals geno/pheno, and report Diff like gcta

    # setup jobs
    runs <- job_setup(snplist = genolist$markers, n_jobs = n_jobs)

    # Run jobs ---------------------------------------------------------------------
    plyr::l_ply(1:n_jobs,
                .fun = split_n_run,
                runs,
                phenotype,
                phenofile,
                pedigree,
                geno,
                dirlist)

    setwd(old_workdir)     # return to old dir
  }

## Run main function
run_asremlParallel(n_jobs, phenotype, region, phenofile, pedigree)
