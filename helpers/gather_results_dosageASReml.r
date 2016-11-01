# gather all results files and make IGV GWAS-file. --------------------------------------

require(tibble)
require(readr)

###### HUSK ######
# Oppdater chrom til aktuelt kromosom.
# Sett setwd() til dir med aktuell egenskap.

chrom <- 13
headers_mapfile <- c("CHR", "BP", "A1", "A2", "AR2_imputation")


# Read mapfile
mapfile <- sprintf(  "/mnt/users/tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/Chr%i_map_info.txt", chrom)

read_mapfile <-  function(x) {
  stopifnot(file.exists(mapfile)) ## This will produce error if map file does not exist.
  map_info <- read_delim(x, "\t", escape_double = FALSE,
                         col_names = headers_mapfile)
  dplyr::mutate(map_info, SNP = paste(CHR, BP, A1, A2, sep = "_"))
}

map_info <- read_mapfile(mapfile)

results <-
  list.files(
    "runfolder",
    recursive = TRUE,
    pattern = "summary_results_.*.\\csv",
    full.names = TRUE
  )

# SNP,1512,46,0,0.777329752766979,-0.005467,0.01977,V2
# SNP,120,72,0,0.583954417361325,-0.007406,0.01341,V3

summary_results <-
  plyr::ldply(results, function(x)
    data.table::fread(x, data.table = FALSE))
names(summary_results) <- c("type", "AA", "AB", "BB", "p", "effect", "se", "SNP")

# Joining in map info. Check all the missing genotypes.

summary_results_full <-
  inner_join(summary_results, map_info, by = "SNP") %>%
  tbl_df %>%
  select(-1)
rm(summary_results)

# write GWAS file for IGV viewing ---------------------------------------------------------
write_tsv(
  select(summary_results_full, p, SNP, CHR, BP) %>% arrange(BP),
  path = sprintf("res_c08_chr%i.GWAS", chrom),
  col_names = T
)







# OLD Manhattan plot code----------------------------

qqman::manhattan(summary_results_full,
                 xlim = c(
                   min(summary_results_full$BP) / 1e+06,
                   max(summary_results_full$BP) / 1e+06
                 ),
                 main = "ASReml random Milk")

temp <- select(summary_results_full,-P) %>% rename(P = LRT)
qqman::manhattan(
  temp,
  logp = F,
  suggestiveline = F,
  genomewideline = F,
  xlim = c(
    min(summary_results_full$BP) / 1e+06,
    max(summary_results_full$BP) / 1e+06
  ),
  main = "ASReml random Milk",
  ylab = "LRT"
)

qplot(data = temp,
      BP / 1e+06,
      P,
      xlab = "Chrom 6 Mbp",
      ylab = "LRT") + ggtitle("ASReml random snp Milk ") +
  theme_light()
# clean up files and memory. Write new results df. ----------------------
rm(list = ls())
gc()
file.remove(c("ainverse.bin", "colnames_summary_results.csv"))
write.table(
  summary_results,
  "summary_results.csv",
  row.names = FALSE,
  sep = ",",
  quote = FALSE,
  col.names = TRUE
)
