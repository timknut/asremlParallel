# gather all results files and plot. --------------------------------------

# WARNING! [1] "Wed Oct 19 13:05:13 2016"
# Script needs to be updated.


# calc explained var ------------------------------------------------------

tot_var = 0.248585E-01 + 2.97818
gen_var = 0.248585E-01

mutate(summary_results,
		 freq_A = (AA + 0.5*AB)/(AA + AB + BB),
		 freq_B = 1 - freq_A,
		 explained_eff = 2*freq_A * freq_B *(effect^2),
		 of_phenovar_explained = explained_eff/tot_var,
		 of_genovar_explained = explained_eff/gen_var) %>% arrange(pval) %>% head

# end ---------------------------------------------------------------------

mapfile <- "/mnt/users/tikn/Projects/Fatty_acids_bovine/GWAS/new_GWAS_mars_2016/genotypes/seqimputed/vcf/final_seqimputed_merged/dosage_format/Chr13_map_info.txt"

## Make manhattan from parallel run.
headers_bim <- c("CHR", "SNP", "cM", "BP", "A1", "A2")

bim <- read.delim("~tikn/dmu/DMU-data/plink-runs/777genos.bim", header=FALSE,
									stringsAsFactors=FALSE, col.names = headers_bim) %>% select(-cM)
bim <- filter(bim, CHR == 6)
results <- list.files("runfolder", recursive = TRUE, pattern = "summary_results_.*.\\csv",
							 full.names = TRUE)
summary_results <- ldply(results, function(x) data.table::fread(x, data.table = FALSE))
names(summary_results) <- c("AA", "AB", "BB", "LogL", "SNP")

summary_results_full <- inner_join(summary_results, bim, by = "SNP") %>%
	tbl_df
# compute LRT and p-value ----------------------------
summary_results_full <- mutate(summary_results_full, LRT = 2* -(min(summary_results$LogL) - summary_results$LogL))
summary_results_full <- mutate(summary_results_full, P = pchisq(q = LRT, df = 1, lower.tail = F))

qqman::manhattan(summary_results_full, xlim = c(min(summary_results_full$BP)/1e+06,
																max(summary_results_full$BP)/1e+06), main = "ASReml random Milk")

temp <- select(summary_results_full, -P) %>% rename(P = LRT)
qqman::manhattan(temp, logp = F, suggestiveline = F, genomewideline = F,
					  xlim = c(min(summary_results_full$BP)/1e+06, max(summary_results_full$BP)/1e+06),
					  main = "ASReml random Milk", ylab = "LRT")

qplot(data = temp, BP/1e+06, P, xlab = "Chrom 6 Mbp", ylab= "LRT") + ggtitle("ASReml random snp Milk ") +
	theme_light()
# clean up files and memory. Write new results df. ----------------------
rm(list = ls())
gc()
file.remove(c("ainverse.bin", "colnames_summary_results.csv"))
write.table(summary_results,"summary_results.csv",
				row.names = FALSE, sep = ",", quote = FALSE, col.names = TRUE)


