#' Process asreml results files
#'
#' MAKE A NICER OUTPUT and get number of individuals per genotype class
#' @param x data.frame data_loop object.
#' @param multicore Locical.
#' @param tempfolder string. Path to temp-folder
#' @keywords asreml
#' @export
parse_results_Tim <-function(x, multicore = FALSE, tempfolder){
  df <- x
  df$SNP <- round(df$SNP) # To count the dosage GTs
  df <- dplyr::count(df, SNP)
  snp_counts <- data.frame(SNP = c(0,1,2))
  snp_counts <- merge(df, snp_counts, all = T)
  class(snp_counts$SNP) <- "character"
  if(any(is.na(snp_counts$SNP))) {
    snp_counts$SNP[is.na(snp_counts$SNP)] <- "missing"
  } else {
    new_df <- dplyr::data_frame(SNP = "missing",
                                count = 0)
    snp_counts <-  dplyr::bind_rows(snp_counts, new_df)
  }
  snp_counts <- as.data.frame(t(snp_counts), stringsAsFactors = F)[2,]
	names(snp_counts)<-c("AA","AB","BB", "missing")
	snp_counts$Source<-"SNP"
	snp_counts[is.na(snp_counts)] <- 0

	ifelse (
	  multicore == TRUE,
	  ANOVA <-
	    readANOVAASREML_multicore(run, SNP, tempfolder = tempfolder),
	  ANOVA <- readANOVAASREML(name)
	)

	# GET SNP EFFECTS
	est <-  data.frame(SNP , read.table(sprintf("%s/%s.sln", tempfolder, SNP),
	                                    header = FALSE))
	colnames(est) <- c("name", "Source", "level", "effect", "se")
	terms <- as.matrix(c('SNP'))
	colnames(terms) <- 'Source'
	EFFECTS <- merge(terms, est, by = 'Source')

	# GET P VALUE
	anova<-NULL
	effects<-NULL
	anova<-subset(ANOVA,Source=="SNP")
	anova$pval<- pf(anova$F_con,anova$NumDF,anova$DenDF_con, lower.tail=F)
	anova<-subset(anova,select=c(Source,pval))
	effects<-as.data.frame(EFFECTS$effect)
	colnames(effects)<-c("effect")
	effects$Source<-"SNP"

	res<-NULL
	# Tim
	res <- merge(snp_counts, EFFECTS)
	res<-merge(res,anova,by="Source")
	res <- res[c("name", "AA", "AB", "BB","missing", "pval", "effect", "se")]
	names(res)[names(res) %in% "name"] <- "marker"
	# tim end
	return(res)
}
