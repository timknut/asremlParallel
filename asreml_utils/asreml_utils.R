#' Process asreml ANOVA output
#'
#' function for process the ANOVA (Don't require changes)
#' @param name Dummy variable.
#' @keywords asreml
#' @export
#
readANOVAASREML <- function(name)
{
	ss <- readLines(paste(SNP,".asr",sep = ""))
	begin_ss <- grep(pattern = "Source of Variation",x = unlist(ss))
	ss_1 <- ss[seq(begin_ss + 1,length(ss))]
	end<- grep(pattern = "Notice:",x = unlist(ss_1))
	end_ss <- ss_1[seq(1, end[1]-1)]
	V_ss <- unlist(strsplit(end_ss,split = " "))
	V_ss2 <- V_ss[V_ss!= ""]
	mat_ss <- as.data.frame(matrix(V_ss2,ncol = 8,byrow = T),stringsAsFactors = FALSE)[,-c(1)]
	names(mat_ss) <- c("Source","NumDF","DenDF_con","F_inc","F_con","M", "P_con")
	mat_ss$NumDF <- as.numeric(mat_ss$NumDF)
	mat_ss$DenDF_con <- as.numeric(mat_ss$DenDF_con)
	mat_ss$F_inc <- as.numeric(mat_ss$F_inc)
	mat_ss$F_con <- as.numeric(mat_ss$F_con)
	ANOVA <- data.frame(SNP,mat_ss)
	return(ANOVA)
}

#' Process asreml results files
#'
#' MAKE A NICER OUTPUT and get number of individuals per genotype class
#' @param x data_loop object.
#' @keywords asreml
#' @export
parse_results <- function(x){
	mat <-NULL
	mat1<-NULL
	mat2<-NULL
	mat3<-NULL
	mat<-as.matrix(table(x$snp))
	mat1<-matrix(NA,nrow(mat),2)
	mat1[,1]<-rownames(mat)
	mat1[,2]<-mat[,1]
	mat1<-as.data.frame(mat1)

	mat2<-as.data.frame(as.matrix(c(0,1,2)))
	names(mat2)<-"V1"

	mat3<-as.matrix(merge(mat2,mat1,by="V1",all.x=T))
	mat3[is.na(mat3)]<-0
	mat3<-as.data.frame(mat3)
	mat3<-mat3[order(mat3$V1),]

	mat3<-as.data.frame(matrix(t(mat3$V2),1,3))
	names(mat3)<-c("AA","AB","BB")
	mat3$Source<-"snp"

	ANOVA <- readANOVAASREML(name)
	#ANOVA[1,2] <- "snp"

	# GET SNP EFFECTS

	est <- data.frame(SNP ,read.table(paste(SNP,".sln",sep =""),head = FALSE))
	colnames(est) <- c("name","Source","level","effect","se")
	terms<-as.matrix(c('snp'))
	colnames(terms)<-'Source'
	EFFECTS<- merge(terms,est,by='Source')

	# GET P VALUE
	anova<-NULL
	effects<-NULL
	anova<-subset(ANOVA,Source=="snp")
	anova$pval<- pf(anova$F_con,anova$NumDF,anova$DenDF_con, lower.tail=F)
	anova<-subset(anova,select=c(Source,pval))

	effects<-as.data.frame(EFFECTS$effect)

	colnames(effects)<-c("effect")
	effects$Source<-"snp"

	res<-NULL
	res<-merge(mat3,anova,by="Source")
	res<-merge(res,effects,by="Source")
	res$snp<-SNP
	return(res)
}

#' Format genotypes data frame.
#'
#' Format genotypes data frame. Make genosubset if fraction = TRUE
#' @param x plink raw format data frame.
#' @param fraction Locical. Will pick a random fraction of markes for testing if TRUE. Default [FALSE]
#' @keywords asreml gwas
#' @export
##
genosubset <- function(x, fraction = FALSE) {
	animals_geno <- select(x, 2)
	## find common samples
	index_common_animals <- match(intersect(pheno$animal, animals_geno$IID), animals_geno$IID)
	## make genosubset with common animals. omit all but animal column
	geno_subset <- x[index_common_animals,] %>% select(2,7:ncol(x))
	if (fraction == TRUE){
		rand_markers <- sample_frac(data_frame(markers = seq(2, ncol(geno_subset))), 0.1)
		x <- select(geno_subset, 1, as.integer(rand_markers$markers))
	} else {
		x  <- geno_subset
	}
	## remove objects
	rm(geno_subset, index_common_animals, animals_geno)
	## and clear memory..
	gc()
	return(x)
}


