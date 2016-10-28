#' Process asreml results files
#'
#' MAKE A NICER OUTPUT and get number of individuals per genotype class
#' @param x data_loop object.
#' @keywords asreml
#' @export
parse_results_Tim <-function(x, multicore = FALSE){
	mat <-NULL
	mat1<-NULL
	mat2<-NULL
	mat3<-NULL
	mat<-as.matrix(table(x$SNP))
	mat1<-matrix(NA,nrow(mat),2)
	mat1[,1]<-rownames(mat)
	mat1[,2]<-mat[,1]
	mat1<-as.data.frame(mat1)

	mat2<-as.data.frame(as.matrix(c(0,1,2)))
	#names(mat2)<-"V1"

	mat3<-as.matrix(merge(mat2,mat1,by="V1",all.x=T))
	mat3[is.na(mat3)]<-0
	mat3<-as.data.frame(mat3)
	mat3<-mat3[order(mat3$V1),]

	mat3<-as.data.frame(matrix(t(mat3$V2),1,3))
	names(mat3)<-c("AA","AB","BB")
	mat3$Source<-"SNP"


	ifelse (multicore == TRUE,
			  ANOVA <- readANOVAASREML_multicore(run, SNP),
			  ANOVA <- readANOVAASREML(name))

	#ANOVA[1,2] <- "snp"

	# GET SNP EFFECTS

	est <- data.frame(SNP ,read.table(sprintf("temp_%s/%s.sln",run,SNP),head = FALSE))
	colnames(est) <- c("name","Source","level","effect","se")
	terms<-as.matrix(c('SNP'))
	colnames(terms)<-'Source'
	EFFECTS<- merge(terms,est,by='Source')

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
	res <- merge(mat3, EFFECTS)
	res<-merge(res,anova,by="Source")
	res <- res[c(1:4,9,7,8,5)]
	names(res)[names(res) %in% "name"] <- "marker"
	# tim end
# 	res<-merge(mat3,anova,by="Source")
# 	res<-merge(res,effects,by="Source")
# 	res$snp<-SNP
	return(res)
}
