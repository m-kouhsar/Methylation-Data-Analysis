
suppressMessages(library(metafor))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(qqman))
suppressMessages(library(data.table))

args<-commandArgs(TRUE)
wd <- args[1]
study1 <- args[2]
study2 <- args[3]
out_pref <- args[4]

print("Arguments: ")
print(paste0("  Working dir= ",wd))
print(paste0("  EWAS file for study 1= ",study1))
print(paste0("  EWAS file for study 2= ",study2))
print(paste0("  Output files prefix= ",out_pref))


setwd(wd)

load(study1)
EWAS_RESULTS.1 <- EWAS_RESULTS.ANNOT
rm(EWAS_RESULTS.ANNOT)
load(study2)
EWAS_RESULTS.2 <- EWAS_RESULTS.ANNOT
rm(EWAS_RESULTS.ANNOT)

index <- rownames(EWAS_RESULTS.1) %in% rownames(EWAS_RESULTS.2)
table(index)
EWAS_RESULTS.1 <- EWAS_RESULTS.1[index,]

index <- rownames(EWAS_RESULTS.2) %in% rownames(EWAS_RESULTS.1)
table(index)
EWAS_RESULTS.2 <- EWAS_RESULTS.2[index,]

index <- match(rownames(EWAS_RESULTS.1), rownames(EWAS_RESULTS.2))
EWAS_RESULTS.2 <- EWAS_RESULTS.2[index,]
identical(rownames(EWAS_RESULTS.1),rownames(EWAS_RESULTS.2))
#[1] TRUE

es<-cbind.data.frame(EWAS_RESULTS.2[,1],EWAS_RESULTS.1[,1])
names(es) <- c('es.2','es.1')
se<-cbind.data.frame(EWAS_RESULTS.2[,2],EWAS_RESULTS.1[,2])
names(se) <- c('se.2','se.1')
data<-cbind.data.frame(es,se)
rownames(data) <- rownames(EWAS_RESULTS.1)
meta <- function(x){
  res1<-metafor::rma(x[c(1,2)],sei=x[c(3,4)], method="DL")
  return(c(res1$zval, res1$pval, res1$QE, res1$QEp, res1$tau2, res1$I2,res1$k,res1$beta,res1$se))
}

cl<-parallel::makeCluster(32)
res.meta<-t(parApply(cl,data,1,meta))
colnames(res.meta)<-c("metafor.zval", "metafor.pval", "metafor.QE", "metafor.QEp", "metafor.tau2", "metafor.I2","metafor.k","metafor.beta","metafor.se")
stopCluster(cl)

results_metafor<-as.data.frame(res.meta)
chisq <- qchisq(1-(results_metafor[,2]),1)
inf_indx <- median(chisq)/qchisq(0.5,1)
print(paste0("meta-analysis inflation index= ",inf_indx))

print("table(results.pval < 1e-4):")
table(results_metafor$metafor.pval < 1e-4)

epicManifest<-fread("/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv", skip = 7,stringsAsFactors = F)
ANNOT<-data.frame(epicManifest)
index1 <- rownames(results_metafor)%in% ANNOT[,1]
index2 <- match(rownames(results_metafor)[index1],ANNOT[,1])
identical(ANNOT[index2,][,1],rownames(results_metafor)[index1])
# True
results_metafor.Annot <- cbind.data.frame(results_metafor[index1,],ANNOT[index2,])

results_metafor.Annot$MAPINFO.1<-(results_metafor.Annot$MAPINFO)+1
results_metafor.Annot <- results_metafor.Annot[order(results_metafor.Annot$metafor.pval,decreasing = F),]
results_metafor.Annot$p.adj.BH <- p.adjust(results_metafor.Annot$metafor.pval,method = "BH")
results_metafor.Annot$p.adj.bonf <- p.adjust(results_metafor.Annot$metafor.pval,method = "bonferroni")

manhattan <- data.frame("SNP"= as.character(rownames(results_metafor.Annot)),
                        "CHR"= 	as.numeric(as.character(results_metafor.Annot$CHR)), 
                        "BP"= 	results_metafor.Annot$MAPINFO,
                        "P"= 	results_metafor.Annot$metafor.pval)

manhattan<-manhattan[order(manhattan$BP),]
manhattan<-manhattan[order(manhattan$CHR),]

tiff(filename = paste0(out_pref,".meta.MANHATTAN.tiff"), units="in", width=7, height=5, res=300)
manhattan(manhattan,annotatePval=1e-6,annotateTop=FALSE,suggestiveline = FALSE, genomewideline=(-log10(1e-6)), cex.axis = 0.6 ,cex = 0.6, na.rm=na.omit,col = c("aquamarine4", "aquamarine2"))
dev.off()
tiff(filename = paste0(out_pref,".meta.QQplot.tiff"), units="in", width=3, height=3, res=300)
qq(manhattan$P, cex = 0.4, cex.axis = 0.6 , na.rm=na.omit)
dev.off()

dmr_meta <- data.frame("chrom" = 	paste0("chr", results_metafor.Annot$CHR), 
                       "start" = 	results_metafor.Annot$MAPINFO,
                       "end" 	= 	results_metafor.Annot$MAPINFO.1,
                       "pvalue" = 	results_metafor.Annot$metafor.pval)

dmr_meta<-dmr_meta[order(dmr_meta[,1],dmr_meta[,2]),]

write.table(dmr_meta, file=paste0(out_pref,".meta.DMR.bed"),sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
save(results_metafor.Annot,file = paste0(out_pref,'.meta.RData'))

#comb-p pipeline -c 4 --dist 500 --seed 1.0e-3 --anno hg19 -p  ./meta.DMR.cnv ./meta.DMR.bed
print("All done!")
#########################################################
