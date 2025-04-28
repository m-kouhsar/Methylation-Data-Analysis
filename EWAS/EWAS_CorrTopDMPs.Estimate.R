
args<-commandArgs(TRUE)
EWAS_file1 <- args[1]
EWAS_file2 <- args[2]
Num_cpg <- as.numeric(args[3])
out_prefix <- args[4]

print("Loading packages...")
suppressMessages(library(ggplot2))


print("Reading inputs...")
var <- load(EWAS_file1)
EWAS_data1 <- get(var[1])
rm(list = var)

var <- load(EWAS_file2)
EWAS_data2 <- get(var[1])
rm(list = var)
rm(var)

EWAS_data1=as.data.frame(EWAS_data1)
EWAS_data2=as.data.frame(EWAS_data2)

print("Selecting top CpGs...")
EWAS_data1 <- EWAS_data1[order(EWAS_data1$p.value,decreasing = F),]
EWAS_data1.selected <- EWAS_data1[1:Num_cpg,]
EWAS_data2.selected <- EWAS_data2[rownames(EWAS_data2) %in% rownames(EWAS_data1.selected),]

print(paste("From",Num_cpg,"selected CpGs in data1",nrow(EWAS_data2.selected),"are presented in data2"))

EWAS_data2.selected <- EWAS_data2.selected[rownames(EWAS_data2.selected) %in% rownames(EWAS_data1.selected),]
index <- match(rownames(EWAS_data2.selected) , rownames(EWAS_data1.selected))
EWAS_data1.selected <- EWAS_data1.selected[index,]

print("Row names are matched?")
identical(rownames(EWAS_data1.selected) , rownames(EWAS_data2.selected))

print("Calculating correlations...")
estimate <- cbind.data.frame(EWAS_data1.selected[,1],EWAS_data2.selected[,1])
names(estimate) <- c("Estimate1","Estimate2")
rownames(estimate) <- rownames(EWAS_data1.selected)
corr_ <- cor.test(estimate$Estimate1,estimate$Estimate2,method="pearson")

print("Saving Plot...")
cor.val <- round(corr_$estimate, 2)
cor.pval <- format(corr_$p.value,digit=4)
cor.label <- paste0("Correlation: ", cor.val," Pvalue: ",cor.pval)

pdf(file = paste0(out_prefix,".Top",Num_cpg,"Cpg.Corr.pdf"))
ggplot(estimate, aes(x = Estimate1, y = Estimate2)) + geom_point() + 
  geom_smooth(method = "lm", se = T) + ggtitle(cor.label)+xlim(+0.03,-0.03) + ylim(+0.03,-0.03)+
  xlab("Estimate values of top 100 probes in PITT-ADRC")+ylab("Estimate values of top 100 probes in BDR")
dev.off()

pdf(file = paste0(out_prefix,".Top",Num_cpg,"Cpg.Corr1.pdf"))
ggplot(estimate, aes(x = Estimate1, y = Estimate2)) + geom_point() + 
   ggtitle(cor.label) + geom_hline(aes(yintercept=0), colour="#BB0000", linetype="dashed") + 
  geom_vline(aes(xintercept=0), colour="#BB0000", linetype="dashed")+xlim(+0.03,-0.03) + ylim(+0.03,-0.03)+
  xlab("Estimate values of top 100 probes in PITT-ADRC")+ylab("Estimate values of top 100 probes in BDR")
dev.off()

