library("funr")
library(stringr)

args <- commandArgs(T)

#wd <- args[1]
input_file <- args[1]
backgr_file <- args[2]
out_pref <- args[3]

#setwd(wd)
print("Reading inputs...")
source(paste0(str_remove(sys.script(),pattern="EWAS_enrichment.R"),"crystalmethEPIC.R"))

sigcpg <- read.csv(file = input_file,header = F,stringsAsFactors = F)[,1]
allcpg <- read.csv(file = backgr_file,header = F,stringsAsFactors = F)[,1]

print("GO enrichment...")
go <- crystalmeth(sig.cpg = sigcpg,
                  all.cpg = allcpg,
                  collection = "GO",array.type = "EPIC" )

print("Saving GO results...")
go<-go[order(go$P.DE),]
#head(go[,1:6], 20)

write.csv(go,file = paste0(out_pref,"_GOEnrichment.csv"))

print("KEGG enrichment...")
kegg <- crystalmeth(sig.cpg = sigcpg,
                    all.cpg = allcpg,
                    collection = "KEGG", array.type = "EPIC" )

print("Saving KEGG results...")                    
kegg<-kegg[order(go$P.DE),]
write.csv(kegg,file = paste0(out_pref,"_KeggEnrichment.csv"))

print("All Done!")