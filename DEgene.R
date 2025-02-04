ls()
rm(list = ls(all.names = TRUE))
library(varhandle)
# varhandle package required
mydatGE = read.csv("prad_entrez.csv",fileEncoding = "UTF-8-BOM")
mydatGE[,1] <- gsub ("[.]" , "-" , mydatGE [,1])

rownames(mydatGE) <- mydatGE[,1]
mydatGE[,1] <- NULL

mydatGE[,1:62] <- log2(mydatGE[,1:62]+1)

mydatGE=unfactor(mydatGE[-1,])

mydatGE_t = as.data.frame(t(mydatGE))
#mydatGE_t = t(mydatGE)

row.names(mydatGE)


group = read.csv("prad_group.csv",row.names=1)
#group <- as.data.frame(group)

row.names(group)

patients = rownames(group)
patients = intersect(patients,rownames(mydatGE))

group = group[patients,]


mydatGE=mydatGE[patients,]
data <- as.data.frame(t(mydatGE))

#FC_threshold = 0.6
#P_threshold = 0.05

controlDat = data[,group$Group %in% 'c']
diseaseDat = data[,group$Group %in% 'd']

#controlDat$na_count <- apply(controlDat, 1, function(x) sum(is.na(x)))
#diseaseDat$na_count <- apply(diseaseDat, 1, function(x) sum(is.na(x)))

controlMean <- apply(controlDat,MARGIN=1,mean, na.rm=TRUE)
diseaseMean <- apply(diseaseDat,MARGIN=1,mean, na.rm=TRUE)
foldChange <- diseaseMean-controlMean

pvalues <- sapply(seq(nrow(controlDat)),function(x) t.test(t(controlDat[x,]),t(diseaseDat[x,]),paired=FALSE)$p.value)
#pvalues <- as.data.frame(pvalues)
names(pvalues) <- names(foldChange)
adj.p <- p.adjust(pvalues, method = "fdr", n = length(pvalues))
#adj.p <- as.data.frame(adj.p)
names(adj.p) <- names(foldChange)

DEGenes_FC <- foldChange
DEGenes_adj.p <- adj.p[names(DEGenes_FC)]
DEGenes_pvalues <- pvalues[names(DEGenes_FC)]
DEGenes_FC <- DEGenes_FC
DEGenes_pvalues <- pvalues[names(DEGenes_FC)]
DEGenes_adj.p <- adj.p[names(DEGenes_FC)]

combined <- cbind(DEGenes_FC, DEGenes_pvalues, DEGenes_adj.p)
View(combined)

