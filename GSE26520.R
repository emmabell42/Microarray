setwd("X:\\Current lab members\\Emma Bell\\Data\\Microarray\\GSE26520")
GSE26520 <- read.table("X:\\Current lab members\\Emma Bell\\Data\\Microarray\\GSE26520\\GSE26520_series_matrix.txt",head=T,skip=112,row.names=1,stringsAsFactors=F,comment.char="",quote="")
GSE26520 <- GSE26520[,c(1:130,147:228,237:ncol(GSE26520))]
GSM  <- read.table("GSM.txt",head=F,sep="\t",stringsAsFactors=F)
GSM <- GSM[c(1:130,147:228,237:nrow(GSM)),]
GSM[grep("Control",GSM[,2]),2] <- "Control"
GSM[,2] <- gsub(" ","",GSM[,2])
GSM[,2] <- gsub("\\-","",GSM[,2])
GSM[3:4,2] <- c("Pou5f1_2_72h_shRNAtrasient","Pou5f1_2_72h_shRNAtrasient")
GSM[,2] <- gsub("4932441K18Rik_72h_shRNAtrasient","A4932441K18Rik_72h_shRNAtrasient",GSM[,2])
boxplot(log(GSE26520,base=2))
plot(hclust(dist(t(GSE26520))), labels=GSM[,2])
# Log10(x)=y 10^y=x
GSE26520.exp <- 10^GSE26520
GSE26520.log <- log2(GSE26520.exp)
GSE26520.design <- array(0,dim=c(nrow(GSM),length(unique(GSM[,2]))))
conditions <- unique(GSM[,2])
colnames(GSE26520.design) <- conditions
for(i in 1:length(conditions)){
GSE26520.design[which(GSM[,2]==conditions[i]),i] <- 1
}
library(limma)
toContrast <- paste0(conditions,"-Control")
toContrast <- gsub("Control-Control",NA,toContrast)
toContrast <- na.omit(toContrast)
contrast.matrix <- makeContrasts(contrasts=toContrast, levels=GSE26520.design)
GSE26520.fit <- lmFit(GSE26520.log,GSE26520.design)
GSE26520.fit2 <- contrasts.fit(GSE26520.fit,contrast.matrix)
GSE26520.fit2 <- eBayes(GSE26520.fit2)
annot <- read.table("annotation.txt",sep="\t",head=T,stringsAsFactors=F,comment.char="",quote="")
geneSymbols <- array(NA,dim=c(nrow(GSE26520),2))
geneSymbols[,1] <- rownames(GSE26520)
colnames(geneSymbols) <- c("ID","SYMBOL")
for(i in 1:nrow(geneSymbols)){
geneSymbols[i,2] <- annot[which(annot[,1]==geneSymbols[i,1]),7]
}
GSE26520.fit2$genes <- geneSymbols
hits <- gsub("-Control","",toContrast)
hitsLists <- as.list(rep(NA,length(hits)))
names(hitsLists) <- hits
for(i in 1:length(hits)){
hitsLists[[i]] <- topTable(GSE26520.fit2,coef=i,number=nrow(GSE26520))
}
############################################################
#
## Plot the gene expression of the shRNA inhibited genes
#
############################################################
genes <- sapply(strsplit(toContrast, "_"),"[",1)
genes[which(!genes %in% geneSymbols[,2])] <- c("Zscan10","4932441K18Rik","Pou5f1","Nanog","Nr0b1","NA","NA","Pou5f1","Pou5f1")
genes <- na.omit(genes)
inhibition <- data.frame(array(NA,dim=c(length(genes),3)))
colnames(inhibition) <- c("Target","LogFC","-10*log(adj. P-value)")
inhibition[,1] <- genes
for(i in 1:length(genes)){
target <- genes[i]
hitsList <- hitsLists[[i]]
hitsList <- hitsList[which(!duplicated(hitsList$SYMBOL)),]
inhibition[i,2] <- as.numeric(hitsList[which(hitsList$SYMBOL==target),3])
inhibition[i,3] <- as.numeric(hitsList[which(hitsList$SYMBOL==target),7])
}
inhibition[,3] <- -10*(log10(inhibition[,3]))
plot(inhibition[,2],inhibition[,3],pch=20,xlab="LogFC",ylab="-10*log(Adj. P-value)",main="shRNA KD")
abline(h=13.0103,lty=2)
