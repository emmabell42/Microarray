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
GSE26520.fit2$genes$Symbol <- getSYMBOL(GSE26520.fit2$genes$ID,"ARRAY")
GSEXXX.hits <- topTable(GSE26520.fit2,coef=2,number=nrow(GSE26520))
write.table(GSE39452LNCaPhits,"GSE39452LNCaPhits.txt",row.names=FALSE,quote=FALSE,sep=" ")
plot(density(topTable(fit2,coef=2,number=nrow(normexprs))$logFC))
qqnorm(topTable(fit2,coef=2,number=nrow(normexprs))$logFC)

