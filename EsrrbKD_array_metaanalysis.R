#Esrrb KD array meta-analysis

#GSE34170 - status: complete

library(lumi)
library(lumiMouseIDMapping)
library(lumiMouseAll.db)
library(limma)

#processing as lumibatch object
data <- read.table(file="C:/Users/Emma/Documents/PhD/Esrrb KD metaanalysis/GSE34170/GSE34170_non-normalized_noheader.txt",head=T)
lumi <- lumiR(file="C:/Users/Emma/Documents/PhD/Esrrb KD metaanalysis/GSE34170/GSE34170_non-normalized_noheader.txt",lib.mapping="lumiMouseIDMapping")

lumi.exprs <- exprs(lumi)
lumi.exprs[which(lumi.exprs<1)] <- 1
log <- log2(lumi.exprs)
norm <- lumiN(log,method="quantile")

plot(hclust(dist(t(norm))),labels=c("CTRL1","CTRL2","CTRL3","CTRL4",))
design <- array(0,dim=c(16,4))
colnames(design) <- c("GFP","EsrrbKD","Sox2KD","EsrrbSox2KD")
design[1:4,1] <- 1
design[5:8,2] <- 1
design[9:12,3] <- 1
design[13:16,4] <- 1
cont <- makeContrasts(EsrrbKD="EsrrbKD-GFP",levels=design)

fit <- lmFit(norm,design)
fit <- contrasts.fit(fit,cont)
fit <- eBayes(fit)

probe.list <- rownames(exprs(lumi))
if (require(lumiMouseAll.db) & require(annotate)) {
gene.symbol <- getSYMBOL(probe.list, 'lumiMouseAll.db')
fit$genes <- data.frame(ID=probe.list, SYMBOL=gene.symbol)
}

fit$genes$SYMBOL <- pData(featureData(data))$SYMBOL
hits <- topTable(fit,number=nrow(norm))
write.table(hits,"GSE34170hits.txt",sep="\t",row.names=F,quote=F)

#GSE4679 - status: ongoing

data <- as.matrix(read.table(file="C:/Users/Emma/Documents/PhD/Esrrb KD metaanalysis/GSE4679/GSE4679-GPL339_series_matrix.txt/GSE4679-GPL339_series_matrix.txt",skip=67,head=T,comment.char="",sep="\t",row.names=1))
data[which(data<1)] <- 1
log <- log2(data)
norm <- normalize.quantiles(log)
normSel <- norm[,1:32]
plot(hclust(dist(t(norm))))

design <- array(0,dim=c(70,16))
colnames(design) <- c("Emptyd0","Emptyd1","Emptyd2","Emptyd3","Emptyd4","Emptyd5","Emptyd6","Emptyd7","EsrrbKDd0","EsrrbKDd1","EsrrbKDd2","EsrrbKDd3","EsrrbKDd4","EsrrbKDd5","EsrrbKDd6","EsrrbKDd7")
design[1:2,1] <- c(1,1)
design[1:2,9] <- c(1,1)
design[3:4,2] <- c(1,1)
design[3:4,10] <- c(1,1)
design[5:6,3] <- c(1,1)
design[5:6,11] <- c(1,1)
design[7:8,4] <- c(1,1)
design[7:8,12] <- c(1,1)
design[9:10,5] <- c(1,1)
design[9:10,13] <- c(1,1)
design[11:12,6] <- c(1,1)
design[11:12,14] <- c(1,1)
design[13:14,7] <- c(1,1)
design[13:14,15] <- c(1,1)
design[15:16,7] <- c(1,1)
design[15:16,16] <- c(1,1)

design <- array(0,dim=c(32,2))
colnames(design) <- c("Empty","EsrrbKD")
design[,1] <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,rep(0,16))
design[,2] <- c(rep(0,16),1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8)

cont <- makeContrasts(EsrrbKD="EsrrbKD-Empty",levels=design)

fit <- lmFit(normSel,design)
fit <- contrasts.fit(fit,cont)
fit <- eBayes(fit)

design <- array(0,dim=c(32,2))
colnames(design) <- c("Empty","EsrrbKD")
design[1:16,1] <- c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7)
design[17:32,2] <- c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7)
fit <- lmFit(normSel,design)
fit <- eBayes(fit)

#this method results in fewer NAs, but more NAs in the top hits
probe.list <- rownames(exprs(lumi))
if (require(lumiMouseAll.db) & require(annotate)) {
gene.symbol <- getSYMBOL(probe.list, 'lumiMouseAll.db')
rownames(norm) <- gene.symbol
}

#best of both
symbol <- getSYMBOL(as.character(data[,1]),"illuminaMousev2.db")
symbols <- cbind(gene.symbol,symbol,rep(NA,length(symbol)))
symbols[,3] <- symbols[,1]
symbols[is.na(symbols[,1]),3] <- symbols[is.na(symbols[,1]),2]
fit$genes$SYMBOL <- symbols[,3]

#GSE26520 - status: complete
setwd("C:/Users/Emma/Documents/PhD/Esrrb KD metaanalysis/GSE26520")
library(limma)
GSE26520.targets <- readTargets("GSE26520targets.txt")
GSE26520.rg <- read.maimages(GSE26520.targets,source="agilent")
GSE26520.bgcorrected <- backgroundCorrect(GSE26520.rg,method="normexp")
GSE26520.ma <- normalizeWithinArrays(GSE26520.bgcorrected,method="loess")
GSE26520.ave <- avereps(GSE26520.ma, ID=GSE26520.ma$genes$ProbeName)
GSE26520.design <- c(0,0,1,1)
GSE26520.design <- modelMatrix(GSE26520.targets, ref="shCTRL")
GSE26520.fit <- lmFit(GSE26520.ave,GSE26520.design)
GSE26520.fit2 <- eBayes(GSE26520.fit)
hits <- topTable(GSE26520.fit2,coef=1,number=nrow(GSE26520.ave))
incl <- hits[which(hits$ProbeName %in% probe.list$OLIGO_NAME),]

Agilent015087 <- read.table("Agilent015087.txt",skip=10,head=T,sep="\t",comment.char="",quote="")
Agilent015087 <- Agilent015087[order(Agilent015087$OLIGO_NAME),]
Agilent015087 <- Agilent015087[9:nrow(Agilent015087),]
probe.list <- Agilent015087
probe.list <- probe.list[,c(2,7)]
probe.list <- probe.list[order(probe.list$OLIGO_NAME),]
probe.list <- probe.list[which(probe.list$OLIGO_NAME %in% incl$ProbeName),]
probe.list <- probe.list[which(!duplicated(probe.list$OLIGO_NAME)),]
identical(as.character(probe.list$OLIGO_NAME),incl$ProbeName)

incl <- Agilent015087[which(Agilent015087$OLIGO_NAME %in% hits$ProbeName),]
incl <- incl[order(incl$OLIGO_NAME),]
incl <- hits[which(hits$ProbeName %in% probe.list$OLIGO_NAME),]
incl <- incl[order(incl$ProbeName),]

annotated <- incl
annotated$GeneName <- probe.list[,2]
annotated <- annotated[order(annotated$B,decreasing=T),]
write.table(annotated,"GSE26520hits.txt",sep="\t",row.names=F,quote=F)
