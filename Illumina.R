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

