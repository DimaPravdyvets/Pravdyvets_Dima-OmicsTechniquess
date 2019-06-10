## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- include=FALSE------------------------------------------------------
if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}

installifnot("pd.mogene.1.0.st.v1")
installifnot("mogene10sttranscriptcluster.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")


## ------------------------------------------------------------------------
targets <- read.csv('./data/targets.csv', header = T,sep = ';')
targets

## ----include=FALSE-------------------------------------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "data")
resultsDir <- file.path(workingDir, "results")
setwd(resultsDir)


## ---- include=F----------------------------------------------------------
CELfiles <- list.celfiles(file.path(dataDir))
rawData <- read.celfiles(file.path(dataDir,CELfiles))


## ------------------------------------------------------------------------
sampleNames <- as.character(targets$ShortNAme)
sampleColor <- as.character(targets$Colors)


## ------------------------------------------------------------------------
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
clust.euclid.averageR <- hclust(dist(t(exprs(rawData))),method="average")
plot(clust.euclid.averageR, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-100000, max(pcX$x[,1])+100000),ylim=c(min(pcX$x[,2])-100000, max(pcX$x[,2])+100000))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

pdf(file.path(resultsDir, "QCPlots_Raw.pdf"))
boxplot(rawData, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
plot(clust.euclid.averageR, labels=sampleNames, main="Hierarchical clustering of samples of RawData", 
     cex=0.7,  hang=-1)
plotPCA(exprs(rawData), labels=sampleNames, dataDesc="raw data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()


## ---- include=F----------------------------------------------------------
eset<-rma(rawData)

write.exprs(eset, file.path(resultsDir, "NormData.txt"))



## ------------------------------------------------------------------------
#BOXPLOT
boxplot(eset, las=2, main="Intensity distribution of Normalized data", cex.axis=0.6, 
        col=sampleColor, names=sampleNames)



## ------------------------------------------------------------------------
clust.euclid.averageN <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.averageN, labels=sampleNames, main="Hierarchical clustering of NormData", 
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
plotPCA <- function ( X, labels=NULL, colors=NULL, dataDesc="", scale=FALSE, formapunts=NULL, myCex=0.8,...)
{
  pcX<-prcomp(t(X), scale=scale) # o prcomp(t(X))
  loads<- round(pcX$sdev^2/sum(pcX$sdev^2)*100,1)
  xlab<-c(paste("PC1",loads[1],"%"))
  ylab<-c(paste("PC2",loads[2],"%"))
  if (is.null(colors)) colors=1
  plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=colors, pch=formapunts, 
       xlim=c(min(pcX$x[,1])-10, max(pcX$x[,1])+10),ylim=c(min(pcX$x[,2])-10, max(pcX$x[,2])+10))
  text(pcX$x[,1],pcX$x[,2], labels, pos=3, cex=myCex)
  title(paste("Plot of first 2 PCs for expressions in", dataDesc, sep=" "), cex=0.8)
}

plotPCA(exprs(eset), labels=sampleNames, dataDesc="NormData", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)

#plots normalized data
pdf(file.path(resultsDir, "QCPlots_Norm.pdf"))
boxplot(eset, which="all",las=2, main="Intensity distribution of Normalized data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)
plot(clust.euclid.averageN, labels=sampleNames, main="Hierarchical clustering of samples of Normalized", 
     cex=0.7,  hang=-1)
plotPCA(exprs(eset), labels=sampleNames, dataDesc="Normalized data", colors=sampleColor,
        formapunts=c(rep(16,4),rep(17,4)), myCex=0.6)
dev.off()

## ------------------------------------------------------------------------
arrayQualityMetrics(rawData,  reporttitle="QualityControl", force=TRUE)
annotation(eset) <- "org.Mm.eg.db"
eset_filtered <- nsFilter(eset, var.func=IQR,
         var.cutoff=0.75, var.filter=TRUE,
         filterByQuantile=TRUE)


## ------------------------------------------------------------------------
print(eset_filtered$filter.log$numLowVar)


## ------------------------------------------------------------------------
print(eset_filtered$eset)


## ------------------------------------------------------------------------
treat <- targets$grupos
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design)<-c("it","ct","ip","cp")
rownames(design) <- sampleNames
print(design)


## ------------------------------------------------------------------------
cont.matrix1 <- makeContrasts (
  itvsct = it-ct,
  itvsip = it-ip,
  itvscp = it-cp,
  ctvsip = ct-ip,
  ctvscp = ct-cp,
  ipvscp = ip-cp,
  levels=design)
comparison1 <- "Infection or control vs polyA or Total"
cont.matrix1


## ------------------------------------------------------------------------
fit1 <- lmFit(eset_filtered$eset, design)
fit.main1 <- contrasts.fit(fit1, cont.matrix1)
fit.main1


## ------------------------------------------------------------------------
fit.main1 <- eBayes(fit.main1)
fit.main1


## ------------------------------------------------------------------------
topTab_itvsct <- topTable (fit.main1, number=nrow(fit.main1), coef="itvsct", adjust="fdr"); 
head(topTab_itvsct)
topTab_itvsip <- topTable (fit.main1, number=nrow(fit.main1), coef="itvsip", adjust="fdr"); head(topTab_itvsip)
topTab_itvscp <- topTable (fit.main1, number=nrow(fit.main1), coef="itvscp", adjust="fdr"); head(topTab_itvsct)
topTab_ctvsip <- topTable (fit.main1, number=nrow(fit.main1), coef="ctvsip", adjust="fdr"); head(topTab_ctvsip)
topTab_ctvscp <- topTable (fit.main1, number=nrow(fit.main1), coef="ctvscp", adjust="fdr"); head(topTab_ctvscp)
topTab_ipvscp <- topTable (fit.main1, number=nrow(fit.main1), coef="ipvscp", adjust="fdr"); head(topTab_ipvscp )


#topTab<-Reduce(merge,list(topTab_itvscp[1:round(NROW(topTab_itvscp)*0.05),],topTab_itvsct[1:round(NROW(topTab_itvsct)*0.05),],topTab_itvsip[1:round(NROW(topTab_itvsip)*0.05),],topTab_ctvscp[1:round(NROW(topTab_ctvscp)*0.05),],topTab_ctvsip[1:round(NROW(topTab_ctvsip)*0.05),],topTab_ipvscp[1:round(NROW(topTab_ipvscp)*0.05),]))
topTab<-do.call('rbind',list(topTab_itvscp[1:round(NROW(topTab_itvscp)*0.05),],topTab_itvsct[1:round(NROW(topTab_itvsct)*0.05),],topTab_itvsip[1:round(NROW(topTab_itvsip)*0.05),],topTab_ctvscp[1:round(NROW(topTab_ctvscp)*0.05),],topTab_ctvsip[1:round(NROW(topTab_ctvsip)*0.05),],topTab_ipvscp[1:round(NROW(topTab_ipvscp)*0.05),]))
topTab


## ------------------------------------------------------------------------
volcanoplot(fit.main1, highlight=10, names=fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep="\n"))
abline(v = c(-3, 3))


pdf(file.path(resultsDir,"Volcanos.pdf"))
volcanoplot(fit.main1, highlight = 10, names = fit.main1$ID, 
            main = paste("Differentially expressed genes", colnames(cont.matrix1), sep = "\n"))
abline(v = c(-3, 3))
dev.off()



## ------------------------------------------------------------------------

my_frame <- data.frame(exprs(eset))
head(my_frame)
HMdata <- merge(my_frame, topTab, by.x = 0, by.y = 0)
rownames(HMdata) <- HMdata$Row.names
HMdata <- HMdata[, -c(1,10:15)]
head(HMdata)
HMdata2 <- data.matrix(HMdata, rownames.force=TRUE)
head(HMdata2)
write.csv2(HMdata2, file = file.path(resultsDir,"Data2HM.csv"))

#HEATMAP PLOT
my_palette <- colorRampPalette(c("blue", "red"))(n = 299)
ncol(my_frame)

heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap raw data",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)

#EXPORT TO PDF FILE
pdf(file.path(resultsDir,"HeatMap.pdf"))
heatmap.2(HMdata2,
          Rowv=TRUE,
          Colv=TRUE,
          main="HeatMap rawdata",
          scale="row",
          col=my_palette,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          cexRow=0.5,
          cexCol=0.9,
          key=TRUE,
          keysize=1.5,
          density.info="histogram",
          ColSideColors=sampleColor,
          tracecol=NULL,
          srtCol=30)
dev.off()


## ------------------------------------------------------------------------
require(pd.hugene.2.0.st)
require(hugene20sttranscriptcluster.db)

columns(hugene20sttranscriptcluster.db)
probes_tot<-rownames(unique(topTab))
write.csv(select(hugene20sttranscriptcluster.db,probes_tot,
                     columns = c("SYMBOL","ENSEMBL","ENTREZID","PROBEID","UNIGENE","UNIPROT","REFSEQ","GENENAME")),"anotations.csv")

