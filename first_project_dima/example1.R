## ----include=FALSE-------------------------------------------------------
require(knitr)
opts_chunk$set(
concordance=FALSE, echo=TRUE, cache=TRUE, warning=FALSE, error=FALSE, message=FALSE)


## ----simulateData--------------------------------------------------------
expressionValues <- matrix (rnorm (300), nrow=30)
colnames(expressionValues) <- paste0("sample",1:10)
head(expressionValues)


## ----simulateCovariates--------------------------------------------------
targets <- data.frame(sampleNames = paste0("sample",1:10), group=c(paste0("CTL",1:5),paste0("TR",1:5)),age = rpois(10, 30), sex=as.factor(sample(c("Male", "Female"),10,replace=TRUE)))
head(targets, n=10)


## ----simulateGeneInfo----------------------------------------------------
myGenes <-  paste0("gene",1:30)


## ----simulateInfo--------------------------------------------------------
myInfo=list(myName="Alex Sanchez", myLab="Bioinformatics Lab", 
          myContact="alex@somemail.com", myTitle="Practical Exercise on ExpressionSets")
show(myInfo)


## ------------------------------------------------------------------------
pcs <- prcomp(expressionValues)
names(pcs)
barplot(pcs$sdev)
plot(pcs$rotation[,1], pcs$rotation[,2], col=targets$group, main="Representation of first two principal components")
text(pcs$rotation[,1], pcs$rotation[,2],targets$sex)


## ------------------------------------------------------------------------
variab <- apply(expressionValues, 1, sd)
orderedGenes <- myGenes[order(variab, decreasing=TRUE)]
head(variab[order(variab, decreasing=TRUE)])
head(orderedGenes)


## ----subsetExpressions---------------------------------------------------
newExpress<- expressionValues[,-9]
newTargets <- targets[-9,]
wrongNewTargets <- targets [-10,]


## ----loadPackage---------------------------------------------------------
require(Biobase)


## ----creaExpressionSet1--------------------------------------------------
myEset <- ExpressionSet(expressionValues)
class(myEset)
show(myEset)


## ----AnnotatedDataFrame2-------------------------------------------------
columnDesc <-  data.frame(labelDescription= c("Sample Names", "Treatment/Control", "Age at disease onset", "Sex of patient (Male/Female"))
myAnnotDF <- new("AnnotatedDataFrame", data=targets, varMetadata= columnDesc)
show(myAnnotDF)


## ------------------------------------------------------------------------
phenoData(myEset) <- myAnnotDF


## ----eval=FALSE----------------------------------------------------------
## # myEset <- ExpressionSet(assayData=expressionValues, phenoData=myAnnotDF)
## # Error in validObject(.Object) :
## #   invalid class ExpressionSet object: 1: sampleNames differ between assayData and phenoData
## # invalid class ExpressionSet object: 2: sampleNames differ between phenoData and protocolData


## ------------------------------------------------------------------------
rownames(pData(myAnnotDF))<-pData(myAnnotDF)$sampleNames
myEset <- ExpressionSet(assayData=expressionValues, phenoData=myAnnotDF)
show(myEset)


## ------------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues, 
                        phenoData=myAnnotDF, 
                        featureNames =myGenes)
# show(myEset)


## ----label=MIAME---------------------------------------------------------
myDesc <- new("MIAME", name= myInfo[["myName"]],
            lab= myInfo[["myLab"]],
            contact= myInfo[["myContact"]] ,
            title=myInfo[["myTitle"]])
print(myDesc)


## ------------------------------------------------------------------------
myEset <- ExpressionSet(assayData=expressionValues, 
                        phenoData=myAnnotDF,
                        fetureNames =myGenes,
                        experimentData = myDesc)
# show(myEset)


## ----usingExpressionSets-------------------------------------------------
dim(exprs(myEset))
class(phenoData(myEset))
class(pData(phenoData(myEset)))
head(pData(phenoData(myEset)))
head(pData(myEset))


## ------------------------------------------------------------------------
smallEset <- myEset[1:15,c(1:3,6:8)]
dim(exprs(smallEset))
dim(pData(smallEset))
head(pData(smallEset))
all(colnames(exprs(smallEset))==rownames(pData(smallEset)))


## ------------------------------------------------------------------------
youngEset <- myEset[,pData(myEset)$age<30]
dim(exprs(youngEset))
head(pData(youngEset))


## ------------------------------------------------------------------------
if (!require(GEOquery)) {
  BiocManager::install("GEOquery")
}
require(GEOquery)
install.packages('GEOquery')
gse <- getGEO("GSE507")
class(gse)
names(gse)
gse[[1]]
esetFromGEO <- gse[[1]]

