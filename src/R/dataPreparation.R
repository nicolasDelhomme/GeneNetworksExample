#' ---
#' title: "Biological QA"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(parallel)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' * Data
load("analysis/gene_network/raw_counts.rda")

#' * Metadata
#' Sample information
samples <- data.frame(
  SampleID=sub("-","_",colnames(counts)),
  Condition=gsub("count\\.|_[1-3]|-[1-3]","",colnames(counts)),
  Dataset=ifelse(grepl("count",colnames(counts)),"D2","D1"))
colnames(counts) <- samples$SampleID

counts[is.na(counts)] <- 0

#' ## Quality Control
#' * Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' * Let us take a look at the sequencing depth, coloring by Dataset
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=Dataset)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by Dataset 
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Condition=samples$Condition[match(ind,samples$SampleID)]) %>% 
  mutate(Dataset=samples$Dataset[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=Dataset)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
dds <- DESeqDataSetFromMatrix(
  countData=counts,
  colData = samples,
  design = ~ Dataset)

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
boxplot(sizeFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation is definitely affected by the datasets variability
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' * Cumulative components effect
#' 
#' We define the number of variable of the model
nvar=1

#' An the number of possible combinations
nlevel=nlevels(dds$Dataset)

#' We plot the percentage explained by the different components, the
#' red line represent the number of variable in the model, the orange line
#' the number of variable combinations.
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)
  
#' ### 2D
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(dds)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Dataset,shape=Condition,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ### Sequencing depth
#' Number of genes expressed per condition at different cutoffs
conds <- factor(dds$Condition)
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' ### Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)
vst.cutoff <- 3

#' * Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
          distfun=pearson.dist,
          hclustfun=function(X){hclust(X,method="ward.D2")},
          labRow = NA,trace = "none",
          labCol = conds,
          col=hpal)

plot(as.hclust(hm$colDendrogram),xlab="",sub="")

#' ### Hierarchical clustering
#' Done to assess the previous dendrogram's reproducibility
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                       method.hclust = "ward.D2", 
                       nboot = 1000, parallel = TRUE)

#' plot the clustering with bp and au
plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

#' bootstrapping results as a table
print(hm.pvclust, digits=3)

#' ## Conclusion QA
#' The difference in sequencing depth is the factor separating the datasets.
#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Normalisation
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked better using the model, but the sequencing depth
#' effect is still visible
#' 
meanSdPlot(vst[rowSums(vst)>0,])

#' # Low expression filtering
#' We will use DESeq2 independence filtering to remove the lowly expressed genes
#' by performing a DE analysis on the Dataset variable
dds <- DESeq(dds)
res <- results(dds)

#' We use the independent filtering cutoff on the expression
#' to select the genes
sel <- res$baseMean > metadata(res)$filterThreshold

sprintf("%s genes are expressed above the independent filtering cutoff",
        sum(sel))

#' # Export
dir.create(here("analysis/gene_network/seidr"),showWarnings=FALSE)

#' * gene by column, without names matrix
write.table(t(vst[sel,]),
            file=here("analysis/gene_network/seidr/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(rownames(vst)[sel]),
            file=here("analysis/gene_network/seidr/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' 
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
