####)1
library(Matrix)
library(Seurat)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(reshape2)
library(ggraph)
load("r_objects/ALLs15.rdata")
source("source/functions.r")
source("source/clustering_tree.R")
load("r_objects/cell.cyclegenes.rdata")

fibcell<-rownames(ALLs15@meta.data)[ALLs15@meta.data$res.1 %in% c(5,9,8,15)]
fibdat<-testInput(ALLs15@raw.data[,fibcell])

# mingenes= 200
# mincells= 3
# total expr= 10000
# [1] "setting up"
# [1] "annotating mito genes"

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767    WT37672  WT4108AC1  WT4108AC2
#         47         47        128         97         14         93         27        138        112         97        104        351        113
#  WT4108BD1  WT4108BD2
#        506        434
# subsetting data
# gene accept hi= 4000

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767  WT4108AC1  WT4108AC2  WT4108BD1
#         47         47        128         97         14         93         27        138        112         97        351        113        506
#  WT4108BD2    WT37672
#        434        104
# nUMI accept hi= 12000

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767  WT4108AC1  WT4108AC2  WT4108BD1
#         43         44        127         97         12         93         23        132        107         94        348        111        491
#  WT4108BD2    WT37672
#        432        102
# mitoCutoff= 0.15

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767  WT4108AC1  WT4108AC2  WT4108BD1
#         43         44        127         97         12         93         23        131        107         94        348        111        491
#  WT4108BD2    WT37672
#        432        102
#  "4743 variable genes"

fibdat<-mkPCA(fibdat)
fibdat<-addmetadata(fibdat)
save(fibdat,file="r_objects/180725_fib_ano.rdata")
pdf(file='plots/180725_fib_qc.pdf')
plotstuff(fibdat)
dev.off()

devtools:: install_version("ggplot2", version = "2.2.1")
install.packages("clustree")

devtools::install_github("lazappi/clustree")
library(clustree)
clustree(fibdat,prefix='res.')
ggsave('180725_fib_clustree.pdf')

plotPie<-function(seuratObj=fibdat,primeFactor="isElfMouse",res=0.7){
    x<-table(seuratObj@meta.data[,primeFactor],seuratObj@meta.data[,paste0('res.',as.character(res))])
    x<-x/rowSums(x)
    apply(x,2,pie)
    plot(1)
return(x)
}
par(mfrow=c(4,6))
plotPie(res=0.1)
plotPie(res=0.4)
dev.copy2pdf(file="180725_fib_pie.pdf")

fibdat<-SetIdent(fibdat,ident.use=paste0(fibdat@meta.data$res.0.1,"_",fibdat@meta.data$isElfMouse))
markers_res.0.1<-FindAllMarkers(fibdat,min.pct = 0.25, thresh.use = 0.25)
table(fibdat@ident)
# 0_Elf5   0_WT 1_Elf5   1_WT 2_Elf5   2_WT
#    189    859    192    415     58    542

library(dplyr)
markers_res.0.1 %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
	DoHeatmap(fibdat,genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
ggsave("180725_fib_top10heatmap.pdf")

###to do

markers_res.0.1 %>% group_by(cluster) %>% top_n(10, abs(pct.1-pct.2)) -> top10


fibdat<-SetIdent(fibdat,ident.use=fibdat@meta.data$isElfMouse)
markers_genotype<-FindAllMarkers(fibdat,min.pct = 0.25, thresh.use = 0.25)
markers_genotype %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
DoHeatmap(fibdat,genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)


fibdat<-SetIdent(fibdat,ident.use=paste0(fibdat@meta.data$res.0.1,"_",fibdat@meta.data$isElfMouse))
VlnPlot(fibdat,features.plot=c("Pdgfra","Pdgfrb"))
ggsave("180725_fib_vln_pdgf.pdf")

VlnPlot(fibdat,features.plot=c("Pdgfra","Pdgfrb","Col1a1","Ctla2a",'Mcam',"Pdpn"))

### Fa wants 
1) plots
2) marker lists (by SSM) WT vs Elf5 and res.0.1 indep of eachother
3) 



##### redoing by cca
load("r_objects/180725_fib_ano.rdata")
genotype<-split(rownames(fibdat@meta.data),fibdat@meta.data$isElfMouse)
remakeObj <- function(seuratObj,cellsUse) {
tmp <- CreateSeuratObject(raw.data = seuratObj@data[,cellsUse])
tmp <- NormalizeData(object = tmp)
tmp <- ScaleData(object = tmp)
tmp <- FindVariableGenes(object = tmp, do.plot = FALSE)
return(tmp)
}
fib_bygeno<-lapply(genotype, function(x) remakeObj(fibdat,x))
vgenes<-lapply(fib_bygeno, function(x) rownames(x@hvg.info)[1:2000])
uvgenes<-union(vgenes[[1]],vgenes[[2]])
save(fib_bygeno,uvgenes,file='r_objects/180801_fibdatByGenotype.rdata')
fib_cca <- RunCCA(fib_bygeno[[1]],fib_bygeno[[2]],genes.use=uvgenes) #crapped out, try on a clean run

fib_cca@meta.data$isElfMouse<-ifelse(grepl("WT",rownames(fib_cca@meta.data)),"WT","Elf5")
fib_cca<-mkPCA(fib_cca)
fib_cca<-addmetadata_cca(fib_cca)
save(fib_cca,file="r_objects/180801_fibCCA_ano.rdata")

clustree(fib_cca,prefix='res.')
ggsave('plots/180801_fib_cca_clusterTree.pdf')
par(mfrow=c(4,4))
plotPie(seuratObj=fib_cca,res=0.1) 
plot(1,1)
plotPie(seuratObj=fib_cca,res=0.4)
dev.copy2pdf(file="plots/180801_fib_cca_pie.pdf")

compPlot(fibdat,fib_cca,groupVar='isElfMouse')
ggsave("plots/180801_fib_comparison_PCA.pdf")
compPlot(fibdat,fib_cca,"tsne")
ggsave("plots/180801_fib_comparison_TSNE.pdf")

compPlot(fibdat,fib_cca,groupVar='res.0.4')
ggsave("plots/180801_fib_comparison_PCA_res.0.4.pdf")
compPlot(fibdat,fib_cca,"tsne",groupVar='res.0.4')
ggsave("plots/180801_fib_comparison_TSNE_res.0.4.pdf")

fib_cca<-SetIdent(fib_cca,ident.use=fib_cca@meta.data$res.0.4)
markers_res.0.4<-FindAllMarkers(fib_cca,min.pct = 0.25, thresh.use = 0.25,reduction.type = "cca.aligned")
save(markers_res.0.4,file="fib_cca_markers0.4.rdata")
library(dplyr)
markers_res.0.4 %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top10

SplitDotPlotGG(fib_cca,genes.plot=unique(top10$gene),grouping.var="isElfMouse", do.return=T) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle('cca markers, cca clusters res 0.45')
ggsave('plots/180801_fib_res.0.4Top5_splitdotplot.pdf')



fib_cca<-SetIdent(fib_cca,ident.use=paste0(fib_cca@meta.data$res.0.4,"_",fib_cca@meta.data$isElfMouse))
markers_res.0.4_geno<-FindAllMarkers(fib_cca,min.pct = 0.25, thresh.use = 0.25,reduction.type = "cca.aligned")
save(markers_res.0.4_geno,file="fib_cca_markers0.4_byGenotype.rdata")
library(dplyr)
markers_res.0.4_geno %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top10

DotPlot(fib_cca,genes.plot=unique(top10$gene), do.return=T) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle('cca markers, cca clusters res 0.45')
ggsave('plots/180801_fib_res.0.4Top5_geno_dotplot.pdf')




##### revised post analysis to exclude cluster 3 from res 0.4

library(Matrix)
library(Seurat)
library(biomaRt)
library(clusterProfiler)
library(ReactomePA)
library(dplyr)
library(reshape2)
library(ggraph)
load("r_objects/ALLs15.rdata")
source("source/functions.r")
source("source/clustering_tree.R")
load("r_objects/cell.cyclegenes.rdata")
load("r_objects/180801_fibCCA_ano.rdata")

fibcell<-rownames(fib_cca@meta.data)[!fib_cca@meta.data$res.0.4 %in% 3]
fibdat<-testInput(fib_cca@raw.data[,fibcell])
# mingenes= 200
# mincells= 3
# total expr= 10000
# [1] "setting up"
# [1] "annotating mito genes"

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767    WT37672
#         42         42        123         96         12         90         23        129        104         91        100
#  WT4108AC1  WT4108AC2  WT4108BD1  WT4108BD2
#        341        107        470        412
# subsetting data
# gene accept hi= 4000

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767    WT37672
#         42         42        123         96         12         90         23        129        104         91        100
#  WT4108AC1  WT4108AC2  WT4108BD1  WT4108BD2
#        341        107        470        412
# nUMI accept hi= 12000

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767    WT37672
#         42         42        123         96         12         90         23        129        104         91        100
#  WT4108AC1  WT4108AC2  WT4108BD1  WT4108BD2
#        341        107        470        412
# mitoCutoff= 0.15

# ELF3759AC4 ELF3759AC7  ELF3759BD   ELF3768A   ELF3768B   ELF3772A   ELF3772B     WT3752     WT3758     WT3767    WT37672
#         42         42        123         96         12         90         23        129        104         91        100
#  WT4108AC1  WT4108AC2  WT4108BD1  WT4108BD2
#        341        107        470        412
# [1] "5081 variable genes"
########################################

fibdat<-mkPCA(fibdat)
fibdat<-addmetadata(fibdat)
save(fibdat,file="r_objects/181017_fib_ano.rdata")
pdf(file='plots/181017_fib_qc.pdf')
plotstuff(fibdat)
dev.off()

devtools:: install_version("ggplot2", version = "2.2.1")
install.packages("clustree")


library(clustree)
clustree(fibdat,prefix='res.')
ggsave('181017_fib_clustree.pdf')


plotPie<-function(seuratObj=fibdat,primeFactor="isElfMouse",res=0.7){
    x<-table(seuratObj@meta.data[,primeFactor],seuratObj@meta.data[,paste0('res.',as.character(res))])
    x<-x/rowSums(x)
    apply(x,2,pie)
    plot(1)
return(x)
}
par(mfrow=c(4,6))
plotPie(res=0.1)
plotPie(res=0.3)
plotPie(res=0.7)
dev.copy2pdf(file="181017_fib_pie.pdf")

##CCA
genotype<-split(rownames(fibdat@meta.data),fibdat@meta.data$isElfMouse)
remakeObj <- function(seuratObj,cellsUse) {
tmp <- CreateSeuratObject(raw.data = seuratObj@data[,cellsUse])
tmp <- NormalizeData(object = tmp)
tmp <- ScaleData(object = tmp)
tmp <- FindVariableGenes(object = tmp, do.plot = FALSE)
return(tmp)
}
fib_bygeno<-lapply(genotype, function(x) remakeObj(fibdat,x))
vgenes<-lapply(fib_bygeno, function(x) rownames(x@hvg.info)[1:2000])
uvgenes<-union(vgenes[[1]],vgenes[[2]])
save(fib_bygeno,uvgenes,file='r_objects/181017_fibdatByGenotype.rdata')
fib_cca <- RunCCA(fib_bygeno[[1]],fib_bygeno[[2]],genes.use=uvgenes) #crapped out, try on a clean run

fib_cca@meta.data$isElfMouse<-ifelse(grepl("WT",rownames(fib_cca@meta.data)),"WT","Elf5")
fib_cca<-mkPCA(fib_cca)
fib_cca<-addmetadata_cca(fib_cca)
save(fib_cca,file="r_objects/181017_fibCCA_ano.rdata")

clustree(fib_cca,prefix='res.')
ggsave('plots/181017_fib_cca_clusterTree.pdf')
par(mfrow=c(4,4))
plotPie(seuratObj=fib_cca,res=0.2) 
plot(1,1)
plotPie(seuratObj=fib_cca,res=0.4)
dev.copy2pdf(file="plots/181017_fib_cca_pie.pdf")

compPlot(fibdat,fib_cca,groupVar='isElfMouse')
ggsave("plots/181017_fib_cca_comparison_PCA.pdf")
compPlot(fibdat,fib_cca,"tsne")
ggsave("plots/181017_fib_cca_comparison_TSNE.pdf")

compPlot(fibdat,fib_cca,groupVar='res.0.4')
ggsave("plots/181017_fib_cca_comparison_PCA_res.0.4.pdf")
compPlot(fibdat,fib_cca,"tsne",groupVar='res.0.4')
ggsave("plots/181017_fib_cca_comparison_TSNE_res.0.4.pdf")

fib_cca<-SetIdent(fib_cca,ident.use=fib_cca@meta.data$res.0.4)
markers_res.0.4<-FindAllMarkers(fib_cca,min.pct = 0.25, thresh.use = 0.25,reduction.type = "cca.aligned")
save(markers_res.0.4,file="r_objects/181017_fib_cca_markers0.4.rdata")
library(dplyr)
markers_res.0.4 %>% group_by(cluster) %>% top_n(5, avg_logFC) -> top10

SplitDotPlotGG(fib_cca,genes.plot=unique(top10$gene),grouping.var="isElfMouse", do.return=T) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +ggtitle('cca markers, cca clusters res 0.4')
ggsave('plots/181017_fib_cca_res.0.4Top5_splitdotplot.pdf')