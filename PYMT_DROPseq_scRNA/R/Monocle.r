# Monocle
## Code used on Monocle object

### Libaries required
library(monocle)
library(Seurat)

### Functions

plotTraj<-function(monoclObj,genevec,filename,metadataDesired,...){
message('depending on your gene list, this might take ages \n you have been warned')
pdf(file=paste0("plots/",filename,".pdf"))
tmp<-monoclObj
tmp <- setOrderingFilter(tmp, genevec)
#print(plot_ordering_genes(tmp))
ncenters <- length(unlist(sapply(metadataDesired, function(x) unique( as(phenoData(HSMM),'data.frame')[,x]))))
tmp <- reduceDimension(tmp, max_components=2,method = 'DDRTree')#,auto_param_selection = F, nCenter = ncenters)
tmp <- orderCells(tmp, reverse=FALSE)
sapply(c(metadataDesired,"Pseudotime"), function(x) print(plot_cell_trajectory(tmp, color_by=x)))
return(tmp)
dev.off()
}

### Translating seurat object for monocle
# for Seurat object called epidat

pd <- AnnotatedDataFrame (epidat@meta.data)

gene.names <- rownames(epidat@raw.data)
fd <- data.frame(gene.names)
rownames(fd) <- gene.names
colnames(fd) <- c("gene_short_name")
fd <- AnnotatedDataFrame(fd)

# To create the monocle object
HSMM <- newCellDataSet(as(as.matrix(epidat@raw.data[rownames(epidat@meta.data)]), "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit=1, expressionFamily=negbinomial.size())
HSMM <- estimateSizeFactors(HSMM)
message('this takes ages, go get a coffee')
HSMM <- estimateDispersions(HSMM) 
disp_table <- dispersionTable(HSMM)

# Plot pseu-ordering of the cells
# First get genes to use in ordering
var_genes <- subset(disp_table,mean_expression >= 0.05 &dispersion_empirical >= 1 * dispersion_fit)$gene_id
all_cc_genes <- cc.genes[cc.genes %in% featureData(HSMM)$gene_short_name]
ordering_genes <- all_cc_genes[all_cc_genes %in% var_genes]


# Using gordon smyth DE genes
gorsmy<-plotTraj(HSMM,gorsmy_markers,"180308_gorsmy_lineage",c("res.2","Phase"))

#downsampling for equal numbers from both genotypes
ano<-pData(gorsmy)
cellkeep<-c(rownames(ano)[ano$isElfMouse=="Elf5"],sample(rownames(ano)[ano$isElfMouse=="WT",length(rownames(ano)[ano$isElfMouse=="Elf5"])])
)
cellkeep<-c(rownames(ano)[ano$isElfMouse=="Elf5"],sample(rownames(ano)[ano$isElfMouse=="WT",length(rownames(ano)[ano$isElfMouse=="Elf5"])])
length(rownames(ano)[ano$isElfMouse=="Elf5"])
cellkeep<-c(rownames(ano)[ano$isElfMouse=="Elf5"],sample(rownames(ano)[ano$isElfMouse=="WT"],length(rownames(ano)[ano$isElfMouse=="Elf5"])])
cellkeep<-c(rownames(ano)[ano$isElfMouse=="Elf5"],sample(rownames(ano)[ano$isElfMouse=="WT"],length(rownames(ano)[ano$isElfMouse=="Elf5"])
)
)

# plots 
x<-gorsmy[,cellkeep]
plot_cell_trajectory(x,markers="Elf5",markers_linear=T,color_by="res.1.2") + facet_grid(~isElfMouse) + ggtitle("SmythMarkerTrajectories, colour by  cluster, point size by Elf5, split by genotype")

