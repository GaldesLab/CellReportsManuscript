##mod in march 2018 to use DDRtree (coz of biorxiv paper)
library(monocle)
load("r_objects/180209_monocleObjs.rdata")
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


varall<-plotTraj(HSMM,var_genes,"180308_epith_var",c("res.2","Phase"))
save(varall,file="r_objects/180308_PseudoVarMonocle.rdata")

gorsmy_downSampleByGenotype <-plotTraj(HSMM,gorsmy_markers,"180612_downSampled_gorsmy_lineage",c("res.2","Phase"))
save(gorsmy_downSampleByGenotype,file="r_objects/180612_PseudoGorsmyMonocle.rdata")
##up to here

ccvar<-plotTraj(HSMM,ordering_genes,"180308_epith_ccVar",c("res.2","Phase"))
ccall<-plotTraj(HSMM,all_cc_genes,"1800308_epith_ccAll",c("res.2","Phase"))






#To Run Monocle
###########################
library(Seurat)
library(monocle)
load("r_objects/new_epi_ano.rdata")
load("r_objects/cell.cyclegenes.rdata")
#Tp remove non-relevant clusters - not part of the  differentiating cells


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

save(HSMM,gorsmy_markers,var_genes,all_cc_genes,ordering_genes, file="r_objects/180209_monocleObjs.rdata")

#Now use these genes to order cells and plot the cell trajectory
plotTraj<-function(monoclObj,genevec,filename,metadataDesired,...){
message('depending on your gene list, this might take ages \n you have been warned')
pdf(file=paste0("plots/",filename,".pdf"))
tmp<-monoclObj
tmp <- setOrderingFilter(tmp, genevec)
#print(plot_ordering_genes(tmp))
ncenters <- length(unlist(sapply(metadataDesired, function(x) unique( as(phenoData(HSMM),'data.frame')[,x]))))
tmp <- reduceDimension(tmp, max_components=2,method = 'DDRTree')#,auto_param_selection = F, nCenter = ncenters)
tmp <- orderCells(tmp, reverse=FALSE)
sapply(metadataDesired, function(x) print(plot_cell_trajectory(tmp, color_by=x)))
return(tmp)
dev.off()
}


varall<-plotTraj(HSMM,var_genes,"180308_epith_var",c("res.2","Phase"))
save(varall,file="r_objects/180308_PseudoVarAllMonocle.rdata")

gorsmy<-plotTraj(HSMM,gorsmy_markers,"180308_gorsmy_lineage",c("res.2","Phase"))
plot_cell_trajectory(gorsmy,color_by="Pseudotime") + facet_wrap(~isElfMouse,nrow=1)
dev.copy2pdf(file="plots/180309_gorsmyTraject_byPseudoAndGenotype.pdf")

gorsmy_BEAM<-lapply(1:3,function(x) BEAM(gorsmy,branch_point=x,cores=6))
gorsmy_trimBeam<-lapply(gorsmy_BEAM, function(x){
    return(x[order(x$qval),c("gene_short_name", "pval", "qval")])
})

plot_genes_branched_heatmap(gorsmy[row.names(subset(gorsmy_trimBeam[[2]],qval<1e-5)),],branch_point=2,cores=2,use_gene_short_name=T,show_rownames=T)
dev.copy2pdf(file="plots/180309_gorsmyTraject_BP1_heatMap.pdf") #rep for 2 and 3, had to set q<1e-50 coz numbers looking at data, it doesnt really look real

#more complicated graphing
plot_cell_trajectory(gorsmy,markers="Elf5",markers_linear=T,color_by="res.1.2") + facet_grid(~isElfMouse) + ggtitle("SmythMarkerTrajectories, colour by  cluster, point size by Elf5, split by genotype")
ggsave("monocle_trajectories_by_genotypeAndElfExpAndPhase.pdf")

save(gorsmy,gorsmy_BEAM,gorsmy_trimBeam,file="r_objects/180308_PseudoGorsmyMonocle.rdata")

###plotting and annotating SeuratObj
pd<-as(pData(gorsmy),'data.frame')
pd$x=y[rownames(pd),1]
pd$y=y[rownames(pd),2]

elf<-data.matrix(t(exprs(gorsmy['Elf5',])))
pd$elf=elf[rownames(pd),1]

cellkeep<-c(rownames(pd)[ano$isElfMouse=="Elf5"],sample(rownames(pd)[pd$isElfMouse=="WT"],length(rownames(pd)[pd$isElfMouse=="Elf5"])
ggplot(pd,aes(x=as.numeric(res.1.2),fill=State)) + geom_bar(position='fill') + facet_wrap(~isElfMouse)
ggsave("plots/180612_normalised_gorsmy_res1.2vsState.pdf")
ggplot(pd[cellkeep,],aes(x=x,y=y,size=elf,col=factor(as.numeric(res.1.2)))) + geom_point() + facet_wrap(~isElfMouse)
ggsave("plots/180612_subsetted_gorsmy_trajectories.pdf")


load("r_objects/180403_epidat_withMonocle.rdata")
epidat@meta.data$gorsmy_dimRed_x=y[rownames(epidat@meta.data),1]
epidat@meta.data$gorsmy_dimRed_y=y[rownames(epidat@meta.data),2]





########original code#######
ordering_genes <- subset(disp_table,mean_expression >= 0.05 &dispersion_empirical >= 1 * dispersion_fit)$gene_id

#Now use these genes to order cells and plot the cell trajectory
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2,auto_param_selection = F)
HSMM <- orderCells(HSMM, reverse=FALSE)
plot_cell_trajectory(HSMM, color_by="res.2")
##############################


## Plot a gene's expression across Pseudo-time
my_genes <- row.names(subset(fData(HSMM),genesymbol %in% c("Pdgfra")))
cds_subset <- HSMM[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="res.2")

### Some other functions

## Plot a gene's expression across multi-branching Pseudo-time
plot_genes_branched_pseudotime(HSMM[c('Pdgfra', 'Pecam1'),], branch_point=21,color_by="cluster",ncol=1)

unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.05 & dispersion_empirical >= 0.5 * dispersion_fit)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

HSMM <- clusterCells(HSMM, num_clusters=7)
plot_cell_trajectory(HSMM, 1, 2, color="cluster")
col.immune <-  c("#99ff66","#00eefd","#aa00ff","#ff8000")
monocle.plot <- plot_cell_trajectory(HSMM, 1, 2, color="cluster")
#if you want to set the color manually
monocle.plot <- plot_cell_trajectory(HSMM, 1, 2, color=col.immune)


