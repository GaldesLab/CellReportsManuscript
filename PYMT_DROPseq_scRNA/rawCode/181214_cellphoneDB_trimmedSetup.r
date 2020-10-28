library(Seurat)
load('r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsva.rdata')

#get 100 cells with highest number of genes from each cluster
epidat<-SetIdent(epidat,ident.use=paste0("epi_",epidat@meta.data$res.0.7))
x<-epidat@meta.data[,c('res.0.7','nGene')]
sx<-split(x,x$res.0.7)
osx<-lapply(sx, function(x) x[order(x$nGene,decreasing=T),])
cells.use_epi<-as.character(unlist(lapply(osx, function(x) head(rownames(x),100))))
ids_epi<-as.character(unlist(lapply(osx, function(x) head(x$res.0.7,100))))

load("r_objects/180725_fib_ano.rdata")
fibdat<-SetIdent(fibdat,ident.use=paste0("fib_",fibdat@meta.data$res.0.1))
x<-fibdat@meta.data[,c('res.0.1','nGene')]
sx<-split(x,x$res.0.1)
osx<-lapply(sx, function(x) x[order(x$nGene,decreasing=T),])
cells.use_fib<-as.character(unlist(lapply(osx, function(x) head(rownames(x),100))))
ids_fib<-as.character(unlist(lapply(osx, function(x) head(x$res.0.7,100))))

all_cells<-c(cells.use_fib,cells.use_epi)
load("r_objects/181129_cellPhoneIntermediates.rdata"))

tda<-da[,colnames(da) %in% all_cells]
tclusterdat<-clusterdat[rownames(clusterdat) %in% all_cells,]
tcd<-data.frame(Cell=tclusterdat[colnames(tda),'cellID'],cell_type=tclusterdat[colnames(tda),'cellType'])

 write.table(data.matrix(tda),file="trimmed_fib_epi_counts.txt",quote=F,sep='\t')
 write.table(tcd,file="trimmed_fib_epi_meta.txt",quote=F,sep='\t',row.names=F)
####### run on cluster analysed with 181214_cellPhoneDBanalysis.r


##NEXT

load("r_objects/ALLs15.rdata")
load("~/Desktop/briglo_backup/scRNA/davgal/paperstuff/r_objects/181128_IDsForCellPhone.rdata")
load("r_objects/181129_cellPhoneIntermediates.rdata"))

ALL_clusterID <-FetchData(ALLs15,c('res.1','nGene'))
ALL_clusterID$res.1<-paste0('orig_',ALL_clusterID$res.1)
ALL_clusterID$newID<-sapply(rownames(ALL_clusterID), function(y) ifelse(y %in% rownames(clusterdat), clusterdat[y,"cellType"],x[y,'res.1']))
save(ALL_clusterID,file='r_objects181214_ALL_reIDforcellphone.rdata')

sx<-split(ALL_clusterID,ALL_clusterID$newID)
ssx<-sx[unlist(lapply(sx,function(X) dim(X)[1])>100)]
osx<-lapply(ssx, function(x) x[order(x$nGene,decreasing=T),])
cells.use<-as.character(unlist(lapply(osx, function(x) head(rownames(x),100))))

tda<-SubsetData(ALLs15,cells.use=cells.use)

getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
		 human<- useMart(biomart='ensembl', dataset = "hsapiens_gene_ensembl")
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=toupper(seuratobj@var.genes),mart=human))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=toupper(rownames(seuratobj@data)),mart=human))
	 }
	 }#p
mappers<-getEntrez(tda,'all') #Elf5 isnt a var.gene(WTF)

 da<-tda@data[toupper(rownames(tda@data)) %in% mappers$hgnc_symbol[!is.na(mappers$ensembl_gene_id)],]
 m<-match(toupper(rownames(da)),mappers$hgnc_symbol)
 rownames(da)<-mappers$ensembl_gene_id[m]
 save(tda,da,mappers,clusterdat,file="r_objects/181214_cellPhoneIntermediates.rdata")

t_da<-da[,colnames(da) %in% cells.use]
tclusterdat<-ALL_clusterID[rownames(ALL_clusterID) %in% cells.use,]
tclusterdat$cellID<-rownames(tclusterdat)
tcd<-data.frame(Cell=tclusterdat[colnames(t_da),'cellID'],cell_type=tclusterdat[colnames(t_da),'newID'])

write.table(data.matrix(t_da),file="trimmed_allish_counts.txt",quote=F,sep='\t')
 write.table(tcd,file="trimmed_allish_meta.txt",quote=F,sep='\t',row.names=F)

#module load briglo/miniconda/3
#source activate cellphone
 #cellphonedb method statistical_analysis trimmed_allish_meta.txt trimmed_allish_counts.txt --project-name all_top100 --threads 15
 ####### run on cluster analysed with 181214_cellPhoneDBanalysis.r




 ##NEXT!!!
library(Seurat)
load("r_objects/180725_fib_ano.rdata")
load('r_objects/181214_ALL_reIDforcellphone.rdata')
old<-FetchData(fibdat,"res.1")
ALL_clusterID$newFibID<-ALL_clusterID$newID
old$res.1<-paste0("fib_",old$res.1)
ALL_clusterID[rownames(old),"newFibID"]=old$res.1
save(ALL_clusterID,file="r_objects/190107_ALL_id_forCellPhone_newFIB.rdata")
rm(fibdat)

nuke=c('fib_6','fib_7','fib_8','orig_4')
tALL_clusterID<-ALL_clusterID[!ALL_clusterID$newFibID %in% nuke,]

sx<-split(tALL_clusterID,tALL_clusterID$newFibID)
ssx<-sx[unlist(lapply(sx,function(X) dim(X)[1])>100)]
osx<-lapply(ssx, function(x) x[order(x$nGene,decreasing=T),])
cells.use<-as.character(unlist(lapply(osx, function(x) head(rownames(x),100))))

load("r_objects/ALLs15.rdata")
tda<-SubsetData(ALLs15,cells.use=cells.use)

getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
		 human<- useMart(biomart='ensembl', dataset = "hsapiens_gene_ensembl")
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=toupper(seuratobj@var.genes),mart=human))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","hgnc_symbol"),filters="hgnc_symbol",values=toupper(rownames(seuratobj@data)),mart=human))
	 }
	 }#p
mappers<-getEntrez(tda,'all') #Elf5 isnt a var.gene(WTF)

 da<-tda@raw.data[toupper(rownames(tda@raw.data)) %in% mappers$hgnc_symbol[!is.na(mappers$ensembl_gene_id)],cells.use]
 m<-match(toupper(rownames(da)),mappers$hgnc_symbol)
 rownames(da)<-mappers$ensembl_gene_id[m]
 save(tda,da,mappers,tALL_clusterID,file="r_objects/190107_cellPhoneIntermediates.rdata")

t_da<-da[,colnames(da) %in% cells.use]
tclusterdat<-tALL_clusterID[rownames(tALL_clusterID) %in% cells.use,]
tclusterdat$cellID<-rownames(tclusterdat)
tcd<-data.frame(Cell=tclusterdat[colnames(t_da),'cellID'],cell_type=tclusterdat[colnames(t_da),'newFibID'])

write.table(data.matrix(t_da),file="190107_trimmed_allish_counts.txt",quote=F,sep='\t')
 write.table(tcd,file="190107_trimmed_allish_meta.txt",quote=F,sep='\t',row.names=F)


####in terminal
#send to cluster
# mv 190107_trimmed_allish_* ~/share/ScratchGeneral/briglo/cellphone/

# #on cluster
# qrsh -pe smp 16
# cd /share/ScratchGeneral/briglo/cellphone
# mkdir out/190107_newFib_allTop100
# module load briglo/miniconda/3
# source activate cellphone
#  cellphonedb method statistical_analysis 190107_trimmed_allish_meta.txt 190107_trimmed_allish_counts.txt --project-name 190107_newFib_allTop100 --threads 15
 ###### run on cluster analysed with 181214_cellPhoneDBanalysis.r
