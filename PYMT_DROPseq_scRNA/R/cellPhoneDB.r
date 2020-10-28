# for Seurat object named epidat and fibdat

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

# build objects for downstream
###EPI
load("r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsva.rdata")
epiID<-data.frame(row.names=rownames(epidat@meta.data),res=paste0("epi_",epidat@meta.data$res.0.7))

###FIB
load("/Users/briglo/Dropbox/PYMT scRNAseq paper 2018/180725_fib_ano.rdata")
fibID<-data.frame(row.names=rownames(fibdat@meta.data),res=paste0("fib_",fibdat@meta.data$res.0.1))

clusterdat<-data.frame(rbind(epiID,fibID))
clusterdat$cellID=rownames(clusterdat)
clusterdat$cellType=as.character(clusterdat$res)
save(clusterdat,file="r_objects/181128_IDsForCellPhone.rdata")

#getData
tda<-SubsetData(ALLs15,cells.use=clusterdat$cellID)

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
 

# write out files for cellphoneDB

tda<-da[,colnames(da) %in% all_cells]
tclusterdat<-clusterdat[rownames(clusterdat) %in% all_cells,]
tcd<-data.frame(Cell=tclusterdat[colnames(tda),'cellID'],cell_type=tclusterdat[colnames(tda),'cellType'])

 write.table(data.matrix(tda),file="trimmed_fib_epi_counts.txt",quote=F,sep='\t')
 write.table(tcd,file="trimmed_fib_epi_meta.txt",quote=F,sep='\t',row.names=F)

#cellphoneDB run according to software instructions

 ### a stepwise analysis of cellphoneDB output
cd("PATH/TO/CELLPHONE/OUT")

#preprocessing: also based loosly on email from teichman lab see previous versions
preprocCellphone<-function(pval=0.01,varval=0.05){
message("you should have run this in an cellphone results directory")
message("pval cutoff=",pval)
message("variation value=",varval)

fnam<-c(pval='pvalues.txt',mean="means.txt")
dat<-lapply(fnam,function(x) read.table(x, header=T, stringsAsFactors = F, sep='\t', comment.char = ''))
dat<-lapply(dat,function(x){
rownames(x) = x$interacting_pair
x<-x[,-c(1:9)]
})
tdat<-lapply(dat,function(x) x[rowSums(dat$pval<pval)>0,colSums(dat$pval<pval)>0])
tmp<-tdat$mean
tmp[tdat$pval>pval]=0
tdat$plotdat<-tmp[matrixStats::rowVars(data.matrix(tmp))>varval,]
#print(pheatmap::pheatmap(data.matrix(tdat$plodat))) #this needs a big display object
return(invisible(tdat))
}
dat<-preprocCellphone()
sdat<-preprocCellphone(varval=0)


#plotting

library(tidyverse)
an=read.table("cellphoneDB/190107_cellAno.txt",header=T,stringsAsFactors=F) #cluster/celltype anntatations
rex<-sdat$mean
target<-unlist(lapply(strsplit(colnames(rex),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(rex),"_"), function(x) paste0(x[1],"_",x[2])))

con<-lapply(seq(0,1,.1), function(x) colSums(rex>x))
y<-data.frame(source,target,do.call(cbind,con))
colnames(y)[3:13]<-paste0('countAboveMean',seq(0,1,.1))
y$targetename<-an$ID[match(y$target,an$name)]
y$sourcename<-an$ID[match(y$source,an$name)]

ggplot(y,aes(x=paste(source,sourcename),y=paste(target,targetename),fill=countAboveMean0.3)) + geom_tile() + scale_fill_gradient(low="white",high='blue') + ggtitle('countAboveMean0.3') + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


#sanityCheck
pdat<-y[,c(1,2,6,14,15)]
uid<-unique(c(pdat$sourcename,pdat$targetename))
x<-matrix(NA,length(uid),length(uid))
rownames(x) <-uid
colnames(x) <-uid
 for (i in uid) for (j in uid) { x[i,j]=ifelse (sum(pdat$sourcename==i & pdat$targetename==j)==1,pdat$countAboveMean0.3[which(pdat$sourcename==i & pdat$targetename==j)],NA)}

#props
corrplot::corrplot(x,is.corr=F,type='lower')
dev.copy2pdf(file="200311_allIntRenamed_pheat.pdf")
$raw
pheatmap::pheatmap(t(x),cluster_rows=F,cluster_cols=F)


gr<-graph_from_data_frame(y, directed = F,vertices=an)
deg=degree(gr, mode ="all")
intgraph<-function(cutoff=20){
   require(ggraph)
   print(ggraph(gr, layout = 'linear') + ylim(-1.5,8) +     geom_edge_arc(aes(alpha =countAboveMean0.3>cutoff)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2,nudge_y=-1) + coord_flip())}



#int graph function
cellphoneDB_data<-y
cluster_anno<-an

intgraph<-function(scoreCut=0.3, numberCut=0, numberSplit=35){
   require(ggraph)
   require(dplyr)
   require(igraph)
   require(cowplot)
   message("use like this: intgraph(scoreCut=0.3, numberCut=0, numberSplit=35)")
 tmp<-cellphoneDB_data[cellphoneDB_data[,paste0("countAboveMean",scoreCut)]>numberCut,c('source','target',paste0("countAboveMean",scoreCut))]
 colnames(tmp)[3]<-"score"
 gr<-graph_from_data_frame(tmp,directed = F,vertices=cluster_anno)
 deg=degree(gr, mode ="all")
 colname<-paste0("countAboveMean",numberCut)
   p1<-ggraph(gr, layout = 'linear') + ylim(-1.5,8) +     geom_edge_arc(aes(width =score,alpha=score>numberSplit)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=count)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2,nudge_y=-1) + coord_flip() + ggtitle(paste0("scoreCut=",scoreCut, ", numberCut=",numberCut,", numberSplit=",numberSplit))
   p2<-ggplot(tmp,aes(x=score)) + geom_bar()
   ggdraw() + draw_plot(p1 + theme(legend.justification="bottom"), 0,0,1,1) + draw_plot(p2 + geom_vline(xintercept=numberSplit) ,0.75, 0.7, 0.2, 0.3)
   }
intgraph(scoreCut=.5, numberCut=10, numberSplit=20)