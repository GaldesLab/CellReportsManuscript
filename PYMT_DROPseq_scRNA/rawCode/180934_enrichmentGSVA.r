mkListToTERMTOGENE<-function(obj){
x<-data.frame(do.call(rbind,lapply(1:length(obj), function(x){
	z<-data.frame(obj[[x]])
	z$term<-names(obj)[x]
	return(z)
	}
	)))
return(data.frame(ont=factor(x[,2]),gene=factor(x[,1])))
	}


load("/Users/briglo/Desktop/briglo_backup/scRNA/davgal/paperstuff/r_objects/mouse_TERM2GENE_c2c6hallmark.rdata")

library(clusterProfiler)
doCluster<-function(markerList,term2geneObj){
	return(GSEA(markerList,TERM2GENE=term2geneObj))
}

load("r_objects/180502_epires1.2posmarkers.rdata")
mouse<- useMart(biomart='ensembl', dataset = "mmusculus_gene_ensembl")
	 getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="mgi_symbol",values=seuratobj@var.genes,mart=mouse))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="mgi_symbol",values=rownames(seuratobj@data),mart=mouse))
	 }
	 }
hasEntrez<-getEntrez(epidat,geneList='all')

 mkEnt<-function(listOids){
	     require(clusterProfiler)
	     require(ReactomePA)
	     return(lapply(listOids,function(x) as.character(na.omit(hasEntrez$entrezgene[match(x,hasEntrez$mgi_symbol)]))))
	  } 
tmark %>% arrange(cluster,-avg_logFC) -> sortedTmark
splMar<-function(marker) lapply(split(marker,marker$avg_logFC>0), function(x) mkEnt(split(x$gene,x$cluster))) # splits a FindAllMarkers object into cluster and up/down


posMarkerEntrez<-mkEnt(smark)

C2_GSEA<-lapply(posMarkerEntrez,function(x) doCluster(x,mouse_C2_TERM2GENE))




expMat<-AverageExpression(epidat)
load("/Users/briglo/Desktop/briglo_backup/scRNA/davgal/paperstuff/r_objects/mouse_TERM2GENE_c2c6hallmark.rdata")

library(biomaRt)
mouse<- useMart(biomart='ensembl', dataset = "mmusculus_gene_ensembl")
	 getEntrez<-function(seuratobj,geneList='var.genes'){
	     require(biomaRt)
	     require(Seurat)
	 if (geneList=="var.genes") {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="mgi_symbol",values=seuratobj@var.genes,mart=mouse))
	 } else {
	 	return(getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="mgi_symbol",values=rownames(seuratobj@data),mart=mouse))
	 }
	 }
hasEntrez<-getEntrez(epidat,geneList='all')
the<-hasEntrez[!is.na(hasEntrez$entrezgene),]
uent<-unique(the$entrezgene)
id<-the[match(uent,the$entrezgene),]
x<-data.frame(do.call(rbind,lapply(split(id$entrezgene,id$mgi_symbol), function(x) x[1])))

library(GSVA)
library(parallel)

doGVSA<-function(markerList,term2geneObj){
	return(gsva(markerList,TERM2GENE=term2geneObj))
}
expMat<-AverageExpression(epidat)



library(biomaRt)
mgi_symbol_match<-getBM(attributes =c("mgi_symbol","entrezgene"), values =  list("*"), mart = mouse)
str(mgi_symbol_match)

mgi_symbol_match_filt<-mgi_symbol_match[!is.na(mgi_symbol_match$entrezgene),] #21420
save(mgi_symbol_match, file="../../extraData/EnrichmentLists/mgi_symbol_match.rdata")




library(GSVA)
library(parallel)
library(Seurat)
load("../../extraData/EnrichmentLists/mgi_symbol_match.rdata")
load("../../extraData/EnrichmentLists/mouse_H_v5p2.rdata")
load("r_objects/180727_epi_ano_monocle_tern_cca.rdata")
gene.set.H <- lapply(Mm.H, function(x){
  # gets the index in mgi_symbol_match
  index <- match(x, mgi_symbol_match$entrezgene)
  # returns the matching gene symbols - may be shorter than
  # original gene ID list if match not found
  unique(mgi_symbol_match[index[which(!is.na(index))],]$mgi_symbol)
})

gene.set.H <- sapply(names(gene.set.H), function(x) unique(gene.set.H[[x]]))
###works fine for GSVA


dsEdat<-SubsetData(epidat,max.cells.per.ident=1000,ident.use=epidat@meta.data$isElfMouse)
dat<-
ann<-data.frame(FetchData(epidat,vars.all=c("res.0.7","isElfMouse")),stringsAsFactors=F)

##########This bit used 21 sept 2018################
expMat<-as.matrix(t(FetchData(epidat,vars.all=epidat@var.genes,use.scaled=T)))
pexpMat<-expMat[,sample(1:length(colnames(expMat)),length(colnames(expMat)),replace=F)] #randomised in case it uses dispersion to estimate

sexpmat<-list(
	a=pexpMat[,1:1000],
	b=pexpMat[,1001:2000],
	c=pexpMat[,2001:3000],
	d=pexpMat[,3001:4000],
	e=pexpMat[,4001:5000],
	f=pexpMat[,5001:6000],
	g=pexpMat[,6001:7000],
	h=pexpMat[,7001:8000],
	i=pexpMat[,8001:9000],
	j=pexpMat[,9001:9708]
)

ssgsva_2<-lapply(sexpmat, function(x) return(gsva(x,gene.set.H,parallel.sz=1)))
save(ssgsva_2,file='r_objects/180928_SingleCellgsva_batch.rdata')


load('r_objects/180921_SingleCellgsva_batch.rdata')
bmat<-data.frame(do.call(rbind,lapply(ssgsva,t)))
obmat<-bmat[rownames(epidat@meta.data),]

bmat2<-data.frame(do.call(rbind,lapply(ssgsva_2,t)))
obmat2<-bmat2[rownames(epidat@meta.data),]


epidat@meta.data<-cbind(epidat@meta.data,obmat)
save(epidat,file='r_objects/180921_epi_ano_monocle_tern_cca_gsva.rdata')

library(ggplot2)
colnames(epidat@meta.data)

ggplot(epidat@meta.data,aes(x=gorsmy_dimRed_x , y=gorsmy_dimRed_y , color=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "wheat2", high = "lightpink3") + theme(legend.position="top")
plo<-function(x="gorsmy_dimRed_x",y="gorsmy_dimRed_y",contVar){
ggplot(epidat@meta.data,aes_string(x=x , y=y, color=contVar)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "blue", high = "red") + theme(legend.position="top",legend.text=element_text(size=5))
}
x<-grep("HALLMARK",colnames(epidat@meta.data),value=T)
pdf(file='180921_hallmarkGSVA_Pseudo.pdf',width=12,height=8)
sapply(x, function(i) print(plo(contVar=i)))
dev.off()

x<-grep("HALLMARK",colnames(epidat@meta.data),value=T)
pheatmap::pheatmap(epidat@meta.data[,x],kmeans_k=100)
dev.copy2pdf(file="180924_100kmeans_heatmap.pdf")
a<-pheatmap(epidat@meta.data[,x],kmeans_k=9)
dev.copy2pdf("180924_6kmeans_heatmap.pdf")


kClusters<-data.frame(a$kmeans$cluster)
epidat@meta.data$HallMark_6kmeans_clusters<-kClusters[rownames(epidat@meta.data),'a.kmeans.cluster']
p1<-ggplot(epidat@meta.data,aes(x=as.factor(as.numeric(as.character(res.0.7))),y=HallMark_6kmeans_clusters,fill=factor(HallMark_6kmeans_clusters))) + geom_bar(stat='identity',position='fill') + facet_wrap(~isElfMouse)
p2<-ggplot(epidat@meta.data,aes(x=isElfMouse,y=HallMark_6kmeans_clusters,fill=factor(HallMark_6kmeans_clusters))) + geom_bar(stat='identity',position='fill')
cowplot::plot_grid(p1,p2,nrow=2)
ggsave(file="190824_k6_breakdown.pdf")

 ################################################

singleCellGSVA<-gsva(expMat,gene.set.H,parallel.sz=4,verbose=T)

pheatmap::pheatmap(singleCellGSVA,main="GSVA by singleCells",show_colnames=F,annotation_col=data.frame(ann$genotype))
dev.copy2pdf(file="plots/180914_epiGSVA_singleCell.pdf")

epidat<-SetIdent(epidat,ident.use=epidat@meta.data$res.0.7)
expMat<-as.matrix(AverageExpression(epidat,genes.use=epidat@var.genes))
clusterGSVA<-gsva(expMat,gene.set.H,parallel.sz=4,verbose=T)
pheatmap::pheatmap(clusterGSVA,main="GSVA by res.0.7",,cluster_cols=F)
dev.copy2pdf(file="plots/180914_epiGSVA_res.0.7.pdf")

epidat<-SetIdent(epidat,ident.use=paste0(epidat@meta.data$res.0.7,'_',epidat@meta.data$isElfMouse))
expMat<-as.matrix(AverageExpression(epidat,genes.use=epidat@var.genes))
cluster_genoGSVA<-gsva(expMat,gene.set.H,parallel.sz=4,verbose=T)
pheatmap::pheatmap(cluster_genoGSVA,main="GSVA by genotype and res.0.7",cluster_cols=F)#,annotation_col=log(an,2))
dev.copy2pdf(file="plots/180914_epiGSVA_geno_res.0.7.pdf")

load("../../extraData/GSE63310_RAW/180227_gorsmy_sigGenes.rdata")
 ngs<-c(gene.set.H[c(1,2,4,5:10,13,14,18,19,20,22,24,25,30,31,34,35,36,37,40,48,49)],signatureGenes)

 cluster_geno_modGSVA<-gsva(expMat,ngs,parallel.sz=4,verbose=T)
 pheatmap::pheatmap(cluster_geno_modGSVA,main="GSVA mod with gorsmy by genotype and res.0.7",cluster_cols=F,annotation_col=log(an,2))
 dev.copy2pdf(file="plots/180914_epi_modGSVA_geno_res.0.7.pdf")


 x=data.frame(rdat[order(paste0(rdat[,c(51)],rdat[,c(52)])),-c(51:52)])
 y=data.frame(rdat[order(paste0(rdat[,c(51)],rdat[,c(52)])),c(51:52)])


x<-t(rdat[order(paste0(rdat[,c(51)],rdat[,c(52)])),-c(51,52)])
 pheatmap(x,annotation_col=y,show_colnames=F,cluster_cols=F)

ggplot(epidat@meta.data,aes(x=HALLMARK_G2M_CHECKPOINT,y=HALLMARK_MITOTIC_SPINDLE, color=isElfMouse)) + geom_point() + facet_wrap(~as.factor(as.numeric(res.0.7)))

 HALLMARK_G2M_CHECKPOINT


x<-grep("HALLMARK",colnames(epidat@meta.data),value=T)
dat<-epidat@meta.data[,x]
ano<-epidat@meta.data[,c('isElfMouse','res.0.7','gorsmy_State')]
sortv<-paste0(ano$res.0.7,"_",ano$isElfMouse)
sortdat<-dat[order(sortv),]
sortano<-ano[order(sortv),]
pheatmap::pheatmap(t(data.matrix(sortdat)),cluster_cols=F,annotation_col=sortano,show_colnames=F)
dev.copy2pdf(file="180927_organisedHeatmap.pdf")

sid<-split(dat,ano$res.0.7)
sano<-split(ano,ano$res.0.7)
for(i in 1:length(sid)) sid[[i]]<-sid[[i]][order(as.character(sano[[i]]$isElfMouse)),]
for(i in 1:length(sid)) sano[[i]]<-sano[[i]][order(as.character(sano[[i]]$isElfMouse)),]

pdf('erm_2.pdf')
for(i in 1:length(sid)) pheatmap::pheatmap(t(sid[[i]]),show_colnames=F,annotation_col=sano[[i]])
dev.off()