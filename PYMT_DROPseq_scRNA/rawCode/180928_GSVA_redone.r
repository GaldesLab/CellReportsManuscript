#########SINGLE CELL###########
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
load("../../extraData/GSE63310_RAW/180227_gorsmy_sigGenes.rdata")
 ngs<-c(gene.set.H,signatureGenes)

save(ngs,file="r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")

library(GSVA)
library(parallel)
library(Seurat)
load("r_objects/180727_epi_ano_monocle_tern_cca.rdata")
load("r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")
expMat<-as.matrix(FetchData(epidat,vars.all=epidat@var.genes))
pexpMat<-expMat[sample(1:length(rownames(expMat)),length(rownames(expMat)),replace=F),] #randomised in case it uses dispersion to estimate
save(pexpMat, file='tmp_randomizedscRNA.rdata')
divvec<-rep(1:50, each=300)
sexpmat<-lapply(split(pexpMat,divvec[1:9708]),t)


library(GSVA)
library(parallel)
library(Seurat)
load("r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")
load('tmp_randomizedscRNA.rdata')
divvec<-rep(1:50, each=400)
sexpmat<-lapply(split(data.frame(pexpMat),divvec[1:9708]), function(x) {
	tmp<-t(x)
	rownames(tmp)<-colnames(pexpMat)
	return(tmp)})
#warn(-1)
sink('gsvaLog_pam50.txt')
scgsva_meta<-lapply(sexpmat, function(x) {
	return(gsva(expr=x,gset.idx.list=pam50,kcdf="Poisson",parallel.sz=4,verbose=T))
	gc()
})
save(scgsva_meta,file='r_objects/181004_SingleCellgsva_pam50_batch.rdata')
sink()

bmat<-data.frame(do.call(rbind,lapply(scgsva_meta,t)))
obmat<-bmat[rownames(epidat@meta.data),]


epidat@meta.data<-cbind(epidat@meta.data,obmat)
plo<-function(x="gorsmy_dimRed_x",y="gorsmy_dimRed_y",contVar){
    #print(grep("_x|_y",colnames(epidat@meta.data),value=T))
    #print(grep("HALLMARK",colnames(epidat@meta.data),value=T))
   require(ggplot2)
ggplot(epidat@meta.data,aes_string(x=x , y=y, color=contVar)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "blue", high = "red") + theme(legend.position="top",legend.text=element_text(size=5))
}
scorelist<-names(ngs)
save(epidat,plo,scorelist,file='r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsva.rdata')

pdf(file='180928_hallmarkgorsmyGSVA_Pseudo.pdf',width=12,height=8)
sapply(scorelist, function(i) print(plo(contVar=i)))
dev.off()

pdf(file='180928_hallmarkgorsmyGSVA_TSNE.pdf',width=12,height=8)
sapply(scorelist, function(i) print(plo(x="seurat_tsne_x",y="seurat_tsne_y",contVar=i)))
dev.off()
####################SINGLE CELL COMPLETE######################
library(GSVA)
library(parallel)
library(Seurat)
load("r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")
load('r_objects/180928_epi_ano_monocle_tern_cca_gsva.rdata')
epidat<-SetIdent(epidat,ident.use=paste0(epidat@meta.data$res.0.7,'_',epidat@meta.data$isElfMouse))
expMat<-as.matrix(AverageExpression(epidat,genes.use=epidat@var.genes))
cluster_genoGSVA<-gsva(expMat,ngs,parallel.sz=4,verbose=T)
pheatmap::pheatmap(cluster_genoGSVA[51:53,],main="GSVA by genotype and res.0.7",cluster_cols=F)#,annotation_col=log(an,2))
dev.copy2pdf(file="plots/181003_epiGSVA_geno_res.0.7.pdf")



epidat<-SetIdent(epidat,ident.use=epidat@meta.data$res.0.7)
expMat<-as.matrix(AverageExpression(epidat,genes.use=epidat@var.genes))
clusterGSVA<-gsva(expMat,ngs,parallel.sz=4,verbose=T)
pheatmap::pheatmap(clusterGSVA,main="GSVA by res.0.7",cluster_cols=F)