#Fetch gene lists
library(GSVA)
library(parallel)
library(Seurat)
load("../../extraData/EnrichmentLists/mgi_symbol_match.rdata") # ID matching object mgi_symbol_match
load("../../extraData/EnrichmentLists/mouse_H_v5p2.rdata") #TERM2GENE from mouse msgidb remapping Mm.H
load("r_objects/epi_ano_monocle_tern_cca.rdata") # Seurat object epidat
gene.set.H <- lapply(Mm.H, function(x){
  # gets the index in mgi_symbol_match
  index <- match(x, mgi_symbol_match$entrezgene)
  # returns the matching gene symbols - may be shorter than
  # original gene ID list if match not found
  unique(mgi_symbol_match[index[which(!is.na(index))],]$mgi_symbol)
})

gene.set.H <- sapply(names(gene.set.H), function(x) unique(gene.set.H[[x]]))

#randomised in case GSVA uses dispersion to estimate
# done twice to check for effect, none observed 
 expMat<-as.matrix(FetchData(epidat,vars.all=epidat@var.genes))
pexpMat<-expMat[sample(1:length(rownames(expMat)),length(rownames(expMat)),replace=F),] 
# split for speed
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

ssgsva<-lapply(sexpmat, function(x) return(gsva(x,gene.set.H,parallel.sz=1)))

#annotate epidat
bmat<-data.frame(do.call(rbind,lapply(ssgsva,t)))
obmat<-bmat[rownames(epidat@meta.data),]


epidat@meta.data<-cbind(epidat@meta.data,obmat)


#plotting
scorelist<-names(ngs)
plo<-function(x="gorsmy_dimRed_x",y="gorsmy_dimRed_y",contVar){
    #print(grep("_x|_y",colnames(epidat@meta.data),value=T))
    #print(grep("HALLMARK",colnames(epidat@meta.data),value=T))
   require(ggplot2)
ggplot(epidat@meta.data,aes_string(x=x , y=y, color=contVar)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "blue", high = "red") + theme(legend.position="top",legend.text=element_text(size=5))
}

sapply(scorelist, function(i) print(plo(x="seurat_tsne_x",y="seurat_tsne_y",contVar=i)))


## bulk is easier
library(GSVA)
library(parallel)
library(Seurat)
load("r_objects/180928_HALLMARKandGORSMYgeneListsforGSVA.rdata")
load('r_objects/180928_epi_ano_monocle_tern_cca_gsva.rdata')
epidat<-SetIdent(epidat,ident.use=paste0(epidat@meta.data$res.0.7,'_',epidat@meta.data$isElfMouse))
expMat<-as.matrix(AverageExpression(epidat,genes.use=epidat@var.genes))
cluster_genoGSVA<-gsva(expMat,ngs,parallel.sz=4,verbose=T)

pheatmap::pheatmap(cluster_genoGSVA[51:53,],main="GSVA by genotype and res.0.7",cluster_cols=F)