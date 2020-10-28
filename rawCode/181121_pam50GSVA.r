library(genefu)
pam50<-apply(pam50.robust$centroids,2, function(x) tools::toTitleCase(tolower(rownames(pam50.robust$centroids)))[x>0])
save(pam50,file="r_objects/181120_pam50offRobustCentroids.rdata")

load("r_objects/tmp_randomizedscRNA.rdata")
 divvec<-rep(1:50, each=400)
sexpmat<-lapply(split(data.frame(pexpMat),divvec[1:9708]), function(x) {
    tmp<-t(x)
    rownames(tmp)<-colnames(pexpMat)
    return(tmp)})

 library(GSVA)
 library(parallel)
 library(Seurat)

sink('gsvaLog_pam50.txt')
scgsva_meta<-lapply(sexpmat, function(x) {
	return(gsva(expr=x,gset.idx.list=pam50,kcdf="Poisson",parallel.sz=4,verbose=T))
	gc()
})
save(scgsva_meta,file='r_objects/181121_SingleCellgsva_pam50_batch.rdata')
sink()

load('r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsva.rdata')
bmat<-data.frame(do.call(rbind,lapply(scgsva_meta,t)))
obmat<-bmat[rownames(epidat@meta.data),]


epidat@meta.data<-cbind(epidat@meta.data,obmat)
plo<-function(x="gorsmy_dimRed_x",y="gorsmy_dimRed_y",contVar){
    #print(grep("_x|_y",colnames(epidat@meta.data),value=T))
    #print(grep("HALLMARK",colnames(epidat@meta.data),value=T))
   require(ggplot2)
ggplot(epidat@meta.data,aes_string(x=x , y=y, color=contVar)) + geom_point(alpha=.5) +facet_wrap(~isElfMouse) + scale_colour_gradient2(low = "blue", high = "red") + theme(legend.position="top",legend.text=element_text(size=5))
}
scorelist<-names(pam50)
save(epidat,plo,scorelist,file='r_objects/181005_epi_ano_monocle_tern_cca_POISSONgsvaHgorsmyPAM50.rdata')


pdf(file='180928_pam50GSVA_Pseudo.pdf',width=12,height=8)
sapply(scorelist, function(i) print(plo(contVar=i)))
dev.off()

