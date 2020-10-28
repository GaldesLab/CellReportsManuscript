### a stepwise analysis of cellphoneDB output


#preprocessing: i had an incorrect pva mod of means so this has fixed that also based loosly on email from teichman lab see previous versions
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
cd("~/share/ScratchGeneral/briglo/cellphone/out/all_top100/")
dat<-preprocCellphone()
sdat<-preprocCellphone(varval=0)
cd("~/Desktop/briglo_backup/scRNA/davgal/paperstuff/")
save(dat,sdat,file="r_objects/181221_cellphoneDB_FIXED.rdata")

### to do plot pairs in ALLs15 where clusters were called
library(reshape2)
library(network)
library(sna)
library(ggplot2)
library(RColorBrewer)
library(ggraph)
library(igraph)
library(pheatmap)
library(cowplot)

#This seems a good framework : STILL NOT SURE ABOUT WHETHER TO USE A VARIANCE CUTOFF MOST PLOTS ARE FOR cutoff....
load("r_objects/181221_cellphoneDB_FIXED.rdata")
x<-sdat$plotdat
target<-unlist(lapply(strsplit(colnames(x),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(x),"_"), function(x) paste0(x[1],"_",x[2])))

#1)test for bidirectional
df<-data.frame(target,source,stringsAsFactors=F)

#make list of genes where forward==reverse (i.e. in epi1_fib1 AND fib1_epi1), SCORE CUTOFF IMPORTANT ALOS HARD CODED WHERE IT DOESNT REPORT BIDIRECTIONAL INTERACTION
z<-vector(mode='list',length=length(df$source))
for (y in c(1:172,174:369,371:length(z))) z[[y]]=intersect(rownames(x)[x[,paste0(df[y,1],"_",df[y,2])]>0],rownames(x)[x[paste0(df[y,2],"_",df[y,1])]>0])

#extract unique pairs 
xr<-sapply(1:length(df$source), function(y) match(paste0(df[y,1],"_",df[y,2]),paste0(df[,2],"_",df[,1])))
k<-sapply(1:length(xr), function(X) ifelse(xr[X]>X,xr[X],NA))

#trim
tdf<-df[!is.na(k),]
tz<-z[!is.na(k)]
names(tz)<-paste0(tdf$source,"_",tdf$target)
ugen<-unique(unlist(tz))

oi<-data.frame(do.call(cbind,lapply(tz,function(x) ugen %in% x)),row.names=ugen)
ss<-data.frame(do.call(rbind,lapply(strsplit(colnames(x),"_"), function(x) c('source'=x[1],'target'=x[3]))),row.names=colnames(x))
tss<-ss[!is.na(k),]
library(pheatmap)
pheatmap(data.matrix(oi),annotation_col=tss[,1:2])


#renetwork
rex<-x[rownames(oi),colnames(oi)]
target<-unlist(lapply(strsplit(colnames(rex),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(rex),"_"), function(x) paste0(x[1],"_",x[2])))

con<-lapply(seq(0,1,.1), function(x) colSums(rex>x))
y<-data.frame(source,target,do.call(cbind,con))
colnames(y)[3:13]<-paste0('countAboveMean',seq(0,1,.1))

p<-lapply(colnames(y)[3:13], function(cutoff){
    xval= "source"
    yval="target"
    return(ggplot(y,aes_string(x=xval,y=yval,fill=cutoff)) + geom_tile() + scale_fill_gradient(low="white",high='blue') + ggtitle(cutoff) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
})
plot_grid(plotlist=p)
ggsave("plots/181221_Heatmaps_numberOfInteractionsByMeanCutoff.pdf")

gr<-graph_from_data_frame(y, directed = F,vertices=an)
deg=degree(gr, mode ="all")
ggraph(gr, layout = 'linear') +    ylim(-1,8) + geom_edge_arc(aes(alpha =countAboveMean0.3>25)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2) + coord_flip()

# this is where i figured i had fixed issue with previous versions, relationships are stable at all levels of cutoffs...

#So get a handle on means, trim out "complex" interactions
toi<-oi[-grep("complex",rownames(oi)),]
trex<-x[rownames(toi),colnames(toi)]
target<-unlist(lapply(strsplit(colnames(trex),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(trex),"_"), function(x) paste0(x[1],"_",x[2])))
ty<-data.frame(source,target,t(trex),id=colnames(trex))
mty<-melt(ty,id.vars=c('source','target','id'))
mty$s<-unlist(lapply(strsplit(as.character(mty$source),"_"), function(x) x[1]))
mty$t<-unlist(lapply(strsplit(as.character(mty$target),"_"), function(x) x[1]))
mty$interaction=paste0(mty$s,"_",mty$t)
ggplot(mty,aes(x=source,y=target,size=value,color=interaction)) + geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + facet_wrap(~variable) + scale_color_brewer(palette="Paired")
ggsave("plots/181221_Dotplots_meansFromCellPhoneNonComplex.pdf")

#redoing graph building, this could probably be wrapped
tcon<-lapply(seq(0,1,.1), function(x) colSums(trex>x))
tty<-data.frame(source,target,do.call(cbind,tcon))
colnames(tty)[3:13]<-paste0('countAboveMean',seq(0,1,.1))
gr<-graph_from_data_frame(tty, directed = F,vertices=an)
deg=degree(gr, mode ="all")

ggraph(gr, layout = 'linear') + ylim(-1.5,8) +     geom_edge_arc(aes(alpha =countAboveMean0.3>12)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2,nudge_y=-1) + coord_flip()
ggsave("forFa_2.pdf")
ggraph(gr, layout = 'kk') +  geom_edge_link(aes(alpha =countAboveMean0.3>8)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2,nudge_x=.1) + coord_flip()
ggsave("forFa.pdf")

# see if reflected in values

#make useful Seurat Object
load('r_objects181214_ALL_reIDforcellphone.rdata') # from when i built the data for cellPhone
load("r_objects/ALLs15.rdata")
load("r_objects/181221_cellphoneDB_FIXED.rdata")
ALLs15@meta.data$cellphoneDB_id=ALL_clusterID[rownames(ALLs15@meta.data),"newID"]

m<-match(ALLs15@meta.data$cellphoneDB_id,an$name)
ALLs15@meta.data$cellphoneDB_designation=an$ID[m]
ALLs15@meta.data$cellphoneDB_designation[is.na(ALLs15@meta.data$cellphoneDB_designation)]="orphan"
ALLs15@meta.data$cellphoneDB_comb<-paste0(ALLs15@meta.data$cellphoneDB_id,"_",ALLs15@meta.data$cellphoneDB_designation)
TSNEPlot(ALLs15,group.by="cellphoneDB_comb",do.label=T, do.return=T,pt.shape="genotype",,pt.size=1.5)
ALLs15<-SetIdent(ALLs15,ident.use=ALLs15@meta.data$cellphoneDB_comb)
TA<-SubsetData(ALLs15,cells.use=rownames(ALLs15@meta.data)[!is.na(m)])

## pull apart gene and plot pairs
sgn<-strsplit(rownames(toi),"_")
df<-as.character(tools::toTitleCase(tolower(unlist(lapply(sgn, function(x) x[1])))))
ip<-as.character(tools::toTitleCase(tolower(unlist(lapply(sgn, function(x) x[2])))))
ddf<-data.frame(g1=df,g2=ip,stringsAsFactors=F)
tddf<-ddf[!(grepl(" ",ddf$g1) | grepl(" ",ddf$g2)),]
tddf<-tddf[tddf$g1 %in% rownames(TA@data) & tddf$g2 %in% rownames(TA@data),]
pp<-apply(tddf,1, function(x) DotPlot(TA,x,do.return=T,x.lab.rot=T))
plot_grid(plotlist=pp[1:20],nrow=2)
ggsave("plots/181221_sig_var_genePairs_inAllcells_1of2.pdf")
plot_grid(plotlist=pp[19:23],nrow=2)
ggsave("plots/181221_sig_var_genePairs_inAllcells_2of2.pdf")


###grouping by pathway gawd i hate this...
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
mappers<-getEntrez(TA,'all') 
map<-mappers[mappers$hgnc_symbol %in% c("TGFBR1",unlist(sgn)),]
m<-match(unlist(sgn),map$hgnc_symbol)

andf<-data.frame(gene=unlist(sgn), entrezgene=map$entrezgene[m],stringsAsFactors=F)

andf$entrezgene[is.na(m)]=7046
andf$gene[is.na(m)]="TGFBR1"
andf$entrezgene<-as.character(andf$entrezgene)
path=bitr_kegg(andf$entrezgene,"ncbi-geneid",c("Path","SYMBOL","hsa")
tandf<-andf[andf$entrezgene %in% path$'ncbi-geneid',]

upath<-unique(path$Path)
library(RCurl)
x<-sapply(upath,function(x) getURL(paste0("http://togows.dbcls.jp/entry/pathway/",x,"/name")))
tx<-data.frame(unlist(lapply(strsplit(x, " - "), function(x) x[1]))

table(path$Path) %>% sort() ->fuck
pathways<-names(fuck)[fuck>1]

tpa<-path[path$Path %in% pathways,]
stpa<-split(tpa$'ncbi-geneid',tpa$Path)
geneLists<-lapply(stpa, function(x) as.character(andf$gene)[andf$entrezgene %in% x])
names(geneLists)<-tx[names(geneLists)]

pathPlots<-lapply(1:length(geneLists), function(x) DotPlot(TA,unique(tools::toTitleCase(tolower(geneLists[[x]]))),do.return=T) + ggtitle(names(geneLists)[x]))

pdf(file='plots/181221_sig_var_genePairs_inAllcell_byPathway.pdf')
pathPlots
dev.off()