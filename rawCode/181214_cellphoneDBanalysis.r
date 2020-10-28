#original from  Mirjana Efremova <me5@sanger.ac.uk>
#Thu, 13 Dec, 21:16
#modified by Brian Gloss

# plot.data is a data frame with columns: colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean'). Here is how I produced it: (it was so mistake riddled i rewrote)

wall_pval = read.table('pvalues.txt', header=T, stringsAsFactors = F, sep='\t', comment.char = '')
all_means = read.table('means.txt', header=T, stringsAsFactors = F, sep='\t', comment.char = '')

rownames(all_pval) = all_pval$interacting_pair
rownames(all_means) = all_means$interacting_pair

all_pval = all_pval[,-c(1:9)]
all_means = all_means[,-c(1:9)]

selected_rows = rowSums(all_pval<0.01)>0
selected_columns = colSums(all_pval<0.01)>0

sel_pval = all_pval[selected_rows, selected_columns]
sel_means = all_means[selected_rows, selected_columns]

plot_data<-sel_means
plot_data[sel_pval>0.01]=0

plo_cols<-grepl("fib",colnames(sel_means))
plo_rows<-matrixStats::rowVars(data.matrix(plot_data))>.05

pd<-plot_data[plo_rows,plo_cols]
pheatmap::pheatmap(pd)
ggsave("181214_cellphonedb_fib_epi_means_heatmap.pdf")

####making a more generic plot including cells not used in cellphone analysis
load("~/Desktop/briglo_backup/scRNA/davgal/paperstuff/r_objects/181129_cellPhoneIntermediates.rdata")
load("~/Desktop/briglo_backup/scRNA/davgal/paperstuff/r_objects/181128_IDsForCellPhone.rdata")
tda@meta.data$cellphoneDB_id<-clusterdat[rownames(tda@meta.data),'cellType']
tda<-SetIdent(tda,ident.use=tda@meta.data$cellphoneDB_id)

x<-strsplit(rownames(pd),"_")
df<-as.character(tools::toTitleCase(tolower(unlist(lapply(x, function(x) x[1])))))
ip<-as.character(tools::toTitleCase(tolower(unlist(lapply(x, function(x) x[2])))))
ddf<-data.frame(g1=df,g2=ip,stringsAsFactors=F)
tddf<-ddf[!(grepl(" ",ddf$g1) | grepl(" ",ddf$g2)),]
pp<-apply(tddf,1, function(x) DotPlot(tda,x,do.return=T))
plot_grid(plotlist=pp,nrow=2)
ggsave("181214_sig_var_genePairs_inAllEpiandFib.pdf")




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
#print(pheatmap::pheatmap(data.matrix(tdat$plodat))) #not working dunno why, image works
return(invisible(tdat))
}
data<-preprocCellphone()
sdat<-preprocCellphone(varval=0)

### to do plot pairs in ALLs15 where clusters were called
library(reshape2)
library(network)
library(sna)
library(ggplot2)
library(RColorBrewer)
library(ggraph)
library(igraph)
library(pheatmap)

dat<-preprocCellphone()
save(dat,an, file="181219_cellphoneDB_allCellTypes.rdata")
x<-dat$plotdat

tx<-data.frame(t(x))
tx$id<-rownames(tx)
mtx<-melt(tx)
mtx$from<-unlist(lapply(strsplit(mtx$id,"_"), function(x) paste0(x[1],"_",x[2])))
mtx$to<-unlist(lapply(strsplit(mtx$id,"_"), function(x) paste0(x[3],"_",x[4])))



x<-dat$plotdat
xx<-x[rowSums(x>0)<150,]

target<-unlist(lapply(strsplit(colnames(xx),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(xx),"_"), function(x) paste0(x[1],"_",x[2])))

con<-lapply(seq(0,1,.1), function(x) colSums(xx>x))
y<-data.frame(source,target,con[[1]],con[[5]])
ty<-y[y$con>0,]

gr<-graph_from_data_frame(ty, directed = T,vertices=an)
deg=degree(gr, mode="all")
ggraph(gr, layout = 'kk') +     geom_edge_link(aes(color= con>10)) + geom_node_point(aes(color=ID,shape=type,size=deg)) + scale_edge_color_manual(values=c("blue","red"))

ggraph(gr, layout = 'kk') +     geom_edge_link(alpha = 0.2) + geom_node_point(aes(color=type,size=deg))

ggraph(gr, layout = 'linear') +
    geom_edge_arc(aes(alpha = ..index..,color=con)) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction') +
scale_edge_colour_gradientn(colours = jet.colors(7)) +
geom_node_point(aes(color=ID,shape=type,size=deg)) 

ggraph(gr, layout = 'linear') +
    geom_edge_arc(aes(alpha = ..index..)) + 
    scale_edge_alpha('Edge direction', guide = 'edge_direction')


#tring another way
load("cellphoneDB_results.rdata")
x<-dat$plotdat
target<-unlist(lapply(strsplit(colnames(x),"_"), function(x) paste0(x[3],"_",x[4])))
source<-unlist(lapply(strsplit(colnames(x),"_"), function(x) paste0(x[1],"_",x[2])))

#1)test for bidirectional
df<-data.frame(target,source,stringsAsFactors=F)

#make list of genes where forward==reverse (i.e. in epi1_fib1 AND fib1_epi1), SCORE CUTOFF IMPORTANT ALOS HARD CODED WHERE IT FUCKED UP
z<-vector(mode='list',length=length(df$source))
for (y in c(1:172,174:369,371:length(z))) z[[y]]=intersect(rownames(x)[x[,paste0(df[y,1],"_",df[y,2])]>0],rownames(x)[x[paste0(df[y,2],"_",df[y,1])]>0])

#extract unique pairs 
xr<-sapply(1:length(df$source), function(y) match(paste0(df[y,1],"_",df[y,2]),paste0(df[,2],"_",df[,1])))
k<-sapply(1:length(xr), function(X) ifelse(xr[X]<X,xr[X],NA))

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
ggsave("plots/181219_Heatmaps_numberOfInteractionsByMeanCutoff.pdf")

gr<-graph_from_data_frame(ty, directed = F,vertices=an)
deg=degree(gr, mode ="all")
ggraph(gr, layout = 'linear') +
    ylim(-1,8) +
    geom_edge_arc(aes(alpha =countAboveMean0.3>50)) + 
    #scale_edge_alpha('Edge direction', guide = 'edge_direction') +
geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) +
geom_node_label(aes(label=paste(name,ID)),label.size=.2) + coord_flip()


#low mean cutoff= more interactions between epi1 and things but higher mean cuttoff more interactions between epi0 and fibs, needss more investigation abouyt numbers and means

#So get a handle on means, trim out complex interactions
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
ggsave("plots/181219_Dotplots_meansFromCellPhoneNonComplex.pdf")

tcon<-lapply(seq(0,1,.1), function(x) colSums(trex>x))
tty<-data.frame(source,target,do.call(cbind,tcon))
colnames(tty)[3:13]<-paste0('countAboveMean',seq(0,1,.1))
gr<-graph_from_data_frame(tty, directed = F,vertices=an)
deg=degree(gr, mode ="all")

ggraph(gr, layout = 'linear') + ylim(-1.5,8) +     geom_edge_arc(aes(alpha =countAboveMean0.3>8)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2,nudge_y=-1) + coord_flip()
ggsave("forFa_2.pdf")
ggraph(gr, layout = 'kk') +  geom_edge_link(aes(alpha =countAboveMean0.3>8)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=deg)) + geom_node_label(aes(label=paste(name,ID)),label.size=.2,nudge_x=.1) + coord_flip()
ggsave("forFa.pdf")

# see if reflected in values

#make useful Seurat Object
load('r_objects181214_ALL_reIDforcellphone.rdata')
load("r_objects/ALLs15.rdata")
load("r_objects/181219_cellphoneDB_allCellTypes.rdata")
ALLs15@meta.data$cellphoneDB_id=ALL_clusterID[rownames(ALLs15@meta.data),"newID"]

m<-match(ALLs15@meta.data$cellphoneDB_id,an$name)
ALLs15@meta.data$cellphoneDB_designation=an$ID[m]
ALLs15@meta.data$cellphoneDB_designation[is.na(ALLs15@meta.data$cellphoneDB_designation)]="orphan"
ALLs15@meta.data$cellphoneDB_comb<-paste0(ALLs15@meta.data$cellphoneDB_id,"_",ALLs15@meta.data$cellphoneDB_designation)
TSNEPlot(ALLs15,group.by="cellphoneDB_comb",do.label=T, do.return=T,pt.shape="genotype",,pt.size=1.5)
ALLs15<-SetIdent(ALLs15,ident.use=ALLs15@meta.data$cellphoneDB_comb)
TA<-SubsetData(ALLs15,cells.use=rownames(ALLs15@meta.data)[!is.na(m)])

## pull apart gene pairs
x<-strsplit(rownames(toi),"_")
df<-as.character(tools::toTitleCase(tolower(unlist(lapply(x, function(x) x[1])))))
ip<-as.character(tools::toTitleCase(tolower(unlist(lapply(x, function(x) x[2])))))
ddf<-data.frame(g1=df,g2=ip,stringsAsFactors=F)
tddf<-ddf[!(grepl(" ",ddf$g1) | grepl(" ",ddf$g2)),]
tddf<-tddf[tddf$g1 %in% rownames(ALLs15@data) & tddf$g2 %in% rownames(ALLs15@data),]
pp<-apply(tddf,1, function(x) DotPlot(TA,x,do.return=T,x.lab.rot=T))
plot_grid(plotlist=pp[1:20],nrow=2)
ggsave("plots/181219_sig_var_genePairs_inAllcells.pdf")