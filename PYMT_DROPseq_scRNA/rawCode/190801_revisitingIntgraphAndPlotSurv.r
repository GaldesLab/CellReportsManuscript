###Modifying the survival analysis
devtools::install_github("https://github.com/briglo/funkyDave.git")
library(funkyDave)
load("~/Desktop/briglo_backup-old/scRNA/davgal/paperstuff/r_objects/181128_survivalObjects.rdata")
plotSurv(genes=sample(rownames(cdat),25),expdat=cdat,df=clindat) #has the pvalue we were working off

load("~/Desktop/briglo_backup-old/scRNA/davgal/paperstuff/r_objects/181206_survivalObjects.rdata")
plotSurv(genes=sample(rownames(cdat),25),expdat=cdat,df=clindat) #has different p value (cant find why, but trimming used attempted to get biult into pltoSurv function)
plotSurv(genes=sample(rownames(cdat),25),expdat=cdat,df=clindat,trimFirst=T) #doesnt seem to make a difference,im probably missing something

#splitting by subtype
sclin<-split(clindat,clindat$CLAUDIN_SUBTYPE)
scdat<-lapply(sclin, function(x) cdat[,rownames(x)])
plotSurv(genes=sample(rownames(cdat),25),expdat=scdat$Basal,df=sclin$Basal)

# getting tcga to work requires dicing clincal data differently (no claudin etc its a PITA)




### revisiting cellphoneDB (to have something pretty JUST IN CASE governer wants some pretty picures)

#final version love it
intgraph<- function(scoreCut=0.3, numberCut=22, numberSplit=60){
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
   p1<-ggraph(gr, layout = 'linear', circular=T) +     geom_edge_arc(aes(width =score,alpha=score>numberSplit,colour=score>numberSplit)) + geom_node_point(aes(color=paste(name,ID),shape=type,size=count)) + geom_node_text(aes(label=paste(name,ID)), repel=F) + ggtitle(paste0("scoreCut=",scoreCut, ", numberCut=",numberCut,", numberSplit=",numberSplit))
   p2<-ggplot(tmp,aes(x=score,colour=score<numberSplit)) + geom_bar() + ggtitle("score dist")
   ggdraw() + draw_plot(p1 + theme(legend.justification="bottom"), 0,0,1,1) + draw_plot(p2  + theme(legend.position = "none") ,.85, 0, 0.15  , 0.15)
   }
intgraph()

#integrating with TSNEplot
load("~/Downloads/190111_newNetworkGraph.rdata")
locs<-data.frame(do.call(rbind,lapply(split(data.frame(TA@dr$tsne@cell.embeddings),TA@meta.data$cellphoneDB_id),colMeans))
colnames(locs)<-c('x','y')
cluster_anno<-cbind(cluster_anno,locs[cluster_anno$name,])

p1<- ggraph(gr,layout = 'manual',node.position=cluster_anno) +     geom_edge_arc(aes(width =score,alpha=score>numberSplit,colour=score)) + geom_node_point(aes(x=tSNE_1,y=tSNE_2,shape=type,size=count)) + geom_node_text(aes(x=tSNE_1,y=tSNE_2,label=paste(name,ID)), repel=T,size=6) + theme(legend.position = "none") + xlim(-50, 50) + ylim(-50, 50) 
p2<-ggplot(cbind(TA@meta.data,TA@dr$tsne@cell.embeddings),aes(x=tSNE_1,y=tSNE_2,colour=cellphoneDB_comb)) + geom_point() + theme(legend.position = "none") + xlim(-50, 50) + ylim(-50, 50) 

comb<- ggdraw() + draw_plot(p2 + theme(legend.position = "none"),0,0,1,1) + draw_plot(p1,0,0,1,1)