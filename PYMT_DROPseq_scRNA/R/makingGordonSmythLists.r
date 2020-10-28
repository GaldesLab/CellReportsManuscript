# got from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63310
# which was published at https://bmccancer.biomedcentral.com/articles/10.1186/s12885-015-1187-z#Sec2
# which was cited by https://www.nature.com/articles/s41467-017-01560-x?elqTrackId=c966e5703dee49788c3a69ee9ba8c969 (PAL)
# it also appears to have top200 lists dor luminal vs basal in supp tab 3
# These are counts, they Voomed and limma-edd 


# So what I did was
library(edgeR)
fnam<-dir(pattern='txt')[1:9]
rdat<-lapply(fnam, function(x) data.frame(read.table(x,header=T,stringsAsFactors=F),row.names=1))
cdat<-data.frame(do.call(cbind,lapply(rdat, function(x) x$Count)))
rownames(cdat)<-row.names(rdat[[1]])
colnames(cdat)<-snam<-gsub(".txt","",fnam)
v<-voom(cdat)
ct<-factor(rep(c('basal','LP','ML'),each=3))


ct<-relevel(ct,ref="ML")
design <- model.matrix(~0+ct)
fit<-lmFit(v,design)
fit2<-fit[,c('ctLP','ctbasal')]
fit2<-eBayes(fit2,trend=T)
results<-decideTests(fit2)
ML.sig<-rowSums(results<0)==2
ML.tab=topTable(fit2[ML.sig,],number=10000)

ct<-relevel(ct,ref="LP")
design <- model.matrix(~ct)
fit<-lmFit(v,design)
fit2<-fit[,c('ctML','ctbasal')]
fit2<-eBayes(fit2,trend=T)
results<-decideTests(fit2)
LP.sig<-rowSums(results<0)==2
LP.tab=topTable(fit2[LP.sig,],number=10000)

ct<-relevel(ct,ref="basal")
design <- model.matrix(~ct)
fit<-lmFit(v,design)
fit2<-fit[,c('ctLP','ctML')]
fit2<-eBayes(fit2,trend=T)
results<-decideTests(fit2)
basal.sig<-rowSums(results<0)==2
basal.tab=topTable(fit2[basal.sig,],number=10000)

save(ct,fit,LP.sig,LP.tab,ML.sig,ML.tab,basal.sig,basal.tab,file="GSE63310_output_allSig.rdata")
sigGene_ent<-list(basal=names(basal.sig)[basal.sig],LP=names(LP.sig)[LP.sig],ML=names(ML.sig)[ML.sig])



#sanity checked with topTable(fit2[basal.sig,])
#redid as appropriate for ML and LP

signatures<-data.frame(do.call(cbind,list(LP=LP.sig,ML=ML.sig,basal=basal.sig)))
library(biomaRt)
mouse<- useMart(biomart='ensembl', dataset = "mmusculus_gene_ensembl")
ano<-getBM(attributes=c("ensembl_gene_id","entrezgene","mgi_symbol"),filters="entrezgene",values=row.names(signatures),mart=mouse)

m<-match(rownames(signatures),ano$entrezgene)
signatures$symbol=ano$mgi_symbol[m]
#sanity checked again at entrezgene

bas<-na.omit(signatures$symbol[signatures$basal])
lp<-na.omit(signatures$symbol[signatures$LP])
ml<-na.omit(signatures$symbol[signatures$ML])
signatureGenes<-list(basal=bas,LP=lp,ML=ml)




save(signatureGenes,sigGene_ent,file="GSE63310_signatureGenes_allSig.rdata") #also done for 2 fold, 4 fold (lfc=1 or 2 in decideTests)

