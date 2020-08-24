#Load necessary packages
library(cgdsr)
library(survminer)
library(survival)
library(BASALT)
library(cowplot)
#Load object containing interesting gene lists
load(file = "C:/Users/Jeron/Dropbox/NatureComs_Rebuttal_DATA2020/SurvivalAnalysisJeron/200120_mygenelist.rdata")

# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
# Test the CGDS endpoint URL using a few simple API tests
test(mycgds)
# Get list of cancer studies at server
getCancerStudies(mycgds)[,c(1,2)]
mycancerstudy="brca_metabric" #Select METABRIC breast cancer cohort

getCaseLists(mycgds,mycancerstudy)[,c(1,2)]
mycaselist="brca_metabric_all" #Get all cases

getGeneticProfiles(mycgds,mycancerstudy)[,c(1,2)]
mygeneticprofile="brca_metabric_mrna_median_Zscores"  #Select which expression data should be downloaded

#make a list of all genes for which to download expression data
genes <- c(unlist(mygenelist), "Elf5")
genes <- unique(genes)

########Option for very large genelists
#Loop over every gene in genes and download expression data individually. Data gets merged into one final object at the end
dat<-lapply(genes, function(x) getProfileData(mycgds,genes=x,mygeneticprofile,mycaselist))
dat <- dat[sapply(dat, function(x) dim(x)[1] != 0)]
dat <- unique(dat)
rowna <- sapply(dat, colnames)
ids<-unique(unlist(lapply(dat,row.names)))
cdat<-data.frame(do.call(rbind,lapply(dat,function(x) t(x[ids,]))))
colnames(cdat) <- rownames(dat[[1]])
rownames(cdat) <- rowna
cdat <- cdat[!is.na(cdat$MB.0000), ]
cdat <- cdat[, !is.na(cdat[1,])]

########Option for smaller genelists
#Download expression data for all genes at once
#dat <- getProfileData(mycgds, genes = genes, mygeneticprofile, mycaselist)
#cdat <- t(dat)

#preparing Object containing all clinical data
clin<-getClinicalData(mycgds,mycaselist)[,]
clindat<-clin[colnames(cdat),]

#censor to 5 disease sepecific survival
clindat$newOS_MONTHS <- ifelse(as.numeric(clindat$OS_MONTHS)>60,60,as.numeric(clindat$OS_MONTHS))
clindat$coded_status<-ifelse(clindat$VITAL_STATUS=="Died of Disease",1,0)
clindat$newcoded_status<-ifelse(as.numeric(clindat$OS_MONTHS)>60,0,clindat$coded_status)
clindat <- clindat[!is.na(clindat$newcoded_status), ]

save(cdat, clindat, file = "200120_METABRIC_survivalObject.rdata")

#Plotting functin that splits data into several bins based on gene expression and then constructs KaplanMeier curves
plotSurv<-function(genes, clindatName = "200120_METABRIC", Bins = "Bins: <33%; >67%"){
  message("Are the binning and Clindat name correct??")
  message("load mygenelist!")
  df=clindat
  expdat=cdat
  genelist=mygenelist

  label_list <- lapply(X = genelist, FUN = function(X) is.element(genes[genes != "Elf5"], X))
  label_list <- lapply(X = label_list, FUN = function(X) all(X == TRUE))
  label_list <- unlist(label_list)
  signatureLabel <- names(mygenelist[label_list])
  SurvObj <- clindatName

  genelist<-toupper(genes[toupper(genes) %in% rownames(expdat)])
  mscore<-colSums(expdat[genelist,rownames(df)],na.rm=T)
  mq<-quantile(mscore,c(.1,.25,.33,.5,.67,.75,.9))
  mb<-ifelse(mscore<mq[3],'low',ifelse(mscore>mq[5],"high","mid"))
  elfExp<-t(expdat['ELF5',rownames(df)])
  eq<-quantile(elfExp,c(.1,.25,.33,.5,.67,.75,.9),na.rm=T)
  eb<-ifelse(elfExp<eq[3],'low',ifelse(elfExp>eq[5],"high","mid"))
  comDF<-data.frame(df,elfExp=as.vector(elfExp),binnedElf5= as.vector(eb), metascore=mscore,metabin=mb)[!is.na(as.vector(eb)),]


  makeScores<-list('elf'=comDF[comDF$binnedElf5 %in% c("high",'low'),],
                   'met'=comDF[comDF$metabin %in% c("high",'low'),],
                   'split'=lapply(split(comDF[comDF$binnedElf5 %in% c("high",'low'),],comDF$binnedElf5[comDF$binnedElf5 %in% c("high",'low')]),function(x) {
                     qs<-quantile(x$metascore,c(.25,.33,.5,.67,.75))
                     binMetascore<-ifelse(x$metascore<qs[2],'low',ifelse(x$metascore>qs[4],"high","mid"))
                     df_temp <- data.frame(x,binScore=binMetascore)
                     df_temp <- df_temp[df_temp$binScore %in% c("high", "low"), ]
                     return(df_temp)
                   }))


  makeFits<-list('combElf'=survfit(Surv(newOS_MONTHS,newcoded_status) ~ binnedElf5, data=makeScores$elf),
                 'combMeta'=survfit(Surv(newOS_MONTHS,newcoded_status) ~ metabin, data=makeScores$met),

                 'spl'=lapply(makeScores$split[1:2], function(X) survfit(Surv(newOS_MONTHS,newcoded_status) ~ binScore, data=X)))



  p1a<-ggsurvplot(makeFits$combElf,data=makeScores$elf, pval=T,conf.int=F) + ggtitle("elf5 bins")
  p1b<-ggsurvplot(makeFits$combMeta,data=makeScores$met, pval=T,conf.int=F) + ggtitle("metascore bins")
  p2a<-ggplot(makeScores$elf,aes(x=binnedElf5,fill=CLAUDIN_SUBTYPE)) + geom_bar(stat='count') + ggtitle("elf5 composition")
  p2b<-ggplot(makeScores$met,aes(x=metabin,fill=CLAUDIN_SUBTYPE)) + geom_bar(stat='count')+ ggtitle("metascore composition")
  p3<-ggsurvplot(makeFits$spl,data=makeScores$split[1:2], pval=T,conf.int=F)
  p4<-lapply(lapply(makeScores[[3]], function(y) data.frame(table(y$binScore,y$CLAUDIN_SUBTYPE))), function(z) ggplot(z,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat='identity'))
  p5<-ggplot(comDF,aes(x=elfExp,y=metascore,color=as.factor(newcoded_status))) + geom_point() + geom_smooth(method = "lm", se=T)
  p_all <- plot_grid(p1a$plot,p2a,p1b$plot,p2b,p3$high$plot,p4$high,p3$low$plot,p4$low,p5,ncol=4)
  p_final <- p_all + draw_figure_label(label = paste0(Bins, "\n", SurvObj, "\n", signatureLabel), position = "bottom.right", size = 30, colour = "red", fontface = "bold")
  plot(p_final)
  return(invisible(comDF))
}
