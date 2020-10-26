library(cgdsr)
library(survminer)
library(survival)
library(cowplot)
library(Matrix)
library(cowplot)
mygenelist <- list("INVOLUTION"=c("Lifr","Il1a","Il1b","Il1r2","Il11","Il13",
"Cxcl1","Lsp1","Cd19","Cd83","Tnfsf7","Tnfsf9","A2m","Serpina1c","Fcer1c",
"Pigr","Sema4b","Sema7a","C1qa","C3ar1","C9","Tnfsf4","Tnfsf6","Tnfsf7",
"Tnfsf10","Tnfsf12","Tnfsf9","Il4i1","Gp5","Il11ra2","Ltbr","Cd80","Tnfrsf7",
"Hpxn","Orm2","Cebpd","Axcam","E4bp4","Inhbb","Pafah1b2","Il12rb1","Ccl9",
"Cd14","Lcn2","Cebpd","Necl1","C2","Tgfb1","Tgfb3","Tnfrsf6","Myd88","Emr1",
"Ecgf1","Pla2g2a","Fadd","Osmr","Il1r1","Il10rb","Csf1r","Ccl6","Ccl7","Ccl8",
"Ccl21a","Cxcl9","Cxcl14","Cd79a","Ly86","Cd86","Tnfrsf4","Lrp1","Hp","Saa1",
"Saa2","Saa3","Saa-ps","Sema4a","Ltf","C3","C4","C2","C1qa","C1qc","Casp1"))
genes <- c(toupper(mygenelist[["INVOLUTION"]]), "ELF5")
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
# Test the CGDS endpoint URL using a few simple API tests
test(mycgds)


################################################################################
#################Get inovlution expression for METABRIC cohort##################
################################################################################

#Select METABRIC cohort from cgda
getCancerStudies(mycgds)[,c(1,2)]
mycancerstudy.meta="brca_metabric"
#Select all patients
getCaseLists(mycgds,mycancerstudy.meta)[,c(1,2)]
mycaselist.meta="brca_metabric_all"
#Select patients with valid mrna data
getGeneticProfiles(mycgds,mycancerstudy.meta)[,c(1,2)]
mygeneticprofile.meta="brca_metabric_mrna"
#Download expression of involution signature genes and ELF5 from METABRIC cohort
dat.meta <- getProfileData(mycgds, genes = genes, mygeneticprofile.meta,
  mycaselist.meta)
cdat.meta <- t(dat.meta)
cdat.meta <- as.data.frame(cdat.meta)
#Create an object with all clinical information from the METABRIC cohort
clin.meta<-getClinicalData(mycgds,mycaselist.meta)[,]
#Select only clinical data from patients with mrna data
clindat.meta<-clin.meta[colnames(cdat.meta),]

#Create the metascore table
metascore.meta <- cdat.meta
metascore.meta <- as.data.frame(metascore.meta)
#Remove all genes without any gene expression
metascore.meta <- metascore.meta[!is.na(metascore.meta$MB.0000), ]
metascore.meta <- t(metascore.meta)
metascore.meta <- as.data.frame(metascore.meta)
#Remove all patients that do not have mrna information for the genes of interest
metascore.meta <- metascore.meta[!is.na(metascore.meta$A2M), ]
metascore.meta <- metascore.meta[rownames(clindat.meta), ]
#Select expression of ELF5 to be used in data.frame for plotting
elf5 <- metascore.meta[, "ELF5"]
metascore.meta$ELF5 <- NULL #Removes ELF5 to determine involution  metascore
#Take the expression sum of all involution related genes per patient.
metascore.meta <- as.matrix(metascore.meta)
temp <- rowSums(metascore.meta)
metascore.meta<-as.data.frame(metascore.meta)
metascore.meta$metascore <- temp

#Create data.frame used for plotting
df.meta <- data.frame("patientID"=rownames(clindat.meta), "elf5"=elf5,
"involution"=metascore.meta$metascore)
df.meta <- df.meta[!is.na(df.meta$elf5),] #Remove patients without data on ELF5
q_e <- quantile(x = df.meta$elf5, c(0.33, 0.67))#Determine ELF5 tertile
df.meta$ELFcategory <- ifelse(df.meta$elf5 < q_e[1], "low",
ifelse(df.meta$elf5 > q_e[2], "high", "mid"))
df.meta$ELFcategory <- factor(df.meta$ELFcategory, levels =
  c("high", "mid", "low"))
#Remove patients in middle ELF5 expression tertile
df.meta <- df.meta[df.meta$ELFcategory != "mid", ]

#Create METABRIC expression plot
involution_signature_expression_METABRIC <- ggplot(data = df.meta,
  aes(y = involution, x = ELFcategory)) + geom_boxplot() +
  theme_minimal(base_size = 20) + theme(legend.position="none") +
  xlab("ELF5 expression category") + ylab("Involution signature metascore")
involution_signature_expression_METABRIC
t.test(df.meta[df.meta$ELFcategory == "high", "involution"],
df.meta[df.meta$ELFcategory == "low", "involution"])


################################################################################
####################Get inovlution expression for TCGA cohort###################
################################################################################

#Select most recent TCGA from cgda
getCancerStudies(mycgds)[,c(1,2)]
mycancerstudy.tcga="brca_tcga_pan_can_atlas_2018"
#Select all patients with mrna data
getCaseLists(mycgds,mycancerstudy.tcga)[,c(1,2)]
mycaselist.tcga="brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
#Select all patients with mrna data
getGeneticProfiles(mycgds,mycancerstudy.tcga)[,c(1,2)]
mygeneticprofile.tcga="brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
#Download gene expression of involution signature genes and ELF5
dat.tcga <- getProfileData(mycgds, genes = genes, mygeneticprofile.tcga,
  mycaselist.tcga)
dat.tcga <- log10(dat.tcga + 1) #log transform
cdat.tcga <- t(dat.tcga)

#Downloading all clinical information from TCGA cohort
clin.tcga<-getClinicalData(mycgds,mycaselist.tcga)[,]
clindat.tcga<-clin.tcga[colnames(cdat.tcga),]
clindat.tcga <- clindat.tcga[clindat.tcga$SEX == "Female",]

#Create the metascore table
cdat.tcga <- as.data.frame(cdat.tcga)
metascore.tcga <- cdat.tcga
metascore.tcga <- as.data.frame(metascore.tcga)
#Remove all genes without any gene expression
metascore.tcga <- metascore.tcga[!is.na(metascore.tcga$TCGA.3C.AAAU.01), ]
metascore.tcga <- t(metascore.tcga)
metascore.tcga <- as.data.frame(metascore.tcga)
#Remove all patients without involution gene expression
metascore.tcga <- metascore.tcga[!is.na(metascore.tcga$A2M), ]
metascore.tcga <- metascore.tcga[rownames(clindat.tcga), ]
elf5 <- metascore.tcga[, "ELF5"]
metascore.tcga$ELF5 <- NULL #Removes ELF5 prior to makeing involution metascore
#Take the expression sum of all involution related genes per patient.
metascore.tcga <- as.matrix(metascore.tcga)
temp <- rowSums(metascore.tcga)
metascore.tcga<-as.data.frame(metascore.tcga)
metascore.tcga$metascore <- temp

#Create data.frame used for plotting
df.tcga <- data.frame("patientID"=rownames(clindat.tcga),"elf5"=elf5,
"involution"=metascore.tcga$metascore)
df.tcga <- df.tcga[!is.na(df.tcga$elf5),] #Remove patients without data on ELF5
q_e <- quantile(x = df.tcga$elf5, c(0.33, 0.67))
df.tcga$ELFcategory <- ifelse(df.tcga$elf5 < q_e[1], "low",
ifelse(df.tcga$elf5 > q_e[2], "high", "mid"))
df.tcga$ELFcategory <- factor(df.tcga$ELFcategory,
  levels = c("high", "mid", "low"))
#Remove patients in middle ELF5 expression tertile
df.tcga <- df.tcga[df.tcga$ELFcategory != "mid", ]
#Create METABRIC expression plot
involution_signature_expression_TCGA <- ggplot(data = df.tcga,
  aes(y = involution, x = ELFcategory)) + geom_boxplot()+
  theme_minimal(base_size = 20) + theme(legend.position="none") +
  xlab("ELF5 expression category") + ylab(NULL)
involution_signature_expression_TCGA
t.test(df.tcga[df.tcga$ELFcategory == "high", "involution"],
df.tcga[df.tcga$ELFcategory == "low", "involution"])


#Combine METABRIC and TCGA expression plot into one plot
expression_plot <- cowplot::plot_grid(plotlist = list(
  involution_signature_expression_METABRIC,
  involution_signature_expression_TCGA),
   labels = c("Metabric", "TCGA"), label_x = 0.5)
expression_plot
