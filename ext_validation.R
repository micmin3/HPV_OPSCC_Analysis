###Validation in external datasets, cell lines and TCGA
###Files and data needed for analysis
# - source functions file
# - merged=Raw UMI matrix, found under GSE182227
# - ex =Collect all single-cell datasets used for validation (see Methods) as well as the OPSCC data as UMI matrices, 
# then score them, without centering, for the G1/S program using the type.score function so you end up with one column of G1/S scores and one
# dataset column for each cell. Split the OPSCC data into HPVon, HPVoff and HPVneg.
# - ma3=CPM matrix from GSE157220, selecting only cellines containing "UPPER_AERO" in their description
# - tcga=TCGA subtype signatures, found in supplement
# - programs=Cancer metaprograms, found in supplement
# - celltypes=cell type signatures found in supplement
# - exp=TCGA expression data for HNSC and CESC, available through GDAC firehose
# - clin=TCGA clinical annotations for HNSC and CESC, available through GDAC firehose
# - hpv_tcga=HPV expression data for TCGA samples, avaliable in supplement of article with PMID 24085110


#Plot fraction of cells in top 3 bins, comparing across multiple scRNAseq datasets
ef<-mclapply(seq(1,100),function(x){
  ep<-unique(ex$Dataset)
  en<-unlist(mclapply(ep,function(x){
    ex %>% filter(Dataset==x) %>% rownames() %>% sample(size=1000,replace = F)
      }))
  
  matrix<-ex[en,]
  bins=10
  matrix$bin<-ntile(matrix$G1S,bins)
  matrix$top<-ifelse(matrix$bin>7,"yes","no")
  
  
  mx<- matrix %>% group_by(Dataset) %>% dplyr::count(top) %>% mutate(exp=n/sum(n)) %>%
    filter(top=="yes") %>% as.data.frame()
  rownames(mx)<-mx$Dataset
  my<-as.data.frame(mx$exp)
  rownames(my)<-rownames(mx)
  my})


ef2<-do.call(cbind.data.frame,ef)
colnames(ef2)<-seq(1,100)
m<-rowMeans(ef2)
s<-apply(ef2,1,sd)
ef3<-cbind.data.frame(m,s)
colnames(ef3)<-c("Mean","SD")
rownames(ef3)<-rownames(ef2)
op<-rownames(ef3)[grep("OPSCC",rownames(ef3))]
ef3$group<-ifelse(rownames(ef3) %in% op,"OPSCC","Others")
rownames(ef3)[1]<-"Chen_2020_NPC"
ef3$Dataset<-rownames(ef3)
ef4<-rbind.data.frame(ef3[op,],ef3[!rownames(ef3) %in% op,])
ef3$Dataset <- factor(ef3$Dataset, levels=rownames(ef4))
ef3$g2<-ifelse(ef3$Dataset=="OPSCC_HPVon","a",ifelse(ef3$Dataset=="OPSCC_HPVoff","b",
                                                     ifelse(ef3$Dataset=="OPSCC_HPVneg","c","d")))

ggplot(ef3,aes(x=factor(Dataset,levels=rownames(ef4)),y=Mean,fill=g2))+geom_col()+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width=0.1)+
  theme_light(base_size = 20)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Fraction cells in top 3 G1S bins")+
  scale_fill_manual(values=c("red","darkgreen","chartreuse","grey"))+guides(fill=FALSE)



#Cell line validation
Celline<-sapply(strsplit(colnames(ma3),split="_",fixed=TRUE),function(x)(x[1]))
HPVgenes<-as.character(c("E1","E2","E5","E6","E7","L1","L2")) #All HPV genes
HPV_score<-apply(ma3[HPVgenes,],2,sum) #Sum of HPV UMIs
HPV<-ifelse(HPV_score>0,"HPV+","HPV-") #Binary score if >0 HPV UMIs
#Log-transform matrix
fm<-filter.umi(ma3,colnames(ma3),start = "tpm")
#Dimension reduction
dtx<-standard.umap(fm)

#Plot CDKN2A in cell lines
ggplot(dtx, aes(x=V1, y=V2, colour=as.numeric(fm["CDKN2A",]))) +
  geom_point(size=1,alpha=0.2) + 
  scale_colour_gradientn(name="CDKN2A",colours=c("seagreen","gold","darkred"))+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=20)# +


#Load programs
px<-programs
px<-px[px$Gene %in% rownames(fm),]
#Score all cells
p2<-unique(px$Program)
prr<-lapply(p2,function(x){
  a<-px[px$Program==x,"Gene"]
  b<-score(fm,a,bin=100)
  b})
names(prr)<-p2
prx<-do.call(cbind.data.frame,prr)
prx$HPV<-HPV[colnames(fm)]
prx$Celline<-Celline<-sapply(strsplit(rownames(prx),split="_",fixed=TRUE),function(x)(x[1]))
prx$p16<-as.matrix(fm)["CDKN2A",]
#Select HPV-positive cell lines
cn<-c("93VU","SCC47","SCC90")
prx$class<-ifelse(prx$Celline %in% cn, "HPVpos","HPVneg")
prx$c2<-ifelse(prx$class=="HPVpos" & prx$HPV=="HPV+","HPVon",ifelse(prx$class=="HPVpos" & prx$HPV=="HPV-","HPVoff",prx$class))
#Plot G1/S distribution for HPVpos cell lines
pp<-mclapply(cn,function(x){
  plot.distribution(prx[prx$Celline==x,],"G1S","HPV",ttest = "HPV+",bins=10,name=x,sample = "HPV",reps = 100,
                    split = "bins")+scale_colour_manual(values=pal1)+theme_light(base_size = 20)
})
do.call(plot_grid,pp)

#TCGA HPV Analysis
library(data.table)
library(survival)
library(survminer)

#Load data
exp<-fread("/home/labs/tirosh/tyler/TCGA_data/tcga_expression_data.csv")
clin<-fread("/home/labs/tirosh/tyler/TCGA_data/tcga_clinical_data.csv")
hpv_tcga<-fread("~/TCGA/hpv_tcga.csv")

#Select OPSCC samples - repeat, exchanging HNSC for CESC to analyse cervical cancer
h1<-as.data.frame(hpv_tcga[hpv_tcga$Cancer=="HNSC" & hpv_tcga$`Virus description`!="HPVneg",])

h1$ppm<-log2(h1$ppm +1)
rownames(h1)<-h1$`Sample barcode`
n1<-sapply(strsplit(rownames(h1), split='.01', fixed=TRUE),function(x)(x[1]))
rownames(h1)<-n1
c1<-as.data.frame(clin[clin$cancer_type=="HNSC",])
cx<-c("c01","c01.9","c09.9","c10.9")
c1$Site<-ifelse(c1$icd_10 %in% cx,"OPSCC","nonOPSCC")
c1<-c1%>%filter(Site=="OPSCC")
rownames(c1)<-c1$id
n2<-intersect(rownames(c1),rownames(h1))
h1<-h1[n2,]
c1<-c1[n2,]
rownames(h1)<-h1$`Sample barcode`
n3<-rownames(h1)
e1<-as.data.frame(exp[exp$id %in% n3,])
rownames(e1)<-e1$id
e1<-e1[n3,]
e1<-e1[,-1]


#Program scores
pr<-programs
pr<-pr[,c("Gene","Program")]
cx<-celltypes[celltypes$Celltype=="Epithelial",1:2]
colnames(cx)[2]<-"Program"
ht<-tcga
ht<-ht[ht$Subtype=="Atypical",]
colnames(ht)[2]<-"Program"
pr<-rbind(pr,cx,ht)
p2<-unique(pr$Program)
progs<-lapply(p2,function(x){
  a<-pr[pr$Program==x,"Gene"]
  b<-score(t(e1),a,bin=NULL)
  b})
names(progs)<-p2
df<-do.call(cbind.data.frame,progs)

#Add HPV data to program scores
df$HPV<-h1$ppm


####Regress out epithelial purity
rn<-colnames(df[,grep("Epithelial",colnames(df),invert = T)])
rf<-lapply(rn,function(x){
  a<-lm(df[,x]~df$Epithelial)
  r<-a$residuals
  names(r)<-rownames(df)
  r
})
df<-do.call(cbind.data.frame,rf)
colnames(df)<-rn
#

#HPV correlation with G1/S
  df2<-df[,c("G1S","HPV")]
  corr.scatter(df2)+theme_light(base_size = 20)


#Survival analysis, splitting into HPVhigh and HPVlow groups
rownames(df)<-sapply(strsplit(rownames(df), split='.01', fixed=TRUE),function(x)(x[1]))

df2<-cbind.data.frame(df,c1)

#Find cutoff using maxstat
x<-"PFI"
y<-paste0(x,".time")

library("survminer")
cc <- surv_cutpoint(
  df2,
  time = y,
  event = x,
  variables = "HPV"
)$cutpoint[1]
df2$HPVclass<-ifelse(df2$HPV>as.numeric(cc),"high","low")

#Plot
hfit<- survfit(Surv(PFI.time,PFI==1)~HPVclass,data = df2)
survdiff(Surv(PFI.time,PFI==1)~HPVclass,rho=0,data = df2)
ggsurvplot(hfit,data=df2,risk.table = T,conf.int = FALSE,pval = TRUE,
           palette = c("red","darkgreen","blue"),tables.height = 0.4)

