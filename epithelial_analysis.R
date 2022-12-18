###Analysis of epithelial cells
###Files and data needed for analysis
# - source functions file
# - merged=Raw UMI matrix, found under GSE182227
# - metadata=Cell metadata table, found in supplement
# - hg19=hg19 genome, found in biomart
# - tcga=TCGA subtype signatures, found in supplement
# - programs=Cancer metaprograms, found in supplement
# - celltypes=cell type signatures found in supplement
# - defmat=data frame with CNA correlation and signal values for all epithelial cells, output from cna function

#Plot all cancer cells
dt5<- metadata %>% filter(CNA=="Cancer") %>% rownames %>% filter.umi(matrix = merged) %>%
  standard.umap(spread=5,n_neighbors = 100,var_genes = 5000)
l1<-c("Patient","HPV")

l2<-lapply(l1,function(x){
ggplot(dt5, aes(x=V1, y=V2, colour=metadata[rownames(dt5),x])) +
  geom_point(size=1,alpha=0.5) + 
  guides(colour = guide_legend(title=x,override.aes = list(size=4,alpha=1)))+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=20) +scale_colour_manual(values=pal1)
})
do.call(plot_grid,l2)

#Plot score cutoffs
fm<- metadata %>% filter(CNA=="Cancer") %>% rownames() %>% filter.umi(matrix=merged)
ss1<-tcga
ssx<-ss1[ss1$Gene %in% rownames(fm),]
tpf<-type.score(fm,ssx,tv="Subtype",bin=100)
tl<-do.call(cbind.data.frame,tpf$typescorelist)
tl$Patient<-metadata[rownames(tl),"Patient"]
ttr<- tl %>% group_by(Patient) %>% summarise_all(mean)
tpr2<-reshape2::melt(ttr)
colnames(tpr2)[2]<-"Subtype"

ggplot(tpr2,
       aes(x=value,colour=Subtype))+
  scale_colour_manual(values=c("dodgerblue","salmon","gold"))+
  geom_density(size=2)+theme_light(base_size = 30)+
  geom_vline(xintercept = c(-0.25,0.15,1),lty=2,size=2,colour=c("salmon","gold","dodgerblue"))+
  ggtitle("Mean scores per tumour")
dev.off()

#Centre TCGA scores by cutoff
tl2<-tl
tl2$Basal<-tl2$Basal-1
tl2$Classical<-tl2$Classical-0.15
tl2$Atypical<-tl2$Atypical+0.25
tl2<-tl2[,1:3]

tz<-numeric(length=nrow(tl2))
for(i in 1:nrow(tl2)){
  
  a<- which.max(tl2[i,])
  b<- max(tl2[i,-a])
  c<- max(tl2[i,])
  tz[i]<- c-b}
head(tz)
tl2$diff<-tz

#Set TCGA type per cell using confidence interval
ci<-0.5
TCGA<-ifelse(tl2$Basal>0 & tl2$Atypical< 0 & tl2$Classical<0 & tl2$diff > ci, "Basal",
             ifelse(tl2$Basal<0 & tl2$Atypical> 0 & tl2$Classical<0 & tl2$diff > ci, "Atypical",
                    ifelse(tl2$Basal<0 & tl2$Atypical< 0 & tl2$Classical>0 & tl2$diff> ci,
                           "Classical","Unresolved")))
names(TCGA)<-rownames(tl2)
metadata[names(TCGA),"TCGA_Type"]<-TCGA
#Plot TCGA type per cell
ggplot(dt5, aes(x=V1, y=V2, colour=TCGA[rownames(dt5)])) +
  geom_point(size=1,alpha=0.5) + 
  guides(colour = guide_legend(title="Subtype",override.aes = list(size=6,alpha=1)))+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=20) +
  scale_colour_manual(values=c("salmon","dodgerblue","gold","darkgrey"))

#Plot TCGA type percentages
names(TCGA)<-rownames(tl2)
hn<-c("OP16","OP19","OP10","OP12","OP8")
tt2<-cbind.data.frame(TCGA,metadata[rownames(tl2),c("Patient","Subclone")])
colnames(tt2)[2]<-"Patient"
tt2<-tt2[tt2$TCGA!="Unresolved",]

tt2<-metadata%>%filter(CNA=="Cancer"&TCGA_Type!="Unresolved")%>%
  select(TCGA_Type,Patient,Subclone)

#Choose either Patient or Subclone to plot percentages by either variable
dd<-as.data.frame(table(tt2$Subclone,tt2$TCGA))
dd<-as.data.frame(table(tt2$Patient,tt2$TCGA))
ppl2<-as.character(unique(dd$Var1[!dd$Var1 %in% hn]))
ppl<-c("OP16","OP19","OP10","OP12","OP8",ppl2)

s1<-unique(tt2$Subclone)
s1<-unique(tt2$Patient)

ppl<-s1[grep("OP19|OP10|OP12|OP8|OP35",s1)]
ppl<-ppl[grep("Nor|Unr",ppl,invert = T)]
ppl<-sort(ppl)

pc<-lapply(ppl,function(z){
  a<-dd[dd$Var1==z,]
  sf<-sum(a$Freq)
  a$Freq<-a$Freq/sf
  
  ggplot(a, aes(x="", y=Freq, fill=Var2)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+theme_void(base_size = 30)+guides(fill=FALSE)+
    scale_fill_manual(values=c("salmon","dodgerblue","gold","darkgrey"))+
    ggtitle(z)
})
p11<-do.call(plot_grid,pc)

p22<-ggplot(dd[dd$Var1=="OP10",], aes(x="", y=Freq, fill=Var2)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+theme_void(base_size = 30)+
  scale_fill_manual(name="Subtype",values=c("salmon","dodgerblue","gold","darkgrey"))

legend <- get_legend(
  p22 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

plot_grid(p11, legend, rel_widths = c(3, .4))

#Significance test for enriched subclones
pl<-c("OP10","OP12","OP8","OP19","OP35")
plf<-lapply(pl,function(x){
  a<-tt2[tt2$TCGA!="Unresolved" & tt2$Patient==x & tt2$Subclone %in% ppl,]
  table(a$TCGA,a$Subclone)
  chisq.test(a$TCGA,a$Subclone)
})
  
  

####4. Plot of most variable genes per patient

#Matrix of mean per patient
hn<-c("OP16","OP19","OP10","OP12","OP8")
meta2<-metadata[metadata$CNA=="Cancer",]
s2<-ifelse(meta2$Sample=="OP34norm","OP34_Invasive",ifelse(meta2$Sample=="OP35LN","OP35_LN",meta2$Patient))
names(s2)<-rownames(meta2)
su<-c("OP16","OP19","OP10","OP12","OP8",as.character(unique(s2[!s2 %in% hn])))
f1<-lapply(su,function(x){
  a<-names(s2[s2==x])
  a1<-fm[,a]
  rowMeans(a1)
  })
f1<-do.call(cbind.data.frame,f1)
colnames(f1)<-su

#Top DE genes
dex<-DE.genes(fm,s2,no=50)
dex2<- dex %>% arrange(factor(Cluster, levels = su))

d44<-lapply(su,function(x){
  dex2 %>% filter(Cluster==x) %>% pull(Gene)
  })
d44<-do.call(cbind.data.frame,d44)
colnames(d44)<-su

#Scores on top
pp2<-programs
pp2<-pp2[pp2$Gene %in% rownames(fm),]
tsx<-type.score(fm,pp2,tv="Program",bin=100)
tss<-do.call(cbind.data.frame,tsx$typescorelist)

sf<-mclapply(su,mc.cores = 8, function(x){
  a<-names(s2[s2==x])
  # b<-fm[,a]
  # c<-type.score(b,pp2,tv="Program",bin=100)
  # d<-do.call(cbind.data.frame,c$typescorelist)
  # e<-colMeans(d)
  e<-colMeans(tss[a,])
  f<-metadata[a,"HPV"]
  g<-length(f[f=="HPV+"])/length(a)
  h<-c(e,g)
  })
sf2<-do.call(rbind.data.frame,sf)
rownames(sf2)<-su
colnames(sf2)<-c(unique(pp2$Program),"HPVperc")
sf2$HPVclin<-ifelse(rownames(sf2) %in% c("OP10","OP12","OP16","OP19"),"HPVneg","HPVpos")

#Annotate and plot
da<-sf2
ca<-columnAnnotation(show_legend=TRUE, Atypical=da$Atypical,
                     Basal=da$Basal,
                     Classical=da$Classical,
                     G1S=da$G1S,
                     G2M=da$G2M,
                     EpiSen=da$EpiSen,
                     pEMT=da$pEMT,
                     Stress=da$Stress,
                     Interferon=da$Interferon,
                     Hypoxia=da$Hypoxia,
                     OxPhos=da$OxPhos,
                     HPVperc=da$HPVperc,
                     HPVclin=da$HPVclin,
                     annotation_name_gp=gpar(fontsize=20), simple_anno_size = unit(1, "cm"),
                     annotation_legend_param=list(grid_width=unit(1, "cm"),
                                                  grid_height=unit(1,"cm"),
                                                  title_gp=gpar(fontsize=20),
                                                  labels_gp=gpar(fontsize=20)),
                     
                     
                     col=list(HPVclin=c("HPVneg"="chartreuse","HPVpos"="red"),
                       Atypical=circlize::colorRamp2(c(min(da$Atypical), median(da$Atypical), max(da$Atypical)), 
                                                     c("white", "salmon", "salmon4")),
                       Basal=circlize::colorRamp2(c(min(da$Basal), median(da$Basal), max(da$Basal)), 
                                                  c("white", "dodgerblue", "dodgerblue4")),
                       Classical=circlize::colorRamp2(c(min(da$Classical), median(da$Classical), max(da$Classical)), 
                                                      c("white", "gold", "gold4")),
                       G1S=circlize::colorRamp2(c(min(da$G1S), median(da$G1S), max(da$G1S)), 
                                                c("white", "magenta", "magenta4")),
                       G2M=circlize::colorRamp2(c(min(da$G2M), median(da$G2M), max(da$G2M)), 
                                                c("white", "deeppink", "deeppink4")),
                       EpiSen=circlize::colorRamp2(c(min(da$EpiSen), median(da$EpiSen), max(da$EpiSen)), 
                                                       c("white", "chartreuse", "chartreuse4")),
                       pEMT=circlize::colorRamp2(c(min(da$pEMT), median(da$pEMT), max(da$pEMT)), 
                                                 c("white", "lightsalmon", "darkred")),
                       HPVperc=circlize::colorRamp2(c(min(da$HPVperc), median(da$HPVperc), max(da$HPVperc)), 
                                                    c("white", "red", "red4")),
                       Stress=circlize::colorRamp2(c(min(da$Stress), median(da$Stress), max(da$Stress)), 
                                                c("white", "seagreen", "seagreen4")),
                       Interferon=circlize::colorRamp2(c(min(da$Interferon), median(da$Interferon), max(da$Interferon)), 
                                                   c("white", "blue", "navy")),
                       Hypoxia=circlize::colorRamp2(c(min(da$Hypoxia), median(da$Hypoxia), max(da$Hypoxia)), 
                                                 c("white", "purple", "purple4")),
                       OxPhos=circlize::colorRamp2(c(min(da$OxPhos), median(da$OxPhos), max(da$OxPhos)), 
                                                    c("white", "cyan", "cyan4"))
                     ))


f2<-f1[dex3$Gene,]
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("steelblue", "white", "darkred"))
Heatmap(as.matrix(f2),show_row_names = F,show_column_names = T,cluster_rows = F,cluster_columns = F,
        col=col_fun,name = "log2(TPM)",right_annotation = hr,
        top_annotation = ca,          column_title_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 20),
        column_names_gp = gpar(fontsize=20),
        row_names_gp = gpar(fontsize=20),
        heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                  grid_height=unit(2,"cm"),
                                  title_gp=gpar(fontsize=20),
                                  labels_gp=gpar(fontsize=20)))
        
#Pseudobulk clustering
f3<-hclustreorder(cor(f1))
f3<-cor(f1)
col_dend = as.dendrogram(hclust(dist(1-f3,method = "euclidean"),method = "average"))
col_fun = circlize::colorRamp2(c(-0.4, 0, 0.4), c("steelblue", "white", "darkred"))
h1<-draw(Heatmap(f3,show_row_names = T,show_column_names = F,cluster_rows = col_dend,
                 cluster_columns = col_dend,col=col_fun,name = "Pearson\nCorrelation",
                 show_column_dend = TRUE,show_row_dend = FALSE,
        top_annotation = ca,column_split = 4,
        row_names_gp = gpar(fontsize=20),
        heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                  grid_height=unit(2,"cm"),
                                  title_gp=gpar(fontsize=20),
                                  labels_gp=gpar(fontsize=20)))
        )

h1

#HPV in nonmalignant cells
#Dimension reduction
dt6<- metadata %>% filter(CNA=="Normal") %>% rownames() %>% filter.umi(matrix=merged) %>%
  standard.umap(spread=5)
#Clustering
cl1<-cluster.coord(dt6,assign_all=T,minpts=50)
#Plot
ggplot(dt6, aes(x=V1, y=V2, colour=cl1)) +
  geom_point(size=1,alpha=0.5) + 
  guides(colour = guide_legend(title="Cluster",override.aes = list(size=4,alpha=1)))+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=20) +scale_colour_manual(values=pal1)
#Select cluster that contains HPV+ cells for comparison within cluster
dd<-cbind.data.frame(cl1,metadata[names(cl1),"HPV"])
colnames(dd)[2]<-"HPV"
cc<-dd%>%group_by(cl1)%>%count(HPV)%>%filter(HPV=="HPV+"&n>1)%>%pull(cl1)

#Select HPV+ nonmalignant cells
n1<-metadata[names(cl1[cl1%in%cc]),]%>%filter(HPV=="HPV+")%>%rownames()

#Add cancer cells to check that nonmalignant HPV+ don't have CNA
a<-table(metadata[n1,"Patient"])
n2<-unlist(lapply(names(a),function(x){
  b<-a[x]
  c<-sample(rownames(metadata[metadata$Patient==x & metadata$CNA=="Cancer",]),b,replace = F)
  c
}))
#Load defmat
#Nonmalignant HPV-negative cells from HPV-positive patients
n3<-metadata[names(cl1[cl1%in%cc]),]%>%filter(CNA=="Normal"&!Patient%in%hn&HPV=="HPV-")%>%
  rownames()
dm2<-defmat[c(n1,n2,n3),]
dm2$HPV<-metadata[c(n1,n2,n3),"HPV"]
dm2$class<-ifelse(rownames(dm2) %in% n2, "Cancer",
                  ifelse(rownames(dm2) %in% n1 & dm2$HPV=="HPV+","HPV+ normal","HPV- normal"))

#Violin plots of CNA signal and score distribution
  ggplot(dm2,aes(x=class,y=cnacor,fill=class))+geom_violin()+ylab("CNA Correlation")
  ggplot(dm2,aes(x=class,y=cnascore,fill=class))+geom_violin()+ylab("CNA Signal")
  

#DE genes HPV+ vs HPV- nomalignant
n1a<- dm2 %>% filter(class!="Cancer")%>%
  rownames()
fn2<-filter.umi(merged,n1a)
h1<-HPV[colnames(fn2)]
de2<-DE.genes(fn2,h1,volcano = TRUE,max_pv = 1,no=50000,min_fc=0)
de2$logFC<- -1*de2$logFC
#Volcano plot
volcano.plot(de2,"HPV+ vs HPV-","Normal epithelial cells",fc=0,pv=0.05)
