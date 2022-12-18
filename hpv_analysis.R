###Specific analysis related to HPV
###Files and data needed for analysis
# - source functions file
# - merged=Raw UMI matrix, found under GSE182227
# - metadata=Cell metadata table, found in supplement
# - tcga=TCGA subtype signatures, found in supplement
# - programs=Cancer metaprograms, found in supplement

#Distribution of HPV genes compared to other gene scores
HPV<-metadata$HPV
names(HPV)<-rownames(metadata)
hn<-c("OP10","OP12","OP16","OP19","OP8")
HPVgenes<-as.character(c("E1","E2","E5","E6","E7")) #All HPV
#Create expression matrix
bm<- metadata %>% filter(!Patient %in% hn & CNA=="Cancer") %>% rownames() %>%
  filter.umi(matrix=merged,centre = "none",sparse = T,whitelist = HPVgenes,log = FALSE)

#Nearest neighbours for each gene by expression
ss<-rowMeans(bm)%>%as.data.frame()
ss$Gene<-rownames(ss)
rownames(ss)<-seq(1:nrow(ss))
bb<-FNN::get.knn(as.matrix(ss[,1]), k = 100)$nn.index
rownames(bb)<-ss$Gene
#Create expression-matched gene scores
sx<-mclapply(seq(1:1000),mc.cores=50,function(x){
  a<-unlist(mclapply(HPVgenes,function(y){
    ss[bb[y,],]%>%filter(!Gene %in% HPVgenes)%>%pull(Gene)%>%sample(.,1)
  }))
  b<-bm[a,]
  c<-colSums(b)
  d<-length(c[c==0])/ncol(bm)
  e<-mean(c)
  c(d,e)
})
ss2<-do.call(rbind.data.frame,sx)
colnames(ss2)<-c("perc0","meanExp")
b<-bm[HPVgenes,]
c<-colSums(b)
d<-length(c[c==0])/ncol(bm)
e<-mean(c)
HPVscore<-c(d,e)
ss2<-rbind.data.frame(HPVscore,ss2)
ss2$name<-ifelse(rownames(ss2)==1,"HPVtotal","")
ss2$Set<-ifelse(ss2$name=="HPVtotal","HPV","Control")

#Z-test on nearest neighbours
bb<-FNN::get.knn(as.matrix(ss2[,2]), k = 50)$nn.index
bb2<-ss2[bb[1,],]
m<-mean(bb2$perc0)
s<-sd(bb2$perc0)
dd<-rep(ss2[1,1],10)
DescTools::ZTest(x=dd,mu=m,sd_pop=s)
dx<-ss2[1,1]
(dx-m)/s

#Plot distribution of zeroes
ggplot(ss2,aes(x=perc0,fill=Set,label=name))+geom_histogram(binwidth = 0.01)+
  geom_text_repel(data=ss2[ss2$name=="HPVtotal",],aes(x=0.38,y=3,label="p<2.2e-16\nHPVtotal"),
                  size=10,nudge_y = 0.1)+
  theme_light(base_size = 30)+
  xlab("Fraction cells with 0 reads")+ylab("Number of gene sets")


#Plot expression of HPV genes per patient - same approach for subclones and cell lines, just change pl for list of subclones/cellines
hn<-c("OP10","OP12","OP16","OP19","OP8")
HPVgenes<-as.character(c("E1","E2","E5","E6","E7","L1","L2")) #All HPV
pl<-unique(metadata$Patient[!metadata$Patient %in% hn])

fh <- metadata %>% filter(CNA=="Cancer") %>% rownames() %>% 
  filter.umi(matrix=merged,whitelist = HPVgenes,centre="none")

hf<-lapply(pl,function(x){
  a<-rownames(metadata[metadata$Patient==x & metadata$CNA=="Cancer",])
  ax<-rownames(metadata[metadata$Patient==x & metadata$CNA=="Cancer" & metadata$HPV=="HPV+",])
  b<-merged[HPVgenes,a]
  bx<-merged[HPVgenes,ax]
  HPVtotal<-apply(b,2,sum)
  b2<-rbind(b,HPVtotal)
  c<-apply(b2,1,function(x)sum(x>0)/ncol(b2))
  d<-rowMeans(b2)

  cx<-matrix(ncol=ncol(bx),nrow=nrow(bx))
  rownames(cx)<-rownames(bx)
  for(i in 1:ncol(bx)){
    for(j in 1:nrow(bx)){
      a1<-bx[j,i]
      a2<-a1/sum(bx[,i])
      cx[j,i]<-a2
      
    }
  }
  
  
  e<-cbind.data.frame(c,d)
  colnames(e)<-c("pExp","MeanReads")
  e$Patient<-x
  e$Gene<-rownames(e)
  e$Gene[8]<-"HPVtotal"
  e$frac<-c(rowMeans(cx),1)
  rownames(e)<-seq(1:nrow(e))
  e
})
h2<-do.call(rbind.data.frame,hf)

z<-cbind.data.frame(0,0,"OP8",c(HPVgenes,"HPVtotal"),0)
colnames(z)<-colnames(h2)
z2<-cbind.data.frame(0,0,"HPVneg",c(HPVgenes,"HPVtotal"),0)
colnames(z2)<-colnames(h2)
h2<-rbind.data.frame(h2,z,z2)
h2$Patient<-factor(h2$Patient,levels=c(pl,"OP8","HPVneg"))
h2$Gene<-factor(h2$Gene,levels=c(HPVgenes,"HPVtotal"))

ggplot(h2,aes(x=Patient,y=Gene,size=pExp,colour=frac))+geom_point()+
  theme_light(base_size = 30)+xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_size_continuous(range=c(0,20),name = "Fraction\nnonzero cells")+
  scale_colour_gradientn(colours=c("lightpink", "salmon", "darkred"),
                         name="Mean fraction of \nHPV reads per cell")


#Metaprogram proportions across the three HPV classes in cancer cells
d33<-metadata[metadata$CNA=="Cancer","Subtype"]
names(d33)<-rownames(metadata[metadata$CNA=="Cancer",])

#Resample 100 cells per patient 100 times
b<- unique(metadata$Patient)
d4<-mclapply(seq(1,100),function(x){
  c<-unlist(lapply(b,function(y){
    c1<-metadata %>% filter(Patient==y & CNA=="Cancer" & Subtype!="Unresolved") %>% rownames()
    c2<-sample(c1,min(100,length(c1)),replace = F)
  }))
  d3<-d33[c]
  
  hc<-ifelse(metadata$Patient %in% hn, "HPVneg",
             ifelse(!metadata$Patient %in% hn & metadata$HPV=="HPV+","HPVon","HPVoff"))
  names(hc)<-rownames(metadata)
  hc<-hc[names(d3)]
  d4<-cbind.data.frame(d3,hc)
  colnames(d4)<-c("Program","HPVclass")
  d4 %>% group_by(Program) %>% count(HPVclass) %>% as.data.frame() %>% acast(Program~HPVclass) %>%
    as.data.frame()
})
d44<-do.call(cbind.data.frame,d4)
d5<-lapply(c("HPVneg","HPVoff","HPVon"),function(x){
  rowMeans(d44[,grep(x,colnames(d44))])
})
d5<-do.call(cbind.data.frame,d5)
colnames(d5)<-c("HPVneg","HPVoff","HPVon")

#Hypergeometric test for each combination
d55<-d5
#Uncomment the two lines below to compare only HPVpos vs HPVneg
# d5$HPVpos<- rowSums(d5[,-1])
# d5<-d5[,-c(2,3)]
d6<-matrix(nrow=nrow(d5),ncol=ncol(d5))
colnames(d6)<-colnames(d5)
rownames(d6)<-rownames(d5)

for(i in 1:nrow(d5)){
  for(j in 1:ncol(d5)){
    
    prog<-sum(d5[i,])
    cl<-sum(d5[,j])
    overlap<-d5[i,j]
    tot=sum(rowSums(d5))
    t<-phyper(q=overlap-1,m=prog,n=tot-prog,k=cl,lower.tail = F)
    d6[i,j]<-t
    
  }
}
d5<-d55
d5
d6
d7<-reshape2::melt(d6)
d7$value<-p.adjust(d7$value)
d8<-acast(d7,Var1~Var2)
d8

#Add summary
d5<-rbind.data.frame(d5,colSums(d5))
rownames(d5)[nrow(d5)]<-"Total"

#Plot
dz<-d5
dz$Program<-rownames(dz)
dz<-melt(dz)
head(dz)
dz$Tumor<-ifelse(dz$variable=="HPVneg","HPV-","HPV+")
z<- d5 %>% filter(!rownames(d5) %in% "Total") %>% mutate(ratio=HPVneg/rowSums(.)) %>% arrange(.,ratio) %>% rownames()
z<-c("Total",z)
zq<-dz%>%filter(Program=="Total")%>%mutate(frac=value/sum(value))%>%filter(variable=="HPVon")%>%pull(frac)
pp1<-ggplot(dz,aes(x=factor(Program,levels=z),fill=variable,y=value))+geom_col(position = "fill")+
  theme_light(base_size = 20)+xlab("")+ylab("")+
  scale_fill_manual(name="Subset",values=c("chartreuse","darkgreen","red"))
pp1


#Scatter plot comparing HPVon/HPVoff/HPVneg gene and program expression

#Expression matrix
fm<- metadata %>% filter(CNA=="Cancer") %>% rownames() %>% 
  filter.umi(matrix = merged)

#Averages between HPVon and HPVoff per patient
pl<-unique(metadata$Patient[!metadata$Patient %in% hn])

f1<-mclapply(pl,function(x){

  a<- metadata %>% filter(CNA=="Cancer" & Patient==x & HPV=="HPV+") %>% rownames()
  b<- metadata %>% filter(CNA=="Cancer" & Patient==x & HPV=="HPV-") %>% rownames()
  
  rowMeans(fm[,a]) - rowMeans(fm[,b])
})  
f2<-do.call(cbind.data.frame,f1)
colnames(f2)<-pl
OnVsOff_log2FC<-rowMeans(f2)
OnVsOff_pval<-apply(f2,1,function(x)t.test(x)$p.value)
OnVsOff_pval<-p.adjust(OnVsOff_pval,method = "BH")

#Patient averages HPVon vs HPVneg
pl<-unique(metadata$Patient)
f1<-lapply(pl,function(x){
  print(x)
  if(x %in% hn){
  a<- metadata %>% filter(CNA=="Cancer" & Patient==x) %>% rownames()
  b<-rowMeans(fm[,a])
  }else{
    a<- metadata %>% filter(CNA=="Cancer" & Patient==x & HPV=="HPV+") %>% rownames()
    b<-rowMeans(fm[,a])
  }

b
})  
f2<-do.call(cbind.data.frame,f1)
colnames(f2)<-pl
v1<-ifelse(colnames(f2) %in% hn,"neg","pos")
names(v1)<-colnames(f2)
OnVsNeg_log2FC<-rowMeans(f2[,!colnames(f2) %in% hn]) - rowMeans(f2[,colnames(f2) %in% hn])
OnVsNeg_pval<-apply(f2,1,function(x)t.test(x[v1=="pos"],x[v1=="neg"])$p.value)
OnVsNeg_pval<-p.adjust(OnVsNeg_pval,method = "BH")

#Data frame with all averages and p-values
t1<-cbind.data.frame(OnVsOff_log2FC,OnVsOff_pval,OnVsNeg_log2FC,OnVsNeg_pval)
pr<-programs
hr<-tcga
colnames(hr)[2]<-"Program"
pq<-rbind.data.frame(hr[hr$Program %in% c("Atypical","Basal","Classical"),],
                     pr[,c("Gene","Program")])

t1$gene<-rownames(t1)
z<-unique(pq$Program)
for(i in 1:length(z)){
  t1[t1$gene %in% pq[pq$Program==z[i],"Gene"],"geneset"]<-z[i]
}
t1[is.na(t1$geneset),"geneset"]<-"rest"

w<- t1 %>% group_by(geneset) %>% summarise_all(mean) %>% as.data.frame()
rownames(w)<-w$geneset
w$geneset<-NULL
w$gene<-"Average"
w$geneset<-rownames(w)
w<-w[w$geneset!="rest",]

#Data frame including programs
t1<-rbind.data.frame(w,t1)

#Plot genes
ggplot(t1[t1$geneset=="rest" & t1$gene!="Average",],aes(x=OnVsOff_log2FC,y=OnVsNeg_log2FC))+
  geom_point()+
  geom_point(data=t1[t1$geneset!="rest",],aes(colour=geneset))+
  scale_colour_manual(values=c("salmon","dodgerblue","gold","magenta","red",
                               "deeppink","seagreen","navy",pal2[3],pal2[4],pal2[2],pal2[10]))+
  theme_light(base_size = 30)+
  guides(size=FALSE,colour=guide_legend(override.aes=list(size=4)))+
  geom_text_repel(data=t1[t1$gene!="Average",],nudge_y = 0.1,box.padding = 0.8,
                  max.overlaps = 40,
                  aes(label=ifelse(abs(OnVsOff_log2FC)>1 |
                                     abs(OnVsNeg_log2FC)>2,
                                   as.character(gene),'')))+
  geom_vline(xintercept = c(-1,1),lty=2)+geom_hline(yintercept = c(-2,2),lty=2)


#Plot just programs with error bars
zz<- t1 %>% filter(gene!="Average" & geneset!="rest") %>% group_by(geneset) %>%
  summarise_all(sd) %>% dplyr::select(OnVsOff_log2FC,OnVsNeg_log2FC) %>%
  as.data.frame()
colnames(zz)<-c("OnOffSD","OnNegSD")
zz2<-cbind.data.frame(t1[1:11,c(1,3,6)],zz)

ggplot(zz2,aes(x=OnVsOff_log2FC,y=OnVsNeg_log2FC,colour=geneset))+geom_point(size=4)+
  geom_errorbar(aes(ymin=OnVsNeg_log2FC-OnNegSD,ymax=OnVsNeg_log2FC+OnNegSD),width=0.3)+
geom_errorbar(aes(xmin=OnVsOff_log2FC-OnOffSD,xmax=OnVsOff_log2FC+OnOffSD),width=0.3)+
  theme_light(base_size = 30)+
  scale_colour_manual(values=c("salmon","dodgerblue","gold","magenta","red",
                               "deeppink","seagreen","navy",pal2[3],pal2[4],pal2[2],pal2[10]))
#Enrichment analysis
library(clusterProfiler)
library(msigdbr)

#Get hallmark msigdb gene seta
mdb <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

#ID conversion and overrepresentation test
fo<-metadata %>% filter(CNA=="Cancer" & HPV=="HPV+") %>% rownames() %>%
  filter.umi(matrix=merged) %>%rownames()
t2<-t1 %>% filter(OnVsOff_pval<0.05 & OnVsOff_log2FC> 1) %>% rownames() %>%
  enricher(TERM2GENE=mdb,universe=fo)
dotplot(t2)


#Distribution of program scores in cancer cells

#Score all programs for all normal/cancer ep cells
fm <- metadata %>% filter(CNA %in% c("Cancer","Normal")) %>% rownames() %>%
  filter.umi(matrix = merged,centre = "none")

pr<-programs
pr<-pr[pr$Gene %in% rownames(fm),]

sm<-type.score(fm,pr,tv="Program",bin=100)
sm<-do.call(cbind.data.frame,sm$typescorelist)
sm$Patient<-metadata[rownames(sm),"Patient"]
sm$HPV<-metadata[rownames(sm),"HPV"]
sm$CNA<-metadata[rownames(sm),"CNA"]
sm$class<-ifelse(sm$CNA=="Cancer" & sm$Patient %in% hn,"HPVneg",
                 ifelse(sm$CNA=="Cancer" & sm$HPV=="HPV+","HPVon",
                        ifelse(sm$CNA=="Normal","Normal","HPVoff")))

#Plot program distributions - sample 100 cells per patient and class
z<-unique(pr$Program)
pf<-mclapply(z,function(x){
p1<-plot.distribution(sm[sm$class!="Normal",],x,"class","Patient",ttest = "chisq",
                  reps=100,split="bins",bins=10,sample_no = 100)+
  scale_colour_manual(values=c("chartreuse","darkgreen","red"))+theme_light(base_size = 20)
p1
})
do.call(plot_grid,pf)

#Heatmap of G1/S distribution per patient

pl<-unique(sm$Patient[!sm$Patient %in% hn])
k=10
#Count cells/bin and perform chi-square test for every patient
ppf<-lapply(pl,function(x){
  
  a<-sm[sm$Patient==x & sm$CNA=="Cancer",]
  a$bins<-ntile(a$G1S,k)
  
  b<-a %>% group_by(class) %>% count(bins) %>% mutate(frac= log2(n/(sum(n)/k))) %>%
    as.data.frame()
  
  print(x)
  print(chisq.test(a$bins,a$class)$p.value)
b$pat<-x
b

})
z2<-do.call(rbind.data.frame,ppf)
#Plot heatmap
p2<-lapply(c("HPVon","HPVoff"),function(x){
  a<- z2[z2$class==x,c("bins","pat","frac")]
  b<-reshape2::acast(a,pat~bins)
  b[is.na(b)]<-0
  col_fun = circlize::colorRamp2(c(-2, 0, 2), c("steelblue", "lightgrey", "red"))
    Heatmap(as.matrix(b),cluster_rows = F,cluster_columns = F,column_names_rot = 0,
          row_title =x,name = "log2(Observed/Expected)",col=col_fun,
          column_title_gp = gpar(fontsize = 20),
          row_title_gp = gpar(fontsize = 20),
          column_names_gp = gpar(fontsize=20),
          row_names_gp = gpar(fontsize=20),
          heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                    grid_height=unit(2,"cm"),
                                    title_gp=gpar(fontsize=20),
                                    labels_gp=gpar(fontsize=20)))

    })
htl<-p2[[1]] %v% p2[[2]]
draw(htl,column_title="G1S bins",column_title_side="bottom")

#Cycling in relation to EpiSen

#Set cutoffs from permuted nonmalignant cells
#Plot nonmalignant epithelial cells
p1<-ggplot(sm[sm$CNA=="Normal" & sm$HPV=="HPV-",],aes(x=G1S,y=G2M))+geom_point()+
  theme_light(base_size = 20)+ggtitle("Noncancer Epithelial Cells")+
  geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  viridis::scale_fill_viridis()
p1

#Plot permuted values of nonmalignant epithelial
fm2 <- metadata %>% filter(CNA %in% c("Cancer","Normal")) %>% rownames() %>%
  filter.umi(matrix = merged)
fp<- sm %>% filter(CNA=="Normal" & HPV=="HPV-") %>% rownames() %>% fm2[,.]
fp2<- t(apply(fp,1,gtools::permute))
fp3<-type.score(fp2,pr[pr$Program %in% c("G1S","G2M"),],tv="Program",bin=NULL,centre = FALSE)$typescorelist
fp3<-do.call(cbind.data.frame,fp3)
p2<-ggplot(fp3,aes(x=G1S,y=G2M))+geom_point()+
  theme_light(base_size = 20)+geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+  ggtitle("Noncancer Epithelial Cells, permuted")+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  viridis::scale_fill_viridis()
p2

#Plot cancer cells with cutoffs
p3<-ggplot(sm[sm$CNA=="Cancer",],aes(G1S, G2M)) + geom_point()+
  stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  viridis::scale_fill_viridis()+
  theme_light(base_size = 20)+geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+ggtitle("Cancer Cells")
plot_grid(p1,p2,p3)

#Set cycling definitions based on cutoffs
sm$cycling<-ifelse(sm$G1S>0|sm$G2M>0,"Cycling","Noncycling")

#Define classes of EpiSen expression
sq<-quantile(sm$EpiSen,probs=c(0.2,0.8))
sm$senclass<-ifelse(sm$EpiSen<sq[1],"Low",ifelse(sm$EpiSen>sq[2],"High","Medium"))

sm22<-sm
sm22$comb<-paste(sm22$cycling,sm22$senclass,sep = "_")

#Repeat 100 resampling runs
rf<-lapply(seq(1:100),function(n){

pn<-unique(sm22$Patient)
pc<-c("HPVon","HPVoff","HPVneg")
pnf<-lapply(pn,function(x){
  pn2<-lapply(pc,function(y){
  a<-sm22 %>%filter(Patient==x & class==y) %>% nrow()
  sm22 %>%filter(Patient==x & class==y) %>% rownames() %>% sample(size=min(100,a))
  
  
})
})

sm2<-sm22[unlist(pnf),]

a<- sm2 %>% group_by(class,cycling) %>% count(senclass) %>% mutate(frac=n/sum(n)) %>%
  filter(cycling=="Noncycling" & senclass=="High")
sx<-cbind.data.frame(a$class,a$frac,n)
colnames(sx)<-c("class","ratio","run")
sx
})
sx<-do.call(rbind.data.frame,rf)
sx<- sx %>% group_by(class) %>% mutate(sd=sd(ratio)) %>% summarise_all(mean) %>% as.data.frame()

#Plot EpiSen-high in relation to cell cycle and HPV class
ggplot(sx[sx$class!="Normal",],aes(x=class,y=ratio,fill=class))+
  geom_col(position = "dodge")+
  theme_light(base_size = 20)+
  guides(fill=guide_legend(title="Subset"))+
  ylab("EpiSen-High noncycling\namong all noncycling cells") +xlab("")+
  scale_fill_manual(values=c("chartreuse","darkgreen","red"))+
  geom_errorbar(aes(ymin=ratio-sd,ymax=ratio+sd),width=0.1)

#Select n neighbours per cell and see how many are EpiSen-high
nf<-c("HPVon","HPVoff","HPVneg")
ff<-mclapply(nf,function(x){
  a<-sm[sm$class==x,]
  b<-FNN::get.knn(as.matrix(a[,c("G1S","G2M")]), k = 20)$nn.index
for(i in 1:nrow(a)){
  a1<- a[b[i,],] %>% count(senclass) %>% mutate(pct=n/sum(n)*100) %>% filter(senclass=="High") %>%
    as.data.frame()
  a[i,"senpct"]<-ifelse(nrow(a1)>=1, as.numeric(a1[a1$senclass=="High","pct"]),0)
                  
}
  a
})
sm3<-do.call(rbind.data.frame,ff)
#Plot
pg<-lapply(nf,function(z){
ggplot(sm3[sm3$class==z,],aes(x=G1S,y=G2M,colour=senpct))+geom_point(alpha=0.5)+
    theme_light(base_size = 20)+
  geom_vline(xintercept = 0,lty=2)+
  geom_hline(yintercept = 0,lty=2)+ggtitle(z)+
  scale_colour_gradientn(name="% EpiSen-High\nneighbours",colours=c("seagreen","gold","darkred"))
})
do.call(plot_grid,pg)

#Plot fractions of cycling and senescent cells per patient and HPV class
sm$comb<-paste(sm$cycling,sm$senclass,sep = "_")
sm$Subclone<-metadata[rownames(sm),"Subclone"]
sm$pat2<-sm$Patient
pl<-unique(sm$Patient)

cc<-lapply(pl,function(x){
  cd<-lapply(c("HPVon","HPVoff","HPVneg"),function(y){
  a<-sm[sm$class==y & sm$pat2==x,]
  
  b1a<-nrow(a[a$comb=="Noncycling_High",])/nrow(a[a$cycling=="Noncycling",])
  b2<-nrow(a[a$cycling=="Cycling",])/nrow(a)
  c(b1a,b2)
})
cd2<-do.call(rbind.data.frame,cd)
colnames(cd2)<-c("pSen","pCyc")
cd2$Patient<-x
cd2$class<-c("HPVon","HPVoff","HPVneg")

cd2
  })
cc2<-do.call(rbind.data.frame,cc)
#Plot
ggplot(cc2,aes(x=pSen,y=pCyc,colour=class,group=Patient))+geom_point()+
  geom_text_repel(aes(label=Patient,size=15),show.legend = FALSE)+
  theme_light(base_size = 30)+
  xlab("EpiSen-High noncycling\namong all noncycling cells")+ylab("Fraction cycling")+
  geom_line(colour="black")+
  scale_colour_manual(values = c("green","darkgreen","red") )+
  guides(colour=guide_legend(title="group",override.aes=list(size=4)))

#Plot fractions of cycling, EpiSen-high and HPV-positive cells per subclone
sm$Subclone<-metadata[rownames(sm),"Subclone"]

q<-sm %>% filter(CNA=="Cancer" & !Patient %in% c("OP5",hn)) %>%
  group_by(Patient,Subclone) %>%count(HPV) %>% mutate(HPVperc=1-(n/sum(n))) %>% 
  filter(HPV=="HPV-") %>% dplyr::select(Subclone,HPVperc) %>% as.data.frame()

q$cycling<-sm %>% filter(CNA=="Cancer" & !Patient %in% c("OP5",hn)) %>% 
  group_by(Patient,Subclone) %>%count(cycling) %>% mutate(frac=1-(n/sum(n))) %>% 
  filter(cycling=="Noncycling") %>% pull(frac)

qs<-sm %>% filter(CNA=="Cancer" & !Patient %in% c("OP5",hn)) %>% 
  group_by(Patient,Subclone) %>%count(senclass) %>% mutate(frac=n/sum(n)) %>%
  filter(senclass=="High") %>% dplyr::select(Subclone,frac)
q[q$Subclone %in% qs$Subclone,"senclass"]<-qs$frac
q[is.na(q$senclass),"senclass"]<-0

q<-q[grep("Unres",q$Subclone,invert = T),]

p<-unique(q$Patient)
qq<-mclapply(p,function(x){
  a<-cor(q[q$Patient==x,"HPVperc"],q[q$Patient==x,"cycling"])
  b<-max(q[q$Patient==x,"HPVperc"])-min(q[q$Patient==x,"HPVperc"])
  c(a,b)
})
qq<-do.call(rbind.data.frame,qq)
colnames(qq)<-c("cycling","range")
rownames(qq)<-p
head(qq)
qq$Patient<-rownames(qq)

colnames(q)[3:5]<-c("Fraction HPVon cells","Fraction cycling cells","Fraction EpiSen-hi cells")
qm<-melt(q)
pl2<-rownames(qq[qq$range>0.15,])
pp<-mclapply(pl2,function(x){
  a<-qm[qm$Patient==x & qm$variable!="G1S",]
  
  aa<-a[a$variable=="Fraction HPVon cells","value"]
  bb<-a[a$variable=="Fraction cycling cells","value"]
  try(print(paste0(x," ",cor.test(aa,bb)$p.value)))
  
  
  z<-arrange(a[a$variable=="Fraction HPVon cells",],desc(value)) %>% pull(Subclone)
  b<-ggplot(a,aes(x=factor(Subclone,levels=z),y=variable,fill=value))+geom_tile()+
    ggtitle(x)+theme_light(base_size = 20)+xlab("")+ylab("")+
    scale_fill_gradientn(name="value",breaks=c(0,0.4,0.8),
                           colours=c("seagreen","gold","darkred"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    guides(fill=FALSE)

  b
  
})
p11<-do.call(plot_grid,pp)

legend <- get_legend(
  b + theme(legend.box.margin = margin(0, 0, 0, 12))
)
p23<-plot_grid(p11, legend, rel_widths = c(3, .4))
p23
