###Preprocessing and cell type annotation
###Files and data needed for analysis
# - source functions file
# - filtered.samplelist=list of cellranger output folders
# - celltypes=cell type signatures found in supplement
#Merge all matrices and name cells after sample of origin
matrices<-lapply(filtered.samplelist,function(x){
  a1<-sub("./","/",x)
  m1<-Matrix::readMM(paste0(getwd(),a1,"/matrix.mtx"))
  b1<-read.table(paste0(getwd(),a1,"/barcodes.tsv"))
  g1<-read.table(paste0(getwd(),a1,"/genes.tsv"))
  colnames(m1)<-b1[,1]
  rownames(m1)<-g1[,2]
  matrix<-m1
  samplename<-unlist(strsplit(x,"/"))[2] #Get sample name from list of samples
  colnames(matrix)<- paste(samplename, colnames(matrix), sep = "_") #Add samplename to cell names
  genespercell<-apply(matrix,2,function(x)sum(x>0)) #Genes expressed per cell
  cellthreshold<-genespercell>=1000 #|genespercell>=quantile(genespercell,0.4)#Local and global complexity threshold
  matrix<-matrix[,cellthreshold] #Filter cells by thresholds
  
  #Doublet identification
  d1<-find.doublets(matrix,method="scdblfinder",scale = 0.006)
  d2<-find.doublets(matrix,method="scds",scale = 0.006)
  d3<-find.doublets(matrix,method="scran",scale=0.006)
  rownames(matrix)<-sub("HPV16","",rownames(matrix))
  if(nrow(matrix)>33701){matrix<-matrix[-c(33702:33727),]}
  doubles<-cbind.data.frame(d1,d2,d3)
  out<-list(matrix=matrix,doubles=doubles)
  out
})
#
print("All single matrices loaded")

#Bind matrices and doublet data
mm<-do.call(rbind,matrices)
dbl<-do.call(rbind,mm[,2])
merged<-do.call(cbind,mm[,1]) #Merge the loaded matrices

#Remove doublets
dbl<-dbl[,c(2,4,6)]
colnames(dbl)<-c("scdblfinder","scds","scran")
dbl<-as.data.frame(ifelse(dbl=="doublet",1,0))
dbl$sum<-apply(dbl,1,sum)
dbl$consensus<-ifelse(dbl$sum>1,"doublet","singlet")

#Remove doublets, filter, centre, logTPM-transform matrix
merged<-merged[,dbl$consensus=="singlet"]
final.matrix<-filter.umi(merged,colnames(merged),cut=8)

#Set metadata
HPVgenes<-as.character(c("E1","E2","E5","E6","E7","L1","L2")) #All HPV genes
HPV_score<-apply(merged[HPVgenes,],2,sum) #Sum of HPV UMIs
HPV<-ifelse(HPV_score>0,"HPV+","HPV-") #Binary score if >0 HPV UMIs
Sample <- sapply(strsplit(colnames(merged), split='_', fixed=TRUE),function(x)(x[1]))
names(Sample)<-colnames(merged)
rm.string <- c("CD45P", "norm", "BOT", "LN", "GEX")
unlist.rmstring <- paste(unlist(rm.string), collapse = "|")
Patient<-gsub(unlist.rmstring, replacement = "", x = Sample)

#Set cell type on a median-centred matrix
f2<-filter.umi(merged,colnames(merged),centre = "median",cut=8)
tt<-type.score(f2,celltypes,conf_int = 0.15,bin=100)
Type<-tt$type

#Save metadata
metadata<-cbind.data.frame(Sample,Patient,HPV,HPV_score,Type)

#Dimension reduction
dt2<-standard.umap(final.matrix,n_neighbors = 300,spread=5,min_dist=0.0001)


im<-c("T-cell","B-cell","Macrophage","Plasma cell","NK-cell","Mast cell","cDC","pDC")
metadata$t2<-ifelse(metadata$Type %in% im,"Immune",metadata$Type)


#Plot type, patient, HPV status
l1<-c("Patient","t2","HPV")
ff<-lapply(l1,function(x){
  ggplot(dt2, aes(x=V1, y=V2, colour=metadata[rownames(dt2),x])) +
    geom_point(size=1,alpha=0.2) + 
    guides(colour = guide_legend(title=as.character(x),override.aes = list(size=4,alpha=1)))+
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_light(base_size=20)
})
do.call(plot_grid,ff)

#Clustering
Cluster<-cluster.coord(dt2,method = "louvain",louvain = 200)
metadata$Cluster<-Cluster
#Plot clusters
dt3<-cbind.data.frame(dt2,Cluster)
clustcent = dt3 %>% group_by(Cluster) %>% dplyr::select(V1,V2) %>% summarize_all(mean) 

ggplot(dt2, aes(x=V1, y=V2, colour=metadata[rownames(dt2),"Cluster"])) +
  geom_point(size=1,alpha=0.5) + 
  guides(colour=FALSE,size=FALSE)+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=20)  +
  geom_text_repel(size=8,aes(x=V1,y=V2,label = Cluster),
                  inherit.aes=F,data = clustcent,max.overlaps = 30)

#DE genes per cluster
dg<-DE.genes(final.matrix,Cluster)

#Assign type by combining clustering and signature scores

t1<-do.call(cbind.data.frame,tt$typescorelist)
t2<-tt$type

#Plot cell type fraction by cluster
tp<-type.perc(final.matrix,cv=Cluster,tv=t2)
tp2<-tp[,c("Epithelial","Fibroblast","Myofibroblast","Endothelial","Lymphovascular",
           "Tcells","NK","Bcells","Plasma cell","Macrophage","DC","Mastcells")]
Heatmap(as.matrix(tp2),cluster_columns = F)

#Set cell type groups
imm<-c("Tcells","NK","Bcells","Plasma cell","Macrophage","DC","Mastcells")
stroma<-c("Fibroblast","Myofibroblast","Endothelial","Lymphovascular")

#Calculate main and secondary type scores + difference for each cell
t3<-matrix(nrow=nrow(t1),ncol=3)
rownames(t3)<-rownames(t1)
colnames(t3)<-c("Main","Secondary","Diff")
for(i in 1:nrow(t3)){
  
  a<- names(which.max(t1[i,]))
  b<- names(which.max(t1[i,grep(a,colnames(t1),invert=TRUE)]))
  c<- max(t1[i,])
  d<-max(t1[i,grep(a,colnames(t1),invert=TRUE)])
  e<- c- (d+0.15*c)
  t3[i,]<-c(a,b,e)
}

#Assign type to each cluster by the highest cell type fraction per cluster
t4<-cbind.data.frame(t1,t3)
t4$HPV<-HPV
t4$cluster<-Cluster
ft<-lapply(unique(Cluster),function(x){
  a<-names(which.max(tp2[x,]))
  b<-cbind.data.frame(names(Cluster[Cluster==x]),x,a)
  rownames(b)<-b[,1]
  b<-b[,2:3]
  colnames(b)<-c("cluster","type")
  b
})
ftt<-do.call(rbind,ft)
ftt<-ftt[rownames(t4),]
t4$clustertype<-ftt$type

#Set cells with discord between individual score and cluster as unresolved
t4$def<-ifelse(t4$HPV=="HPV+" & t4$clustertype !="Epithelial","Unresolved",
               ifelse(t4$clustertype=="Epithelial" & !t4$Main %in% c("Epithelial","Fibroblast","Myofibroblast"),
                      "Unresolved",ifelse(t4$cluster=="cluster0","Unresolved",
                                          ifelse(t4$clustertype=="Epithelial" & !t4$Secondary %in% c("Fibroblast","Myofibroblast") &
                                                   t4$Main=="Epithelial" &t4$Diff<=0,"Unresolved",
                                                 ifelse(t4$clustertype %in% stroma & 
                                                          !t4$Main %in% stroma,"Unresolved", 
                                                        ifelse(t4$clustertype %in% stroma & 
                                                                 t4$Diff<=0 &
                                                                 !t4$Secondary %in% stroma & 
                                                                 t4$Main %in% stroma,"Unresolved",
                                                               ifelse(t4$clustertype %in% imm & !t4$Main %in% imm,"Unresolved",
                                                                      ifelse(t4$clustertype %in% imm & !t4$Secondary %in% imm & t4$Main %in% imm
                                                                             & t4$Diff <=0,"Unresolved",t4$Main))))))))
t4$d2<-ifelse(t4$def!="Unresolved",t4$clustertype,"Unresolved")

#Plot reassigned cell types
Type<-t4$d2
metadata$Type<-Type
names(Type)<-rownames(t4)
n1<-length(unique(Type))-1

ggplot(dt2, aes(x=V1, y=V2, colour=Type)) +
  geom_point(size=1,alpha=0.2) + 
  guides(colour = guide_legend(override.aes = list(size=4,alpha=1)))+
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=30)

# Separately cluster and assign immune cells
Type<-metadata$Type #metadata vector not done yet at this point
names(Type)<-rownames(metadata)
Patient<-metadata$Patient
names(Patient)<-rownames(metadata)
immune<-c("Tcells","Bcells","NK","Plasma cell","pDC","cDC","Macrophage","Mastcells","DC")

#New matrix for just immune
n1<- names(Type[Type %in% immune])
mi<-filter.umi(merged,n1)
di<-standard.umap(mi)

#New score
t2<-type.score(mi,
               celltypes[celltypes$Gene %in% rownames(mi) & celltypes$Celltype %in% immune,],
               conf_int = 0.15,bin=100)
t3<-t2$type

mi2<-batch.correct(mi,c("OP4","OP6","OP14","OP9"),Patient,t3)
di2<-standard.umap(mi2,var_genes = 2000)

cl<-cluster.coord(di2,method = "louvain",louvain=100)
#Plot cluster purity
tp<-type.perc(mi2,cv=cl,tv=t3)
Heatmap(tp[,1:8],cluster_columns = F)
#DE genes
dq<-DE.genes(mi2,cl)
dq<-arrange(dq,Cluster)
c1<-rownames(tp)
for(i in 1:length (c1)){
  b<-c1[i]
  dq[dq$Cluster==b,"Assignment"]<-
    names(which.max(tp[b,1:8]))
}

#
clt<-cl
clu<-unique(clt)
for(i in 1:length(clu)){
  clt[clt==clu[i]]<-dq[dq$Cluster==clu[i],"Assignment"][1]
}
#Plot cluster type
l1<-list(Patient=Patient,Type=clt)
ff<-lapply(names(l1)[1:2],function(x){
  
  ggplot(di2, aes(x=V1, y=V2, colour=l1[[x]][rownames(di2)])) +
    geom_point(size=1,alpha=0.2) + 
    guides(colour = guide_legend(title=as.character(x),override.aes = list(size=4,alpha=1)))+
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_light(base_size=30)
  
})
do.call(plot_grid,ff)

#Add new immune cell types to metadata
metadata[names(clt),"Type"]<-clt





#Plot marker genes
#Gene order
l1<-c("CDKN2A","KRT5","KRT19","COL1A1","MMP2","ACTA2","MYL9",
      "VWF","ENG","CCL21","TFF3","PTPRC")
gene.plot(dt2,l1,final.matrix)
#Plot immune marker genes
l2<-c("CD3D","CD7","KLRD1","NKG7","CD19","CD79A","JCHAIN","MZB1","CD14","AIF1",
      "LY75","LILRA4","CLEC4C","MS4A2","TPSB2")
gene.plot(di2,l2,final.matrix[,rownames(di2)])

##Plot accuracy of cell type assignment
#Split epithelial into HPV+ and HPV-
ttr<-ifelse(metadata$Type=="Epithelial" & metadata$HPV=="HPV+","EpHPV+",
            ifelse(metadata$Type=="Epithelial" & metadata$HPV=="HPV-",
                   "EpHPV-",metadata$Type))
names(ttr)<-rownames(metadata)

#Set order for plotting
l3<-c("EpHPV+","EpHPV-","Fibroblast","Myofibroblast","Endothelial","Lymphovascular",
      "T-cell","NK-cell","B-cell","Plasma cell","Macrophage","cDC","pDC","Mast cell")

dot.plot(c(l1,l2),final.matrix,merged,ttr,l3)

#Create tumor score

#DE genes for each patient with core and margin samples
pp<-c("OP33","OP34","OP35")
sf<-mclapply(pp,mc.cores=3,function(x){
  a<- metadata[metadata$Patient==x & metadata$Type=="Epithelial",]
  b<- a[grep("norm",rownames(a)),]
  a1<-setdiff(rownames(a),rownames(b))
  b1<-rownames(b[b$HPV=="HPV-",])
  c<-c(a1,b1)
  d<-filter.umi(merged,c)
  e<-ifelse(colnames(d) %in% a1,"Tumour","Normal")
  f<-DE.genes(d,e,100)
  f$Patient<-x
  f
})
sf2<-do.call(rbind,sf)

#Take common DE genes and export
mf<-lapply(c("Tumour","Normal"),function(x){
  a<-sf2[sf2$Cluster==x,]
  b<-extract.metaprogs(a,a$Patient,x,TRUE,100)
  b$topgenes
})
mf22<-do.call(rbind,mf)
mf2<-mf22[mf22$Number>1,]

#Score cells by substracting normal gene score from tumour gene score
tx<-score(final.matrix,mf2[mf2$Cluster=="Tumour","Gene"],bin = 100) - 
  score(final.matrix,mf2[mf2$Cluster=="Normal","Gene"],bin = 100)
#Set negative values to 0 and plot
tx2<-ifelse(tx<=0,0,tx)
metadata[names(tx2),"tumourscore"]<-tx2


ggplot(dt2, aes(x=V1, y=V2, colour=metadata[rownames(dt2),"tumourscore"])) +
  geom_point(size=1,alpha=0.2) + 
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size=30)+ 
  scale_colour_gradientn(name="HPVpos\nTumor score",colours=c("seagreen","gold","darkred"),
                         labels=c(0,2.5,5),breaks=c(0,2.5,5))
