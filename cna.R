###Assignment of malignant cells and subclones
###Files and data needed for analysis
# - source functions file
# - merged=Raw UMI matrix, found under GSE182227
# - metadata=Cell metadata table, found in supplement
# - hg19=hg19 genome, found in biomart
# - tcga=TCGA subtype signatures, found in supplement
# - programs=Cancer metaprograms, found in supplement
# - celltypes=cell type signatures found in supplement


#CNA calling
#Prepare CNA function inputs
HPV<-metadata$HPV
names(HPV)<-rownames(metadata)
HPV_score<-metadata$HPV_score
names(HPV_score)<-rownames(metadata)
tumour_score<-metadata$tumourscore
names(tumour_score)<-rownames(metadata)
gg<-NULL #Use patient-specific rather than global thresholds
acut=0.15 #Cutoff parameter for chromosomal event calling


pl<-unique(metadata$Patient)
sf<-lapply(pl,function(pat){
  #Select epithelial and reference cells among all patients
  ep<-metadata %>% filter(Type=="Epithelial" & Patient==pat) %>% rownames()
  refs<-c("Fibroblast","Myofibroblast","Lymphovascular","Endothelial")
  refcells<-lapply(unique(refs),function(x){
    metadata %>% filter(Type==x & Patient==pat) %>% rownames()
  })
  rn<-c(ep,unlist(refcells))
  #Create uncentred matrix with only epithelial and stromal cells
  m1<-filter.umi(matrix=merged,cells=rc(ep,rn),centre="none",log=TRUE,cut=8)
  
  #Create CNA matrix
  mat<-calc.cna(query=query,genome=hg19,matrix=m1,ref=refcells,startlog=TRUE,per_chr = TRUE,
                noise=0.15,mav_n=100,scale=NULL,top_genes=NULL,range=c(-3,3))
  
  #Split into subclones
  a2<-cluster.clone(mat,query,by_cell = FALSE,knn=TRUE,cluster_method = "louvain",louvain_k = 15,
                    cell_matrix = "numeric",merge = TRUE,armcut=10,name=pat,merge_method = "both",
                    expcut = 0.1,sdcut = 5,top_region = 2/3,clonesize = 15,adregion = 2/3,
                    adgenes = "top",dimred = "uwot",genome=hg19,
                    adcut=acut,dcut=acut)
  
#Assign malignant cells per subclone
  if(length(a2)>1){
    ad<-a2$ad
    epmat<-a2$epmat
    tcl<-sort(a2$tcl)
    dm<-a2$distmat
    
    a3<-reassign.clones(mat,query,ad,tcl,distmat=dm,em=epmat,tscore=tumour_score,
                        adcut=acut,adgenes="top",adregion=2/3,order_clones = TRUE,genome=hg19,
                        name=pat,refilt=FALSE,gt=NULL,goldstandard=gg,clonesize = 15)
    #
    ad<-a3$ad
    epmat<-a3$epmat
    tcl<-sort(a3$clones)
    thr<-a3$thresholds
    defmat<-a3$defmat
    
  }
  
  #Assign status if there is only one subclone
  if(length(a2)<2){
    ax<-define.cna(mat,query,unlist(refcells),tumour_score=tumour_score,top_cells = 1/4,
                   gs=gg,no_ts=TRUE,name=pat,print_ref = TRUE,extra_score = NULL)
    defmat<-ax$def_mat
    defmat<-defmat[,!colnames(defmat) %in% c("knn","prob")]
    thr<-ax$thresholds
    dmx<-defmat[defmat$Origin!="Reference",]
    tcl<-ifelse(dmx$cna3=="Cancer",paste0(pat,"_subclone_A"),
                ifelse(dmx$cna3=="Normal",paste0(pat,"_Normal"),
                       paste0(pat,"_Unresolved")))
    names(tcl)<-rownames(dmx)
    tcl<-sort(tcl)
    
    br1<-br.fun(genes=rownames(mat),separate = "arm")
    br_vec<-br1$br_vec
    labels<-br1$bdf$labels
    rm<-br1$rm
    
    if(length(rm)!=0){
      mat2<-mat[-rm,]}else{mat2<-mat}
    #
    
    ax2<-admat.fun(mat2,br_vec,tcl[tcl==paste0(pat,"_subclone_A")],
                   labels,cut=acut,genes="top",region=2/3)
    ad<-ax2$ad
    epmat<-ax2$epmat
  }
  
  #Plot CNA matrix
  plot<-invisible(plot.cna2(matrix=mat, cells=names(tcl),separate="chr",order=FALSE,
                      cluster_rows=FALSE,name=pat,row_split=tcl,genome = hg19,
                      column_title_gp = gpar(fontsize = 20),column_title_rot=90,
                      row_title_gp = gpar(fontsize = 20),row_title_rot=0,
                      heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                                grid_height=unit(2,"cm"),
                                                title_gp=gpar(fontsize=20),
                                                labels_gp=gpar(fontsize=20))))
  
  print(paste0(pat," done!"))
  #Outputs
  out<-list(clones=tcl,epmat=epmat,defmat=defmat,admat=ad,thresholds=thr,mat=mat,plot=plot)
  out
})

#Save CNA data
names(sf)<-pl
saveRDS(sf,"clonedata.RDS")

#Prepare clone data
sf<-readRDS("clonedata.RDS")
clonedata<-do.call(rbind,sf)
clones<-unlist(clonedata[,1])
epmat<-do.call(cbind,clonedata[,2])
defmat<-do.call(rbind,clonedata[,3])


#Remove naming artifacts
names(clones)<-sapply(strsplit(names(clones),split = ".",fixed = TRUE),function(x)(x[2]))
rownames(defmat)<-sapply(strsplit(rownames(defmat),split = ".",fixed = TRUE),function(x)(x[2]))
colnames(epmat)<-sapply(strsplit(colnames(epmat),split = ".",fixed = TRUE),function(x)(x[2]))

#Set CNA metadata
CNA<-defmat$cna3
names(CNA)<-rownames(defmat)
metadata[names(CNA),"CNA"]<-ifelse(CNA!="Reference",CNA,"Nonepithelial")
metadata$CNA[is.na(metadata$CNA)]<-"Nonepithelial"
 
#Set subclone metadata
metadata[names(clones),"Subclone"]<-clones
metadata$Subclone[is.na(metadata$Subclone)]<-"Nonepithelial"

##########5. Summary CNA plot
fm<- metadata %>% filter(CNA=="Cancer") %>% rownames() %>% filter.umi(matrix=merged)
ssx<-ss1[ss1$Gene %in% rownames(fm),]
ttx<-type.score(fm,ssx,tv="Subtype",bin=100)
Basal<-ttx$typescorelist$Basal
Classical<-ttx$typescorelist$Classical
Atypical<-ttx$typescorelist$Atypical

#Tumour and HPV score per clone
clones<-metadata$Subclone
names(clones)<-rownames(metadata)
clones<-clones[clones %in% colnames(epmat)]
scores<-lapply(colnames(epmat),function(x){
  a<-names(clones[clones==x])
  b<-mean(tumour_score[a])
  c<-HPV[a]
  d<-length(c[c=="HPV+"])/length(a)
  e<-mean(Atypical[a])
  f<-mean(Basal[a])
  g<-mean(Classical[a])
  c(b,d,e,f,g)
  
})
scores<-do.call(rbind.data.frame,scores)
rownames(scores)<-colnames(epmat)
colnames(scores)<-c("HPVpos Tumour score","%HPVpos","Atypical","Basal","Classical")

#Scores and categories for ordering
sx<-scores
sx$clone<-rownames(sx)
sx$pat<-sapply(strsplit(rownames(sx), split='_', fixed=TRUE),function(x)(x[1]))
sx$hpv<-ifelse(sx$pat %in% hn, "HPV-","HPV+")

pt<-table(sx$pat)
for(i in 1:length(unique(names(pt)))){
  sx[sx$pat==names(pt[i]),"no"]<-pt[i]
}

#Bypass bug where 10 is sorted after 1
sx[sx$no==10,"no"]<-9

#Order
sx$code<-paste0(sx$no,sx$pat)
sx2<- sx %>% arrange(desc(hpv),desc(code),desc(sx$`%HPVpos`))

#Order plotting matrix
em2<-epmat[,rownames(sx2)]
pat1<-sx2$pat
names(pat1)<-rownames(sx2)
scores<-scores[rownames(sx2),]

#Plot
col_fun = circlize::colorRamp2(c(-1, 0, 1), c("steelblue", "white", "darkred"))
Heatmap(t(em2),cluster_columns = F,cluster_rows = F,col=col_fun,show_row_names = F,
        right_annotation = rowAnnotation(df=scores,show_legend=FALSE,
                                         annotation_name_gp=gpar(fontsize=20), simple_anno_size = unit(1, "cm")),
        name = "Inferred Copy Number\nLog-ratio",row_title_rot = 0,
        row_split = factor(pat1, levels = unique(pat1)), border = TRUE,row_gap = unit(0,"mm"),
        row_title_gp = gpar(fontsize=12),column_names_gp = gpar(fontsize=20),
        heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                  grid_height=unit(1,"cm"),
                                  title_gp=gpar(fontsize=15),
                                  labels_gp=gpar(fontsize=15)))



#DE genes for malignant cells in margin of OP34
#Set data
meta2<-metadata[metadata$Patient=="OP34" & metadata$Type=="Epithelial"& metadata$CNA!="Unresolved",]
class<-ifelse(meta2$CNA=="Normal","Normal Epithelial Cells",
              ifelse(meta2$CNA=="Cancer" & meta2$Sample=="OP34norm","Invasive Cancer Cells",
                     "Core Cancer Cells"))
names(class)<-rownames(meta2)

fg<-filter.umi(merged,names(class))
#Find DE genes
deg2<-uniqueDEG(merged,class,30000,fc=1)


cu<-c("Normal Epithelial Cells","Core Cancer Cells","Invasive Cancer Cells")
dd<-unlist(lapply(cu,function(x){
  deg2 %>% filter(Cluster==x) %>% pull(Gene) %>% head(n=50)
}))
dd<-c(dd[1:90],"CDKN2A",dd[91:140])



#Arrange cells based on their signature scores
n4<-unlist(lapply(cu,function(x){
  a<- deg2 %>% filter(Cluster==x) %>% pull(Gene)
  b<-fg[,names(class[class==x])] %>% score(genes=a) %>%sort(decreasing = T) %>%head(n=29) %>%names()
  c<-hclustreorder(cor(as.matrix(fg[dd,b]))) %>% colnames()
  
})
)
cl2<-class[n4]
#OP34 heatmap
col_fun = circlize::colorRamp2(c(-5, 0, 5), c("steelblue", "white", "darkred"))
rl<-c("SOX2","SOX4","CDKN2A","SGK1","E5","E7",
      "APOBEC3B","COL17A1","MARCKSL1","KRT8","KRT14","IL1RN","CD74",
      "IL1RN","KRT4","SPRR3","RHCG","S100A9","AQP3")

rl2<-intersect(dd,rl)
length(rl2)
hr = rowAnnotation(q=anno_mark(at = which(dd %in% rl2), labels =rl2),
                   annotation_legend_param=list(labels_gp=gpar(fontsize=20)))
ff<-filter.umi(merged,n4,whitelist=dd)
h1<-Heatmap(as.matrix(fg[dd,n4]),show_row_names = F,show_column_names = F,
            top_annotation = columnAnnotation(Origin=cl2,
                                              col=list(
                                                Origin=c("Core Cancer Cells"="dodgerblue",
                                                         "Invasive Cancer Cells"="orange",
                                                         "Normal Epithelial Cells"="darkgreen")),
                                              annotation_legend_param=list(
                                                grid_width=unit(1, "cm"),
                                                grid_height=unit(2,"cm"),
                                                title_gp=gpar(fontsize=20),
                                                labels_gp=gpar(fontsize=20))
            ),
            right_annotation = hr,
            heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                      grid_height=unit(2,"cm"),
                                      title_gp=gpar(fontsize=20),
                                      labels_gp=gpar(fontsize=20)),
            cluster_rows = F,cluster_columns = F,col=col_fun,name = "log2(CPM)")
h1


#OP34 CNA plot
pat<-"OP34"
meta2<-metadata[metadata$Patient==pat & metadata$CNA %in% c("Normal","Cancer"),]
su<-unique(meta2$Subclone)
sff<-unlist(mclapply(su,function(x){
  a<-rownames(meta2[meta2$Subclone==x,])
  sample(a,23)
}))

meta2<-meta2[sff,]


mat<-sf[[pat]]$mat[,rownames(meta2)]
tcl<-meta2$Subclone
names(tcl)<-rownames(meta2)
tcl<-sort(tcl)
norm<-tcl[grep("norm",names(tcl))]
Origin<-ifelse(names(tcl) %in% names(norm),"Margin","Core")
names(Origin)<-names(tcl)


invisible(plot.cna2(matrix=mat,cells=names(tcl),separate="chr",order=FALSE,
                    cluster_rows=FALSE,name=pat,row_split=tcl,genome=hg19,
                    column_title_gp = gpar(fontsize = 20),column_title_rot=90,
                    row_title_gp = gpar(fontsize = 20),row_title_rot=0,
                    right_annotation=rowAnnotation(Origin=Origin[names(tcl)],
                                                   col=list(
                                                     Origin=c("Core"="dodgerblue","Margin"="orange")),
                                                   
                                                   annotation_name_gp=gpar(fontsize=20), 
                                                   simple_anno_size = unit(1, "cm"),
                                                   annotation_legend_param=list(
                                                     grid_width=unit(1, "cm"),
                                                     grid_height=unit(2,"cm"),
                                                     title_gp=gpar(fontsize=20),
                                                     labels_gp=gpar(fontsize=20))),
                    heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                              grid_height=unit(2,"cm"),
                                              title_gp=gpar(fontsize=20),
                                              labels_gp=gpar(fontsize=20))))

#OP9 mutations - change plot.cna2 to plot both sides for op9

#Load data
pat<-"OP9"
mat<-sf[[pat]]$mat

tcl<-metadata[metadata$Patient==pat & metadata$CNA=="Cancer","Subclone"]
names(tcl)<-rownames(metadata[metadata$Patient==pat & metadata$CNA=="Cancer",])
#Keep only cancer
tc2<-tcl[grep("Unres|Norm",tcl,invert = T)]
#Compare major subclones only
tc2[grep("A|B",tc2)]<-"OP9_subclone_A/B"
tc2[grep("A|B",tc2,invert = TRUE)]<-"OP9_subclone_C/D/E/F"

#HPV% in subclones
z<-metadata[names(tc2),]
z$tc2<-tc2
z %>% group_by(tc2) %>%count(HPV) %>% mutate(f=n/sum(n)) %>% filter(HPV=="HPV+")
chisq.test(z$tc2,z$HPV)

#Cell names for subsetting
#Reorder cells per clone randomly
tn<-unique(tc2)
tf<-lapply(tn,function(x){
  a<-tc2[tc2==x]
  b<-sample(names(a),length(a),replace = F)
  c<-a[b]
})
tc2<-unlist(tf)

query<-names(tc2)

#Function for creating variant matrices - VARTRIX OUTPUT NEEDED FOR THIS
vf<-lapply(c("mutect","germline"),function(y){
  a<-list.files(pattern = paste0(y,".mtx"))
  a<-a[grep(pat,a)]
  a<-a[grep("comb",a,invert = TRUE)]
  b<-sub(paste0(".",y,".mtx"),"",a)
  c<-lapply(b,function(x){
    m1<-readMM(paste0(x,".",y,".mtx"))
    v1<-read.table(paste0(x,".",y,".txt"))
    b1<-read.table(paste0("~/processed/OPSCC/",x,"/filtered_feature_bc_matrix/barcodes.tsv"))
    colnames(m1)<-paste0(x,"_",b1[,1])
    rownames(m1)<-v1[,1]
    m1
  })
  vm<-do.call(cbind,c)
  
  #Keep only query cells
  mdf<-as.matrix(vm)
  n2<-intersect(colnames(mdf),query)
  md2<-mdf[,n2]
  
  #Keep only variants present in at 5% OR 10 cells
  s1<-apply(md2,1,function(x)sum(x>1))
  if(sum(s1>10)>50){varcut<- 0.05*ncol(md2)}else{varcut<-10}
  md3<-md2[s1>=varcut,]
  
  #Set value to 0 if no evidence of alt allele, 1 if any evidence
  md3[md3==1]<- 0
  md3[md3>1]<- 1
  md3
})

#Unpack matrix
v1<-vf[[1]]
d1<-as.data.frame(t(v1[rowSums(v1)>19,]))
v2<-vf[[2]]
rn<-sample(rownames(v2),ncol(d1),replace = F)
d2<-as.data.frame(t(v2[rn,]))

#Plot variants

invisible(plot.cna2(matrix=mat,cells=names(tc2),separate="chr",order=FALSE,
                    cluster_rows=FALSE,name=pat,row_split=tc2,genome=hg19,
                    column_title_gp = gpar(fontsize = 20),
                    left_annotation = rowAnnotation(df=d1[names(tc2),],
                                                    show_legend=FALSE, annotation_name_gp=gpar(fontsize=20), 
                                                    simple_anno_size = unit(1, "cm")),
                    right_annotation=rowAnnotation(df=d2[names(tc2),],
                                                   show_legend=FALSE, annotation_name_gp=gpar(fontsize=20), 
                                                   simple_anno_size = unit(1, "cm")),
                    row_title_gp = gpar(fontsize = 20),row_title_rot=0,column_title_rot=90,
                    heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                              grid_height=unit(2,"cm"),
                                              title_gp=gpar(fontsize=20),
                                              labels_gp=gpar(fontsize=20))))


#Compare pEMT, OP9-Mes and Fibroblast scores

#Cells to use
n1<-rownames(metadata[metadata$Type=="Fibroblast" | metadata$Type=="Myofibroblast",])
n2<-rownames(metadata[metadata$Patient!="OP9" & metadata$CNA=="Cancer",])
n3<-names(tc2[grep("A/B",tc2,invert = TRUE)])
n4<-names(tc2[grep("A/B",tc2,invert = FALSE)])



#Subsample non-OP9
n1<-sample(n1,300,replace = F)
n2<-sample(n2,300,replace = F)
n3<-sample(n3,300,replace = F)

#Genes to use
ss1<-tcga
ss1<-ss1[ss1$Subtype=="Mesenchymal","Gene"]
ss2<-programs
ss2<-ss2[ss2$Program=="pEMT","Gene"]
ss3<-celltypes
ss3<-ss3[ss3$Celltype=="Fibroblast","Gene"]

#New matrix
fq<-filter.umi(merged,c(n1,n2,n3,n4),whitelist = c(ss1,ss2,ss3))
f9<- rownames(filter.umi(matrix=merged,cells=names(tc2[tc2=="OP9_subclone_C/D/E/F"])))
o9m<-intersect(ss1,f9)
#Scores
o9<-score(fq,o9m)
mes<-score(fq,ss1)
pemt<-score(fq,ss2)
fibro<-score(fq,ss3)

#Plot score matrix
d1<-cbind.data.frame(o9,mes,pemt,fibro)
d1<-d1[c(n1,n2,n3,n4),]
d1$class<-ifelse(rownames(d1) %in% n1, "Fibroblasts",ifelse(rownames(d1) %in% n2,
                                                            "Other cancer cells",
                                                            ifelse(rownames(d1) %in% n3, "OP9 C/D/E/F", "OP9 A/B")))
d1<- arrange(d1,class,o9,mes,fibro,pemt)
d1<-arrange(d1,class)
colnames(d1)<-c("OP9_Mes","Mesenchymal","pEMT","Fibroblast","Cells")
class<-d1$Cells

#T-test
cn<-unique(d1$Cells[d1$Cells!="OP9 C/D/E/F"])
tt<-unlist(lapply(cn,function(x){
  a<-t.test(d1[d1$Cells=="OP9 C/D/E/F",1],d1[d1$Cells==x,1])$p.value
  b<-t.test(d1[d1$Cells=="OP9 C/D/E/F",2],d1[d1$Cells==x,2])$p.value
  c(a,b)
}))
p.adjust(tt)
#Plot
dx<-melt(d1)
dx<-dx[grep("Mes",dx$variable),]
ggplot(dx,aes(x=variable,fill=Cells,y=value))+geom_violin()+theme_light(base_size = 30)+xlab("")
