###Files and data needed for analysis
# - source functions file
# - merged=Raw UMI matrix, found under GSE182227
# - metadata=Cell metadata table, found in supplement
# It is suggested to, if possible, run the NMF function for each patient in parallel
# Here, example code is provided for running everything within an R environment
# The example uses cancer cells. For any other cell type, replace CNA=="Cancer" with Type==x,
# where x is the cell type

#List of patients
patlist<-metadata%>%filter(CNA=="Cancer")%>%count(Patient)%>%filter(n>30)%>%pull(Patient)
library(parallel)
library(NMF)
#Apply function to all
meta<-mclapply(patlist,function(pat){
  metadata %>% filter(CNA=="Cancer" & Patient==pat) %>% rownames() %>%
    filter.umi(matrix=merged)  %>%    
    meta.nmf(k=10,name=pat,plot=TRUE,size=5,method="snmf/r")
})
qq<-do.call(rbind,meta) #Save all data from NMF runs
s<-unlist(qq[,1]) #Character vector assigning every cell to an NMF cluster
m2<-do.call(rbind.data.frame,qq[,2]) #DE genes for each cluster

zq<-jac.mat(m2,vector = m2$Cluster,add_diagonal = TRUE,jacmin=0.2) #Jaccard similarity between clusters
jm<-jac.heatmap(zq) #Transform jaccard values to jaccard correlation heatmap

#Plot metaclusters
col_fun = circlize::colorRamp2(c(0, 0.05, 0.5), c("steelblue", "white", "darkred"))
col_dend = as.dendrogram(hclust(dist(1-jm,method = "euclidean"),
                                method = "average"))
n=8 #Number of metaclusters, change manually depending on how the clustering looks
h2<-Heatmap(jm,col=col_fun,show_row_names = F,show_column_names = F,
            name="Jaccard\nCorrelation",cluster_columns = col_dend,
            cluster_rows=col_dend,column_split=n,row_split =n,
            column_dend_height = unit(3,"mm"),row_dend_width = unit(3,"mm"),
            heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                      grid_height=unit(2,"cm"),
                                      title_gp=gpar(fontsize=20),
                                      labels_gp=gpar(fontsize=20)))
h2
#Extract top recurrent genes/metacluster
co<-column_order(h2)
for(i in 1:length(co)){
  names(co)[i]<-paste0("Metacluster_",i)}

mp<-mclapply(names(co),mc.cores=length(co),function(x){
  nq<-colnames(jm)[co[[x]]]
  mq<-m2[m2$Cluster %in% nq,]
  zzd<-extract.metaprogs(mq,mq$Cluster,as.character(x),clusters = T,no=100)
  zzd$topgenes
})
mp2<-do.call(rbind,mp) #Output table

#This will be the metaprogram file, print it and annotate manually
mp2<-mp2[mp2$Patfreq>0.5,]


#Score cells for expression of metaprograms
a<-mp2
a<-a[grep("LowQual",a$Program,invert = T),]
a<-a[!is.na(a$Program) & a$Patfreq>0.5,2:8]
a<- arrange(a,desc(Patfreq,desc(logFC))) %>% filter(!duplicated(Gene)) %>%
  arrange(Program)
table(a$Program)

a<- arrange(a,desc(Patfreq,desc(logFC))) %>% filter(!duplicated(Gene)) %>%
  arrange(Program)

d<-filter.umi(merged,rownames(metadata[metadata$CNA=="Cancer",]),whitelist = a$Gene,centre = "none",sparse = T)

a<-a[a$Gene %in% rownames(d),]
a2<-a$Program
names(a2)<-a$Gene
f<-d[names(a2),]
#Annotate genes
an<-unique(a$Program)
tn<-unlist(mclapply(an,function(x){
  a %>% filter(Program==x) %>% arrange(desc(Patfreq),desc(logFC)) %>%
    head(n=round(50/length(an))) %>%pull(Gene)
}))

hr<- rowAnnotation(q=anno_mark(at = which(rownames(f) %in% tn), labels =tn),
                   annotation_legend_param=list(
                     labels_gp=gpar(fontsize=20)))

#Score all cells for programs
z<-type.score(d,a,tv="Program",conf_int = 0.15)
z1<-do.call(cbind.data.frame,z$typescorelist)
#Arrange cells by expression of scores
z1 <- z1[sample(nrow(z1)),]
for(i in 1:nrow(z1)){z1[i,"max"]<-names(which.max(z1[i,]))}
h2<-arrange(z1,max)
#Add subtype to metadata
d2<-z$type
d3<-d2[d2!="Unresolved"]
metadata[names(d3),"Subtype"]<-d3


p2 <- metadata[rownames(h2),"Patient"]
ha<-columnAnnotation(annotation_name_gp=gpar(fontsize=20), simple_anno_size = unit(1, "cm"),
                     annotation_legend_param=list(grid_width=unit(1, "cm"),
                                                  grid_height=unit(1,"cm"),
                                                  title_gp=gpar(fontsize=20),
                                                  labels_gp=gpar(fontsize=20)),
                     
                     Patient=p2,col=list(Patient=c("OP13"="gold",
                                                   "OP14"=pal2[2],
                                                   "OP5"=pal2[3],
                                                   "OP6"=pal2[4],
                                                   "OP8"=pal2[25],
                                                   "OP9"="lightgreen",
                                                   "OP4"=pal2[7],
                                                   "OP33"=pal2[8],
                                                   "OP34"=pal2[9],
                                                   "OP35"=pal2[10],
                                                   "OP10"=pal2[20],
                                                   "OP12"=pal2[17],
                                                   "OP16"=pal2[33],
                                                   "OP17"=pal2[28],
                                                   "OP19"=pal2[15],
                                                   "OP20"=pal2[19])))


#Plot of all cells
col_fun = circlize::colorRamp2(c(-4, 0, 4), c("steelblue", "white", "darkred"))
p1<-Heatmap(as.matrix(f[,rownames(h2)]),col=col_fun,show_row_names=F,show_column_names=F,
            cluster_rows=F,cluster_columns=F,right_annotation = hr,
            row_split = a2,row_title_gp = gpar(fontsize = 20),row_title_rot=0,
            name="log2(TPM)",top_annotation = ha,
            column_title_gp = gpar(fontsize = 30),
            heatmap_legend_param=list(grid_width=unit(1, "cm"),
                                      grid_height=unit(2,"cm"),
                                      title_gp=gpar(fontsize=20),
                                      labels_gp=gpar(fontsize=20)))
          
p1
