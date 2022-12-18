library("ggplot2")
library("ggrepel")
library("reshape2")
library("cowplot")
library("dplyr")
#library("magick")
library("ComplexHeatmap")
library("Matrix")
library("parallel")
#

##########Take genes x cells UMI matrix, filter genes, logtransform and centre
filter.umi<-function(matrix,
                     cells,
                     min_umi=5,
                     min_cells=20,
                     min_avlog=4,
                     centre="mean",
                     ref=NULL,
                     log=TRUE,
                     cut=NULL,
                     start="umi",
                     whitelist=NULL,
                     sparse=FALSE){
  m2<-matrix[,colnames(matrix) %in% cells]
  #Merge doublet genes
  m2<-Matrix.utils::aggregate.Matrix(m2,rownames(m2),fun="sum")

  #UMI to TPM
  if(start=="umi"){
  count_sum <- colSums(m2)  #apply(m2, 2, sum)
  tpmmat <- (t(t(m2)/count_sum))*1000000}
  if(start=="tpm"){tpmmat<-m2}
  #Average logtransformed tpm for filtering
  tpmav<- rowMeans(tpmmat) #apply(tpmmat,1,mean)
  logfilt<-log2((tpmav)+1)
  #UMI threshold for genes highly expressed in few cells
  if(isFALSE(sparse)){
  if(is.null(cut)){umi_threshold<-apply(m2,1,function(x)sum(x>=min_umi))}else{
    z1<-ntile(cells,cut)
    zn<-unique(z1)
    l1<-parallel::mclapply(zn,mc.cores=cut,function(x){
      apply(m2[,z1==x],1,function(y)sum(y>=min_umi))
    })
    d1<-do.call(cbind,l1)
    umi_threshold<-rowSums(d1)
      }
  #Filter genes with mean log expression > 4 OR passing UMI threshold
  genestokeep <- logfilt>=min_avlog|umi_threshold>=min_cells|names(logfilt) %in% whitelist}else{
  genestokeep <- logfilt>=min_avlog|names(logfilt) %in% whitelist
  }
  filtered.matrix <- tpmmat[genestokeep,]
  #Logtransform and centre
  if(isFALSE(log)){uncmat<-filtered.matrix}
  if(isTRUE(log)){
    if(isTRUE(sparse)){
    uncmat<-filtered.matrix
    uncmat[uncmat>0]<-log2( (uncmat[uncmat>0]/10) +1)
    }else{uncmat<-log2( (filtered.matrix/10) +1)}
  }
  if(centre=="mean"){avg <- rowMeans(uncmat)} #apply(uncmat, 1, mean)}
  if(centre=="median"){avg <- rowMedians(uncmat)}
  if(centre=="median" | centre=="mean"){centmat <- sweep(uncmat,1,avg)}
  if(centre=="ref"){centmat<-sweep(uncmat,1,ref[names(ref) %in% rownames(uncmat)])}
  if(centre!="none"){return(centmat)}
  if(centre=="none"){return(uncmat)}
}

TPMtoLog<-function(x,centre=FALSE){
  
  logged<-log2((x/10) +1)
  if(isTRUE(centre)){
    av<-rowMeans(logged)
    logged<-sweep(logged,1,av)
      }
  return(logged)
}

LogtoTPM<-function(x){
  
  tpmed<-10*(2^x-1)
  return(tpmed)
}
#Perform hierarchical clustering with average distance and reorder correlation matrix
hclustreorder <- function(x){
  d <- as.dist((1-x)/2)
  hc <- hclust(d, method = "average")
  x <-x[hc$order, hc$order]
}



####################Functions to fiddle with chromosomes

####chr.order - with gene list as input, reorder genes by chromosome  position and get breaks for chromosomes and arms

chr.order<-function(genes,local=NULL,
                    attributes=c("hgnc_symbol", "start_position", "chromosome_name","band"),h="uswest.ensembl.org"){
  
  
  #Choose which species to use and server to download from
  if(is.null(local)){
  mart <- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host=h, port = 80)
  #Get attributes for gene list
  cna_attr <- getBM(attributes=attributes,  filters="hgnc_symbol",values=genes, mart=mart, uniqueRows=T)
  }
  if(!is.null(local)){
  cna_attr<-as.data.frame(local)
  rownames(cna_attr)<-cna_attr[,1]
  cna_attr<-cna_attr[rownames(cna_attr) %in% genes,]
  colnames(cna_attr)[1]<-"hgnc_symbol"
  
  }
  #Remove genes not mapped to proper chromosomes
  real<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")
  rmlist<-setdiff(unique(cna_attr$chromosome_name),real)
  
  rmgenes<-lapply(rmlist,function(x){
    
    rm1<-cna_attr[grep(x,cna_attr$chromosome_name),]
  })
  allrm<-do.call(rbind,rmgenes)
  cna_attr<-cna_attr[!is.element(rownames(cna_attr),rownames(allrm)),]
  #Remove doubles
  cna_attr<-cna_attr[!duplicated(cna_attr$hgnc_symbol),]
  #Remove NAs
  cna_attr<- na.omit(cna_attr)
  #Change X and Y chromosomes to numbers for ordering
  cna_attr$chromosome_name<-gsub("X",23,cna_attr$chromosome_name)
  cna_attr$chromosome_name<-gsub("Y",24,cna_attr$chromosome_name)
  cna_attr$chromosome_name<-as.numeric(cna_attr$chromosome_name)
  #Order
  cna_attr<-cna_attr[order(cna_attr$chromosome_name,cna_attr$start_position, decreasing = F),]
  #Chromosome number as sorting vector
  chr_names<-as.character(unique(cna_attr$chromosome_name))
  ##Chromosome length
  chr_length1<- lapply(chr_names, function(x){
    
    l<-nrow(subset(cna_attr, cna_attr$chromosome_name ==x))
    
  })
  chr_length<-do.call(rbind,chr_length1)
  #Break vector at each chromosome end
  chr_breaks<-cumsum(chr_length)
  ######Split into chromosome arms
  #Set all genes as p or q
  cna_attr$band = gsub("[[:digit:]]", "", cna_attr$band)
  cna_attr$band = gsub("p.", "p", cna_attr$band)
  cna_attr$band = gsub("q.", "q", cna_attr$band)
  #Chromosome arm lenghts
  arm_breaks<-integer(length = length(chr_breaks))
  for(i in 1:length(chr_breaks)){
    
    a<-chr_breaks[i]
    b<-nrow(cna_attr[cna_attr$chromosome_name==i&cna_attr$band=="q",])
    c<-a-b
    arm_breaks[i]<-round(c,digits = 0)
  }
  #Breaks and labels for arms and full chromosomes
  full_breaks<-sort(c(1,chr_breaks,arm_breaks))
  #Add p and q labels
  qq<-paste(seq(1,length(chr_breaks)),"q",sep = "")
  pp<-paste(seq(1,length(chr_breaks)),"p",sep = "")
  full_labels<-sort(c(seq(1,length(chr_breaks)),seq(1,length(chr_breaks))))
  full_labels[seq(2,length(full_labels), 2)]<-qq
  full_labels[seq(1,length(full_labels), 2)]<-pp
  #Empty label at end of Y chromosome
  full_labels<-c(full_labels," ")
  #Name X and Y chromosomes
  full_labels[45:48]<-c("Xp","Xq","Yp","Yq")
  chr_breaks<-c(1,chr_breaks)
  chr_names[23:25]<-c("X","Y"," ")
  
  output<-list(cna_genes=cna_attr$hgnc_symbol,
               chr_breaks=chr_breaks,
               chr_names=chr_names,
               arm_breaks=arm_breaks,
               full_breaks=full_breaks,
               full_labels=full_labels,
               all=cna_attr)
  return(output)
}

############calc.cna - calculate CNA scores per cell given a  gene  x cell matrix, query cells and two reference cell lists
calc.cna<-function(matrix,
                   query,
                   ref,
                   top_genes=NULL,
                   mav_n=100,
                   range=c(-3,3),
                   per_chr=TRUE,
                   scale=0.05,
                   noise=NULL,
                   startlog=FALSE,
                   centred=FALSE,
                   genome=NULL,min_rl=10){
  
  
  #Remove small references
  rl<-lapply(ref,length)
  ref<-ref[rl>min_rl]
  
  
  #Matrix to work on
  cna_mat<-matrix[,c(query,unlist(ref))]
  
  #Order chromosomes
  a1<-chr.order(rownames(cna_mat),local=genome)
  genes<-a1$cna_genes
  chr_breaks<-a1$chr_breaks
 
  
  #Optional list of top expressed genes
  if(!is.null(top_genes)){
    
    # if(isTRUE(startlog)){
    #   cna_mat<-LogtoTPM(cna_mat)}
    
    m1<-cna_mat[genes,]
    m2<-apply(cna_mat,1,mean)
    m3<-names(tail(sort(m2),n=top_genes))
    a1<-chr.order(m3,local=genome)
    genes<-a1$cna_genes
    chr_breaks<-a1$chr_breaks
    
    
    # expmat<-cna_mat[genes,]
    # gex<-apply(expmat,1,mean)
    # gsort<-sort(gex,decreasing = T)
    # top<-gsort[1:top_genes]
    # g2<-intersect(genes,names(top))
    # genes<-g2
    # temp<-chr.order(genes)
    # cna_genes<-temp$cna_genes
    # chr_breaks<-temp$chr_breaks
    # chr_names<-temp$chr_names
    # full_labels<-temp$full_labels
    # full_breaks<-temp$full_breaks
    # genes<-cna_genes
  }  
  
  
  
  #Reorder
  ordered_mat<-cna_mat[genes,]
  #Log before first centering
  if(isFALSE(startlog)){
    ordered_mat<-TPMtoLog(ordered_mat)}
  #First row centering step
  if(isFALSE(centred)){
  avg<-apply(ordered_mat,1,mean)
  ordered_mat<-sweep(ordered_mat,1,avg)
  }
  #Set 3 and -3 as extreme values
  ordered_mat<- apply(ordered_mat,2,function(x)pmax(x,range[1]))
  ordered_mat<- apply(ordered_mat,2,function(x)pmin(x,range[2]))
  #TPM again for moving average
  ordered_mat<-LogtoTPM(ordered_mat)
  #Calculate mav per chromosome
  if(isTRUE(per_chr))
  {
    num<-seq(1:(length(chr_breaks)-1))
    perchr<-lapply(num,function(y){
      
      
      if(y==length(num)) {end=nrow(ordered_mat)}
      if(y!=length(num)){end<-chr_breaks[y+1]-1}
      chr<-ordered_mat[chr_breaks[y]:end,]
      #
      chr_mat<-apply(chr,2,function(x)caTools::runmean(x,k=mav_n,endrule = "mean"))
    })
    #
    calc_mat<-do.call(rbind,perchr)  
  }
  #Calculate moving average for all genes as one chromosome
  if(isFALSE(per_chr))
  {calc_mat<-apply(ordered_mat,2,function(x)caTools::runmean(x,k=mav_n,endrule = "mean"))}
  #Log before second centering
  calc_mat<-TPMtoLog(calc_mat)
  #Substract median per cell
  cmed<-apply(calc_mat,2,median)
  calc_mat <- sweep(calc_mat,2,cmed)
  #Back to TPM for reference removal
  calc_mat<-LogtoTPM(calc_mat)
  #Create max/minvalues per gene from reference cells
  refl<-seq(1,length(ref))
  meanmat<-lapply(refl,function(x){
    
    
    r1<-ref[[x]]
    m1<-apply(calc_mat[,r1],1,mean)
    m1})
  refmat<-do.call(cbind,meanmat)
  #Log references
  refmat<-TPMtoLog(refmat)
  refmax<-apply(refmat,1,max)
  refmin<-apply(refmat,1,min)
  #Expand reference boundaries by scaling percentage
  if(!is.null(scale)){
    rmax<-refmax+scale*abs(refmax)
    rmin<-refmin-scale*abs(refmin)
  }
  #Or expand by fixed noise factor
  if(!is.null(noise)){
    rmax<-refmax+noise
    rmin<-refmin-noise 
  }
  #Log cna matrix
  calc_mat<-TPMtoLog(calc_mat)
  #Centre by reference
  score_mat<- ifelse(calc_mat>rmax,calc_mat-rmax,
                     ifelse(calc_mat<rmin,calc_mat-rmin,0))
  rownames(score_mat)<-rownames(ordered_mat)
  
  return(score_mat)
}


####Adapt breaks and labels to remove genes without moving average
cut.flank<-function(genes,separate="chr",n=50){
  
  if(separate=="chr") {breaks = chr_breaks
  labels=chr_names}
  if(separate=="arm") {breaks = full_breaks
  labels=full_labels}
  
  #Shift breaks
  breaks<-pmax(1,breaks-n)
  #Remove labels outside new boundaries
  bdf<-cbind.data.frame(breaks,labels)
  a1<-length(genes)-2*n
  bdf2<-bdf[bdf$breaks<a1,]
  a2<-nrow(bdf2)
  bb1<-c(bdf2$breaks,a1)
  bb2<-bdf$labels[0:a2+1]
  g2<-genes[-c(1:n)]
  g3<-length(g2)-n
  g4<-g2[1:g3]
  outputs<-list(breaks=bb1,labels=bb2,genes=g4)
  return(outputs)
}


##########Plot CNA with score matrix and cells of interest as input, separate by chromosome or arm, cut by half of moving av.
plot.cna<-function(matrix,cells,name=NULL,separate="chr",order=TRUE){
  
  
  if(separate=="chr") {breaks = chr_breaks
  labels=chr_names}
  if(separate=="arm") {breaks = full_breaks
  labels=full_labels}
  
  plot_matrix<-matrix[,cells]
  if(isTRUE(order)){
    #Correlation matrix
    cormat<-cor(as.matrix(plot_matrix), method = "pearson")
    #Reorder matrix
    ordmat<-hclustreorder(cormat)
    ordnames<-colnames(ordmat)
    plot_matrix<-plot_matrix[,ordnames]
  }
  
  
  rownames(plot_matrix)<-seq(1,nrow(plot_matrix))
  colnames(plot_matrix)<-seq(1,ncol(plot_matrix))
  melted_matrix<-reshape2::melt(plot_matrix)
  #Create plot
  p1<-ggplot(melted_matrix, aes(Var1, Var2, fill = value))+
    geom_tile()+
    scale_fill_gradient2(low="steelblue", mid = "white", high="darkred", midpoint=0,
                         limits=c(-1,1), space = "Lab", 
                         name="Gene\nExpression") + 
    theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
    scale_x_continuous(breaks=breaks,labels=labels, expand = c(0,0)) + scale_y_continuous(expand = c(0,0))+
    geom_vline(xintercept=breaks) +ggtitle(paste(name, "CNA plot",sep = " "))
  return(p1)
}
#

#######Define CNA from score matrix
define.cna<-function(cna_matrix,
                     epi_cells,
                     ref_cells,
                     name=NULL,
                     tumour_score=hpvpos_score,
                     gs=HPV_score,
                     gscut=1,
                     min_gs=5,
                     top_region=1/3,
                     top_cells=1/3,
                     spec_cut=0.99,
                     uppercut=2,
                     lowercut=2,
                     extra_score=ifelse(HPV[epi_cells]=="HPV+",1,0),
                     knnprob=0.9,
                     reass=FALSE,
                     no_ts=FALSE,
                     nd=NULL,
                     glob_thr=NULL,
                     cnasig="abs",
                     print_ref=TRUE,
                     genome=NULL){
  
  
  if(!is.null(intersect(ref_cells,epi_cells))){rmcells<-intersect(ref_cells,epi_cells)
  ref_cells<-setdiff(ref_cells,rmcells)
  }    
  cna_matrix<-cna_matrix[,c(epi_cells,ref_cells)]
  all_cells<-colnames(cna_matrix)
  epi_matrix<-cna_matrix[,epi_cells]
  #Absolute CNA value per gene for all epithelial cells
  absvals<-apply(epi_matrix,1,function(x)mean(abs(x)))
  if(top_region<1){
  #Top third genes with highest absolute value
  abs2<-ifelse(absvals>quantile(absvals, prob = 1-top_region),absvals,0)
  #Define matrices where to calculate cell averages - only the regions with most CNA
  score_cna<-cna_matrix[abs2>0,]
  }else{
  #Select arms with top absolute value  
    br1<-br.fun(genes=rownames(cna_matrix),separate = "arm",
                genecut=5,genome=genome)
    br_vec<-br1$br_vec
    labels<-br1$bdf$labels
    rm<-br1$rm
    
    if(length(rm)!=0){
      cna_matrix<-cna_matrix[-rm,]
      epi_matrix<-epi_matrix[-rm,]
      absvals<-absvals[-rm]}
    bn<-unique(br_vec)
    
    bf1<-unlist(mclapply(bn,function(x){
mean(absvals[br_vec==x])
    }))
    
    #Create mixture model
    mix<-nor1mix::norMixEM(bf1, m = 2, maxiter = 2500)
    #Get probabilities
    post<-nor1mix::estep.nm(bf1, mix)
    s1<-mean(bf1[which(post[,1]>top_region/100)])-mean(bf1[which(post[,2]>top_region/100)])
    if(s1>0){
      a<-which(post[,1]>top_region/100)
      sel<-bn[a]}else{
        a<-which(post[,2]>top_region/100)
        sel<-bn[a]
      }
    score_cna<-cna_matrix[br_vec %in% sel,]
  }
  #CNA score for each tumour cell across relevant region
  if(cnasig=="abs"){
  cnascore<-apply(score_cna,2,function(x)mean(abs(x)))}
  if(cnasig=="sqrt"){cnascore<-apply(score_cna,2,function(x)mean(sqrt(x^2)))}
  #Calculate correlation only with epithelial cells with strongest signal
  cellcut<- quantile(cnascore, probs = 1-top_cells)
  if(top_region<1){epi_cna<-epi_matrix[abs2>0,cnascore[epi_cells]>cellcut]
  }else{epi_cna<-epi_matrix[br_vec %in% sel,cnascore[epi_cells]>cellcut]}
  #Correlation vector
  cor1<-cor(score_cna,epi_cna)
  cnacor<-apply(cor1,1,mean)
  cnacor[is.na(cnacor)]<-0
  #substract normcells in subclone reassignments
  if(isTRUE(reass)){
    normdef<-nd
    normref<-mean(normdef[normdef$cna3=="Normal","cnacor"])
    cnacor[epi_cells]<-cnacor[epi_cells]-normref
  }
  #Plot by type
  r<-as.data.frame(t(rbind.data.frame(cnacor,cnascore)))
  rownames(r)<-all_cells
  colnames(r)<-c("cnacor","cnascore")
  
  Origin<-with(r, ifelse(rownames(r) %in% epi_cells,"Epithelial","Reference"))
  #
  p1<-ggplot(r,aes(x=cnacor, y=cnascore, colour = tumour_score[all_cells]))+
    theme_bw()+
    scale_colour_gradientn(name="Tumour Score",colours=c("#00AFBB", "#E7B800", "#FC4E07"))+
    geom_point(aes(size=Origin))+
    scale_size_manual(name="Origin",values=c(2,0))+
    ggtitle(paste(name,"CNA scores",sep = " "))
  ###########################Define gold standard cells
  #If no tumour score is used, create empty placeholder
  if(is.null(tumour_score)){
    #tumour_score<-runif(length(all_cells),0,1)
    tumour_score<-numeric(length = length(all_cells))
    for(i in 1:length(tumour_score)){tumour_score[i]<-0}
    names(tumour_score)<-all_cells}
    
    r$tscore<-tumour_score[all_cells]
  if(!is.null(gs)){
    r$cna1<-ifelse(Origin=="Epithelial"&gs[all_cells]>gscut,"Cancer",ifelse(Origin=="Reference","Normal","Unresolved"))
  }
  if(is.null(gs)){g1=0
  r$cna1<-ifelse(Origin=="Epithelial","Unresolved",ifelse(Origin=="Reference","Normal","Other"))}
  #
  r$Origin<-Origin
  p2<-ggplot(r,aes(x=cnacor, y=cnascore, colour = cna1)) + theme_bw()+
    geom_point(aes(size=Origin))+
    scale_size_manual(name="Origin",values=c(2,0))+
    ggtitle(paste(name,"Gold standard definition",sep = " "))
  ###########Check whether there are enough gold standard cells to use GS definition
  if(!is.null(gs)){
    
    gscore<-gs[epi_cells]
    g1<-length(gscore[gscore>gscut])}
 
  if(spec_cut<1){ 
  if(g1<min_gs | is.null(gs)){
    ##########Local thresholds for cells without gold standard
    cor_thr<-quantile(r[r$Origin=="Reference","cnacor"],probs = spec_cut)
    score_thr<-quantile(r[r$Origin=="Reference","cnascore"],probs = spec_cut)
    ts_thr<-quantile(r[r$Origin=="Reference","tscore"],probs = spec_cut)
  }
  }else{
    cor_thr<-median(r[r$Origin=="Reference","cnacor"])+spec_cut*sd(r[r$Origin=="Reference","cnacor"])
    score_thr<-median(r[r$Origin=="Reference","cnascore"])+spec_cut*sd(r[r$Origin=="Reference","cnascore"])
    ts_thr<-median(r[r$Origin=="Reference","tscore"])+spec_cut*sd(r[r$Origin=="Reference","tscore"])
      }
  
  ###########Set local thresholds based on GS
  if(!is.null(gs) & g1>=min_gs){
    
    fix_class<-r$cna1[r$cna1!="Unresolved"]
    
    #
    varlist<-c("cnacor","cnascore","tscore")
    
    thresholds<-lapply(varlist,function(var){
      
      temp_var<-r[r$cna1!="Unresolved",var]
      temp_roc<-roc(fix_class ~ temp_var, levels=c("Normal","Cancer"), plot=FALSE)
      temp_coords<-coords(temp_roc, "local maximas", ret=c("threshold","sens", "spec", "ppv", "npv"),transpose=TRUE)
      temp_coords<-as.data.frame(t(temp_coords))
      filt_coords<-temp_coords[temp_coords$specificity>spec_cut,]
      thr<-filt_coords[1,"threshold"]
      thr1<-round(filt_coords[1,"threshold"],digits = 2)
      sp1<-round(filt_coords[1,"specificity"],digits=2)
      sn1<-round(filt_coords[1,"sensitivity"],digits=2)
      
      df<-cbind.data.frame(temp_roc$specificities,temp_roc$sensitivities)
      colnames(df)<-c("spec","sens")
      
      p3<-ggplot(df,aes(x=spec,y=sens))+geom_line()+theme_bw()+
        geom_point(x=sp1, y=sn1,colour="red")+
        geom_text(x=sp1-0.2,y=sn1-0.1,label=paste(thr1,sp1,sn1,sep = ","))+
        ggtitle(paste(name,var,"ROC",sep=" "))
      #
      outs<-list(thr,p3)
      outs
    })
    
    
    thr.out<-do.call(rbind,thresholds)
    rocplots<-thr.out[,2]
    p3<-do.call(plot_grid,rocplots)
    thresholds<-thr.out[,1]
    cor_thr<-thresholds[[1]]
    score_thr<-thresholds[[2]]
    ts_thr<-thresholds[[3]]
  }
  ############Create cell classification score
  if(!is.null(glob_thr)){
    score_thr<-max(score_thr,glob_thr$score)
    cor_thr<-max(cor_thr,glob_thr$cor)
    ts_thr<-max(ts_thr,glob_thr$ts)
  }
  
  
  score_bin<-ifelse(r[epi_cells,"cnascore"] > score_thr,1,0)
  cor_bin<-ifelse(r[epi_cells,"cnacor"] > cor_thr,1,0)
  ts_bin<-ifelse(r[epi_cells,"tscore"]>ts_thr,1,0)
  if(isTRUE(no_ts)){ts_bin<-0
  uppercut<-uppercut-1
  lowercut<-lowercut-1}
  if(!is.null(extra_score)){
    cna_bin<-score_bin+cor_bin+ts_bin+extra_score}
  if(is.null(extra_score)){cna_bin<-score_bin+cor_bin+ts_bin}
  
  
  #
  r[epi_cells,"cna2"]<-ifelse(cna_bin>uppercut,"Cancer",ifelse(cna_bin<lowercut,"Normal","Unresolved"))
  
  #
  r$totscore<-with(r,ifelse(rownames(r) %in% epi_cells,cna_bin,-1))
  r$cna3<-with(r,ifelse(rownames(r) %in% epi_cells,r$cna2,"Reference"))
  #Plot definition with cutoffs
  p4<-ggplot(r,aes(x=cnacor, y=cnascore, colour = cna3)) + theme_bw()+
    geom_point(aes(size=Origin))+
    scale_size_manual(name="Origin",values=c(2,0))+
    geom_vline(xintercept=cor_thr) +
    geom_hline(yintercept = score_thr)  + 
    ggtitle(paste(name,"CNA, 1st definition with cutoffs",sep = " "))
  
  ##########Add KNN as parameter
  unres_cells<-rownames(r[r$cna3=="Unresolved",])
  cancer_cells<-rownames(r[r$cna3=="Cancer",])
  normal_cells<-rownames(r[r$cna3=="Normal",])
  if(length(unres_cells)>1 &length(cancer_cells)>1&length(normal_cells)>1){
    knn_cna<-class::knn(train = r[r$cna3!="Unresolved"&r$cna3!="Reference",1:2],
                        test = r[r$cna3=="Unresolved",1:2],k=max(1,round(log(length(epi_cells)))),
                        cl=r$cna3[r$cna3!="Unresolved"&r$cna3!="Reference"],prob=T)
    r[unres_cells,c("knn","prob")]<-c(as.character(knn_cna),attr(knn_cna,"prob"))
    r$cna4<-ifelse(!rownames(r) %in% unres_cells, r$cna3, 
                   ifelse(rownames(r) %in% unres_cells & r$knn=="Cancer"& r$prob>=knnprob,"Cancer",
                          ifelse(rownames(r) %in% unres_cells & r$knn=="Normal"& r$prob>=knnprob,"Normal","Unresolved")))
  }else{r$cna4<-r$cna3}
  #Plot final definition
  p5<-ggplot(r,aes(x=cnacor, y=cnascore, colour = cna4)) + theme_bw()+
    geom_point(aes(size=Origin))+
    scale_size_manual(name="Origin",values=c(2,0))+
    geom_vline(xintercept=cor_thr) +
    geom_hline(yintercept = score_thr)  + 
    ggtitle(paste(name,"CNA, final definition",sep = " "))
  #
  plots<-plot_grid(p1,p2,p4,p5,ncol = 2,nrow = 2)
  if(isFALSE(print_ref)){
  def_mat<-r[epi_cells,]}else{def_mat<-r}
  output<-list(plots=plots,def_mat=def_mat,if(!is.null(gs) & g1>=min_gs){roc=p3},thresholds=c(cor_thr,score_thr,ts_thr))
  return(output)
}
#
br.fun<-function(genes,genecut=10,separate="chr",genome=NULL){
  
  a1<-chr.order(genes,local=genome)
  genes<-a1$cna_genes
  chr_breaks<-a1$chr_breaks
  chr_names<-a1$chr_names
  full_breaks<-a1$full_breaks
  full_labels<-a1$full_labels
  
  
  if(separate=="chr") {breaks = chr_breaks
  labels=chr_names}
  if(separate=="arm") {breaks = full_breaks
  labels=full_labels}
  if(length(breaks)!=length(labels)){labels<-labels[1:length(breaks)]}
  #Breaks as factor
  br_vec<-character(length = length(genes))
  
  for(i in 1:(length(breaks)-1)){
    
    if(i==(length(breaks)-1)) {end=length(genes)}
    if(i!=(length(breaks)-1)){end<-breaks[i+1]-1}
    v<- end - breaks[i]
    
    if(v<=genecut & v>0){br_vec[breaks[i]:end]<-NA}
    if(v==0){br_vec[breaks[i]]<-NA}
    if(v > genecut){
      br_vec[breaks[i]:end]<-i}
  }
  #Remove genes on chromosomes with too few genes
  rm<-which(is.na(br_vec))
  br_vec<-na.omit(br_vec)
  br_vec<-sort(as.numeric(br_vec))
  #Keep only labels with chromosomes that have enough genes
  bdf<-cbind.data.frame(breaks,labels)
  rownames(bdf)<-seq(1,nrow(bdf))
  b2<-bdf[unique(br_vec),]
  labels<-b2[,"labels"]
  
  out<-list(br_vec=br_vec,bdf=b2,rm=rm)
  return(out)
}
################Complex heatmap of chromosomes

plot.cna2<-function(matrix,
                    cells,
                    separate="chr",
                    order=TRUE,
                    scores=NULL,
                    scorenames=NULL,
                    name=NULL,typeplot=FALSE,numscale=NULL,genome=NULL,...){
  
  
  br1<-br.fun(genes=rownames(matrix),separate = separate,genome=genome)
  br_vec<-br1$br_vec
  labels<-br1$bdf$labels
  rm<-br1$rm
  
  #Filter matrix if genes were removed
  if(length(rm)!=0){
    matrix<-matrix[-rm,]}
  
  #Subset matrix
  plot_matrix<-matrix[,cells]
  if(isTRUE(order)){
    #Correlation matrix
    cormat<-cor(as.matrix(plot_matrix), method = "pearson")
    #Reorder matrix
    ordmat<-hclustreorder(cormat)
    ordnames<-colnames(ordmat)
    plot_matrix<-plot_matrix[,ordnames]
    cells<-ordnames
  }
  
  #Function to get main colours
  col_fun = circlize::colorRamp2(c(-1, 0, 1), c("steelblue", "white", "darkred"))
  #Use arbitrary number of scores for annotating
  
  if(is.list(scores)){
    nn<-seq(1,length(scores))
    sfun<-lapply(nn,function(x){
      
      s1<-scores[[x]][cells]
    })
    scoremat<-do.call(cbind.data.frame,sfun)
    colnames(scoremat)<-scorenames
    row_ann = rowAnnotation(df=scoremat)
  }
  if(isTRUE(typeplot)){
    cf2<-lapply(scorenames,function(x){
      
      a1<-as.character(x)
      a<-circlize::colorRamp2(c(min(scoremat[,a1]),median(scoremat[,a1]),max(scoremat[,a1])),c("#00AFBB", "white", "#FC4E07"))
      a})
    names(cf2)<-scorenames
    
    row_ann = rowAnnotation(df=scoremat,col=cf2)
  }
  if(!is.null(numscale)){
    
    q<-numscale
    cf2<-lapply(scorenames,function(x){
      
      a1<-as.character(x)
      a<-circlize::colorRamp2(c(-1*q,0,q),c("#00AFBB", "white", "#FC4E07"))
      a})
    names(cf2)<-scorenames
    
    row_ann = rowAnnotation(df=scoremat,col=cf2)
  }
  if(is.null(scores)){row_ann=NULL}
  #Row annotation
  #TS = tumour_score[cells]
  #cf2<-circlize::colorRamp2(c(min(TS),median(TS),max(TS)),c("#00AFBB", "#E7B800", "#FC4E07"))
  
  
  pp2<-invisible(Heatmap(t(plot_matrix), name = "Inferred Copy Number\nLog-ratio",
                         cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,col = col_fun,
                         column_title_side = "bottom", left_annotation = row_ann,cluster_rows = FALSE,
                         column_split = br_vec, column_title = labels,border = TRUE,column_gap = unit(0,"mm"),...))
  
  pp22<-invisible(draw(pp2, padding = unit(c(10, 10, 10, 10), "mm"),   column_title=paste(name,"CNA plot", sep = " "), column_title_side = "top",
                       column_title_gp=gpar(fontsize=20)))
  
  return(pp22)
  
}
##################Check for bimodality
check.bim<-function(matrix,genes,prob=0.95,freq=0.8,size=10){
  m1<-matrix[genes,]
  arm_av<-colMeans(m1)
  #Create mixture model
  mix<-nor1mix::norMixEM(arm_av, m = 2, maxiter = 2500)
  #Get probabilities
  post<-nor1mix::estep.nm(arm_av, mix)
  #Show which cells pass cutoff
  passed = post >=prob
  #
  c1<-colSums(passed)[1]
  c2<-colSums(passed)[2]
  cfreq<- sum(colSums(passed))/length(arm_av)
  if(c1>size & c2>size &cfreq>freq){res=TRUE}else{res=FALSE}
  
  #Assign cells to subclones
  rownames(passed)<-colnames(matrix)
  d1<-rownames(passed[passed[,1]==TRUE,])
  d2<-rownames(passed[passed[,2]==TRUE,])
  
  cf<-ifelse(rownames(passed) %in% d1, "a",
             ifelse(rownames(passed) %in% d2, "b","none"))
  #
  
  outdf<-c(c1,c2,cfreq,res)
  names(outdf)<-c("c1","c2","cfreq","res")
  out<-list(outdf,cf)
  return(out)
}
#############Find subclones
clone.fun<-function(matrix,cells,separate="arm",clonesize=10,...){
  
  #Break vector
  br1<-br.fun(rownames(matrix),separate=separate)
  br_vec<-br1$br_vec
  bdf<-br1$bdf
  rm<-br1$rm
  #Filter matrix if genes were removed
  if(length(rm)!=0){
    matrix<-matrix[-rm,]}
  #Subset matrix
  matrix<-matrix[,cells]
  #Check for bimodality per split
  bn<-unique(br_vec)
  fit.bim<-lapply(bn,function(x){
    s<-check.bim(matrix,genes=rownames(matrix[br_vec==x,]),...)
    s
  })
  #
  fitout<-do.call(rbind,fit.bim)
  fitres<-do.call(rbind,fitout[,1])
  fitres<-cbind.data.frame(bdf,fitres)
  #Data frame of clone combinations
  clonemat<-do.call(cbind.data.frame,fitout[,2])
  colnames(clonemat)<-bdf$labels
  rownames(clonemat)<-colnames(matrix)
  #Arms to keep
  keep<-as.character(fitres[fitres$res==1,"labels"])
  #Subset clone matrix to arms that are bimodal
  clonemat<-clonemat[,keep]
  #Count combinations
  combmat<-plyr::count(clonemat)
  #Select clones with minimum cell number
  combmat<-combmat[combmat$freq>clonesize,]
  #Tidy combination matrix
  colnames(combmat)<-gsub("X","",colnames(combmat))
  combmat2<-combmat
  combmat[,"freq"]<-NULL
  #Create strings of combinations
  combs<-do.call(paste0,combmat)
  clonecomb<-do.call(paste0,clonemat)
  
  #Select cells in subclones
  cl1<-seq(1,length(combs))
  clonefun<-lapply(cl1,function(x){
    
    a<-paste0("subclone",x)
    b<-rownames(clonemat[clonecomb==combs[x],])
    c<-cbind.data.frame(b,a)
    colnames(c)<-c("cells","clone")
    c
  })
  #
  ###Bind subclone cells to rest
  cl2<-do.call(rbind,clonefun)
  s<-setdiff(colnames(matrix),cl2$cells)
  cl3<-cbind.data.frame(s,"rest")
  colnames(cl3)<-colnames(cl2)
  cl4<-rbind.data.frame(cl2,cl3)
  clones<-cl4
  rownames(clones)<-clones$cells
  #
  out<-list(clones=clones,combs=combmat2)
  return(out)
}
#
###Extract clusters by DBScan, finding the knee point of knn distance graph
dbscan.cluster<-function(matrix,minpts=log(nrow(matrix)),probs=3/4,diffcut=10,assign_all=FALSE){
  #Calculate knn distance
  dist <- dbscan::kNNdist(matrix, minpts)
  #Find knee point through derivative of knn graph
  d2<-quantile(sort(dist),probs = probs)
  d3<-dist[dist>d2]
  dist<-d3
  # order result
  dist <- dist[order(dist)]
  # scale
  dist <- dist / max(dist)
  # derivative
  ddist <- diff(dist) / (1/length(dist))
  # get first point where derivative is higher than 1
  knee <- dist[length(ddist)- length(ddist[ddist > diffcut])] +d2
  #Use knee point in plot as eps for dbscan algorithm
  db2<-dbscan::dbscan(matrix, eps = knee, minPts = minpts)
  #Extract clusters
  out<-as.factor(db2$cluster)
  
  if(isTRUE(assign_all)){
    matrix$cluster<-out
    c1<-class::knn(train = matrix[matrix$cluster!=0,1:2],
                        test = matrix[matrix$cluster==0,1:2],k=minpts,
                        cl=matrix$cluster[matrix$cluster!=0],prob=F)
    names(c1)<-rownames(matrix[matrix$cluster==0,])
    matrix[names(c1),"cluster"]<-c1
    out<-as.factor(matrix$cluster)
      }
  
  return(out)
}

cluster.clone<-function(matrix,
                        cells,
                        adgenes="all",
                        adregion=1,
                        genecut=10,
                        separate="arm",
                        top_method="abs",
                        dimred="uwot",
                        umap_pca=NULL,
                        pca_dims=10,
                        ica_comp=10,
                        cluster_method="dbscan",
                        louvain_k=10,
                        merge_method="ad",
                        top_region=1/3,
                        diffcut=10,
                        minpts=NULL,
                        ptscale=NULL,
                        knn=TRUE,
                        clonesize=10,
                        name=NULL,
                        merge=TRUE,
                        corcut=0.95,
                        adcut=0.2,
                        dcut=0.2,
                        armfilt=FALSE,
                        arm_n=10,
                        sdcut=0.15,
                        expcut=0.15,
                        hclust_k=15,
                        by_cell=FALSE,
                        armcut=0,
                        cell_matrix="numeric",
                        perm_bins=50,
                        perm_cut=0.01,
                        perm_trials=100,genome=NULL,...){
  
  
  
  ##Subset by chromosome
  br1<-br.fun(genes=rownames(matrix),separate = separate,
              genecut=genecut,genome=genome)
  br_vec<-br1$br_vec
  labels<-br1$bdf$labels
  rm<-br1$rm
  
  if(length(rm)!=0){
    matrix<-matrix[-rm,]}
  bn<-unique(br_vec)
  
  ###Set working and reference matrices
  refmat<-matrix
  matrix<-matrix[,cells] 
  names(br_vec)<-rownames(matrix)
  
  #Select most variant gene regions  
  if(top_method=="sd"){
    sd1<-apply(matrix,1,sd)
    m2<-matrix[sd1>quantile(sd1, prob = 1-top_region),]
  }
  
  ###Select most CNA-rich regions
  if(top_method=="abs"){
    absvals<-apply(matrix,1,function(x)mean(abs(x)))
    m2<-matrix[absvals>quantile(absvals, prob = 1-top_region),]
  }
  
  tm1<-as.matrix(t(m2))
  
  #Cluster without dimred
  if(isFALSE(by_cell)){
  if(dimred=="none"){tm3<-tm1}
  
  if(dimred=="pca"){
    p1<-prcomp(tm1, center = FALSE, scale. = FALSE)
    tm3<-p1$x[,1:pca_dims]}
  if(dimred=="svd"){  
    p1<-svd::propack.svd(tm1,neig = pca_dims)
    tm3<-as.data.frame(p1$u)
  }
   
  if(dimred=="ica"){
    i1<-fastICA::fastICA(tm1,n.comp=ica_comp,method="C")
    tm3<-i1$S
  }
  
  
  
  #Run tSNE
  if(dimred=="tsne"){
    tm2 = Rtsne::Rtsne(tm1, check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2,
                       pca_center=FALSE,pca_scale=FALSE,...)
    tm3 = as.data.frame(tm2$Y) #data frame for ggplotting
  }
  #Run Umap
  if(dimred=="umap"){
    
    if(!is.null(umap_pca)){
      p1<-prcomp(tm1, center = FALSE, scale. = FALSE)
      tm1<-p1$x[,1:umap_pca]}
    
    
    umap.object<-umap::umap(tm1, metric="pearson", spread=5,min_dist=0.01,n_neighbors=clonesize)
    tm3<-cbind.data.frame(umap.object$layout[,1],umap.object$layout[,2])
    colnames(tm3)<-c("V1","V2")  
  }
  if(dimred=="uwot"){
    tm3<-standard.umap(t(tm1),n_neighbors = clonesize,spread=5,min_dist = 0.01,pca=umap_pca)
      }
    
    
  
  rownames(tm3)<-colnames(matrix)
  
  #########Louvain clustering
  if(cluster_method=="louvain"){
    tcl<-cluster.coord(tm3,method = "louvain",louvain = louvain_k)
  }
 
  #####DBScan
  if(cluster_method=="dbscan"){
    #Select minpts for knn graph
    if(is.null(minpts) & is.null(ptscale)){minpts<-log(nrow(tm3))}
    if(!is.null(ptscale)){minpts<-ptscale*log(nrow(tm3))}
    if(!is.null(minpts)){minpts=minpts}
    #Perform dbscan with set minpts
    tcl<-dbscan.cluster(tm3,minpts=minpts,diffcut=diffcut)}
  
  #####Hierachical clustering
  if(cluster_method=="hclust"){
  hcc<-cor(t(tm3))
  dfc <- as.dist((1-hcc)/2)
  hc1 <- hclust(dfc, method = "average")
  tcl<-cutree(hc1,k=hclust_k)
  }
  
  }

  ###Methods that don't rely on dimred - cluster bycell  
  if(isTRUE(by_cell)){
    
    #Subset break vector by filtered matrix
    br2<-br_vec[rownames(m2)]
    bn2<-unique(br2)
    bl<-lapply(bn2,function(x)length(br2[br2==x]))
    bnx<-bn2[bl>genecut]
    br3<-br2[br2 %in% bnx]
    bnx<-unique(br3)
    m2<-m2[names(br3),]
    ll1<-cbind.data.frame(bn,labels)
    ll2<-ll1[ll1$bn %in% bnx,]
    
    #Select arms with highest expression and sd
    if(armcut!=0){
      bf1<-lapply(bnx,function(x){
        a<-m2[br3==x,]
        b<-mean(abs(colMeans(a)))
        c<-mean(apply(a,2,sd))
        c(b,c)
      })
      bf1<-do.call(rbind.data.frame,bf1)
      colnames(bf1)<-c("abs","sd")
      rownames(bf1)<-bnx
      bf1$chr<-ll2$labels
      bf2<-bf1[bf1$abs %in% tail(sort(bf1$abs),n=armcut) | bf1$sd %in% tail(sort(bf1$sd),n=armcut),]
      
      br3<-br3[br3 %in% rownames(bf2)]
      for(i in 1:nrow(bf2)){
        br3[br3==rownames(bf2)[i]]<-bf2$chr[i]
      }
      bnx<-unique(br3)
      m2<-m2[names(br3),]
    }
    
    
    ##########Chromosome expression matrix
    bf<-lapply(bnx,function(x){
      a<-m2[br3==x,]
      b<-colMeans(a)
      c<-apply(a,2,sd)
      d<-ifelse(b>expcut & c < sdcut,"amp",
                ifelse(b< -1*expcut & c<sdcut,"del","none"))
      d
    })
    bfm<-do.call(cbind.data.frame,bf)
    colnames(bfm)<-bnx
    
    nf<-lapply(bnx,function(x){
      a<-m2[br3==x,]
      b<-colMeans(a)
      b})
    nfm<-do.call(cbind.data.frame,nf)
    colnames(nfm)<-bnx
    
    if(cluster_method=="hclust"){
      if(cell_matrix=="binary"){
        bfm[bfm=="amp"]<-1
        bfm[bfm=="del"]<- -1
        bfm[bfm=="none"]<-0
        bfm<-apply(bfm,1,as.numeric)
        hc<-hclust(dist(t(bfm)), method = "average")
      }
      if(cell_matrix=="numeric"){
        cm1<-cor(t(nfm))
        d <- as.dist((1-cm1)/2)
        hc <- hclust(d, method = "average")
        #hc<-hclust(dist(nfm), method = "average")
        }
      tcl<-cutree(hc,k=hclust_k)
      }
    
    if(cluster_method=="paste"){
      zf3<-do.call(paste0,bfm)
      names(zf3)<-rownames(bfm)
      tcl<-zf3
      }
    
    if(cluster_method=="perm"){
      zf<-parallel::mclapply(bnx,mc.cores=min(5,length(bnx)),function(x){
        
        a<-mean(colMeans(m2[names(br3[br3==x]),]))
        if(a<0){nn<-"del";mm<-"negative"}else{nn<-"amp";mm<-"positive"}
        nn2<-paste0(x,"_",nn,"_")
        z3<-permutation.score(refmat[rownames(m2),],names(br3[br3==x]),mode=mm,
                              name = nn2,cutoff=perm_cut,trials=perm_trials,method = "sd")
        z3
      })
      zf2<-do.call(cbind.data.frame,zf)
      colnames(zf2)<-bnx
      zf3<-do.call(paste0,zf2)
      names(zf3)<-rownames(zf2)
      tcl<-zf3[names(zf3) %in% cells]
    }
    tm3<-nfm
  }
  
  #Rename clusters
  ren.fun2<-function(x){paste("subclone",x,sep = "")}
  tcl<-ren.fun2(tcl)
  names(tcl)<-colnames(matrix)
  
  ##Assign cells in small clusters by knn method if no clusters are larger than min size
  if(isTRUE(knn) & max(table(tcl))<clonesize){
    t2<-table(tcl)<=2
    t3<-t2[t2==TRUE]
    tm4<-tm3
  tm4$subclone<-tcl
  
  unres_cells<-rownames(tm4[tm4$subclone %in% names(t3),])
  tm4<-tm4[,-ncol(tm4)]
  #
  knn_sc<-class::knn(train = tm4[!rownames(tm4) %in% unres_cells,],
                     test = tm4[rownames(tm4) %in% unres_cells,],k=max(1,round(log(nrow(tm4)))),
                     cl=tcl[!names(tcl) %in% unres_cells],prob=T)
  #
  tm4[unres_cells,c("knn","prob")]<-c(as.character(knn_sc),attr(knn_sc,"prob"))
  #
  tcl<-ifelse(names(tcl) %in% unres_cells, tm4$knn,tcl)
  names(tcl)<-rownames(tm4)
  }
  
  ##Assign cells in small clusters by knn method
  t2<-table(tcl)<clonesize
  t3<-t2[t2==TRUE]
  
  if(isTRUE(knn) && length(t3)>0){
    
    tm3$subclone<-tcl
    unres_cells<-rownames(tm3[tm3$subclone %in% names(t3),])
    tm3<-tm3[,-ncol(tm3)]
    #
    knn_sc<-class::knn(train = tm3[!rownames(tm3) %in% unres_cells,],
                       test = tm3[rownames(tm3) %in% unres_cells,],k=max(1,round(log(nrow(tm3)))),
                       cl=tcl[!names(tcl) %in% unres_cells],prob=T)
    #
    tm3[unres_cells,c("knn","prob")]<-c(as.character(knn_sc),attr(knn_sc,"prob"))
    #
    tcl<-ifelse(names(tcl) %in% unres_cells, tm3$knn,tcl)
    names(tcl)<-rownames(tm3)
    
    
  }
  
  tn1<-sort(unique(tcl))
  let<-letters[1:length(tn1)]
  for(i in 1:length(tn1)){
    
    tcl[tcl==tn1[i]]<-paste(let[i],tn1[i],sep="_")
    
  }
  #
  
  
  tcl2<-tcl
  
  ###Merge clusters based on correlation of mean expression across chromosome arms
  
  if(isTRUE(merge) && length(unique(tcl)) >1){
    
    if(merge_method=="cor"){
      cormax=1
      
      while(corcut<cormax && length(unique(tcl))>1){
        
        comb<-data.frame(t(combn(unique(tcl), 2)))
        ad1<-admat.fun(matrix,br_vec,tcl,labels,cut=adcut,genes=adgenes,region=adregion)  
        epmat<-ad1$epmat
        ad<-ad1$ad
        if(isTRUE(armfilt)){
          ad<-varfilt(epmat,ad,cut=NULL,n=arm_n)$ad
          epmat<-varfilt(epmat,ad,cut=NULL,n=arm_n)$epmat}
        for(i in 1:nrow(comb)){
          
          a<-as.character(comb[i,1])
          b<-as.character(comb[i,2])
          comb[i,"cor"]<-cor(epmat[,a],epmat[,b])
          cormax<-max(comb$cor)
          
          tcl[tcl==a]<-ifelse(comb[i,"cor"] > corcut,b,a)
          
          
          
          
        }
      }
    }
    ###Merge subclones where the maximal distance between two chromosome arms is lower than cut
    if(merge_method=="dist" | merge_method=="both"){
      dmin=0
      
      while(dmin<dcut && length(unique(tcl))>1){
        
        comb<-data.frame(t(combn(unique(tcl), 2)))
        ad1<-admat.fun(matrix,br_vec,tcl,labels,cut=adcut,genes=adgenes,region=adregion)  
        epmat<-ad1$epmat
        ad<-ad1$ad
        if(isTRUE(armfilt)){
          ad<-varfilt(epmat,ad,cut=NULL,n=arm_n)$ad
          epmat<-varfilt(epmat,ad,cut=NULL,n=arm_n)$epmat}
        
        for(i in 1:nrow(comb)){
          
          a<-as.character(comb[i,1])
          b<-as.character(comb[i,2])
          #c<-stringdist::stringdist(ad[,a],ad[,b])
          #if(sum(c)>0){
          comb[i,"dist"]<-max(abs(epmat[,a]-epmat[,b]))
          #if(sum(c)==0){comb[i,"dist"]<-0}
          dmin<-min(comb$dist)
          
          tcl[tcl==a]<-ifelse(comb[i,"dist"] < dcut,b,a)
          
          
          
        }
      }
      #
    }
    #
    
    
    ####Amplification/deletion matrix
    #Matrix of difference to reference per subclone and arm
    
    if(merge_method=="ad"){tcl<-tcl2}
    if(merge_method=="ad" | merge_method=="both"){
      if(length(unique(tcl))>1){
        stringmin=0
        while(stringmin<1 && length(unique(tcl))>1){
          
          comb<-data.frame(t(combn(unique(tcl), 2)))
          ad1<-admat.fun(matrix,br_vec,tcl,labels,cut=adcut,genes=adgenes,region=adregion)  
          epmat<-ad1$epmat
          ad<-ad1$ad
          if(isTRUE(armfilt)){
            ad<-varfilt(epmat,ad,cut=NULL,n=arm_n)$ad
            epmat<-varfilt(epmat,ad,cut=NULL,n=arm_n)$epmat}
          for(i in 1:nrow(comb)){
            
            a<-as.character(comb[i,1])
            b<-as.character(comb[i,2])
            comb[i,"stringdist"]<-sum(stringdist::stringdist(ad[,a],ad[,b]))
            tcl[tcl==a]<-ifelse(comb[i,"stringdist"]==0,b,a)
            stringmin<-min(comb$stringdist)
            
          }
        }
      }
    }
  }
  
  
  
  #Get back to full ampdel matrix if using only top arms
  if(isTRUE(armfilt) || isFALSE(merge) & length(unique(tcl))>1){
    ad1<-admat.fun(matrix,br_vec,tcl,labels,cut=adcut,genes=adgenes,region=adregion)  
    epmat<-ad1$epmat
    ad<-ad1$ad}
  
  
  #Rename merged clusters
  if(length(unique(tcl))>1){
    
    sn<-sort(unique(tcl))
    
    ad2<-ad[,sn]
    em2<-epmat[,sn]
    l2<-LETTERS[1:length(sn)]
    for(i in 1:length(sn)){
      tcl<-ifelse(tcl==sn[i],paste("subclone",l2[i],sep = "_"),tcl)
    }
    tcl<-paste(name,tcl,sep = "_")
    names(tcl)<-rownames(tm1)
    #Print ampdel matrix
    colnames(ad2)<-sort(unique(tcl))
    colnames(em2)<-sort(unique(tcl))
    
    
    out<-list(ad=ad2,tcl=tcl,epmat=em2,distmat=tm3)}
  
  if(length(unique(tcl))<2){tcl<-paste(name,"subclone1",sep = "_")
  out<-tcl}
  
  
  return(out)
}
#

##Create amplification/deletion matrix
admat.fun<-function(matrix,breaks,vector,lab,cut=0.15,genes="top",region=1/2,refmat=NULL){
  bn<-unique(breaks)
  tn<-sort(unique(vector))
  #Split matrix per subclone
  cut.mat<-lapply(tn,function(x){
    n1<-names(vector[vector==x])
    m1<-matrix[,n1]
  })
  #Mean arm expression per subclone
  epmat<-matrix(nrow = length(bn),ncol = length(tn))
  for(i in 1:length(tn)){
    for(j in 1:length(bn)){
      m1<-cut.mat[[i]]
      m2<-m1[breaks==bn[j],]
      absvals<-apply(m2,1,function(x)mean(abs(x)))
      m3<-m2[absvals>quantile(absvals, prob = 1-region),]
      if(!is.null(dim(m3))){c<-mean(colMeans(m3))}else{c=0}
      d<-mean(colMeans(m2))
      if(genes=="top"){epmat[j,i]<-c}
      if(genes=="all"){epmat[j,i]<-d}
      if(genes=="ref" && !is.null(refmat)){
        mr<-refmat[breaks==bn[j],!colnames(refmat) %in% cells]
        mc<-mean(colMeans(mr))
        epmat[j,i]<-c-mc}
      
    }
  }
  rownames(epmat)<-lab
  colnames(epmat)<-tn
  #Amplification/deletion matrix
  ad<-matrix(ncol = length(tn),nrow = length(bn))
  
  for(i in 1:length(bn)){
    for(j in 1:length(tn)){
      ad[i,j]<-ifelse(epmat[i,j]>cut,"amp",ifelse(epmat[i,j]< -1*cut,"del","none"))
    }}
  rownames(ad)<-lab
  colnames(ad)<-tn
  
  ad<-as.data.frame(ad)
  ad[is.na(ad)]<-"none"
  epmat<-as.data.frame(epmat)
  epmat[is.na(epmat)]<-0
  out<-list(ad=ad,epmat=epmat)
  return(out)
}
##Filter arms by variance cutoff or top variant arms
varfilt<-function(epmat,ad,cut=0.01,n=10){
  armvar<-apply(epmat,1,var)
  if(is.null(cut)){
    z<-names(tail(sort(armvar),n))
    epmat<-epmat[z,]
    ad<-ad[z,]}
  if(!is.null(cut)){
    epmat<-epmat[armvar>cut,]
    ad<-ad[armvar>cut,] 
  }
  out<-list(ad=ad,epmat=epmat)
  return(out)
}

####Reassign subclones
reassign.clones<-function(matrix,
                          cells,
                          admat,
                          clones,
                          distmat,
                          em,
                          tscore,
                          define_region=1/3,
                          define_cells=1/4,
                          goldstandard=NULL,
                          clonesize=10,
                          knnprob=0.9,
                          hpvcut=0.1,
                          separate="arm",
                          genecut=10,
                          adgenes="all",
                          adcut=0.2,
                          adregion=1,
                          name="x",
                          refilt=FALSE,
                          gt=NULL,
                          noHPV=FALSE,
                          cnasig="abs",
                          order_clones=FALSE,
                          genome=NULL,
                          spec=0.99
){
  ##0. Subset by chromosome, needed for new admat and epmat
  br1<-br.fun(genes=rownames(matrix),separate = separate,
              genecut=genecut,genome=genome)
  br_vec<-br1$br_vec
  labels<-br1$bdf$labels
  rm<-br1$rm
  
  if(length(rm)!=0){
    matrix<-matrix[-rm,]}
  bn<-unique(br_vec)
  
  #####1. Find normal subclone
  ad2<-admat
  ad2<-ifelse(ad2=="none",0,1)
  ad2<-data.frame(ad2)
  ad2$w<-0
  ad2$q<-1
  
  d<-colSums(ad2)
  ld<-length(d[d==0])
  if(ld==2){
    normclone<-colnames(ad2[,colSums(ad2)==0])[1]
  }
  
  #Select weakest signal clone if several dont pass ampdel cut
  if(ld>2){
    ab1<-sort(colSums(abs(em)),decreasing = F)
    normclone<-names(ab1)[1]
  }
  #If no normal clones exist set normclone as null
  if(ld==1){normclone=NULL}
  
  #Reference cells
  reference= colnames(matrix[,!colnames(matrix) %in% cells])
  rx<-0
  #Name for potential new cluster
  new<-paste(name,"subclone",LETTERS[length(unique(clones))+1],sep = "_")
  #If applying to tumour type without HPV classification
  if(isTRUE(noHPV)){
    HPV<-character(length = length(colnames(matrix)))
    names(HPV)<-colnames(matrix)
    for(i in 1:length(HPV)){HPV[i]<-"HPV-"}
  }
  if(!is.null(normclone)){normclone<-gsub(".","-",normclone,fixed = T)}
  #####Possibly remove below functions, to keep normal HPVpos cells
  
  #Set aside HPV+ cells as true cancer
  # if(!is.null(normclone)){
  # h1<-HPV[names(clones[clones==normclone])]
  # hp<-names(h1[h1=="HPV+"])
  # hn<-names(h1[h1=="HPV-"])
  # }
  # #Dont set HPV-rich clones as normal
  # if(!is.null(normclone) && length(hp)/length(c(hp,hn)) > hpvcut){normclone=NULL}
  ####2. Set true normal as normal, true cancer as cancer, rest unresolved
  
  
  #Define CNA with all cells, to set proper correlation cutoffs for the normal clone.
  #Also print values for the reference cells
  normdef<-define.cna(matrix,cells,reference,tumour_score=tscore,cnasig=cnasig,
                      top_cells=define_cells,top_region=define_region,extra_score=NULL,gs=NULL,
                      uppercut=2,lowercut=2,reass=FALSE,no_ts=TRUE,name="All cells",glob_thr = gt,
                      print_ref = TRUE,spec_cut = spec)
  
  normplots<-normdef$plots
  nts<-normdef$thresholds
  nd1<-normdef$def_mat
  refdef<-nd1[reference,]
  
  if(!is.null(normclone)){
    normdef<-nd1[names(clones[clones==normclone]),]
    
    
    #Reomve - keep HPVpos normals
    #can1<-c(hp,rownames(normdef[normdef$cna3=="Cancer",]))
    #cancer<-unique(can1)
    
    normal<-rownames(normdef[normdef$cna3=="Normal",])
    cancer<-rownames(normdef[normdef$cna3=="Cancer",])
    unres<-rownames(normdef[normdef$cna3=="Unresolved",])
    clones[unres]<-paste(name,"Unresolved",sep = "_")
    
    ####3. Reassign newly defined cancer cells to subclones
    if(length(cancer)>5){
      
      unres_cells<-cancer
      tm1<-distmat[,!colnames(distmat) %in% c("subclone","knn","prob")]
      #
      knn_sc<-class::knn(train = tm1[!rownames(tm1) %in% unres_cells,],
                         test = tm1[rownames(tm1) %in% unres_cells,],k=max(1,round(log(nrow(tm1)))),
                         cl=clones[!names(clones) %in% unres_cells],prob=T)
      #
      normdef[unres_cells,"sc1"]<-clones[unres_cells]
      normdef[unres_cells,c("sc2","scprob")]<-c(as.character(knn_sc),attr(knn_sc,"prob"))
      
      normdef$sc3<-ifelse(rownames(normdef) %in% unres_cells & normdef$scprob>=knnprob,normdef$sc2,
                          ifelse(rownames(normdef) %in% unres_cells & normdef$scprob<knnprob,new,normdef$sc1))
      #
      clones[unres_cells]<-normdef[unres_cells,"sc3"]
      #If new clone too small set as unresolved
      if(length(clones[clones==new]) < clonesize){clones[clones==new]<-paste(name,"Unresolved",sep = "_")}
    }
    
    normdef<-normdef[names(clones[clones==normclone | 
                    clones==paste(name,"Unresolved",sep = "_")]),]
    rx<-length(rownames(normdef[normdef$cna3=="Normal",]))
  }
  
  
  ####4. Check for normal cells in other subclones
  cn<-colnames(ad2[,!colnames(ad2) %in% c("q","w")])
  cn<-gsub(".","-",cn,fixed = T)
  if(!is.null(normclone)){cn<-cn[-grep(normclone,cn)]}
  
  
  if(length(clones[clones==new])>=clonesize){cn<-c(cn,new)}
  #Use gold standard in HPV+ patients
  
  #Run define CNA, skip tumour score, substract cnacor from normal cells
  normcheck<-lapply(cn,function(x){
    a<-define.cna(matrix,names(clones[clones==x]),reference,tumour_score=tscore,
                  top_cells=define_cells,top_region=define_region,uppercut=2,
                  gs=goldstandard,extra_score=NULL,glob_thr = gt,
                  lowercut=2,reass=ifelse(is.null(normclone)|rx<clonesize,FALSE,TRUE),
                  no_ts=TRUE,nd=normdef,name=x,spec_cut=spec)
    
    a$def_mat<-a$def_mat[,!colnames(a$def_mat) %in% c("knn","prob")]
    out<-list(defmat=a$def_mat,plots=a$plots,thresholds=a$thresholds)
    out})
  #
  normcheck<-do.call(rbind,normcheck)
  if(length(cn)>1){
    defmat<-do.call(rbind.data.frame,normcheck[,1])}else{defmat<-normcheck[[1,1]]}
  plots<-normcheck[,2]
  thresholds<-do.call(rbind,normcheck[,3])
  rownames(thresholds)<-cn
  
  #Skip this, keep HPVpos normals as normals
  #defmat$cna3<-ifelse(rownames(defmat) %in% names(HPV[HPV=="HPV+"]) & defmat$cna3 == "Normal", "Unresolved",defmat$cna3)
  
  
  #Rerun normcheck if there wasn't any first normal clone to substract from CNAcor
  if(is.null(normclone) & nrow(defmat[defmat$cna3 == "Normal",]) >= clonesize & isTRUE(refilt)){
    n2<-rownames(defmat[defmat$cna3=="Normal",])
    defmat<-defmat[defmat$cna3 != "Normal",]
    clones[n2]<-paste(name,"Normal",sep = "_")
    normclone<-unique(clones[n2])
    normdef<-define.cna(matrix,cells,reference,tumour_score=tscore,glob_thr = gt,
                        top_cells=define_cells,top_region=define_region,extra_score=NULL,gs=goldstandard,uppercut=2,
                        lowercut=2,reass=FALSE,no_ts=TRUE,name="All cells",spec_cut=spec)
    normplots<-normdef$plots
    nts<-normdef$thresholds
    normdef<-normdef$def_mat
    normdef<-normdef[n2,]
    normdef$cna3<-"Normal"
    
    normcheck<-lapply(cn,function(x){
      a<-define.cna(matrix,names(clones[clones==x]),reference,tumour_score=tscore,
                    top_cells=define_cells,top_region=define_region,uppercut=2,
                    gs=goldstandard,extra_score=NULL,glob_thr = gt,
                    lowercut=2,reass=ifelse(is.null(normclone),FALSE,TRUE),no_ts=TRUE,nd=normdef,name=x,spec_cut = spec)
      
      a$def_mat<-a$def_mat[,!colnames(a$def_mat) %in% c("knn","prob")]
      out<-list(defmat=a$def_mat,plots=a$plots,thresholds=a$thresholds)
      out})
    
    #
    normcheck<-do.call(rbind,normcheck)
    if(length(cn)>1){
      defmat<-do.call(rbind.data.frame,normcheck[,1])}else{defmat<-normcheck[[1,1]]}
    plots<-normcheck[,2]
    thresholds<-do.call(rbind,normcheck[,3])
    rownames(thresholds)<-cn
    
    
    #defmat$cna3<-ifelse(rownames(defmat) %in% names(HPV[HPV=="HPV+"]) & defmat$cna3 == "Normal", "Unresolved",defmat$cna3)
  }
  
  ur2<-rownames(defmat[defmat$cna3=="Unresolved",])
  n2<-rownames(defmat[defmat$cna3=="Normal",])
  clones[ur2]<-paste(name,"Unresolved",sep = "_")
  if(!is.null(normclone)){clones[n2]<-normclone
  clones[clones==normclone]<-paste(name,"Normal",sep = "_")
  normdef<-normdef[,colnames(defmat)]
  refdef<-refdef[,colnames(defmat)]
  defmat<-rbind.data.frame(defmat,normdef,refdef)
  }
  if(is.null(normclone)){
    clones[n2]<-paste(name,"Normal",sep = "_")
    normplots<-plot(sort(runif(100)))
    nts<-c(10,10,10)
  }
  
  #Remove clones that were completely filtered out
  r1<-setdiff(cn,clones)
  thresholds<-thresholds[!rownames(thresholds) %in% r1,]
  
  ####Set small clones as unresolved and remove small clone thresholds
  c1<-unique(clones)
  c2<-c1[grep("Normal|Unresolved",c1,invert = T)]
  l1<-lapply(c2,function(x)length(clones[clones==x]))
  l2<-unlist(l1)
  names(l2)<-c2
  l3<-names(l2[l2<clonesize])
  if(length(c2)>=2){
    thresholds<-thresholds[!rownames(thresholds) %in% l3,]
  }
  clones[clones %in% l3]<-paste(name,"Unresolved",sep = "_")
  
  
  
  #Rename clones
  c1<-unique(clones)
  c2<-c1[grep("Normal|Unresolved",c1,invert = T)]
  
  for(i in 1:length(c2)){
    if(isTRUE(order_clones) & length(unique(clones))>4){
    clones[clones %in% c2[i]] <- paste(name,"subclone",letters[i],sep = "_")
    }else{
      clones[clones %in% c2[i]] <- paste(name,"subclone",LETTERS[i],sep = "_")
    }
  }
  
  #Save thresholds
  thresholds<-rbind.data.frame(nts,thresholds)
  colnames(thresholds)<-c("cor","sig","tscore")
  c1<-unique(clones)
  c2<-c1[grep("Normal|Unresolved",c1,invert = T)]
  rownames(thresholds)<-c(paste0(name,"_Normal"),sort(c2))
  #
  
  ###New epmat and admat
  if(adgenes=="ref"){refmat<-matrix[,reference]
  cells<-cells}
  
  ad1<-admat.fun(matrix,br_vec,clones[clones!=paste0(name,"_Normal") & clones!=paste0(name,"_Unresolved")],labels,cut=adcut,genes=adgenes,region=adregion)  
  epmat<-ad1$epmat
  ad<-ad1$ad
  
  if(isTRUE(order_clones) & length(unique(clones))>4){
    c1<-Heatmap(t(epmat),cluster_columns = FALSE)
    c2<-colnames(epmat[,row_order(c1)])
      for(i in 1:length(c2)){
        clones[clones == c2[i]] <- paste(name,"subclone",LETTERS[i],sep = "_")
      }
    # c1<-unique(clones)
    # c2<-c1[grep("Normal|Unresolved",c1,invert = T)]
    # 
    # for(i in 1:length(c2)){
    #   clones[clones %in% c2[i]] <- paste(name,"subclone",LETTERS[i],sep = "_")
    # }
    ad1<-admat.fun(matrix,br_vec,clones[clones!=paste0(name,"_Normal") & clones!=paste0(name,"_Unresolved")],labels,cut=adcut,genes=adgenes,region=adregion)  
    epmat<-ad1$epmat
    ad<-ad1$ad
          }
    
    
  defmat<-defmat[intersect(rownames(defmat),colnames(matrix)),]
  
  #
  ###Print plots, admat, epmat, matrix  
  out<-list(plots=plots,normplots=normplots,defmat=defmat,ad=ad,epmat=epmat,clones=clones,thresholds=thresholds)
  
  return(out)
}
#######TCGA scores
tcga.types<-function(matrix=merged,cells,ref,centre=FALSE,log=FALSE,nofilt=FALSE,centref=NULL,min_umi=5,
                     min_cells=20,min_avlog=4,...){
  
  #Select only cancer cells
  c2<-cells
  #Put reference in order
  ref<-ref[!duplicated(ref$Gene),]
  rownames(ref)<-ref$Gene
  m1<-matrix[rownames(matrix) %in% rownames(ref),]
  if(isFALSE(nofilt)){
    #Transform values to log or TPM, optional centering
    mat<-filter.umi(matrix=m1,cells=c2,ref=centref,...)
  }
  #Option where matrix isnt filtered
  if(isTRUE(nofilt)){mat<-m1[,c2]}
  #Calculate values per subtype
  s1<-unique(ref$Subtype)
  sf<-lapply(s1,function(x){
    
    r1<-rownames(ref[ref$Subtype==x,])
    m1<-mat[rownames(mat) %in% r1,]
    a1<-mean(colMeans(m1))
  })
  s2<-unlist(sf)
  names(s2)<-s1
  return(s2)
}
####TCGA Plot
tcga.plot<-function(matrix,scores=NULL,scorenames,...){
  
  mat<-matrix
  #Function to get main colours
  col_fun = circlize::colorRamp2(c(min(mat), median(mat), max(mat)), c("steelblue", "white", "darkred"))
  #Use arbitrary number of scores for annotating
  
  if(is.list(scores)){
    nn<-seq(1,length(scores))
    sfun<-lapply(nn,function(x){
      
      s1<-scores[[x]]
    })
    scoremat<-do.call(cbind.data.frame,sfun)
    colnames(scoremat)<-scorenames
    row_ann = rowAnnotation(df=scoremat)
  }
  if(is.null(scores)){row_ann=NULL}
  
  pp2<-invisible(Heatmap(mat, name = "Relative\nProgram Scores",cluster_columns=FALSE,
                         col = col_fun,left_annotation = row_ann,...))
  return(pp2)
}

####Correlation plot
corr.plot<-function(matrix,cormet="pearson",scores=NULL,scorenames=NULL,name=NULL, n=50,typeplot=TRUE,numscale=NULL){
  
  
  
  m2<-hclustreorder(cor(matrix,method = cormet))
  
  #Function to get main colours
  col_fun = circlize::colorRamp2(c(-0.4, 0, 0.4), c("steelblue", "white", "darkred"))
  #Use arbitrary number of scores for annotating
  
  cells=colnames(m2)
  #Make column index with breaks
  nam <- rep(" ", ncol(m2))
  k<- plyr::round_any( (ncol(m2)/10),n)
  nam[seq(n,ncol(m2), k)] <- seq(n, ncol(m2), k)
  colnames(m2) <- nam
  #Make heatmap
  
  if(is.list(scores)){
    nn<-seq(1,length(scores))
    sfun<-lapply(nn,function(x){
      
      s1<-scores[[x]][cells]
    })
    scoremat<-do.call(cbind.data.frame,sfun)
    colnames(scoremat)<-scorenames
    row_ann = rowAnnotation(df=scoremat)
    
    if(isTRUE(typeplot)){
      cf2<-lapply(scorenames,function(x){
        
        a1<-as.character(x)
        a<-circlize::colorRamp2(c(min(scoremat[,a1]),median(scoremat[,a1]),max(scoremat[,a1])),c("#00AFBB", "white", "#FC4E07"))
        a})
      names(cf2)<-scorenames
      
      row_ann = rowAnnotation(df=scoremat,col=cf2)
    }
    if(!is.null(numscale)){
      
      q<-numscale
      cf2<-lapply(scorenames,function(x){
        
        a1<-as.character(x)
        a<-circlize::colorRamp2(c(-1*q,0,q),c("#00AFBB", "white", "#FC4E07"))
        a})
      names(cf2)<-scorenames
      
      row_ann = rowAnnotation(df=scoremat,col=cf2)
    }
    
  }
  if(is.null(scores)){row_ann=NULL}
  #Row annotation
  #TS = tumour_score[cells]
  
  pp2<-invisible(Heatmap(m2, name = "Pearson\nCorrelation",
                         cluster_rows=FALSE,cluster_columns = FALSE, show_row_names = FALSE, show_column_names = TRUE,col = col_fun,
                         left_annotation = row_ann,column_title = name))
  return(pp2)
}
#

#############Jaccard barcode similarity
jac.barcode<-function(vector){

w<-vector  

comb<-data.frame(t(combn(unique(w), 2)))
for(i in 1:nrow(comb)){
  s1<-as.character(comb[i,1])
  s2<-as.character(comb[i,2])
  
  n1<-names(w[w==s1])
  n2<-names(w[w==s2])
  n11<-sapply(strsplit(n1, split='_', fixed=TRUE),function(x)(x[2]))
  n12<-sapply(strsplit(n2, split='_', fixed=TRUE),function(x)(x[2]))
  comb[i,"jac"]<-(length(intersect(n11,n12))/length(union(n11,n12)))
}
return(comb[order(comb$jac,decreasing = T),])
}


#Differentially expressed genes
DE.genes<-function(matrix,vector,no=50,just="BH",max_pv=0.05,min_fc=2,volcano=FALSE,
                   pairwise=FALSE){
  
  clusterlist<-unique(vector)  
  comparisons <- mclapply(clusterlist, function(x){
    
    if(isFALSE(pairwise)){
      logfc<- rowMeans(matrix[,vector==x]) - rowMeans(matrix[,vector!=x])
      b<-names(logfc)
      if(!is.null(max_pv)){
        pval  <- apply(matrix,1,function(y) t.test(y[vector==x],y[vector!=x])$p.value)
        padj<- p.adjust(pval, method = just)
        a<-cbind.data.frame(logfc,padj,x)
        c<-cbind.data.frame(b,a)
        d<-subset(c, c$logfc > min_fc & c$padj < max_pv)
        
        colnames(d)<-c("Gene","logFC","adj.pval","Cluster")
      }else{
        d<-cbind.data.frame(b,logfc,x)
        d<-d[d$logfc> min_fc,]
        colnames(d)<-c("Gene","logFC","Cluster")
      }
    }else{
      cl2<-clusterlist[clusterlist!=x]
      fx<-mclapply(cl2,function(y){
        logfc<- rowMeans(matrix[,vector==x]) - rowMeans(matrix[,vector==y])
        b<-names(logfc)
        if(!is.null(max_pv)){
          pval  <- apply(matrix,1,function(z) t.test(z[vector==x],z[vector==y])$p.value)
          padj<- p.adjust(pval, method = just)
          a<-cbind.data.frame(logfc,padj,x,y)
          c<-cbind.data.frame(b,a)
          d<-subset(c, c$logfc > min_fc & c$padj < max_pv)
          
          colnames(d)<-c("Gene","logFC","adj.pval","Cluster","Comparison")
        }else{
          d<-cbind.data.frame(b,logfc,x,y)
          d<-d[d$logfc> min_fc,]
          colnames(d)<-c("Gene","logFC","Cluster","Comparison")
        }
        d})
      fx2<-do.call(rbind.data.frame,fx)
      gl<-fx2 %>% count(Gene) %>% filter(n==length(cl2)) %>% pull(Gene)
      d<-fx2[fx2$Gene %in% gl,-ncol(fx2)] %>% group_by(Gene,Cluster) %>% summarise_all(mean) %>%
        as.data.frame()
    }
    
    if(!is.null(no)){
      e<-head(d[order(d$logFC, decreasing = T),], n=no)}else{
        e<-d[order(d$logFC, decreasing = T),]
      }
    
    e
    
  })
  
  deg<- do.call(rbind,comparisons)
  if(nrow(deg)>0){
    rownames(deg)<-seq(1,nrow(deg))}
  if(isTRUE(volcano)){
    deg[deg$Cluster==clusterlist[2],"logFC"]<- -1*deg[deg$Cluster==clusterlist[2],"logFC"]
  }
  
  return(deg)
}
#
#Calculate Jaccard index for each pairwise cluster comparison
jac.mat<-function(cluster_matrix,list=NULL,matrix=NULL,vector,add_diagonal=FALSE,
                  jacmin=NULL){
  
  jac_mat <- data.frame(matrix( ncol = 6)) #Create empty matrix
  colnames(jac_mat)<-c("C1","C2","Jaccard","L1","L2","Diff")
  v1<-as.character(cluster_matrix$Cluster)
  combs <- data.frame(t(combn(unique(v1), 2)), stringsAsFactors = F) #Al
  
  if(is.null(list) & is.null(matrix)){jac_mat<-jac_mat[,1:3]}
  
  for(i in 1:nrow(combs)){
    
    
    a1<-cluster_matrix[v1==combs[i,1],]
    a2<-cluster_matrix[v1==combs[i,2],]
    a1$Gene<-as.character(a1$Gene)
    a2$Gene<-as.character(a2$Gene)
    
    j1<-length(intersect(a1$Gene,a2$Gene))/length(union(a1$Gene,a2$Gene)) #Jaccard similarity
    
    if(!is.null(list)){
      b1<-unlist(strsplit(combs[i,1],split="_"))[3]
      b2<-unlist(strsplit(combs[i,2],split="_"))[3]
      jac_mat[i,4]<-length(list[[b1]] )
      jac_mat[i,5]<-length(list[[b2]]) #Cluster 2 size
      jac_mat[i,6]<-jac_mat[i,4]-jac_mat[i,5] 
    }
    if(!is.null(matrix)){
      jac_mat[i,4]<-ncol(matrix[,vector==combs[i,1]]) #Cluster 1 size
      jac_mat[i,5]<-ncol(matrix[,vector==combs[i,2]]) #Cluster 2 size
      jac_mat[i,6]<-jac_mat[i,4]-jac_mat[i,5] 
    }
    
               #Size difference
    jac_mat[i,3]<-j1
    jac_mat[i,2]<-combs[i,2]
    jac_mat[i,1]<-combs[i,1]
  }
  
  if(!is.null(jacmin)){
    zn<-unique(v1)
    jf<-lapply(zn,function(x){
      jac_mat %>% filter(C1==x | C2==x) %>% pull(Jaccard) %>% max()
          })
    zn2<-cbind.data.frame(zn,unlist(jf))
    colnames(zn2)<-c("Cluster","Jaccard")
    l1<- zn2 %>% filter(Jaccard<jacmin) %>% pull(Cluster)
    
  }
  
  
  if(isTRUE(add_diagonal)){
    zn<-unique(v1)
    zqq<-matrix(ncol=ncol(jac_mat),nrow=length(zn))
    colnames(zqq)<-colnames(jac_mat)
    for(i in 1:length(zn)){
      zqq[i,1]<-zn[i]
      zqq[i,2]<-zn[i]
      zqq[i,3]<-1
      if(ncol(jac_mat)==6){
        zqq[i,4]<-5000
        zqq[i,5]<-5000
        zqq[i,6]<-0
      }
      
      
      }
    jac_mat<-rbind(zqq,jac_mat)
    jac_mat$Jaccard<-as.numeric(jac_mat$Jaccard)
  }
  if(!is.null(jacmin)){
  jac_mat<-jac_mat %>%filter(!C1 %in% l1 & !C2 %in% l1)
  }
  
  jac_mat<-jac_mat[order(jac_mat$Jaccard,decreasing = T),]
  
  return(jac_mat)
}
####Clean clusters from cutree
hclean<-function(cut,min=10,max=0.5){
  
  k <- lapply(1:ncol(cut), function(i) {
    tapply(names(cut[,i]), cut[,i], c)
  })
  # Flatten list. Cell names to top level.
  k <- unlist(k, recursive=FALSE)
  
  # remove duplicate clusters
  k <- k[!duplicated(k)]
  
  # remove clusters smaller than min val
  k <- k[lapply(k, length) >= min]
  
  # remove clusters bigger than max val
  k <- k[lapply(k, function(i) length(i)/ncol(cut)) <= max]
  
  names(k)<-1:length(k)
  return(k)
}
#
#####Extract metaprograms from clusters




extract.metaprogs<-function(matrix,vector,name,clusters=FALSE,no=50){
  
  
  c1<-unique(vector)
  
  if(length(grep("_",vector))>0){
  v1 <- sapply(strsplit(vector, split='_', fixed=TRUE),function(x)(x[1]))
  }else{v1<-vector}
  ###Create matrix of top50 genes per cluster in metaprogram
  genes<-lapply(c1,function(x){
    a1<- matrix[vector == x,]
    a2<-as.character(I(a1$Gene))
  })
  
  for(i in 1:length(genes)){
    a<- no-length(genes[[i]])
    if(a>0){
      b<-character(length = a)
      for(j in 1:length(b)){
        b[j]<-"Temp"
      }
      genes[[i]]<-c(genes[[i]],b)
    }
    
  }
  names(genes)<-c1

  genemat<-do.call(cbind.data.frame,genes)
  
  colnames(genemat)<-names(genes)
  
  
  ###Find common genes
  c2<-as.character(unique(matrix[vector %in% colnames(genemat),"Gene"]))
  genefreq<-lapply(c2,function(y){
    
    
    ax<-paste("\\b",y,"\\b",sep = "")
    a1<-length(grep(ax,as.matrix(genemat)))
    a2<-length(grep(ax,as.matrix(genemat)))/length(colnames(genemat))
    a4<-cbind.data.frame(y,a1,a2)
    a4$Cluster<-name
    colnames(a4)<-c("Gene","Number","Frequency","Cluster")
    if(isTRUE(clusters)){
      a3<-colnames(genemat)[apply(genemat, 2, function(col) any(grepl(ax, col)))]
      a4$clusterlist<-list(as.character(I(a3)))
    }
    a4
  })
  freqmat<-do.call(rbind,genefreq)
  f2<-freqmat[order(freqmat$Frequency,decreasing = T),]
  f3<-f2[f2$Number>0,]
  
  #Add mean logFC
  lf <- matrix %>% filter(Gene %in% f3$Gene) %>% group_by(Gene) %>% summarise(mean(logFC))
  for(i in 1:nrow(lf)){
    f3[f3$Gene==as.character(lf[i,1]),"logFC"]<-as.numeric(lf[i,2])
      }
  
  
  if(isTRUE(clusters) ){
    #Frequency in patients, not clusters
    qw<-f3$clusterlist
    for(i in 1:length(qw)){
      c<-strsplit(qw[[i]],split = "_")
      a<-do.call(rbind,c)
      f3[i,"Patfreq"]<-length(unique(a[,1]))/length(unique(v1))
      f3[i,"Pat.No"]<-length(unique(a[,1]))
    }
    f3<-f3[order(f3$Patfreq,decreasing = T),c(1,2,3,8,7,6,4,5)]
    f3$clusterlist <- vapply(f3$clusterlist, paste, collapse = ", ", character(1L))
  }
  out<-list(genemat=genemat,topgenes=f3)
  return(out)
}
###Gene expression clustering
#Cluster coordinate matrix with cells as rows
cluster.coord<-function(matrix,method="dbscan",louvain=25,...){
  
  if(method=="dbscan"){
    cluster<-dbscan.cluster(matrix,...)
  }
  if(method=="louvain"){
    k=louvain
    knn.norm = FNN::get.knn(as.matrix(matrix), k = k)
    knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),k), to = as.vector(knn.norm$nn.index), 
                          weight = 1/(1 +as.vector(knn.norm$nn.dist)))
    nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
    nw.norm = igraph::simplify(nw.norm)
    lc.norm = igraph::cluster_louvain(nw.norm)
    cluster<-igraph::membership(lc.norm)}
  
  ren.fun<-function(x){paste("cluster",x,sep = "")}
  cluster<-ren.fun(cluster)
  names(cluster)<-rownames(matrix)
  
  return(cluster)
}

####Merge clusters by jaccard similarity
jac.merge<-function(matrix,cv,cutoff=0.3,size=1000,rename=TRUE,mpv=0.05,geneno=50,lfc=1){
  
  s1<-Sys.time()
  jacmax=1
  sizemin=1
  
  while(jacmax>cutoff & sizemin < size){
    
    print("Starting merging round")
    
    c1<-unique(cv)
    
    deg<-parallel::mclapply(c1,mc.cores=length(c1),function(x){
      
      
      v<-ifelse(colnames(matrix) %in% names(cv[cv==x]),"Y","N")
      
      d<-DE.genes(matrix=matrix,vector=v,max_pv = mpv,no=geneno,min_fc = lfc)
      d<-d[d$Cluster=="Y",]
      d$Cluster<-x
      d})
    
    print("Differential Expression done")
    
    d3<-do.call(rbind,deg)
    rownames(d3)<-1:nrow(d3)
    
    jac_mat<-jac.mat(cluster_matrix=d3,list=NULL,matrix=matrix,vector=cv)
    
    print("Jaccard matrix formed")
    
    jacmax<-max(jac_mat$Jaccard)
    if(jacmax>cutoff){
    sizemin<-min(jac_mat[jac_mat$Jaccard > cutoff,c(4:5)])
    }else{sizemin<-1}
    
    for (i in 1:nrow(jac_mat)) {
      
      
      if (jac_mat[i, 3] > cutoff & jac_mat[i,5]<size & jac_mat[i,6]>=0) {
        cv[cv==jac_mat[i,2]] <-jac_mat[i,1]
        
      }
      
      if (jac_mat[i, 3] > cutoff & jac_mat[i,4]<size  & jac_mat[i,6]<0) {
        cv[cv==jac_mat[i,1]] <-jac_mat[i,2]
        
      }
    }
    s2<-Sys.time()
    if(jacmax>cutoff & sizemin < size){
    print( paste(length(c1),"clusters merged to",length(unique(cv)),sep=" "))
    print(s2-s1)
    }
    print(paste0("Max similarity ",jacmax))
    print(paste0("Min size of cluster over cutoff is ",sizemin))
  }
  
  if(isTRUE(rename)){
  for(i in 1:length(c1)){
    d3[d3$Cluster==c1[i],"Cluster"]<-paste0("Cluster_",i)
    cv[cv==c1[i]]<-paste0("Cluster_",i)
    
  }
  }
  print("Finally done!")
  print(s2-s1)
 
  
  
  
  out<-list(cluster=cv,clustermat=d3)
  return(out)
}

#####Percentage of cell types in cluster
type.perc<-function(matrix,cv=Cluster,tv=Type){
  
  c1<-unique(cv)
  t1<-unique(tv)
  c2<-lapply(c1,function(x){
    
    a1<-matrix[,cv==x]
    
    
    t2<-parallel::mclapply(t1,mc.cores=length(t1),function(x){
      b1<-a1[,colnames(a1) %in% names(tv[tv==x])]
      b2<-length(colnames(b1))/length(colnames(a1))
      b2})
    t3<-unlist(t2)
    t3<-c(t3,length(colnames(a1)))
  })
  o1<-do.call(rbind.data.frame,c2)
  rownames(o1)<-c1
  colnames(o1)<-c(t1,"Size")
  return(o1)
}

####Add scores
score<-function(matrix,genes,bin=100,centre=TRUE,no_bins=30){
  
  
  genes<-genes[genes %in% rownames(matrix)]
  a1<-matrix[rownames(matrix) %in% genes,]
  if(is.null(bin)){score<-colMeans(a1)}
  if(!is.null(bin)){
    
    aa<-rowMeans(matrix)
    bb<-ntile(aa,no_bins)
    names(bb)<-names(aa)
    
    fc<-mclapply(genes,function(x){
      a2<-bb[x]
      b2<-sample(names(bb[bb==a2]),size = bin,replace = F)
      a3<-matrix[x,]
      b3<-colMeans(matrix[b2,])
      a3 -b3
    })
    fc2<-do.call(rbind,fc)
    score<-colMeans(fc2)
  }
  
  
  if(isTRUE(centre)){score<-score - mean(score)}
  return(score)
}

#####Consensus clustering of genes x cells matrix
consclust<-function(matrix,maxk=10,r=100,setk=NULL,dist="pearson",size=10,name=NULL,
                    extraplot=FALSE,plot=TRUE,geneno=50){
  
  if(is.null(name)){name="x"}  
  
  m1<-as.matrix(matrix)
  ####Compute consensus clustering for a range of k
  cc2<-invisible(ConsensusClusterPlus::ConsensusClusterPlus(m1,maxK = maxk,reps=r,distance = dist))
  if(is.null(setk)){
  #Compute item consensus
  icl<-invisible(ConsensusClusterPlus::calcICL(cc2))
  
  #Table of cluster and item consensus values for range of k
  ccs <- left_join(as_tibble(icl$clusterConsensus) %>% group_by(k) %>% summarise(CC = mean(clusterConsensus, na.rm = TRUE)),
                   
                   as_tibble(icl$itemConsensus) %>% group_by(item, k) %>% summarise(MIC = max(itemConsensus)) 
                   %>% group_by(k) %>% summarise(IC = mean(MIC, na.rm = TRUE)),by="k")
  
  ccx<-as.data.frame(ccs)
  ccx<-ccx[ccx$k!=2,]
  #Take number of clusters that maximises cluster*item consensus
  ccx$score<-ccx$CC*ccx$IC
  mk<-ccx$k[ccx$score==max(ccx$score)]
  }
  if(!is.null(setk)){mk<-setk}
  #Cluster vector with selected k
  clusters <- cc2[[mk]]$consensusClass
  #Remove small clusters
  cu<-unique(clusters)
  for(i in 1:length(cu)){
    
    if(length(clusters[clusters==cu[i]])<size){
      clusters[clusters==unique(clusters)[i]]<-"X"}
    
  }
  #
  cu2<-unique(clusters[clusters!="X"])
  for(i in 1:length(cu2)){
    clusters[clusters==cu2[i]]<-paste0(name,"_Cluster_",LETTERS[i])
  }
  cu2<-unique(clusters[clusters!="X"])
  
  #DE test to get top overexpressed genes
  deg<-lapply(cu2,function(x){
    
    a<-names(clusters[clusters==x])
    v<-ifelse(colnames(m1) %in% a,"Y","N")
    
    d<-DE.genes(matrix=m1,vector=v,no=geneno)
    d<-d[d$Cluster=="Y",]
    d$Cluster<-x
    d})
  
  d3<-do.call(rbind,deg)
  rownames(d3)<-1:nrow(d3)
  
  if(isTRUE(plot)){
  p1<-invisible(meta.plot(m1,d3,clusters,name,extra = extraplot))
  }
  
  out<-list(clusters=clusters,clustermat=d3,plot=p1)
  
}
#######Find doublets by a bunch of algorithms
find.doublets<-function(matrix,method="scdblfinder",scransubset=rownames(matrix),sctrans=F,scpc=10,doubleprop=NULL,
                        scale=0.004){
  
  
  if(is.null(doubleprop)){doubleprop<-scale*ncol(matrix)/500}
  
  
  
  if(method != "seurat"){
    sce<-SingleCellExperiment::SingleCellExperiment(assays=list(counts=matrix))
  }
  if(method=="scdblfinder"){
    sce <- scDblFinder::scDblFinder(sce,dbr = doubleprop)
    
    out<-cbind.data.frame(sce$scDblFinder.score,sce$scDblFinder.class)
  }
  if(method=="scds"){
    sce<- scds::cxds_bcds_hybrid(sce)
    sce$bin<-ifelse(sce$hybrid_score>quantile(sce$hybrid_score,probs=1-doubleprop),"doublet","singlet")
    out<-cbind.data.frame(sce$hybrid_score,sce$bin)
  }
  if(method=="scran"){
    
    scd<-scDblFinder::computeDoubletDensity(sce,subset.row=intersect(scransubset,rownames(matrix)))
    scd<-log2(scd+1)
    sce$scran_dblscore <- scd
    sce$bin<-ifelse(sce$scran_dblscore>quantile(sce$scran_dblscore,probs=1-doubleprop),"doublet","singlet")
    out<-cbind.data.frame(sce$scran_dblscore,sce$bin)
  }
  
  if(method=="seurat"){
    library("Seurat")
    library("DoubletFinder")
    library("sctransform")
    seu <- CreateSeuratObject(counts = matrix, project = "x", min.cells = 20, min.features = 1000)
    if(isFALSE(sctrans)){
      seu <- NormalizeData(seu)
      seu <- ScaleData(seu)
      seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
    }
    if(isTRUE(sctrans)){
      seu<-SCTransform(seu)
    }
    seu<-RunPCA(seu)
    seu <- FindNeighbors(seu, reduction = "pca", dims = 1:scpc)
    seu <- FindClusters(seu, resolution = 0.5)
    
    sweeplist <- paramSweep_v3(seu, PCs = 1:scpc, sct = sctrans)
    sweepstats <- summarizeSweep(sweeplist, GT = FALSE)
    bcmvn <- find.pK(sweepstats)
    
    pk<-as.character(bcmvn[which.max(bcmvn$BCmetric),"pK"])
    pk<-as.numeric(pk)
    
    homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)  
    nExp_poi <- round(doubleprop*ncol(seu))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    seu <- doubletFinder_v3(seu, PCs = 1:scpc, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = sctrans)
    
    n1<-ncol(seu@meta.data)
    
    out<-seu@meta.data[,(n1-1):n1]
  }
  
  colnames(out)<-c("DoubletScore","DoubletClass")
  rownames(out)<-colnames(matrix)
  return(out)
}
####Score cells for the percentage of UMIs coming from a gene pattern
pattern.score<-function(matrix,pattern){
  
  tot<- colSums(matrix)
  mit<-matrix[grep(pattern,rownames(matrix)),]
  ms<-colSums(mit)
  frac<-ms/tot
  return(frac)
}

#####Classify cells by using a celltype matrix
type.score<-function(matrix,typemat,gv='Gene',tv='Celltype',size=2,conf_int=0,bin=100,centre=TRUE,
                     conf="scale",no_bins=30){
  
  #Only genes existing in tested matrix
  typemat<-typemat[typemat[,gv] %in% rownames(matrix),]
  #Filter types with few genes
  tl<- typemat %>% group_by(typemat[,tv]) %>% summarise_all(length)
  tl<-tl[tl$Gene>size,]
  tl2<-as.data.frame(tl[,1])
  c1<-as.character(tl2[,1])
  typemat<-typemat[typemat[,tv] %in% c1,]
  
  #Score cells
  t1<-unique(typemat[,tv])
  typescorelist<-parallel::mclapply(t1,mc.cores=length(t1),function(x){
      typescore<-score(matrix,typemat[typemat[,tv]==x,gv],bin=bin,centre=centre,no_bins = no_bins)

  })
  typescores<-do.call(cbind,typescorelist)
  colnames(typescores)<-t1
  names(typescorelist)<-t1
  #Classify cells by type with the highest expression
  Type<-character(length=nrow(typescores))
  names(Type)<-rownames(typescores)
  for(i in 1:nrow(typescores)){
    
    a<- which.max(typescores[i,])
    b<- max(typescores[i,-a])
    c<- max(typescores[i,])
    #
    if(conf_int!=0){
      d<-ifelse(conf=="scale",conf_int*c,conf_int)
    Type[i]<-ifelse( c-b > d, names(a),"Unresolved")}else{
    Type[i]<-names(a)}
  }
  
  out<-list(type=Type,typescorelist=typescorelist)
  return(out)
}
###############Plot gene programs during collection of metaprograms
meta.plot<-function(matrix,clustermat,clusters,name="x",extra=FALSE){
  
  matrix<-as.matrix(matrix)
  d3<-clustermat
  cc<-as.character(unique(d3$Gene))
  
  #For overclustering method
  if(is.list(clusters)){
    c1<-cor(matrix)
    m3<-hclustreorder(c1)
    n1<-colnames(m3)
    m22<-matrix[cc,n1]
    
    #Annotate clusters
    sfun<-lapply(names(clusters),function(x){
      
      a<-clusters[[x]]
      s1<-ifelse(colnames(m22) %in% a,"Y",NA)
      
    })
    scoremat<-do.call(cbind.data.frame,sfun)
    colnames(scoremat)<-names(clusters)
    hl = columnAnnotation(df=scoremat,na_col="white",show_legend=FALSE)
  }
  if(!is.list(clusters) &isTRUE(extra)){
    m22<-matrix[cc,names(sort(clusters))]
    hl=columnAnnotation(Cluster=clusters[colnames(m22)],
                        HPV=HPV[colnames(m22)],
                                Subclone=tcl[colnames(m22)])
  }
  if(!is.list(clusters) &isFALSE(extra)){
    m22<-matrix[cc,names(sort(clusters))]
    hl=columnAnnotation(Cluster=clusters[colnames(m22)])
  }
  
  
  #Annotate genes
  tn<- d3[d3$Gene %in% cc,] %>% group_by(Cluster) %>%  top_n(60/length(unique(clusters)),logFC) 
  rl<-as.character(unique(tn$Gene))
  hr = rowAnnotation(q=anno_mark(at = which(cc %in% rl), labels =rl))
  #Number columns
  nam <- rep(" ", ncol(m22))
  cn<- plyr::round_any( (ncol(m22)/10),25)
  if(cn !=0){
    nam[seq(25,ncol(m22), cn)] <- seq(25, ncol(m22), cn)}
  colnames(m22) <- nam
  #Make heatmap
  col_fun = circlize::colorRamp2(c(min(m22), 0, max(m22)), c("steelblue", "white", "darkred"))
  p1<-Heatmap(m22,col=col_fun,show_row_names=F,show_column_names=T,cluster_rows=F,cluster_columns=F,top_annotation=hl,
                        right_annotation = hr,name="log2(CPM)",column_title = paste0(name," Clusters"))
  return(p1)
}

###########Find metaprograms by taking cutting hclust tree at every height and merging similar clusters
meta.hclust<-function(matrix,cormet="pearson",distmet="average",minsize=10,maxfrac=0.5,
                      jaccut=0.3,jacsize=1000,name="x",extraplot=FALSE,geneno=50){
  
  
  ###Extract all clusters by hclust
  m1<-as.matrix(matrix)
  c1<-cor(m1,method = cormet)
  d<-as.dist((1-c1)/2)
  hc <- hclust(d, method = distmet)
  
  hcut<-cutree(hc,h=hc$height)
  cl1<-hclean(hcut,min=minsize,max=maxfrac)
  ####diffexp
  jacmax=1
  cutoff=jaccut
  size=jacsize
  while(jacmax>cutoff & length(cl1) >1){
    
    
    deg<-lapply(names(cl1),function(x){
      
      a<-cl1[[x]]
      v<-ifelse(colnames(m1) %in% a,"Y","N")
      
      d<-DE.genes(matrix=m1,vector=v,no=geneno)
      d<-d[d$Cluster=="Y",]
      d$Cluster<-paste0(name,"_Cluster_",x)
      d})
    
    d3<-do.call(rbind,deg)
    rownames(d3)<-1:nrow(d3)
    
    jac_mat<-jac.mat(cluster_matrix=d3,list=cl1,vector=d3$Cluster)
    
    b1<-lapply(strsplit(jac_mat$C1,split="_"),function(x)x[3])
    b2<-lapply(strsplit(jac_mat$C2,split="_"),function(x)x[3])
    
    jacmax<-max(jac_mat$Jaccard)
    
    
    for (i in 1:nrow(jac_mat)) {
      
      
      if (jac_mat[i, 3] > cutoff & jac_mat[i,5]<size & jac_mat[i,6]>=0) {
        cl1[[b1[[i]]]]<-union(cl1[[b1[[i]]]],cl1[[b2[[i]]]])
        cl1[[b2[[i]]]]<-NULL
        
      }
      
      if (jac_mat[i, 3] > cutoff & jac_mat[i,4]<size  & jac_mat[i,6]<0) {
        cl1[[b2[[i]]]]<-union(cl1[[b1[[i]]]],cl1[[b2[[i]]]])
        cl1[[b1[[i]]]]<-NULL
        
        
      }
    }
    names(cl1)<-1:length(cl1)
  }
  #
  #Filter clusters with few DE genes
  cu<-unique(d3$Cluster)
  for(i in 1:length(cu)){
    
    a<-nrow(d3[d3$Cluster==cu[i],])
    d3[d3$Cluster==cu[i],"sc"]<-ifelse(a<geneno,"rm","keep")
  }
  #
  d3<-d3[d3$sc=="keep",]
  d3<-d3[,-5]
  cu<-unique(d3$Cluster)
  s2 <- sapply(strsplit(cu, split='_', fixed=TRUE),function(x)(x[3]))
  cl1<-cl1[names(cl1) %in% s2]
  
  # for(i in 1:length(cu)){
  #   d3[d3$Cluster==cu[i],"Cluster"]<-paste(name,"Cluster",LETTERS[i],sep = "_")
  #   names(cl1)[i]<-paste(name,"Cluster",LETTERS[i],sep = "_")
  # }
  
  for(i in 1:length(cu)){
    d3[d3$Cluster==cu[i],"Cluster"]<-paste(name,"Cluster",i,sep = "_")
    names(cl1)[i]<-paste(name,"Cluster",i,sep = "_")
  }
  
  
  p1<-invisible(meta.plot(m1,d3,cl1,name,extra = extraplot))
  
  out<-list(clusters=cl1,clustermat=d3,plot=p1)
  return(out)
}
batch.correct<-function(matrix,pats,pv=Patient,tv=Type[Type!="Epithelial"]){
  
  tn<-unique(tv)
  
  for(i in 1:length(pats)){
    f1<-mclapply(tn,function(x){
      n1<-intersect(names(pv[!pv %in% pats]),names(tv[tv==x]))
      n2<-intersect(names(pv[pv == pats[i]]),names(tv[tv==x]))
      a<-rowMeans(matrix[,n1])
      b<-rowMeans(matrix[,n2])
      c<- b-a
      c
      
    })
    
    for(j in 1:length(tn)){
      n3<-intersect(names(tv[tv==tn[j]]),names(pv[pv == pats[i]]))
      a<- matrix[,n3 ]
      matrix[,n3 ] <- a - f1[[j]]
    }
  }
  return(matrix)
}
#########Read hdf5 files and calculate moving average of values at allelic positions
#######Need chr_names etc loaded
#library("rhdf5")  
wes.cna<-function(wd,pat,mav=5000,noise=0.15,per_chr=FALSE){
  setwd(wd)
  
  l1<-list.files(pattern = "hdf5")
  a<-l1[grep(pat,l1)]  #Files from patient
  tc<- log2(rhdf5::h5read(a[grep("tumor",a)],"counts")$values +1)  #Tumour counts
  nc<- log2(rhdf5::h5read(a[grep("normal",a)],"counts")$values +1) #normal conts
  d1<-tc-nc
  
  tx<-rhdf5::h5read(a[grep("tumor",a)],"intervals")
  txn<-tx$indexed_contig_names
  txm<-as.data.frame(tx$transposed_index_start_end)
  colnames(txm)<-c("chr","start","end")
  for(i in 1:length(unique(txn))){
    txm[txm$chr==unique(txm$chr)[i],"chr"]<-unique(txn)[i]}
  
  rownames(d1)<-rownames(txm)
  txm<-txm[grep("_",txm$chr,invert = TRUE),]
  d1<-d1[rownames(txm),]
  
  #Split by chromosome
  chrn<-as.character(unique(txm$chr))
  ##Chromosome length
  chrl1<- lapply(chrn, function(x){
    
    l<-nrow(subset(txm, chr ==x))
    
  })
  chrl<-do.call(rbind,chrl1)
  #Break vector at each chromosome end
  chrb<-cumsum(chrl)
  if(isTRUE(per_chr)){
    chrb<-c(1,chrb)
    chrn<-c(chrn,"")
    
    num<-seq(1:(length(chrb)-1))
    perchr<-lapply(num,function(y){
      
      
      if(y==length(num)) {end=length(d1)}
      if(y!=length(num)){end<-chrb[y+1]-1}
      chr<-d1[chrb[y]:end]
      #
      chr_mat<-caTools::runmean(chr,k=mav,endrule = "mean")
    })
    #
    m1<-unlist(perchr)
  }
  #Calculate moving average for all gen
  if(isFALSE(per_chr)){ m1<-caTools::runmean(d1,k=mav,endrule = "mean")}
  m2<-m1-median(m1)
  m3<- ifelse(m2>noise,m2-noise,
              ifelse(m2< noise* -1 ,m2- noise* -1,0))
  
  
  #Plot
  plot(m3,xaxt="n",xlab = "",ylab = paste0(mav,"-gene MA"))
  axis(side=1,at=chrb,labels = chrn,las=2)
  abline(v=chrb)
  abline(h=median(m3))
  title(main = paste0(pat," WES CNA"))
  pq<-recordPlot()
  plot.new()
  #
  txm$diff<-d1
  txm$ma<-m3
  
  cg2<-cna_genes[1:chr_breaks[23]]
  b1<-seq(1,length(d1),by=length(d1)/length(cg2))
  
  d2<-numeric(length=length(cg2))
  for(i in 1:(length(cg2)-1)){
    d2[i]<-mean(txm[b1[i]:b1[i+1],"ma"])}
  
  
  
  
  return(list(plot=pq,wesmat=txm,short=d2))
}
###############Compare CNAs between those inferred from RNASeq and WES
comp.cna<-function(wes,matrix=mat,clones=tcl,noise=0.15,name){
  
  cg2<-rownames(matrix)[1:chr_breaks[23]]
  brx<-br.fun(cg2,separate="arm")
  l2<-length(brx$bdf$breaks[brx$bdf$breaks<length(cg2)])
  breaks<-brx$bdf$breaks[1:l2]
  labels<-brx$bdf$labels[1:l2]
  
  ca<-unique(clones)[grep("Normal|Unresolved",unique(clones),invert = TRUE)]
  cm<-matrix(nrow=length(cg2),ncol=length(ca))
  for(i in 1:length(ca)){
    a<-matrix[cg2,names(clones[clones==ca[i]])]
    t1<-apply(a,1,mean)
    t2<-t1-median(t1)
    t3<- ifelse(t2>noise,t2-noise,
                ifelse(t2< noise* -1 ,t2- noise* -1,0))
    cm[,i]<-t3
  }
  colnames(cm)<-ca
  #rownames(cm)<-cg2
  
  z<-cbind(cm,wes)
  
  rm<-brx$rm
  if(length(rm)!=0){
    z<-z[-rm,]}
  
  colnames(z)[ncol(z)]<-"WES"
  zz<-reshape2::melt(z)
  
  p1<-ggplot(data=zz, aes(x=Var1,y=value,colour=Var2)) + geom_point(alpha=0.3)+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
    scale_x_continuous(breaks=breaks,labels=labels, expand = c(0,0))+
    geom_vline(xintercept=breaks)+theme_light()+
    xlab("")+ylab("Inferred Copy Number Log-Ratio")+
    guides(colour = guide_legend(title="Subclone",override.aes = list(size=4,alpha=1)))+
    ggtitle(paste0(name," Inferred CNA per Subclone"))
  
  #
  a1<-as.character("z[,ncol(z)] ~ ")
  a2<-character(length = length(ca))
  for(i in 1:length(a2)){ a2[i]<-paste0("z[,",i,"]")}
  a3<-paste(a2,collapse = " + ")
  
  lq<-lm(as.formula(paste0(a1,a3)))
  names(lq$coefficients)[2:ncol(z)]<-ca
  
  return(list(plot=p1,compmat=z,linreg=lq))
}
#Volcano plot
volcano.plot<-function(de_mat,title=NULL,main="",fc=2,pv=0.05){
  
  if(is.null(title)){
    a<-unique(de_mat$Cluster)
    title<-paste0(a[1]," vs ",a[2])
    
  }
  if(length(de_mat$Gene)==0){de_mat$Gene<-rownames(de_mat)}
  
  ggplot(de_mat,aes(y=-log10(adj.pval) , x=logFC)) +geom_point()+
    theme_light(base_size = 30)+guides(colour="none")+
    geom_text_repel(size=10,aes(label=ifelse(abs(logFC)>fc&adj.pval<pv,
                                     as.character(Gene),''),
                        colour=ifelse(logFC>0,"green","red")))+
    geom_vline(xintercept = c(-fc,fc),lty=2)+
    geom_hline(yintercept = -log10(pv),lty=2)+ggtitle(main)+
    xlab(paste0("log2FC " ,title))+
    ylab("-log10(p-value)")}

##XY plot
xyplot<-function(matrix,diff=NULL,common=NULL,no=20,shift=0){
  matrix[,2]<-matrix[,2]+shift
  matrix$diff<-matrix[,1]-matrix[,2]
  matrix$gene<-rownames(matrix)
  if(is.null(diff)){
    da<-tail(sort(abs(matrix$diff)),n=no)[1]
  }else{da<-diff}
  if(is.null(common)){
    xa<-tail(sort(matrix[,1]),n=no)[1]
    xb<-head(sort(matrix[,1]),n=no)[no]
    ya<-tail(sort(matrix[,2]),n=no)[1]
    yb<-head(sort(matrix[,2]),n=no)[no]
  }else{
    xa<-common
    xb<- -1*common
    ya<-common
    yb<- -1*common
  }
  
  pa<-ggplot(matrix,aes(x=matrix[,1],y=matrix[,2]))+geom_point()+theme_light()+
    geom_text_repel(max.overlaps = 50, aes(label=ifelse(abs(diff)>da| matrix[,1]>xa &matrix[,2] >ya|
                                       matrix[,1]< xb & matrix[,2] < yb, 
                                     as.character(gene),''),
                        colour=ifelse(diff>da,"green",
                                      ifelse(diff< -1*da,"red","black"))))+guides(colour=FALSE)+
    xlab(colnames(matrix)[1])+ylab(colnames(matrix)[2])
  return(pa)
}
####Annotate DE gene matrix with attributes from biomart and uniprot
annotate.genes<-function(matrix,host="useast.ensembl.org",
                         mart=biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",
                                               host=host),filter="hgnc_symbol",
                         biomart_attributes=c("hgnc_symbol", "uniprotswissprot", "wikigene_description"),
                         uniprot_tax=9606,
                         uniprot_attributes=c("SUBCELLULAR-LOCATIONS","FUNCTION"),
                         uniprot=TRUE,
                         select="hgnc_symbol",
                         order="uniprotswissprot"){
f2<-as.character(filter)
f3<-as.character(order)
f4<-as.character(select)

if(class(matrix)!="data.frame"){matrix<-as.data.frame(matrix)}

geneconv <- biomaRt::getBM(attributes=biomart_attributes,
                  filters=filter,values=matrix[,f2], mart=mart, uniqueRows=T)
#Merge downloaded ids with your matrix. Youll need to manually remove some doubles as one uniprot sometimes maps to several genes and you manually have to choose which hgnc symbol to keep
newmatrix = merge(x = matrix, y = geneconv, by=f2,all.x=TRUE)
ordnew<-newmatrix[order(newmatrix[,order], decreasing = TRUE),]
nm<-ordnew[!duplicated(ordnew[,select]),]
#Make blanks to NAs, biomart output is pretty buggy
nm<-as.data.frame(apply(nm, 2, function(x) gsub("^$|^ $", NA, x)))

##Get descriptions from uniprot
if(isTRUE(uniprot)){
  up <- UniProt.ws::UniProt.ws(taxId=uniprot_tax)
  #List of whatever genes you're interested in - needs uniprot id
  lx<-na.omit(nm$uniprotswissprot) #List of uniprot IDs
  
  ## Get descriptions for your gene list - do it in chunks of ca 80, otherwise process times out
  genedesc<-UniProt.ws::select(up, lx, uniprot_attributes, "UNIPROTKB")
  genedesc<-genedesc[!duplicated(genedesc$UNIPROTKB),]
  #Merge descriptions with original file
  names(genedesc)[names(genedesc) == 'UNIPROTKB'] <- 'uniprotswissprot'
  desc2<-merge(x = nm, y = genedesc, by="uniprotswissprot",all.x=TRUE)
  nm<-desc2
}
return(nm)
}

#####Standard UMAP
standard.umap<-function(matrix,metric="correlation",
                 n_neighbors=50,spread=10,min_dist=0.01,var_genes=NULL,...){
  if(!is.null(var_genes)){
    s1<-apply(matrix,1,sd)
    n1<-names(tail(sort(s1),n=var_genes))
    matrix<-matrix[n1,]
      }
  
  
  tm2<-t(as.matrix(matrix))
  um2<-uwot::umap(tm2,metric=metric,n_neighbors=n_neighbors,
                  spread=spread,min_dist=min_dist,
                  pca_center = FALSE, fast_sgd = TRUE,n_threads = 100,...)
  colnames(um2)<-c("V1","V2")
  rownames(um2)<-rownames(tm2)
  dtsne22<-as.data.frame(um2)
  return(dtsne22)
  }

##########Separate groups by permutation testing
permutation.score<-function(matrix,signature,trials=100,correction="fdr",name="x",method="max",
                            mode="positive",cutoff=0.01,cna=NULL,cna_cells=NULL){
  


  g1<-signature
  if(is.null(cna)){s1 <- score(matrix,g1,bin=NULL)}else{
  if(cna=="abs"){s1 <- apply(matrix[g1,],2,function(x)mean(abs(x)))}
  if(cna=="sqrt"){s1 <- apply(matrix[g1,],2,function(x)mean(sqrt(x^2)))}
  if(cna=="cor"){s1<-rowMeans(cor(matrix[g1,],matrix[g1,cna_cells]))}
}
  N <- trials
  
  perm_test <- mclapply(1:N, function(i) {
    permuted <- t(apply(matrix, 1, gtools::permute))
    
    colnames(permuted) <- colnames(matrix); rownames(permuted) <- rownames(matrix)
    
    if(is.null(cna)){scx <- score(permuted,g1,bin=NULL)}else{
    if(cna=="abs"){scx <- apply(permuted[g1,],2,function(x)mean(abs(x)))}
    if(cna=="sqrt"){scx <- apply(permuted[g1,],2,function(x)mean(sqrt(x^2)))}
    if(cna=="cor"){scx<-rowMeans(cor(permuted[g1,],permuted[g1,cna_cells]))}
    }
    
    if(mode=="positive" & method=="max"){
      res = s1 > max(scx)}
    if(mode=="negative" &method=="max"){
      res = s1 < min(scx)
    }
    if(mode=="positive" & method=="sd"){
      res = s1 > (mean(scx) + 2*sd(scx))
    }
    if(mode=="negative" & method=="sd"){
      res = s1 < (mean(scx) - 2*sd(scx))
    }
    res
    
  })
  perm_test<-do.call(cbind,perm_test)
  # Count the number of times each cell was an outlier
  perm_test_res <- apply(perm_test, 1, function(x) length(which(x == TRUE)))
  prob_success <- sum(perm_test_res) / (N * length(perm_test_res))
  pv1 <- sapply(perm_test_res,
                function(x) binom.test(x = x, n = N, p = prob_success, 
                                       alternative = "greater")$p.value)
  pv2<-p.adjust(p = pv1, method = correction)
  out<-ifelse(pv2<cutoff,name,"rest")
  return(out)
}

#######Plot gene expression on reduced dimensions
gene.plot<-function(dimred,genes,matrix){
  
  final.matrix<-matrix  
  dt2<-dimred
  l1<-genes
  
  pf<-lapply(l1,function(x){
    a<-as.character(x)
    
    ggplot(dt2,aes(x=V1, y=V2, colour = final.matrix[a,]))+
      geom_point(size=1,alpha=0.5) + 
      scale_colour_gradientn(name=x,colours=c("lightgrey", "lightsalmon", "darkred"))+
      guides(colour=FALSE)+
      xlab("") + ylab("") +
      theme_light(base_size=20)+ggtitle(a)
    
  })
  p1<-do.call(plot_grid,pf)
  
  p22<-ggplot(dt2,aes(x=V1, y=V2, colour = final.matrix[l1[1],]))+
    geom_point(size=1,alpha=0.5) + 
    scale_colour_gradientn(name="log2(CPM)",colours=c("lightgrey", "lightsalmon", "darkred"))+
    theme_light(base_size=20)
  
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    p22 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  p3<-plot_grid(p1, legend, rel_widths = c(3, .4))
  return(p3)
}


######Dot plot of genes per cell type
dot.plot<-function(genes,logmatrix,umimatrix,types,order){
  

  
  df<-t(logmatrix[genes,])
  df<-as.data.frame(as.matrix(df))
  
  df$Type<-types
  
  d3<- df %>% group_by(Type) %>% summarise_all(mean)
  d4<-as.matrix(d3[1:length(order),2:ncol(df)])
  rownames(d4)<-d3$Type[1:length(order)]
  d4<-as.data.frame(d4)
  d4<-d4[order,]
  d4<-na.omit(d4)
  d4$Type<-rownames(d4)
  d5<-reshape2::melt(d4)
  mx<-umimatrix[genes,colnames(logmatrix)]
  mx[mx>0]<-1
  mx<-t(as.data.frame(as.matrix(mx)))
  mx<-as.data.frame(mx)
  mx$Type<-types
  mf<- mx %>% group_by(Type) %>% summarise_all(function(x)sum(x)/length(x))
  mf2<-as.matrix(mf[1:length(order),2:ncol(df)])
  rownames(mf2)<-mf$Type[1:length(order)]
  mf2<-as.data.frame(mf2)
  mf2<-mf2[order,]
  mf2<-na.omit(mf2)
  mf2$Type<-rownames(mf2)
  mm<-reshape2::melt(mf2)
  d5$frac<-mm$value
  d5$value<-ifelse(d5$value<=0,0,d5$value)
  
  ggplot(as.data.frame(d5[d5$frac>0.5,]),aes(x=factor(Type,levels=order),
                  y=factor(variable,levels=genes),colour=frac,size=value))+geom_point()+
    theme_light()+xlab("")+ylab("")+
    scale_size_continuous(range=c(1,10),name = "Mean\nlog2(CPM)")+
    scale_colour_gradientn(colours=c("lightpink", "salmon", "darkred"),name="Fraction Positive\nCells")+ 
    theme(text = element_text(size = 20))+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
}
#Correlation scatter plot
corr.scatter<-function(df,title=NULL,label=''){
  df$l2<-label
  a1<-round(cor(df[,1],df[,2]),digits = 2)
  a2<-round(cor.test(df[,1],df[,2])$p.value,digits=2)
  if(a2==0){pv2<- " < 2.2e-16"}else{pv2<-paste0(" = ",a2)}
  
  ggplot(df,aes(x=df[,1],y=df[,2]))+geom_point(size=8)+geom_smooth(method="lm")+
    theme_light(base_size = 30)+ggtitle(paste0(title,"  r=",a1," p",pv2))+
    xlab(colnames(df)[1])+ylab(colnames(df)[2])+
    geom_text_repel(size=20,aes(label=l2))
}

############Cluster for metaprogs by NNMF
meta.nmf<-function(matrix,k=10,size=10,geneno=100,name="x",plot=TRUE,extraplot=FALSE,
                   jaccut=NULL,...){
  
  
  matrix<-as.matrix(matrix)
  #Set negatives to 0
  ma<-ifelse(matrix<=0,0,matrix)
  #Remove zero rows (side-effect of batch correction)
  dd<-apply(ma,1,function(x)sum(x==0))
  dn<- names(dd[dd==ncol(matrix)])
  ma<-ma[!rownames(ma) %in% dn,]
  #Run NMF
  n1<-NMF::nmf(ma,k,...)
  
  #Extract top genes per factor
  if(!is.null(geneno)){
    fs <- NMF::extractFeatures(n1, method = geneno)
    names(fs)<-seq(1,k)}else{
      fs<-NMF::extractFeatures(n1)
      names(fs)<-seq(1,k)
      f2<-lapply(fs,length)
      fs<-fs[f2>=size]
    }
  
  #Assign each cell to a factor
  cs<- apply(NMF::coefficients(n1), 2, which.max)
  cs<-cs[cs %in% names(fs)]
  
  #Merge by jaccard similarity
  if(!is.null(jaccut)){
    
    s1<-Sys.time()
    jacmax=1
    
    while(jacmax>jaccut){
      
      print("Starting merging round")
      cv<-cs
      c1<-unique(cv)
      genes<-unique(unlist(lapply(fs, function(f) rownames(NMF::basis(n1))[f]), recursive = FALSE))
      
      #Assign genes to cells
      jf<-lapply(names(fs),function(x){
        Gene<-rownames(NMF::basis(n1))[fs[[x]]]
        Cluster<-x
        cbind.data.frame(Gene,Cluster)
      })
      jf<-do.call(rbind.data.frame,jf)
      
      jac_mat<-jac.mat(cluster_matrix=jf,list=NULL,matrix=matrix[genes,],vector=cs)
      
      print("Jaccard matrix formed")
      
      jacmax<-max(jac_mat$Jaccard)
      
      for (i in 1:nrow(jac_mat)) {
        
        
        if (jac_mat[i, 3] > jaccut & jac_mat[i,6]>=0) {
          cv[cv==jac_mat[i,2]] <-jac_mat[i,1]
          
        }
        
        if (jac_mat[i, 3] > jaccut  & jac_mat[i,6]<0) {
          cv[cv==jac_mat[i,1]] <-jac_mat[i,2]
          
        }
      }
      s2<-Sys.time()
      if(jacmax>jaccut){
        print( paste(length(c1),"clusters merged to",length(unique(cv)),sep=" "))
        print(s2-s1)
      }
      print(paste0("Max similarity ",jacmax))
      cs<-cv
      cx<-unique(cs)
      fs<-fs[cx]
    }
  }
  
  #Remove factors with few cells
  cx<-table(cs)>=size
  cx<-names(cx[cx==TRUE])
  cs<-cs[cs %in% cx]
  fs<-fs[cx]
  
  #Assign genes to cells
  jf<-lapply(names(fs),function(x){
    Gene<-rownames(NMF::basis(n1))[fs[[x]]]
    Cluster<-x
    cbind.data.frame(Gene,Cluster)
  })
  jf<-do.call(rbind.data.frame,jf)
  
  #Rename clusters
  cu2<-sort(unique(cs))
  for(i in 1:length(cu2)){
    cs[cs==cu2[i]]<-paste0(name,"_Cluster_",LETTERS[i])
    jf[jf$Cluster==cu2[i],"Cluster"]<-paste0(name,"_Cluster_",LETTERS[i])
  }
  
  #Set logFC
  jf$logFC<-unlist(lapply(unique(jf$Cluster),function(x){
    m2<-matrix[jf[jf$Cluster==x,"Gene"],]
    g1 <- names(cs[cs==x])
    g2 <- names(cs[cs!=x])
    l1<- rowMeans(m2[, g1]) - rowMeans(m2[, g2])
    
  }))
  
  
  #Plot
  if(isTRUE(plot)){
    p1<-invisible(meta.plot(matrix,jf,cs,name,extra = extraplot))
  }else{p1="no_plot"}
  
  out<-list(clusters=cs,clustermat=jf,plot=p1)
}

########Plot distribution of variables across groups
plot.distribution<-function(matrix,column,group,sample,bins=20,ttest=NULL,name=NULL,split="bins",reps=1,
                            sample_no=200){
  
  f1<-parallel::mclapply(1:reps,function(x){
    
    matrix$group<-matrix[,group]
    matrix$column<-matrix[,column]
    matrix$sample<-matrix[,sample]
    
    d1<-matrix %>% group_by(sample,group) %>% dplyr::count() %>% as.data.frame()
    d2<-unlist(lapply(seq(1:nrow(d1)),function(y){
      sample(rownames(matrix[matrix$sample==as.character(d1[y,1]) & 
                               matrix$group==as.character(d1[y,2]),]),min(sample_no,d1[y,3]),replace = F)
    })
    )
    
    matrix<-matrix[d2,]
    
    
    
    
    if(split=="bins"){
      matrix$bin<-ntile(matrix$column,bins)
      sx<- matrix %>% group_by(group) %>% dplyr::count(bin) %>% mutate(exp=sum(n)/length(unique(bin))) %>%
        mutate(frac=log2(n/exp)) %>% as.data.frame()
      
    }
    if(split=="value"){

      matrix$bin<- cut(matrix$column,bins)
      
      ub<-as.character(unique(sort(matrix$bin)))
      
      matrix$bin<- as.character(cut(matrix$column,bins))
      
      for(i in 1:bins){
        matrix[matrix$bin==ub[i],"bin"]<-i
      }
      
      sx<- matrix %>% group_by(group) %>%  dplyr::count(bin) %>% mutate(e1=sum(n)) %>% 
        ungroup %>% group_by(bin)  %>% mutate(exp=sum(n)*e1/sum(unique(e1))) %>% 
        mutate(frac=log2(n/exp)) %>% as.data.frame()
      
    }
    
    
    
    
    #Get G1S values per bin
    mm<- matrix %>% group_by(bin) %>%  summarise(mean(column)) %>% as.data.frame()
    for(i in 1:nrow(mm)){
      sx[sx$bin==mm[i,1],column]<-mm[i,2]
    }
    
    if(!is.null(ttest) & ttest!="chisq"){
      pv<-round(t.test(matrix[matrix$group==ttest,column],matrix[matrix$group!=ttest,column])$p.value,digits=2)

      
      
    }else{pv=1}
    if(ttest=="chisq"){
      pv<-round(chisq.test(matrix$group,matrix$bin)$p.value,digits=2)
      
      
    }
    
    sx$rep=x
    list(sx,pv)
  })
  f2<-do.call(rbind,f1)
  ff<-do.call(rbind,f2[,1])
  sx<- ff %>% group_by(group,bin) %>% mutate(sd=sd(frac)) %>% summarise_all(mean) %>% as.data.frame()
  pv<-mean(unlist(f2[,2]))
  pv2<-character(length=1)
  if(pv==0){pv2<- " < 2.2e-16"}else{pv2<-paste0(" = ",pv)}
  
  
  ggplot(sx, aes(colour=group,y=frac, x=sx[,ifelse(split=="value",column,"bin")],group=group)) +
    ylab("log2(Obs/Exp cells/bin)")+xlab(ifelse(split=="value",column,paste0(column," bins")))+ 
    geom_line()+theme_light()+geom_point()+
    ggtitle(ifelse(is.null(ttest),"",paste0(name," p",pv2))) +
    geom_errorbar(aes(ymin=frac-sd,ymax=frac+sd),width=0.1)+
    scale_x_continuous(breaks=seq(1,bins, 1))
  
  
}
#########Create heatmap from jaccard data frame
jac.heatmap<-function(jacmat){
zq<-jacmat  
zr<-unique(union(zq$C1,zq$C2))
zm<-matrix(ncol=length(zr),nrow=length(zr))
colnames(zm)<-zr
rownames(zm)<-zr
for(i in 1:length(zr)){
  for(j in 1:length(zr)){
    
    a<-zq[zq$C1==zr[i] & zq$C2==zr[j],]
    b<-zq[zq$C2==zr[i] & zq$C1==zr[j],]
    if(length(a$Jaccard)==1){zm[i,j]<-a[1,3]
    zm[j,i]<-a[1,3]}else{zm[i,j]<-b[1,3]
    zm[j,i]<-b[1,3]}
    
  }
}
return(zm)
}

##Multiple DE gene comparisons - select only genes significant in all comparisons
uniqueDEG<-function(umimat,vector,no=50,fc=1,...){
  
  vn<-unique(vector)
  v1<-mclapply(vn,function(x){
    vn2<-vn[vn!=x]
    v2<-mclapply(vn2,function(y){
      v3<-vector[vector==x | vector ==y]
      a<-filter.umi(umimat,names(v3))
      b<-DE.genes(a,v3[colnames(a)],no=no,...)
      b$comp<-y
      b<-b[b$Cluster==x,]
      b
    })
    v4<-do.call(rbind,v2)
    v4<-v4[v4$logFC>fc,]
    v5<-extract.metaprogs(v4,v4$comp,x,no=no)$topgenes
    v5<-v5[v5$Number==length(vn2),]
    v5[,c("Gene","Cluster","logFC")]
  })
  v6<-do.call(rbind,v1)
  vd<-v6$Gene[duplicated(v6$Gene)]
  v6<-v6[!v6$Gene %in% vd,]
  v6<-arrange(v6,Cluster,desc(logFC))
  return(v6)
}

#Complexity of large matrices
complexity<-function(matrix,cut=6){
  
  z1<-ntile(colnames(matrix),cut)
  zn<-unique(z1)
  l1<-unlist(parallel::mclapply(zn,function(y){
    apply(matrix[,z1==y],2,function(x)sum(x>0))
  }))
  return(l1)
}
#Enrichment wrapper
enrich<-function(genes,universe,category="C8",species="Homo sapiens"){
  
  if(length(category)>1){
    cn<-unique(category)
    cn2<-mclapply(cn,function(x){
      msigdbr::msigdbr(species = species, category = x) %>% 
        dplyr::select(gs_name, gene_symbol)
      })
    mdb<-do.call(rbind,cn2)
    
  }else{
  mdb <- msigdbr::msigdbr(species = species, category = category) %>% 
    dplyr::select(gs_name, gene_symbol)
  }
  r<-clusterProfiler::enricher(genes,TERM2GENE=mdb,universe=universe)
  return(list(result=r@result,plot=barplot(r)))
}
#