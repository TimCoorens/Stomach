options(stringsAsFactors = F)
library(readxl)
library(lmerTest)
cnvs=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=4))
svs=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=5))
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=2))
met_patient=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=1))

rownames(wgs_data)=wgs_data$Sample
rownames(met_patient)=met_patient$Donor
nsamples=table(wgs_data$Donor[wgs_data$Feature!="Tumour"])

svs_unique=svs[!svs$in_CNV,]
sv_cnv_type=c(cnvs$CNV, svs$Type)
sv_cnv_type[grepl('whole_gain',sv_cnv_type)]="Trisomy"
sv_cnv_type[grepl('LOH',sv_cnv_type)]="Arm cnn-LOH"
sv_cnv_list_full=data.frame(Type=sv_cnv_type,
                       Sample=c(cnvs$sample,svs$Sample))
sv_cnv_list=unique(sv_cnv_list_full)
sv_cnv_list$Donor=substr(sv_cnv_list$Sample,1,7)
sv_cnv_list$Age=met_patient[sv_cnv_list$Donor,"Age"]

SV_CNV_freq=table(sv_cnv_list$Type,substr(sv_cnv_list$Sample,1,7))
SV_CNV_burden=table(sv_cnv_list_full$Type,sv_cnv_list_full$Sample)


SV_CNV_freq_full=matrix(0,nrow=nrow(SV_CNV_freq),ncol=length(nsamples))
rownames(SV_CNV_freq_full)=rownames(SV_CNV_freq)
colnames(SV_CNV_freq_full)=names(nsamples)
SV_CNV_freq_full[,colnames(SV_CNV_freq)]=SV_CNV_freq

select=wgs_data$Feature=="Normal gland"
SV_CNV_burden_all=data.frame(Sample=wgs_data$Sample[select],
                             Small_burden=0,
                             Del_burden=0,
                             Dup_burden=0,
                             Inv_burden=0,
                             Trisomy_burden=0,
                             LOH_burden=0,
                             Total_burden=0,
                             Donor=wgs_data$Donor[select],
                             Age=wgs_data$Age[select],
                             IM=wgs_data$Intestinal_Metaplasia_in_Sample[select],
                             Cohort=wgs_data$Cohort[select],
                             CI=wgs_data$Chronic_inflammation[select],
                             has_CI=wgs_data$Chronic_inflammation[select]%in%c("Severe"))

rownames(SV_CNV_burden_all)=SV_CNV_burden_all$Sample
SV_CNV_burden_all[colnames(SV_CNV_burden),"Total_burden"]=colSums(SV_CNV_burden)
SV_CNV_burden_all[colnames(SV_CNV_burden),"Small_burden"]=colSums(SV_CNV_burden[c("Deletion","Duplication","Inversion"),])
SV_CNV_burden_all[colnames(SV_CNV_burden),"Del_burden"]=SV_CNV_burden[c("Deletion"),]
SV_CNV_burden_all[colnames(SV_CNV_burden),"Dup_burden"]=SV_CNV_burden[c("Duplication"),]
SV_CNV_burden_all[colnames(SV_CNV_burden),"Inv_burden"]=SV_CNV_burden[c("Inversion"),]
SV_CNV_burden_all[colnames(SV_CNV_burden),"LOH_burden"]=SV_CNV_burden[c("Arm cnn-LOH"),]
SV_CNV_burden_all[colnames(SV_CNV_burden),"Trisomy_burden"]=SV_CNV_burden[c("Trisomy"),]

SV_CNV_burden_all$has_trisomy=SV_CNV_burden_all$Trisomy_burden>0
test1=lmer(Total_burden~Age+(1|Donor),data=SV_CNV_burden_all) #p=0.12

test1=lm(Trisomy_burden~Age,data=SV_CNV_burden_all) #p=0.12

test1=lmer(Trisomy_burden~Age+(1|Donor),data=SV_CNV_burden_all) #p=0.12
test2=lmer(Trisomy_burden~IM+(1|Donor),data=SV_CNV_burden_all) 
test2=lmer(Trisomy_burden~CI+(1|Donor),data=SV_CNV_burden_all) #p=0.12

test2=lmer(Trisomy_burden~has_CI+(1|Donor),data=SV_CNV_burden_all) #p=0.12
test2=lmer(Trisomy_burden~Cohort+(1|Donor),data=SV_CNV_burden_all) #p=0.12

fisher.test(table(SV_CNV_burden_all$CI,SV_CNV_burden_all$has_trisomy))
fisher.test(table(SV_CNV_burden_all$CI,SV_CNV_burden_all$has_trisomy))

anova(test1,test2)

test1=lmer(Small_burden~Age+(1|Donor),data=SV_CNV_burden_all) #p=0.12
test2=lmer(Small_burden~Age+IM+(1|Donor),data=SV_CNV_burden_all) #p=0.12

test2=lmer(Small_burden~Age+has_trisomy+(1|Donor),data=SV_CNV_burden_all) #p=0.12


anova(test1,test2)
library(vioplot)

vioplot(SV_CNV_burden_older$Small_burden~SV_CNV_burden_older$IM)
vioplot(SV_CNV_burden_older$LOH_burden~SV_CNV_burden_older$IM)
vioplot(SV_CNV_burden_older$Trisomy_burden~SV_CNV_burden_older$IM)
vioplot(SV_CNV_burden_older$Del_burden~SV_CNV_burden_older$IM)
vioplot(SV_CNV_burden_older$Dup_burden~SV_CNV_burden_older$IM)

fisher.test(rbind(c(sum(SV_CNV_burden_older$Del_burden[SV_CNV_burden_older$IM=="Absent"]),sum(SV_CNV_burden_older$Del_burden[SV_CNV_burden_older$IM=="Present"])),
                  table(SV_CNV_burden_older$IM)))
fisher.test(rbind(c(sum(SV_CNV_burden_older$Dup_burden[SV_CNV_burden_older$IM=="Absent"]),sum(SV_CNV_burden_older$Dup_burden[SV_CNV_burden_older$IM=="Present"])),
                  table(SV_CNV_burden_older$IM)))

fisher.test(table(SV_CNV_burden_all$Trisomy_burden>0,SV_CNV_burden_all$IM))

sample_cols=rep("grey60",nrow(SV_CNV_burden_all))
sample_cols[SV_CNV_burden_all$IM=="Present"]="firebrick"


plot(SV_CNV_burden_all$Small_burden,x=SV_CNV_burden_all$Age,pch=21,bg=sample_cols)

test1=lm(colSums(SV_CNV_freq_full)/nsamples~met_patient[colnames(SV_CNV_freq_full),"Age"]) #p=0.12

SV_CNV_burden_older=SV_CNV_burden_all[SV_CNV_burden_all$Age>49,]


fisher.test(rbind(c(sum(SV_CNV_burden_all$Small_burden[SV_CNV_burden_all$IM=="Absent"]),sum(SV_CNV_burden_all$Small_burden[SV_CNV_burden_all$IM=="Present"])),
      table(SV_CNV_burden_all$IM)))
fisher.test(rbind(c(sum(SV_CNV_burden_older$Small_burden[SV_CNV_burden_older$IM=="Absent"]),sum(SV_CNV_burden_older$Small_burden[SV_CNV_burden_older$IM=="Present"])),
                  table(SV_CNV_burden_older$IM)))

test2=lm(colSums(SV_CNV_freq_full[c("Deletion","Duplication"),])/nsamples~met_patient[colnames(SV_CNV_freq_full),"Age"]) #p-value: 0.00951
test3=lm(SV_CNV_freq_full["Trisomy",]/nsamples~met_patient[colnames(SV_CNV_freq_full),"Age"]) #p=0.4


cnvs$CNV_id=paste0(cnvs$CNV,"_",cnvs$Fragile_site,"_",cnvs$Donor,"_",cnvs$Event)
CNV_uniq=unique(cnvs$CNV_id)
svs_unique$id=paste0(svs_unique$Type,"_",svs_unique$Fragile_site,"_",svs_unique$patient,"_",svs_unique$Event)
SV_uniq=unique(svs_unique$id)
CNV_SV_uniq=c(CNV_uniq,SV_uniq)

fragile_sites=unique(c(svs_unique$Fragile_site,cnvs$Fragile_site))

c(sum(cnv_type%in%c("Deletion","Duplication")),)

CNV_SV_uniq_types=gsub("_PD[4-5].*","",CNV_SV_uniq)
CNV_SV_uniq_types=CNV_SV_uniq_types[grepl("Deletion|Duplication",CNV_SV_uniq_types)]
freqs=table(CNV_SV_uniq_types)
freq_order_vec=c("Deletion",
                 names(sort(freqs[grepl("Del",names(freqs))],decreasing=T)),
                 "Duplication",
                 names(sort(freqs[grepl("Dup",names(freqs))],decreasing=T)),
                 "Inversion")
short=rep(0,length(freq_order_vec))
for(n in 1:length(freq_order_vec)){
  short[n]=sum(grepl(freq_order_vec[n],CNV_SV_uniq))
}

LOH=c(sum(grepl("cnnLOH",CNV_SV_uniq)),
      sum(grepl("chr11q_cnnLOH",CNV_SV_uniq)),
      sum(grepl("chr15q_cnnLOH",CNV_SV_uniq)),
      sum(grepl("chr17q_cnnLOH",CNV_SV_uniq)))
tri=c(sum(grepl("whole_gain",CNV_SV_uniq)),
      sum(grepl("chr20_whole_gain",CNV_SV_uniq)),
      sum(grepl("chr13_whole_gain",CNV_SV_uniq)),
      sum(grepl("chr7_whole_gain",CNV_SV_uniq)),
      sum(grepl("chr18_whole_gain",CNV_SV_uniq)),
      sum(grepl("chr10_whole_gain",CNV_SV_uniq)),
      sum(grepl("chr16_whole_gain",CNV_SV_uniq)))

#Fig. 3a
pdf("~/Desktop/Gastric/SV_CNV_barplot.pdf",width=10,height=5)
barplot(c(short,0,0,0,
          LOH,0,0,0,
          tri),
        col=c('dodgerblue3',rep('lightskyblue',9),'firebrick',rep('salmon',5),"orange",rep('white',3),
              'black',rep("grey70",3),rep('white',3),
              'darkgreen',rep("olivedrab1",6)))
dev.off()

pdf("~/Desktop/Gastric/SV_CNV_intra_barplot.pdf",width=10,height=5)
barplot(short,col=c('dodgerblue3',rep('lightskyblue',9),'firebrick',rep('salmon',5),"orange"))
dev.off()

pdf("~/Desktop/Gastric/SV_CNV_loh_barplot.pdf",width=3,height=5)
barplot(LOH,col=c('black',rep("grey70",3)),ylim=c(0,20))
dev.off()

pdf("~/Desktop/Gastric/SV_CNV_trisomy_barplot.pdf",width=8,height=5)
barplot(tri,col=c('darkgreen',rep("olivedrab1",6)))
dev.off()


patients=met_patient$Donor[order(met_patient$Sex,met_patient$Age,decreasing = F)]

#Fig. 3b
SV_CNV_freq_full=SV_CNV_freq_full[c('Trisomy','Arm cnn-LOH','Deletion',"Duplication","Inversion"),patients[patients%in%colnames(SV_CNV_freq_full)]]
grid = expand.grid(x=rownames(SV_CNV_freq_full), y=colnames(SV_CNV_freq_full))
grid$z = c(SV_CNV_freq_full)
#dev.new(width=6.5, height=6)
library(lattice)
pdf("~/Desktop/Gastric/SV_CNV_heatmap.pdf",width=10,height=3.5)
color.palette = colorRampPalette(c("white", "darkorange", "darkorchid2","darkorchid4"))
levelplot(z~y*x, grid, col.regions=color.palette, scales = list(tck = c(0,0), y = list(cex=1), x = list(rot=90)), ylab="CNV/SV types", xlab="Samples", colorkey=list(space="bottom"), 
          panel=function(...) { arg <- list(...)
          panel.levelplot(...)
          panel.text(arg$x, arg$y, arg$z)})
dev.off()


meta = rbind(nsamples,met_patient[names(nsamples),"Age"])

grid = expand.grid(x=rownames(meta), y=colnames(meta))
grid$z = c(meta)
#dev.new(width=6.5, height=2.5)
library(lattice)
color.palette = colorRampPalette(c("white","grey50"))
levelplot(z~y*x, grid, col.regions=color.palette, scales = list(tck = c(0,0), y = list(cex=1), x = list(rot=90)), ylab="", xlab="Samples", colorkey=list(space="bottom"), 
          panel=function(...) { arg <- list(...)
          panel.levelplot(...)
          panel.text(arg$x, arg$y, arg$z)})

#Enrichment of trisomies
library(ape)
tree=read.tree("~/Desktop/Gastric/PD40293.tree")
SV_CNV_burden_all$plot_col="grey60"
SV_CNV_burden_all$plot_col[SV_CNV_burden_all$IM=="Present"]="firebrick"
pdf("~/Desktop/Gastric/PD40293_SVs.pdf",width=5,height=2.5)
plot(SV_CNV_burden_all[tree$tip.label,"Small_burden"],pch=21,bg=SV_CNV_burden_all[tree$tip.label,"plot_col"],ylab="SVs",xlab="",cex=1.5,lwd=0.5)
dev.off()

tree=read.tree("~/Desktop/Gastric/PD41767.tree")
pdf("~/Desktop/Gastric/PD41767_SVs.pdf",width=5,height=2.5)
plot(SV_CNV_burden_all[tree$tip.label[tree$tip.label%in%rownames(SV_CNV_burden_all)],"Small_burden"],pch=21,bg="grey60",ylab="SVs",xlab="",cex=1.5,lwd=0.5)
dev.off()



wgs_data_sub=wgs_data[wgs_data$Feature!="Tumour"&wgs_data$Site!="Unknown",]
has_trisomy=wgs_data_sub$Sample%in%cnvs$sample[grepl("whole_gain",cnvs$CNV)]

fisher.test(table(has_trisomy,wgs_data_sub$Site)) #p-value = 0.2525

wgs_data_sub=wgs_data[wgs_data$Feature!="Tumour",]
has_trisomy=wgs_data_sub$Sample%in%cnvs$sample[grepl("whole_gain",cnvs$CNV)]
fisher.test(table(has_trisomy,wgs_data_sub$Cohort)) #p-value = 0.4844
fisher.test(table(has_trisomy,wgs_data_sub$Intestinal_Metaplasia_in_Sample)) #p-value = 0.41

drivers=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=5))
has_driver=wgs_data_sub$Sample%in%drivers$Sample
#-------

#TIMING
options(stringsAsFactors = FALSE)
library(data.table)
library(ape)
library(ggtree)
library(readxl)

#Functions

find_children = function(node, tree=tree_df){
  child_nodes = tree$node[tree$parent==node]
  tips = c()
  for (k in 1:length(child_nodes)){
    if (tree$isTip[tree$node==child_nodes[k]]){
      tips=c(tips,child_nodes[k])
    }
    else{
      tips=c(tips,find_children(child_nodes[k]))
    }
  }
  return(tips)
}

find_parent=function(sample,tree_df){
  child=which(tree_df$label==sample)
  parent=tree_df$parent[child]
  all_nodes=c(child,parent)
  while(child!=parent){
    child=parent
    parent=tree_df$parent[child]
    all_nodes=c(all_nodes,parent)
  }
  return(unique(all_nodes))
}

dbinomtrunc = function(x, size, prob, minx=3) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode,minx){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i],minx=minx)
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i],minx=minx)
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode,minx){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode,minx=minx)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode,minx=minx)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6, mode="Full", minx=3){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode,minx=minx)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}

cnvs=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=4))
meta=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=1))
rownames(meta)=meta$Donor
gains=grepl("cnnLOH|whole",cnvs$CNV)
cnvs$CNV_event=paste0(cnvs$Donor,"_",cnvs$CNV,"_",cnvs$Event)
cnvs$num_muts_dup=cnvs$num_muts_nondup=cnvs$Timing_confint_high=cnvs$Timing_confint_low=cnvs$Timing=cnvs$Timing_branch=NA
donors_w_gains=unique(cnvs$Donor[gains])

for(patient in donors_w_gains){
  Muts_to_branches=read.table(paste0("~/Desktop/Gastric/Sequoia/",patient,"_snv_assigned_to_branches.txt"),sep="\t",header=T)
  rownames(Muts_to_branches)=paste(Muts_to_branches$Chr,Muts_to_branches$Pos,Muts_to_branches$Ref,Muts_to_branches$Alt,sep="_")
  
  NV=read.table(paste0("~/Desktop/Gastric/Sequoia/",patient,"_snv_NV_filtered_all.txt"))
  NR=read.table(paste0("~/Desktop/Gastric/Sequoia/",patient,"_snv_NR_filtered_all.txt"))
  tree=read.tree(paste0("~/Desktop/Gastric/Sequoia/",patient,"_snv_tree_with_branch_length.tree"))
  tree_df=as.data.frame(fortify(tree))
  for(gain in unique(cnvs$CNV_event[gains&cnvs$Donor==patient])){
    i=cnvs$CNV_event==gain
    tips=cnvs$sample[i]
    chrom=unique(cnvs$chr[i])
    start=unique(cnvs$start[i])
    end=unique(cnvs$end[i])
    if(length(tips)==1){
      branch=paste0(patient,"_",tree_df$node[tree_df$label==tips&!is.na(tree_df$label)])
    }else{
      branch=paste0(patient,"_",getMRCA(tree,which(tree_df$label%in%tips)))
    }
    if(length(branch)>1) print(patient)
    cnvs$Timing_branch[i]=branch
    
    Muts_present=rownames(Muts_to_branches)[Muts_to_branches$SampleID==branch&
                                              Muts_to_branches$Chr==chrom&
                                              Muts_to_branches$Pos>start&
                                              Muts_to_branches$Pos<end]
    Muts_outside_cnv=rownames(Muts_to_branches)[Muts_to_branches$SampleID==branch&Muts_to_branches$Chr=="chr1"]
    if(length(tips)==1){
      res = binom_mix(x=NV[Muts_present,tips],size=NR[Muts_present,tips],nrange=1:3,mode='Truncated')
      purity=2*mean(NV[Muts_outside_cnv,tips]/NR[Muts_outside_cnv,tips])
    }else{
      res = binom_mix(x=rowSums(NV[Muts_present,tips]),size=rowSums(NR[Muts_present,tips]),nrange=1:3,mode='Truncated')
      purity=2*mean(rowSums(NV[Muts_outside_cnv,tips])/rowSums(NR[Muts_outside_cnv,tips]))
    }
    if(grepl("whole",gain)){
      totCN=3
      dupCN=2
    }
    if(grepl("cnnLOH",gain)){
      totCN=2
      dupCN=2
    }
      cutoff=0.1
      Major_cluster=which(res$p>(dupCN/totCN*purity-cutoff))
      Minor_cluster=which(abs(res$p-1/totCN*purity)<cutoff)
      
      if(length(Major_cluster)){
        if(length(Major_cluster)>1) Major_cluster=which.max(res$p[Major_cluster])
        cnvs$Prop_Duplicated[i]=res$prop[Major_cluster]
      }else{
        cnvs$Prop_Duplicated[i]=0
      }
      if(length(Minor_cluster)){
        cnvs$Prop_NonDuplicated[i]=sum(res$prop[Minor_cluster])
      }else{
        cnvs$Prop_NonDuplicated[i]=0
      }
      cnvs$Timing[i]=totCN/(dupCN+unique(cnvs$Prop_NonDuplicated[i])/unique(cnvs$Prop_Duplicated[i]))
        conf.intv=poisson.test(c(round(length(Muts_present)*unique(cnvs$Prop_NonDuplicated[i])),
                                 round(length(Muts_present)*unique(cnvs$Prop_Duplicated[i]))))$conf.int
        cnvs$Timing_confint_high[i]=totCN/(dupCN+conf.intv[1])
        cnvs$Timing_confint_low[i]=totCN/(dupCN+conf.intv[2])
        cnvs$num_muts_dup[i]=round(length(Muts_present)*unique(cnvs$Prop_Duplicated[i]))
        cnvs$num_muts_nondup[i]=round(length(Muts_present)*unique(cnvs$Prop_NonDuplicated[i]))
  }
}

cnvs$Timing_confint_high=pmin(cnvs$Timing_confint_high,1)


patient="PD40293"
for(patient in unique(cnvs$Donor[!is.na(cnvs$Timing)])){
tree=read.tree(paste0("~/Desktop/Gastric/Sequoia/",patient,"_snv_tree_with_branch_length.tree"))
tree_df=as.data.frame(fortify(tree))

pdf(paste0("~/Desktop/Gastric/Sequoia/",patient,"_CNV_tree.pdf"),useDingbats = F)
plot(tree,label.offset=0.01*max(tree_df$x))
cnvs_sub=cnvs[cnvs$Donor==patient&!is.na(cnvs$Timing),]
for (k in 1:nrow(cnvs_sub)){
  n=as.numeric(strsplit(cnvs_sub$Timing_branch[k],split="_")[[1]][2])
  x_end=tree_df$x[n]
  x_start=tree_df$x[tree_df$parent[n]]
  x_intv=x_end-x_start
  y=node.height(tree)[n]
  tipnum=sum(tree_df$isTip)
  y_adjust=0.2
  col='black'
  if(cnvs_sub$CNV[k]=="chr20_whole_gain") y_adjust=-0.2
  if(cnvs_sub$CNV[k]=="chr20_whole_gain") col='firebrick'
  if(cnvs_sub$CNV[k]=="chr13_whole_gain") col='steelblue'
  if(cnvs_sub$CNV[k]=="chr18_whole_gain") col='forestgreen'
  if(grepl("cnnLOH",cnvs_sub$CNV[k])) col='grey60'
  
  points(x=x_start+min(1,cnvs_sub$Timing[k])*x_intv,y=y+y_adjust,col=col,pch=16,cex=1.5)
  segments(x0=x_start+cnvs_sub$Timing_confint_low[k]*x_intv,x1=x_start+cnvs_sub$Timing_confint_high[k]*x_intv,lwd=2,y0=y+y_adjust,col=col)
}
  axisPhylo(side = 1,backward=F)
  dev.off()
}
points(x=rep(6000,3),cex=2,y=7:9,col=c("black","steelblue","firebrick"),pch=16)
text(x=rep(6100,3),pos=4,y=7:9,labels=c("chr17q LOH","Gain of chr13","Gain of chr20"),pch=16)
dev.off()

cnvs$Age_gain=NA
cnvs$Age_gain_min=NA
cnvs$Age_gain_max=NA

for (k in 1:nrow(cnvs)){
  if(!is.na(cnvs$Timing[k])){
  patient=cnvs$Donor[k]
  tree=read.tree(paste0("~/Desktop/Gastric/Sequoia/",patient,"_snv_tree_with_branch_length.tree"))
  tree_df=as.data.frame(fortify(tree))
  n=as.numeric(strsplit(cnvs$Timing_branch[k],split="_")[[1]][2])
  x_end=tree_df$x[n]
  x_start=tree_df$x[tree_df$parent[n]]
  x_intv=x_end-x_start
  burden_gain=x_start+min(1,cnvs$Timing[k])*x_intv
  burden_min=x_start+cnvs$Timing_confint_low[k]*x_intv
  burden_max=x_start+cnvs$Timing_confint_high[k]*x_intv
  if(tree_df$isTip[n]){
    tot_burden=x_end
  }else{
    tot_burden=mean(tree_df$x[find_children(n)])
  }
  cnvs$Age_gain[k]=burden_gain/tot_burden*meta$Age[meta$Donor==patient]
  cnvs$Age_gain_min[k]=burden_min/tot_burden*meta$Age[meta$Donor==patient]
  cnvs$Age_gain_max[k]=burden_max/tot_burden*meta$Age[meta$Donor==patient]
  }

}
cnvs_sub=unique(cnvs[!is.na(cnvs$Timing),10:24])
cnvs_sub$Age=meta[cnvs_sub$Donor,"Age"]
cnvs_sub=cnvs_sub[order(cnvs_sub$Age,cnvs_sub$Donor,cnvs_sub$CNV),]
col=rep('black',nrow(cnvs_sub))
col[cnvs_sub$CNV=="chr20_whole_gain"]='firebrick'
col[cnvs_sub$CNV=="chr13_whole_gain"]='forestgreen'
col[cnvs_sub$CNV=="chr18_whole_gain"]='dodgerblue'
col[grepl("cnnLOH",cnvs_sub$CNV)]='grey60'

donor_switch=cnvs_sub$Donor[1:(nrow(cnvs_sub)-1)]!=cnvs_sub$Donor[2:nrow(cnvs_sub)]

pdf("~/Desktop/Gastric/gain_timing.pdf",width=6,height=4)
plot(cnvs_sub$Age_gain,pch=21,bg=col,ylim=c(0,80),ylab="Age at which gain occured",xaxt="n",xlab="")
segments(x0=1:nrow(cnvs_sub),y0=cnvs_sub$Age_gain_min,y1=cnvs_sub$Age_gain_max,col=col)
abline(v=which(donor_switch)+0.5,lty="dashed")
segments(x0=c(0,which(donor_switch)+0.5),x1=c(which(donor_switch)+0.5,nrow(cnvs_sub)+1),y0=meta[unique(cnvs_sub$Donor),"Age"],lwd=2)
axis(1, at=(c(0,which(donor_switch)+0.5)+c(which(donor_switch)+0.5,nrow(cnvs_sub)+1))/2,labels=unique(cnvs_sub$Donor), las=2)
dev.off()


tgs=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1.xlsx",sheet=3))

meta=read.csv("~/Desktop/Gastric/master_metadata_gastric_lcm_31OCT2019.csv")

meta_TGS=meta[meta$TGSseq,c("PD_ID","PD_STEM","COLLABORATOR","COLLABORATOR_DONOR_ID","SLIDE_NAME","WELL","SIZE_UM2")]

unique(meta_TGS$SLIDE_NAME)%in%unique(meta$SLIDE_NAME[meta$WGSseq])

slides=unique(meta_TGS$SLIDE_NAME)
slides[!slides%in%unique(meta$SLIDE_NAME[meta$WGSseq])]
meta_TGS=meta_TGS[meta_TGS$PD_ID%in%tgs$Sample,]
rownames(tgs)=tgs$Sample
meta_TGS[,c("Age","Sex","Site","Number of glands")]=tgs[meta_TGS$PD_ID,c("Age","Sex","Site","Number of glands")]
colnames(meta_TGS)=c("Sample","Donor","Collaborator","Collaborator_ID","Slide_name","Well","Size","Age","Sex","Site","Number of glands")
write.csv(meta_TGS,"~/Desktop/Gastric/meta_TGS.csv",row.names = F)



#----
#At least 5 reads

mut_df_all=read.table("~/Desktop/Gastric/mut_df_all_chr20_PD40293.txt")
samples=unique(mut_df_all$sample)
samples=samples[samples%in%tgs_data$Sample]
samples_tri=pval_vec=vaf_diff=c()
for(sample in samples){
  mut_df=mut_df_all[mut_df_all$sample==sample&mut_df_all$NR>5,]
  ll0=sum(dbinom(x=mut_df$NV,size=mut_df$NR,p=sum(mut_df$NV)/sum(mut_df$NR),log=T))
  ll1=sum(dbinom(x=mut_df$NV[mut_df$hap=="hap_1"],size=mut_df$NR[mut_df$hap=="hap_1"],
                 p=sum(mut_df$NV[mut_df$hap=="hap_1"])/sum(mut_df$NR[mut_df$hap=="hap_1"]),log=T))+
    sum(dbinom(x=mut_df$NV[mut_df$hap=="hap_2"],size=mut_df$NR[mut_df$hap=="hap_2"],
               p=sum(mut_df$NV[mut_df$hap=="hap_2"])/sum(mut_df$NR[mut_df$hap=="hap_2"]),log=T))
  pval = 1-pchisq(2*(ll1-ll0),1)
  pval_vec=c(pval_vec,pval)
  
  meanVAF_hap1=sum(mut_df$NV[mut_df$hap=="hap_1"])/sum(mut_df$NR[mut_df$hap=="hap_1"])
  meanVAF_hap2=sum(mut_df$NV[mut_df$hap=="hap_2"])/sum(mut_df$NR[mut_df$hap=="hap_2"])
  vaf_diff=c(vaf_diff,meanVAF_hap1-meanVAF_hap2)
  #pdf(paste0(samples[n],"_VAF_chr2_boxplot.pdf"))
  #boxplot(VAF~hap,data=mut_df[mut_df$NR>5,],main=paste0(samples[n],"; diff: ",round(abs(mean_VAF2-mean_VAF1),digits=2),"; pval: ",pval))
  #dev.off()
  #mut_df_all=rbind(mut_df_all,mut_df)
}
#write.table(mut_df_all,"mut_df_all_chr20_PD40293.txt")
names(pval_vec)=names(vaf_diff)=samples

qval_vec=p.adjust(pval_vec,method="BH")
colvec=rep("grey60",length(qval_vec))
colvec[qval_vec<0.1&vaf_diff>0.1]="dodgerblue"
colvec[qval_vec<0.1&vaf_diff<(-0.1)]="orange"


pdf("~/Desktop/Gastric/PD40293_chr20_volcano.pdf",height=4,width=4)
plot(y=-log10(qval_vec),vaf_diff,xlim=c(-0.33,0.33),pch=21,bg=colvec,cex=1.5,lwd=0.5)
abline(h=-log10(0.1),lwd=2,lty="dashed",col="firebrick")
abline(v=c(-0.1,0.1),lwd=2,lty="dashed",col="firebrick")
dev.off()

samples_tri=names(qval_vec)[qval_vec<0.1&abs(vaf_diff)>0.1]


chr20_tgs=sy_tgs[sy_tgs$Donor=="PD40293",c("Sample","Slide_name","Well","Number of glands","Site","Chronic_inflammation","Intestinal_Metaplasia_in_Microdissection_region")]
chr20_tgs$Site[grepl("body",chr20_tgs$Slide_name)]="Body"
chr20_tgs$Site[grepl("antrum",chr20_tgs$Slide_name)]="Antrum"
chr20_tgs$VAF_diff=vaf_diff[chr20_tgs$Sample]
chr20_tgs$Tri=chr20_tgs$Sample%in%samples_tri



fisher.test(rbind(c(24,11,13),c(2,6,21)))
