options(stringsAsFactors = F)
library(readxl)
drivers=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=6))
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=2))
tgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=3))
tgs_data=tgs_data[!is.na(tgs_data$Chronic_inflammation),]
met_patient=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=1))
cnvs=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=4))

wgs_data_all_normal=wgs_data[wgs_data$Feature=="Normal gland",]

fisher.test(table(wgs_data_all_normal$Sample%in%drivers$Sample[drivers$Seq=="WGS"],wgs_data_all_normal$Intestinal_Metaplasia_in_Sample))
fisher.test(table(wgs_data_all_normal$Sample%in%drivers$Sample[drivers$Seq=="WGS"],wgs_data_all_normal$Chronic_inflammation%in%"Severe"))

saveRDS(dndsout,"dndsout.RDS")

sel_cv = dndsout$sel_cv
rownames(sel_cv)=sel_cv$gene_name
#sel_cv=sel_cv[!grepl("HLA",sel_cv$gene_name),]
#sel_cv=sel_cv[!sel_cv$gene_name%in%c("TTN","EEF1A1"),]

print(head(sel_cv,n=10), digits = 3)

signif_genes = sel_cv[sel_cv$qglobal_cv<0.1,]
print(dndsout_rest$globaldnds, digits = 3)

plot_mat=sel_cv[rownames(signif_genes),2:6]
plot_mat=plot_mat[order(rowSums(plot_mat),decreasing = T),]
pdf("dnds_counts_R1.pdf")
barplot(t(plot_mat),las=2, col=c("grey80","dodgerblue2","firebrick","salmon","goldenrod1"))
dev.off()

plot_mat2=sel_cv[rownames(signif_genes),c(7,9,10)]
plot_mat2=plot_mat2[order(rowSums(plot_mat),decreasing = T),]

pdf("dnds_ratios_R1.pdf")
barplot(t(plot_mat2),beside=T,las=2, col=c("dodgerblue2","firebrick","goldenrod1"))
dev.off()


annot$Patient=substr(annot$sampleID,1,7)
annot_nonsym=annot[annot$impact!="Synonymous"&annot$gene%in%rownames(signif_genes),]

driver_freq=table(annot_nonsym$gene,annot_nonsym$Patient)

meta_data=read.csv("master_metadata_gastric_lcm_31OCT2019.csv")

met_patient=unique(meta_data[,c("PD_STEM","AGE","SEX")])
rownames(met_patient)=met_patient$PD_STEM
patients=met_patient$PD_STEM[order(met_patient$SEX,met_patient$AGE,decreasing = F)]
patients=patients[patients%in%unique(muts$sampleID)]

driver_freq_full=matrix(0,nrow=nrow(driver_freq),ncol=length(patients))
rownames(driver_freq_full)=rownames(driver_freq)
colnames(driver_freq_full)=patients
driver_freq_full[,colnames(driver_freq)]=driver_freq
driver_freq_full=driver_freq_full[order(rowSums(driver_freq_full)),]

grid = expand.grid(x=rownames(driver_freq_full), y=colnames(driver_freq_full))
grid$z = c(driver_freq_full)
pdf("driver_heatmap_R1.pdf",width=10,height=4,useDingbats = F)
#dev.new(width=6.5, height=6)
library(lattice)
color.palette = colorRampPalette(c("white", "darkorange", "darkorchid2","darkorchid4"))
levelplot(z~y*x, grid, col.regions=color.palette, scales = list(tck = c(0,0), y = list(cex=1), x = list(rot=90)), ylab="Genes", xlab="Samples", colorkey=list(space="bottom"), 
          panel=function(...) { arg <- list(...)
          panel.levelplot(...)
          panel.text(arg$x, arg$y, arg$z)})
#dev.copy(pdf, file="Driver_heatmap_bladder_finalcalls.pdf", width=6.5, height=6, useDingbats=F); dev.off(); dev.off()
dev.off()

normal_TGS_samples=read.table("TGS/normal_samples.txt")[,1]

nsamples = table(meta_data$PD_STEM[(meta_data$TGSseq|meta_data$WGSseq)&meta_data$GASTRIC_SITE!="Small intestine"&meta_data$TUM_NORM_ADJNORM!="Tumour"])[colnames(drivers_mat_combined_select)]
#sample_rm=colnames(drivers_mat_combined_select)[is.na(nsamples)]
#drivers_mat_combined_select=drivers_mat_combined_select[,colnames(drivers_mat_combined_select)!=sample_rm]
#met_patient=met_patient[rownames(met_patient)!=sample_rm,]
#nsamples=nsamples[!is.na(nsamples)]
meta = rbind(nsamples,met_patient$AGE)
rownames(meta) = c("samples","age")

meta_data=read.csv("master_metadata_gastric_lcm_31OCT2019.csv")
rownames(meta_data)=meta_data$PD_ID
meta_data_norm=meta_data[meta_data$PD_ID%in%c(normal_TGS_samples,mut_burden_select$Sample),]
nsamples=table(meta_data_norm$PD_STEM,meta_data_norm$WGSseq)

meta_data_norm$CLONE_ESTIMATE[is.na(meta_data_norm$CLONE_ESTIMATE)]=10

patient_glands_surveyed=rep(0,length(patients))
names(patient_glands_surveyed)=patients
for(p in patients){
  patient_glands_surveyed[p]=sum(meta_data_norm$CLONE_ESTIMATE[meta_data_norm$PD_STEM==p],na.rm=T)
}


met_patient=met_patient[patients,]
meta = rbind(patient_glands_surveyed,met_patient$AGE)

pdf("drivers_meta_2024.pdf",width=10,height=3,useDingbats = F)
grid = expand.grid(x=rownames(meta), y=colnames(meta))
grid$z = c(meta)
#dev.new(width=6.5, height=2.5)
library(lattice)
color.palette = colorRampPalette(c("white","grey50"))
levelplot(z~y*x, grid, col.regions=color.palette, scales = list(tck = c(0,0), y = list(cex=1), x = list(rot=90)), ylab="", xlab="Samples", colorkey=list(space="bottom"), 
          panel=function(...) { arg <- list(...)
          panel.levelplot(...)
          panel.text(arg$x, arg$y, arg$z)})
dev.off()


#----
#EDF 9

wgs_data_all_normal$CI_score=tgs_data$CI_score=0
wgs_data_all_normal$CI_score[wgs_data_all_normal$Chronic_inflammation=="Mild"]=1
wgs_data_all_normal$CI_score[wgs_data_all_normal$Chronic_inflammation=="Moderate"]=2
wgs_data_all_normal$CI_score[wgs_data_all_normal$Chronic_inflammation=="Severe"]=3
tgs_data$CI_score[tgs_data$Chronic_inflammation=="Mild"]=1
tgs_data$CI_score[tgs_data$Chronic_inflammation=="Moderate"]=2
tgs_data$CI_score[tgs_data$Chronic_inflammation=="Severe"]=3

wgs_data_all_normal$IM_score=tgs_data$IM_score=0
wgs_data_all_normal$IM_score[wgs_data_all_normal$Intestinal_Metaplasia_in_Sample=="Present"]=1
tgs_data$IM_score[tgs_data$Intestinal_Metaplasia_in_Microdissection_region=="Present"]=1

size_df=as.data.frame(rbind(wgs_data_all_normal[,c("Sample","Donor","Size","CI_score","IM_score")],
                            tgs_data[,c("Sample","Donor","Size","CI_score","IM_score")]))
rownames(size_df)=size_df$Sample
patients=met_patient$Donor
male_patients=met_patient$Donor[met_patient$Sex=="Male"]
mutant_epithelium=patient_epithelium_surveyed=rep(0,length(patients))
names(mutant_epithelium)=names(patient_epithelium_surveyed)=patients

for(n in 1:nrow(drivers)){
  if(drivers$Chr[n]=="chrX"&drivers$Donor[n]%in%male_patients){
    mutant_size=drivers$VAF[n]*size_df$Size[size_df$Sample==drivers$Sample[n]]
  }else{
    mutant_size=drivers$VAF[n]*size_df$Size[size_df$Sample==drivers$Sample[n]]*2
  }
  mutant_epithelium[drivers$Donor[n]]=mutant_epithelium[drivers$Donor[n]]+mutant_size
}
for(d in patients){
  patient_epithelium_surveyed[d]=sum(size_df$Size[size_df$Donor==d])
}

plot(y=as.numeric(mutant_epithelium)/as.numeric(total_epithelium),x=donors$Age)

met_patient$mean_snv_burden=NA
met_patient$CI=NA
met_patient$IM=NA

for(n in 1:nrow(met_patient)){
  patient=met_patient$Donor[n]
  met_patient$mean_snv_burden[n]=mean(wgs_data_all_normal$SNV_Burden_adjusted[wgs_data_all_normal$Donor==patient])+mean(wgs_data_all_normal$Indel_burden_adjusted[wgs_data_all_normal$Donor==patient])
  
  # ci_status=unique(wgs_data_all_normal$Chronic_inflammation[wgs_data_all_normal$Donor==met_patient$PD_STEM[n]])
  # im_status=unique(wgs_data_all_normal$Intestinal_Metaplasia_extent_surrounding_region[wgs_data_all_normal$Donor==met_patient$PD_STEM[n]])
  # 
  # if(length(ci_status)>1){
  #   if("Severe"%in%ci_status)ci_status="Severe"
  #   if("Moderate"%in%ci_status)ci_status="Moderate"
  #   if("Mild"%in%ci_status)ci_status="Mild"
  # }
  # if(length(im_status)>1){
  #   if("Severe"%in%im_status)im_status="Severe"
  #   if("Moderate"%in%im_status)im_status="Moderate"
  #   if("Mild"%in%im_status)im_status="Mild"
  # }
  met_patient$CI[n]=sum(size_df$CI_score[size_df$Donor==patient]*size_df$Size[size_df$Donor==patient])/sum(size_df$Size[size_df$Donor==patient])
  met_patient$IM[n]=sum(size_df$IM_score[size_df$Donor==patient]*size_df$Size[size_df$Donor==patient])/sum(size_df$Size[size_df$Donor==patient])
  
}

met_patient$col="grey65"
met_patient$col[met_patient$CI>0.5]="peachpuff"
met_patient$col[met_patient$CI>1.5]="salmon"
met_patient$col[met_patient$CI>2.5]="firebrick"

met_patient$Avg_CI="Absent"
met_patient$Avg_CI[met_patient$CI>0.5]="Mild"
met_patient$Avg_CI[met_patient$CI>1.5]="Moderate"
met_patient$Avg_CI[met_patient$CI>2.5]="Severe"


select=patient_epithelium_surveyed>500000

plot(y=mutant_epithelium[select]/patient_epithelium_surveyed[select],x=met_patient$Age[select],xlim=c(0,90),
     pch=21,bg='steelblue',xlab = "Age", ylab="")


test=lm(mutant_epithelium[select]/patient_epithelium_surveyed[select]~met_patient$Age[select])

test_CI=lm(mutant_epithelium[select]/patient_epithelium_surveyed[select]~met_patient$Age[select]+met_patient$Avg_CI[select])
confint_lm=confint(test_CI)

met_patient$col="grey65"
met_patient$col[met_patient$CI>0.5]="peachpuff"
met_patient$col[met_patient$CI>1.5]="salmon"
met_patient$col[met_patient$CI>2.5]="firebrick"
anova(test,test_CI) #p=0.002

pdf("~/Desktop/Gastric/gastric_driver_age_CI_v2.pdf",width=6,height=4,useDingbats = F)
plot(x=met_patient$Age,y=mutant_epithelium/patient_epithelium_surveyed,
     bg='white',col='white',pch=21,xlab="Age",ylab="Proportion of epithelium with driver",cex=0.8,xlim=c(0,90),ylim=c(-0.01,0.35),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[1,2],confint_lm[1,2]+90*confint_lm[2,2],confint_lm[1,1]+90*confint_lm[2,1],confint_lm[1,1]),border = NA,col = "grey90")
points(x=met_patient$Age[select],y=mutant_epithelium[select]/patient_epithelium_surveyed[select],bg=met_patient$col[select],pch=21,cex=1.6,lwd=0.6)
axis(1)
abline(a=coef(test_CI)[1],b=coef(test_CI)[2],lwd=2,lty='dashed',col='firebrick')
dev.off()

met_patient$col="grey65"
met_patient$col[met_patient$IM>0]="peachpuff"
met_patient$col[met_patient$IM>0.1]="salmon"
met_patient$col[met_patient$IM>0.25]="firebrick"

test_IM=lm(mutant_epithelium[select]/patient_epithelium_surveyed[select]~met_patient$Age[select]+met_patient$IM[select])
anova(test,test_IM) #p=0.86
confint_lm=confint(test_CI)

pdf("~/Desktop/Gastric/gastric_driver_age_IM_v2.pdf",width=6,height=4,useDingbats = F)
plot(x=met_patient$Age,y=mutant_epithelium/patient_epithelium_surveyed,
     bg='white',col='white',pch=21,xlab="Age",ylab="Proportion of epithelium with driver",cex=0.8,xlim=c(0,90),ylim=c(-0.01,0.35),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[1,2],confint_lm[1,2]+90*confint_lm[2,2],confint_lm[1,1]+90*confint_lm[2,1],confint_lm[1,1]),border = NA,col = "grey90")
points(x=met_patient$Age[select],y=mutant_epithelium[select]/patient_epithelium_surveyed[select],bg=met_patient$col[select],pch=21,cex=1.6,lwd=0.6)
axis(1)
abline(a=coef(test)[1],b=coef(test)[2],lwd=2,lty='dashed',col='firebrick')
dev.off()
plot(mutant_epithelium/patient_epithelium_surveyed,log10(patient_epithelium_surveyed))

test_burden=lm(mutant_epithelium[select]/patient_epithelium_surveyed[select]/met_patient$mean_snv_burden[select]~met_patient$Age[select]+met_patient$Avg_CI[select])
confint_lm=confint(test_burden)
pdf("~/Desktop/Gastric/gastric_driver_age_CI_SNV_norm_v2.pdf",width=6,height=4,useDingbats = F)
plot(x=met_patient$mean_snv_burden,y=mutant_epithelium/patient_epithelium_surveyed,
     bg='white',col='white',pch=21,xlab="Age",ylab="Proportion of epithelium with driver",cex=0.8,xlim=c(0,90),ylim=c(-0.000005,0.00015),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[1,2],confint_lm[1,2]+90*confint_lm[2,2],confint_lm[1,1]+90*confint_lm[2,1],confint_lm[1,1]),border = NA,col = "grey90")
points(x=met_patient$Age[select],y=mutant_epithelium[select]/patient_epithelium_surveyed[select]/met_patient$mean_snv_burden[select],bg=met_patient$col[select],pch=21,cex=1.6,lwd=0.6)
axis(1)
abline(a=coef(test_burden)[1],b=coef(test_burden)[2],lwd=2,lty='dashed',col='firebrick')
dev.off()

test_burden=lm(100*mutant_epithelium[select]/patient_epithelium_surveyed[select]~met_patient$mean_snv_burden[select])
test_burden_CI=lm(100*mutant_epithelium[select]/patient_epithelium_surveyed[select]~met_patient$mean_snv_burden[select]+met_patient$Avg_CI[select])

confint_lm=confint(test_burden_CI)
pdf("~/Desktop/Gastric/gastric_driver_burden_CI_v2.pdf",width=6,height=4,useDingbats = F)
plot(x=met_patient$mean_snv_burden,y=100*mutant_epithelium/patient_epithelium_surveyed,
     bg='white',col='white',pch=21,xlab="Mean SNV+Indel burden",ylab="Percentage of epithelium with driver",cex=0.8,xlim=c(0,4900),ylim=c(-1,35),xaxs = "i",yaxs = "i")
polygon(c(0.1,5000,5000,0.1),c(confint_lm[1,2],confint_lm[1,2]+5000*confint_lm[2,2],confint_lm[1,1]+5000*confint_lm[2,1],confint_lm[1,1]),border = NA,col = "grey90")
points(x=met_patient$mean_snv_burden[select],y=100*mutant_epithelium[select]/patient_epithelium_surveyed[select],bg=met_patient$col[select],pch=21,cex=1.6,lwd=0.6)
axis(1)
abline(a=coef(test_burden_CI)[1],b=coef(test_burden_CI)[2],lwd=2,lty='dashed',col='firebrick')
dev.off()

test_burden=lm(mutant_epithelium[select]/patient_epithelium_surveyed[select]~met_patient$Age[select]+met_patient$mean_snv_burden[select])
test_burden=lm(mutant_epithelium[select]/patient_epithelium_surveyed[select]/met_patient$mean_snv_burden[select]~met_patient$Age[select])
plot(mutant_epithelium[select]/patient_epithelium_surveyed[select],x=met_patient$mean_snv_burden[select],bg=met_patient$col[select],pch=21)

anova(test,test_burden)

pdf("~/Desktop/Gastric/gastric_driver_epithelium_surveyed.pdf",width=6,height=4,useDingbats = F)
plot(y=100*mutant_epithelium/patient_epithelium_surveyed,x=patient_epithelium_surveyed,log="x",pch=21,bg="steelblue",xlab="Total area epithelium surveyed (um2)",ylab="Percentage of epithelium with driver",cex=1.6,lwd=0.6)
abline(v=500000,lty="dashed",col="firebrick")
dev.off()
