# Script for the analysis of mutation burdens and rates in normal stomach
# Produces elements of Figure 1, Figure 2 and Extended Data Figure 2, and all statistical results presented in the section "Mutation rates of normal gastric epithelium"
# Tim Coorens
# November 2024

options(stringsAsFactors = F)
library(readxl)
library(lmerTest)
library(ape)
library(ggtree)
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=2))
wgs_data_all_normal=wgs_data[wgs_data$Feature=="Normal gland",]

#--------------------------------
# Figure 1
#--------------------------------

#Fig. 1c
#pdf("median_VAF_hist.pdf",width=3.5,height=4)
hist(as.numeric(wgs_data$median_VAF[wgs_data$Feature=="Normal gland"]),
     xlab="median VAF",col="steelblue",ylab="Number of microdissected single glands",main="",breaks=13)
abline(v=0.25,col="firebrick",lty="dashed",lwd=2)
#dev.off()

#Fig. 1d: SNV and indel burden vs age in normal glands from non-cancer donors
wgs_data_noncancer=wgs_data[wgs_data$Cohort=="Non-cancer donor",]

test_SNV_burden=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_noncancer)
summary(test_SNV_burden) #Intercept: 63.232; Age: 28.124; Age p-val: 0.00122
confint_lm=confint(test_SNV_burden) #CI: 16.8-39.4

pdf("~/Desktop/Gastric/SNV_burden_normal_only_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg="grey60",cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()

test_indel=lmer(Indel_burden_adjusted~Age+(1|Donor),data=wgs_data_noncancer)
summary(test_indel) #Intercept: -3.439; Age: 2.000; Age p-val: 0.0117
confint_lm_indel=confint(test_indel) #CI: 0.74-3.28

pdf("~/Desktop/Gastric/indel_burden_normal_only_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$Indel_burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of indels",cex=0.8,xlim=c(0,90),ylim=c(0,1000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm_indel[3,2],confint_lm_indel[3,2]+90*confint_lm_indel["Age",2],confint_lm_indel[3,1]+90*confint_lm_indel["Age",1],confint_lm_indel[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$Indel_burden_adjusted,bg="grey60",pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=-3.439,b=2.000,lwd=2,lty='dashed',col='firebrick')
dev.off()

#--------------------------------
#Extended Data Figure 1
#--------------------------------

#Extended Data Figure 1a
#Does median VAF correlate corrected mutation burden?
cor.test(as.numeric(wgs_data_noncancer$median_VAF),wgs_data_noncancer$SNV_Burden_adjusted)
# Correlation=0.03, p-value=0.76
pdf("~/Desktop/Gastric/medVAF_normal.pdf",width=4,height=4)
plot(wgs_data_all_normal$median_VAF[select],wgs_data_all_normal$SNV_Burden_adjusted[select],
     xlab="median VAF",ylab="SNV burden, adjusted",pch=21,bg="grey60",cex=1.5,ylim=c(0,5000))
dev.off()

#Extended Data Figure 1b
#Median VAF and intestinal metaplasia
wgs_data_all_normal$median_VAF=as.numeric(wgs_data_all_normal$median_VAF)
wilcox.test(wgs_data_all_normal$median_VAF[wgs_data_all_normal$Intestinal_Metaplasia_in_Sample=="Present"],
            wgs_data_all_normal$median_VAF[wgs_data_all_normal$Intestinal_Metaplasia_in_Sample=="Absent"])
#p-value=0.005664
pdf("~/Desktop/Gastric/medVAF_IM_boxplot.pdf",width=4,height=4)
boxplot(wgs_data_all_normal$median_VAF~wgs_data_all_normal$Intestinal_Metaplasia_in_Sample,ylab="Median VAF",xlab="Intestinal metaplasia")
dev.off()

#Extended Data Figure 1c
#SNV burden and HP status
donor_meta=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_R1_final.xlsx",sheet=1))
rownames(donor_meta)=donor_meta$Donor
wgs_data$HP_status=donor_meta[wgs_data$Donor,"H pylori status"]
wgs_data$HP_status[wgs_data$HP_status=="Previous"]="Positive" #Include "previous" in the known positive group

wgs_data_subset=wgs_data[wgs_data$HP_status%in%c("Positive","Negative")&wgs_data$Feature!="Tumour",]
test_SNV_sub=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_subset)
test_SNV_sub_HP=lmer(SNV_Burden_adjusted~Age+(1|Donor)+HP_status,data=wgs_data_subset)
anova(test_SNV_sub,test_SNV_sub_HP) #p=0.744

col=rep("steelblue",nrow(wgs_data_subset))
col[wgs_data_subset$HP_status=="Positive"]="firebrick"
pdf("~/Desktop/Gastric/pylori_scatter.pdf",height=4,width=4)
plot(wgs_data_subset$SNV_Burden_adjusted,x=wgs_data_subset$Age,bg=col,pch=21, xlab="Age",ylab="SNV burden",ylim=c(0,8000),xlim=c(0,90),cex=1.2,lwd=0.6)
dev.off()

#Extended Data Figure 1d
#Tumour SNV burden (clonal/subclonal)
wgs_data_cancer=wgs_data[wgs_data$Cohort=="Cancer patient"&wgs_data$Feature=="Tumour",]

tree_PD41759=as.data.frame(fortify(read.tree("~/Desktop/Gastric/Sequoia/PD41759_snv_tree_with_branch_length.tree")))
tree_PD41762=as.data.frame(fortify(read.tree("~/Desktop/Gastric/Sequoia/PD41762_snv_tree_with_branch_length.tree")))

pdf("~/Desktop/Gastric/SNV_burden_tumour_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_cancer$Age,y=wgs_data_cancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,75000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer$Age,y=wgs_data_cancer$SNV_Burden_adjusted,bg='steelblue',cex=1.2,lwd=0.6,pch=21)
points(x=c(59,60),y=c(max(tree_PD41759$branch.length),max(tree_PD41762$branch.length)),bg='black',cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=56.336,b=27.830,lwd=2,lty='dashed',col='firebrick')
dev.off()

#Extended Data Figure 1e
#Tumour indel burden (clonal/subclonal)
tree_PD41759=as.data.frame(fortify(read.tree("~/Desktop/Gastric/Sequoia/PD41759_indel_on_snv_tree_with_branch_length.tree")))
tree_PD41762=as.data.frame(fortify(read.tree("~/Desktop/Gastric/Sequoia/PD41762_indel_on_snv_tree_with_branch_length.tree")))

pdf("~/Desktop/Gastric/Indel_burden_tumour_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_cancer$Age,y=wgs_data_cancer$Indel_burden_adjusted,bg='white',pch=21,xlab="Age",ylab="Number of indel",cex=0.8,xlim=c(0,90),ylim=c(0,2000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm_indel[3,2],confint_lm_indel[3,2]+90*confint_lm_indel["Age",2],confint_lm_indel[3,1]+90*confint_lm_indel["Age",1],confint_lm_indel[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer$Age,y=wgs_data_cancer$Indel_burden_adjusted,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
points(x=c(59,60),y=c(max(tree_PD41759$branch.length),max(tree_PD41762$branch.length)),bg='black',cex=1.2,lwd=0.6,pch=21)

axis(1)
abline(a=-3.439,b=2,lwd=2,lty='dashed',col='firebrick')
dev.off()

#Extended Data Figure 1f
#Telomere length

wgs_data_all_glands=wgs_data[wgs_data$Feature=="Normal gland"&!wgs_data$Sample%in%novaseq&wgs_data$median_VAF>0.25,]
wgs_data_all_glands$has_trisomy=wgs_data_all_glands$Sample%in%SV_CNV_burden_all$Sample[SV_CNV_burden_all$Trisomy_burden>0]
wgs_data_all_glands$has_sv_cnv=wgs_data_all_glands$Sample%in%SV_CNV_burden_all$Sample[rowSums(SV_CNV_burden_all[,grepl("burden",colnames(SV_CNV_burden_all))])>0]
wgs_data_all_glands$has_CI=wgs_data_all_glands$Chronic_inflammation%in%c("Moderate","Severe")

test1=lmer(telomere_length~Age+(1|Donor),data=wgs_data_all_glands)
test2=lmer(telomere_length~Age+Intestinal_Metaplasia_in_Sample+(1|Donor),data=wgs_data_all_glands)


test2=lmer(telomere_length~Age+has_trisomy+(1|Donor),data=wgs_data_all_glands)
test2=lmer(telomere_length~Age+has_sv_cnv+(1|Donor),data=wgs_data_all_glands)
test2=lmer(telomere_length~Age+has_CI+(1|Donor),data=wgs_data_all_glands)
test2=lmer(telomere_length~Age+Cohort+(1|Donor),data=wgs_data_all_glands)

test3=lmer(telomere_length~Age+has_CI+median_VAF+(1|Donor),data=wgs_data_all_glands)

wgs_data_all_glands$plot_col="grey60"
wgs_data_all_glands$plot_col[wgs_data_all_glands$has_CI]="firebrick"

anova(test1,test2) #
anova(test2,test3) #

summary(test1) #Intercept: 56.336; Age: 27.830; Age p-val: 0.055
summary(test2)
confint_lm=confint(test2)
plot(y=wgs_data_all_glands$telomere_length,wgs_data_all_glands$SNV_Burden_adjusted,bg=wgs_data_all_glands$plot_col,cex=1.2,lwd=0.6,pch=21)


plot(y=wgs_data_all_glands$telomere_length,wgs_data_all_glands$Age,
     bg=wgs_data_all_glands$plot_col,cex=1.2,lwd=0.6,pch=21,ylab="Telomere length",xlab="Age")

pdf("~/Desktop/Gastric/telomere_length_combined_R1_CI.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_all_glands$Age,y=wgs_data_all_glands$telomere_length,bg='white',col='white',pch=21,xlab="Age",ylab="Telomere Length",cex=0.8,xlim=c(0,90),ylim=c(0,8000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_glands$Age,y=wgs_data_all_glands$telomere_length,bg=wgs_data_all_glands$plot_col,cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=5300.101,b=-38.557,lwd=2,lty='dashed',col='firebrick')
dev.off()
#--------------------------------
# Figure 2
#--------------------------------
wgs_data$site_cols=rep("goldenrod1",nrow(wgs_data))
wgs_data$site_cols[wgs_data$Site=="Cardia"]="firebrick"
wgs_data$site_cols[wgs_data$Site=="Antrum"]="steelblue"
wgs_data$site_cols[wgs_data$Site=="Body"]="chartreuse3"
wgs_data$site_cols[wgs_data$Site=="Pylorus"]="mediumorchid3"
wgs_data$site_cols[wgs_data$Site=="Incisuria"]="lightblue1"

### PLOTS OF MUTATION BURDEN VS AGE
wgs_data$IM_col="grey90"
wgs_data$IM_col[wgs_data$IM_type_crypt=="Incomplete"]="grey60"
wgs_data$IM_col[wgs_data$IM_type_crypt=="Complete"]="black"

wgs_data$CI_col="grey90"
wgs_data$CI_col[wgs_data$Chronic_inflammation=="Mild"]="peachpuff"
wgs_data$CI_col[wgs_data$Chronic_inflammation=="Moderate"]="salmon"
wgs_data$CI_col[wgs_data$Chronic_inflammation=="Severe"]="firebrick"

wgs_data_cancer_normal=wgs_data[wgs_data$Cohort=="Cancer patient"&wgs_data$Feature=="Normal gland",]
wgs_data_noncancer=wgs_data[wgs_data$Cohort=="Non-cancer donor",]

pdf("~/Desktop/Gastric/mutburden_normal_site_cols_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg=alpha(wgs_data_noncancer$site_cols,1),cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/mutburden_adj_site_cols_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg=alpha(wgs_data_cancer_normal$site_cols,1),cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/mutburden_normal_IM_cols_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg=alpha(wgs_data_noncancer$IM_col,1),cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/mutburden_adj_IM_cols_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg=alpha(wgs_data_cancer_normal$IM_col,1),cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/mutburden_normal_CI_cols_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg=alpha(wgs_data_noncancer$CI_col,1),cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/mutburden_adj_CI_cols_R1.pdf",useDingbats = F,width=4,height=4)
plot(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg=alpha(wgs_data_cancer_normal$CI_col,1),cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=63.232,b=28.124,lwd=2,lty='dashed',col='firebrick')
dev.off()
