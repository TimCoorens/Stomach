# Script for the analysis of mutation burdens and rates in normal stomach
# Produces elements of Figure 3, and all statistical results presented in the section "Mutational signatures and processes in normal gastric epithelium"
# Tim Coorens
# November 2024

options(stringsAsFactors = F)
library(readxl)
library(lmerTest)
library(ape)
library(ggtree)
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/Extended_Data_Table_1_7_final.xlsx",sheet=2))

#Extended Data Figure 5
wgs_data$SBS1_burden=wgs_data$SNV_Burden_adjusted*wgs_data$SBS1
wgs_data$SBS5_40_burden=wgs_data$SNV_Burden_adjusted*wgs_data$SBS5_40
wgs_data$SBS18_burden=wgs_data$SNV_Burden_adjusted*wgs_data$SBS18

wgs_data_noncancer=wgs_data[wgs_data$Cohort=="Non-cancer donor",]
wgs_data_all_normal_not_im=wgs_data[wgs_data$Feature=="Normal gland"&wgs_data$Intestinal_Metaplasia_in_Sample=="Absent",]

library(lmerTest)
test=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_CI=lmer(SNV_Burden_adjusted~Age+(1|Donor)+has_CI,data=wgs_data_all_normal_not_im)

test_SBS1=lmer(SBS1_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_SBS1_CI=lmer(SBS1_burden~Age+(1|Donor)+has_CI,data=wgs_data_all_normal_not_im)
anova(test_SBS1,test_SBS1_CI)

test_SBS5=lmer(SBS5_40_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_SBS5_CI=lmer(SBS5_40_burden~Age+(1|Donor)+has_CI,data=wgs_data_all_normal_not_im)
anova(test_SBS5,test_SBS5_CI)

test_SBS18=lmer(SBS18_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_SBS18_CI=lmer(SBS18_burden~Age+(1|Donor)+has_CI,data=wgs_data_all_normal_not_im)
anova(test_SBS18,test_SBS18_CI)

pdf("SBS1_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_SBS1)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$SBS1_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS1",cex=0.8,xlim=c(0,90),ylim=c(0,1500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$SBS1_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS1)$Donor[,1]),b=coef(test_SBS1)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("SBS5_40_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_SBS5)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$SBS5_40_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS5/40",cex=0.8,xlim=c(0,90),ylim=c(0,3000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$SBS5_40_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS5)$Donor[,1]),b=coef(test_SBS5)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("SBS18_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_SBS18)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$SBS18_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS18",cex=0.8,xlim=c(0,90),ylim=c(0,500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$SBS18_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS18)$Donor[,1]),b=coef(test_SBS18)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

test_SBS1=lmer(SBS1_burden~Age+(1|Donor),data=wgs_data_noncancer)
test_SBS5=lmer(SBS5_40_burden~Age+(1|Donor),data=wgs_data_noncancer)
test_SBS18=lmer(SBS18_burden~Age+Chronic_inflammation+(1|Donor),data=wgs_data_all_normal_not_im)

pdf("~/Desktop/Gastric/SBS1_burden_normal_2024_R1.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_SBS1)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SBS1_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS1",cex=0.8,xlim=c(0,90),ylim=c(0,1500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SBS1_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS1)$Donor[,1]),b=coef(test_SBS1)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/SBS5_40_burden_normal_2024_R1.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_SBS5)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SBS5_40_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS5/40",cex=0.8,xlim=c(0,90),ylim=c(0,3000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SBS5_40_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS5)$Donor[,1]),b=coef(test_SBS5)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("~/Desktop/Gastric/SBS18_burden_normal_2024_R1.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_SBS18)
plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SBS18_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS18",cex=0.8,xlim=c(0,90),ylim=c(0,500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SBS18_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS18)$Donor[,1]),b=coef(test_SBS18)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

#Fig. 3e
diff_mat=exp_mat=c()
for(sample in wgs_data$Sample[wgs_data$Intestinal_Metaplasia_in_Sample=="Present"&!is.na(wgs_data$Intestinal_Metaplasia_in_Sample)]){
  obs=wgs_data[sample,"SNV_Burden_adjusted"]*wgs_data[sample,grepl("SBS",colnames(wgs_data))]
  exp=wgs_data[sample,"Expected_SNV_burden"]*colMeans(wgs_data[wgs_data$Donor==substr(sample,1,7)&!is.na(wgs_data$Intestinal_Metaplasia_in_Sample=="Absent")&wgs_data$Intestinal_Metaplasia_in_Sample=="Absent",grepl("SBS",colnames(wgs_data))])
  diff=obs-exp
  diff=diff[,c("SBS1","SBS5_40","SBS18")]
  mat=rbind(round(diff),round(exp))
  mat=mat[,c("SBS1","SBS5_40","SBS18")]
  #chisq.test(mat)
  diff_mat=rbind(diff_mat,diff)
  exp_mat=rbind(exp_mat,exp[c("SBS1","SBS5_40","SBS18")])
  
}

pdf("~/Desktop/Gastric/fold_increase_hypermut_R1.pdf",width=4,height=6)
boxplot(diff_mat/exp_mat,col='grey60',ylab="Fold Increase Signature")
dev.off()

fc_mat=diff_mat/exp_mat
wilcox.test(fc_mat[,1],fc_mat[,2]) #p-value = 8.206e-07
wilcox.test(fc_mat[,2],fc_mat[,3]) #p-value = 9.037e-08

#Fig 3. f and g - indels
plot_cols=rep('steelblue',nrow(wgs_data))
plot_cols[wgs_data$Intestinal_Metaplasia_in_Sample=="Present"]='firebrick'
plot_cols[wgs_data$Feature=="Tumour"]='grey60'

ID2_dels=wgs_data$ID2*wgs_data$Indel_burden_adjusted
ID1_ins=wgs_data$ID1*wgs_data$Indel_burden_adjusted
id_ratio=data.frame(Cat="Normal",
                    Ratio=ID2_dels/ID1_ins)

id_ratio$Cat[wgs_data$Cohort=="Cancer patient"]="Normal, cancer donor"
id_ratio$Cat[wgs_data$Intestinal_Metaplasia_in_Sample=="Present"]="Intestinal Metaplasia"
id_ratio$Cat[wgs_data$Feature=="Tumour"&wgs_data$Donor=="PD41759"]="Tumour, PD41759"
id_ratio$Cat[wgs_data$Feature=="Tumour"&wgs_data$Donor=="PD41762"]="Tumour, PD41762"

id_ratio$Cat=factor(id_ratio$Cat , levels=c("Normal", "Normal, cancer donor","Intestinal Metaplasia", "Tumour, PD41762", "Tumour, PD41759"))

wilcox.test(id_ratio$Ratio[id_ratio$Cat=="Normal"],
            id_ratio$Ratio[id_ratio$Cat=="Normal, cancer donor"])
#p-value = 0.002398

wilcox.test(id_ratio$Ratio[id_ratio$Cat=="Intestinal Metaplasia"],
            id_ratio$Ratio[id_ratio$Cat=="Normal"])
#p-value = 7.917e-11

wilcox.test(id_ratio$Ratio[id_ratio$Cat=="Intestinal Metaplasia"],
            id_ratio$Ratio[id_ratio$Cat=="Normal, cancer donor"])
#p-value = 6.363e-12

pdf("~/Desktop/Gastric/id1_id2_burden_R1.pdf",width=4,height=4,useDingbats = F)
plot(ID2_dels,ID1_ins,pch=21,bg=plot_cols,cex=1.2,lwd=0.6,xlab="ID2",ylab="ID1")
dev.off()

pdf("~/Desktop/Gastric/id1_id2_ratio_R1.pdf",width=4,height=4,useDingbats = F)
boxplot(Ratio~Cat,data=id_ratio[incl,],col=c("steelblue","steelblue","firebrick","grey60","grey60"))
dev.off()

#Extended Data Figure 6
wgs_data$ID1_burden=wgs_data$Indel_burden_adjusted*wgs_data$ID1
wgs_data$ID2_burden=wgs_data$Indel_burden_adjusted*wgs_data$ID2
wgs_data$ID5_burden=wgs_data$Indel_burden_adjusted*wgs_data$ID5
wgs_data$ID9_burden=wgs_data$Indel_burden_adjusted*wgs_data$ID9


wgs_data_noncancer=wgs_data[wgs_data$Cohort=="Non-cancer donor",]
wgs_data_all_normal_not_im=wgs_data[wgs_data$Feature=="Normal gland"&wgs_data$Intestinal_Metaplasia_in_Sample=="Absent",]

library(lmerTest)
test=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_CI=lmer(SNV_Burden_adjusted~Age+(1|Donor)+has_CI,data=wgs_data_all_normal_not_im)

test_ID1=lmer(ID1_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_ID2=lmer(ID2_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_ID5=lmer(ID5_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)
test_ID9=lmer(ID9_burden~Age+(1|Donor),data=wgs_data_all_normal_not_im)

pdf("ID1_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_ID1)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID1_burden,bg='white',col='white',pch=21,xlab="Age",ylab="Indel burden due to ID1",cex=0.8,xlim=c(0,90),ylim=c(0,200),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID1_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_ID1)$Donor[,1]),b=coef(test_ID1)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("ID2_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_ID2)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID2_burden,bg='white',col='white',pch=21,xlab="Age",ylab="Indel burden due to ID2",cex=0.8,xlim=c(0,90),ylim=c(0,150),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID2_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_ID2)$Donor[,1]),b=coef(test_ID2)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("ID5_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_ID5)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID5_burden,bg='white',col='white',pch=21,xlab="Age",ylab="Indel burden due to ID5",cex=0.8,xlim=c(0,90),ylim=c(0,150),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID5_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_ID5)$Donor[,1]),b=coef(test_ID5)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

pdf("ID9_burden_non_IM_2024.pdf",useDingbats = F,width=4,height=4)
confint_lm=confint(test_ID9)
plot(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID9_burden,bg='white',col='white',pch=21,xlab="Age",ylab="Indel burden due to ID9",cex=0.8,xlim=c(0,90),ylim=c(0,30),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_all_normal_not_im$Age,y=wgs_data_all_normal_not_im$ID9_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_ID9)$Donor[,1]),b=coef(test_ID9)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()
