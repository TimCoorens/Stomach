# Script for the analysis of mutational signatures in normal stomach
# Produces elements of Figure 3, Extended Data Figure 5 and all statistical results presented in the section "Mutational signatures and processes in normal gastric epithelium"
# Tim Coorens
# March 2024

options(stringsAsFactors = F)
library(readxl)
library(lmerTest)
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=2))

avg_normal_sig_dist=colMeans(wgs_data[!wgs_data$Hypermut,grepl("SBS",colnames(wgs_data))])
sigs=colnames(wgs_data)[grepl("SBS",colnames(wgs_data))]

wgs_data$SBS5_40=wgs_data$SBS5+wgs_data$SBS40 #Combined SBS5 and SBS40
wgs_data=wgs_data[,!colnames(wgs_data)%in%c("SBS5","SBS40")]
wgs_data_sub=wgs_data[wgs_data$Feature!="Tumour",]
rownames(wgs_data_sub)=wgs_data_sub$Sample
diff_mat=exp_mat=c()

for(sample in wgs_data_sub$Sample[wgs_data_sub$Hypermut]){
  obs=wgs_data_sub[sample,"SNV_Burden_adjusted"]*wgs_data_sub[sample,grepl("SBS",colnames(wgs_data_sub))]
  exp=wgs_data_sub[sample,"Expected_SNV_burden"]*colMeans(wgs_data_sub[wgs_data_sub$Donor==substr(sample,1,7)&!wgs_data_sub$Hypermut,grepl("SBS",colnames(wgs_data_sub))])
  diff=obs-exp
  diff=diff[,c("SBS1","SBS5_40","SBS18")]
  mat=rbind(round(diff),round(exp))
  mat=mat[,c("SBS1","SBS5_40","SBS18")]
  #chisq.test(mat)
  diff_mat=rbind(diff_mat,diff)
  exp_mat=rbind(exp_mat,exp[c("SBS1","SBS5_40","SBS18")])
  
}

#Fig. 2e
boxplot(diff_mat/exp_mat,col='grey60',ylab="Fold Increase Signature")
dev.off()

fold_increase=as.data.frame(diff_mat/exp_mat)
wilcox.test(fold_increase$SBS1,fold_increase$SBS5_40) #p-value = 1.549e-08
wilcox.test(fold_increase$SBS5_40,fold_increase$SBS18) #p-value = 4.483e-06

#-----
#Extended Data Figure 5

wgs_data_sub$SBS1_burden=wgs_data_sub$SNV_Burden_adjusted*wgs_data_sub$SBS1
wgs_data_sub$SBS5_40_burden=wgs_data_sub$SNV_Burden_adjusted*wgs_data_sub$SBS5_40
wgs_data_sub$SBS18_burden=wgs_data_sub$SNV_Burden_adjusted*wgs_data_sub$SBS18

library(lmerTest)

select=wgs_data_sub$Cohort=="Non-cancer donor"&!wgs_data_sub$Hypermut
test_SBS1=lmer(SBS1_burden~Age+(1|Donor),data=wgs_data_sub[select,])
test_SBS5=lmer(SBS5_40_burden~Age+(1|Donor),data=wgs_data_sub[select,])
test_SBS18=lmer(SBS18_burden~Age+(1|Donor),data=wgs_data_sub[select,])

confint_lm=confint(test_SBS1)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS1_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS1",cex=0.8,xlim=c(0,90),ylim=c(0,1500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS1_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS1)$Donor[,1]),b=coef(test_SBS1)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

confint_lm=confint(test_SBS5)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS5_40_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS5/40",cex=0.8,xlim=c(0,90),ylim=c(0,3000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS5_40_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS5)$Donor[,1]),b=coef(test_SBS5)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

confint_lm=confint(test_SBS18)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS18_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS18",cex=0.8,xlim=c(0,90),ylim=c(0,500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS18_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS18)$Donor[,1]),b=coef(test_SBS18)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

# Repeat for all glands, except hypermutants
select=!wgs_data_sub$Hypermut
test_SBS1=lmer(SBS1_burden~Age+(1|Donor),data=wgs_data_sub[select,])
test_SBS5=lmer(SBS5_40_burden~Age+(1|Donor),data=wgs_data_sub[select,])
test_SBS18=lmer(SBS18_burden~Age+(1|Donor),data=wgs_data_sub[select,])

confint_lm=confint(test_SBS1)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS1_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS1",cex=0.8,xlim=c(0,90),ylim=c(0,1500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS1_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS1)$Donor[,1]),b=coef(test_SBS1)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

confint_lm=confint(test_SBS5)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS5_40_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS5/40",cex=0.8,xlim=c(0,90),ylim=c(0,3000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS5_40_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS5)$Donor[,1]),b=coef(test_SBS5)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

confint_lm=confint(test_SBS18)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS18_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS18",cex=0.8,xlim=c(0,90),ylim=c(0,500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS18_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS18)$Donor[,1]),b=coef(test_SBS18)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

#-----

select=wgs_data_sub$Cohort=="Non-cancer donor"&!wgs_data_sub$Hypermut
test_SBS1=lmer(SBS1_burden~Age+(1|Donor),data=wgs_data_sub[select,])
test_SBS5=lmer(SBS5_40_burden~Age+(1|Donor),data=wgs_data_sub[select,])
test_SBS18=lmer(SBS18_burden~Age+(1|Donor),data=wgs_data_sub[select,])

confint_lm=confint(test_SBS1)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS1_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS1",cex=0.8,xlim=c(0,90),ylim=c(0,1500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS1_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS1)$Donor[,1]),b=coef(test_SBS1)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

confint_lm=confint(test_SBS5)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS5_40_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS5/40",cex=0.8,xlim=c(0,90),ylim=c(0,3000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS5_40_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS5)$Donor[,1]),b=coef(test_SBS5)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')
dev.off()

confint_lm=confint(test_SBS18)
plot(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS18_burden,bg='white',col='white',pch=21,xlab="Age",ylab="SNV burden due to SBS18",cex=0.8,xlim=c(0,90),ylim=c(0,500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_sub[select,]$Age,y=wgs_data_sub[select,]$SBS18_burden,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=mean(coef(test_SBS18)$Donor[,1]),b=coef(test_SBS18)$Donor[1,2],lwd=2,lty='dashed',col='firebrick')

#------------
# Fig. 2f

cols=rep('steelblue',nrow(wgs_data))
cols[wgs_data$Hypermut]='firebrick'
cols[wgs_data$Feature=="Tumour"]='grey60'
plot(wgs_data$ID2_dels,wgs_data$ID1_ins,pch=21,bg=cols)

incl=(wgs_data$ID1_ins+wgs_data$ID2_dels)>25

id_ratio=data.frame(Cat="Normal",
                    Ratio=wgs_data$ID2_dels/wgs_data$ID1_ins)
id_ratio$Cat[wgs_data$Cohort=="Cancer patient"]="Normal, cancer donor"
id_ratio$Cat[wgs_data$Hypermut]="Hypermutant"
id_ratio$Cat[wgs_data$Feature=="Tumour"&wgs_data$Donor=="PD41759"]="Tumour, PD41759"
id_ratio$Cat[wgs_data$Feature=="Tumour"&wgs_data$Donor=="PD41762"]="Tumour, PD41762"

id_ratio$Cat=factor(id_ratio$Cat , levels=c("Normal", "Normal, cancer donor","Hypermutant", "Tumour, PD41759", "Tumour, PD41762"))

# Fig. 2g
boxplot(Ratio~Cat,data=id_ratio[incl,],col=c("steelblue","steelblue","firebrick","grey60","grey60"))

wilcox.test(id_ratio$Ratio[id_ratio$Cat=="Normal"&incl],
            id_ratio$Ratio[id_ratio$Cat=="Normal, cancer donor"&incl]) #p-value = 0.002742

wilcox.test(id_ratio$Ratio[id_ratio$Cat=="Hypermutant"&incl],
            id_ratio$Ratio[id_ratio$Cat=="Normal"&incl]) #p-value = 2.598e-06

wilcox.test(id_ratio$Ratio[id_ratio$Cat=="Hypermutant"&incl],
            id_ratio$Ratio[id_ratio$Cat=="Normal, cancer donor"&incl]) #p-value = 4.146e-05

