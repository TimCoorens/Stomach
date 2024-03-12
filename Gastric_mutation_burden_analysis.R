# Script for the analysis of mutation burdens and rates in normal stomach
# Produces elements of Figure 1, Figure 2 and Extended Data Figure 2, and all statistical results presented in the section "Mutation rates of normal gastric epithelium"
# Tim Coorens
# March 2024

options(stringsAsFactors = F)
library(readxl)
library(lmerTest)
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=2))

#Median coverage
median(wgs_data$Mean_coverage)

### PLOTS OF MUTATION BURDEN VS AGE

#Figure 1c: SNV burden vs age in normal glands from non-cancer donors
wgs_data_noncancer=wgs_data[wgs_data$Cohort=="Non-cancer donor",]

test=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_noncancer)
summary(test) #Intercept: 56.336; Age: 27.830; Age p-val: 0.00147
confint_lm=confint(test) #CI: 16.2-39.4

plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='steelblue',cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=56.336,b=27.830,lwd=2,lty='dashed',col='firebrick')

#Figure 1d: SNV burden vs age in normal glands from cancer patients
wgs_data_cancer_normal=wgs_data[wgs_data$Cohort=="Cancer patient"&wgs_data$Feature=="Normal gland",]

plot(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg='steelblue',cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=56.336,b=27.830,lwd=2,lty='dashed',col='firebrick')

#Extended Data Figure 2b: SNV burden vs age in tumours
wgs_data_cancer=wgs_data[wgs_data$Cohort=="Cancer patient"&wgs_data$Feature=="Tumour",]

plot(x=wgs_data_cancer$Age,y=wgs_data_cancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,75000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm[3,2],confint_lm[3,2]+90*confint_lm["Age",2],confint_lm[3,1]+90*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer$Age,y=wgs_data_cancer$SNV_Burden_adjusted,bg='steelblue',cex=1.2,lwd=0.6,pch=21)
axis(1)
abline(a=56.336,b=27.830,lwd=2,lty='dashed',col='firebrick')

#Figure 1e:indel burden vs age in normal glands from non-cancer donors
test_indel=lmer(Indel_burden_adjusted~Age+(1|Donor),data=wgs_data_noncancer)
summary(test_indel) #Intercept: -5.059; Age: 1.997; Age p-val: 0.0162
confint_lm_indel=confint(test_indel) #CI: 0.67-3.35

plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$Indel_burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of indels",cex=0.8,xlim=c(0,90),ylim=c(0,1000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm_indel[3,2],confint_lm_indel[3,2]+90*confint_lm_indel["Age",2],confint_lm_indel[3,1]+90*confint_lm_indel["Age",1],confint_lm_indel[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$Indel_burden_adjusted,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=-5.059,b=1.997,lwd=2,lty='dashed',col='firebrick')
dev.off()

#Extended Data Figure 2c: indel burden vs age in tumours
plot(x=wgs_data_cancer$Age,y=wgs_data_cancer$Indel_burden_adjusted,bg='white',pch=21,xlab="Age",ylab="Number of indel",cex=0.8,xlim=c(0,90),ylim=c(0,2000),xaxs = "i",yaxs = "i")
polygon(c(0.1,90,90,0.1),c(confint_lm_indel[3,2],confint_lm_indel[3,2]+90*confint_lm_indel["Age",2],confint_lm_indel[3,1]+90*confint_lm_indel["Age",1],confint_lm_indel[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer$Age,y=wgs_data_cancer$Indel_burden_adjusted,bg='steelblue',pch=21,cex=1.2,lwd=0.6)
axis(1)
abline(a=-5.059,b=1.997,lwd=2,lty='dashed',col='firebrick')

### HYPERMUTANT GLANDS

#Note these columns are already present in Extended Data Table 3
wgs_data$Expected_burden=56.336+27.830*wgs_data$Age
wgs_data$Upper_burden=confint_lm[3,2]+wgs_data$Age*confint_lm["Age",2]
wgs_data$Hypermut=wgs_data$SNV_Burden_adjusted>wgs_data$Upper_burden #Define the hypermutant glands
select=wgs_data$Feature=="Normal gland" #Exclude tumours from analysis

fisher.test(table(wgs_data$Hypermut[select],wgs_data$Cohort[select])) #p-val: 5.9x10-5

wgs_data$indel_snv_ratio=wgs_data$Indel_burden_adjusted/wgs_data$SNV_Burden_adjusted

Comp=wgs_data$Cohort
Comp[wgs_data$Hypermut]=paste0(wgs_data$Cohort[wgs_data$Hypermut],"_hypermut")

#Fig 1g                                                          
boxplot(wgs_data$indel_snv_ratio[select]~Comp[select],ylim=c(0,0.15),col='steelblue',ylab="indel/SNV ratio")

wilcox.test(wgs_data$indel_snv_ratio[Comp=="Non-cancer donor"&select],
            wgs_data$indel_snv_ratio[Comp=="Cancer patient"&select]) #p=0.026
wilcox.test(wgs_data$indel_snv_ratio[Comp=="Cancer patient_hypermut"&select],
            wgs_data$indel_snv_ratio[Comp=="Cancer patient"&select]) #p=0.0009
wilcox.test(wgs_data$indel_snv_ratio[Comp=="Cancer patient_hypermut"&select],
            wgs_data$indel_snv_ratio[Comp=="Non-cancer donor"&select]) #p=6.403e-05

#-----
#Test whether fitting stomach section-specific mutation rate improves the model
test_wo_site=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_noncancer[wgs_data_noncancer$Site!="Unknown",])
test_w_site=lmer(SNV_Burden_adjusted~Site:Age+(1|Donor),data=wgs_data_noncancer[wgs_data_noncancer$Site!="Unknown",])
anova(test_wo_site,test_w_site)#p=0.1471

#Test enrichment for hypermutants in antrum
fisher.test(table(wgs_data$Hypermut[select],wgs_data$Site[select]=="Antrum")) #p-val: 6.218e-15

#Fig. 2
site_cols=rep("goldenrod1",nrow(wgs_data))
site_cols[wgs_data$Site=="Cardia"]="firebrick"
site_cols[wgs_data$Site=="Antrum"]="steelblue"
site_cols[wgs_data$Site=="Body"]="chartreuse3"
site_cols[wgs_data$Site=="Pylorus"]="mediumorchid3"
site_cols[wgs_data$Site=="Incisuria"]="lightblue1"
site_cols[wgs_data$Site=="Unknown"]="grey80"
names(site_cols)=wgs_data$Sample

plot(x=wgs_data_noncancer$Age,y=wgs_data_noncancer$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_noncancer$Age[select],y=wgs_data_noncancer$SNV_Burden_adjusted[select],bg=alpha(site_cols[wgs_data_noncancer$Sample],1),cex=1.5,lwd=0.6,pch=21)
axis(1)
abline(a=56.336,b=27.830,lwd=2,lty='dashed',col='firebrick')

plot(x=wgs_data_cancer_normal$Age,y=wgs_data_cancer_normal$SNV_Burden_adjusted,bg='white',col='white',pch=21,xlab="Age",ylab="Number of SNVs",cex=0.8,xlim=c(0,90),ylim=c(0,7500),xaxs = "i",yaxs = "i")
polygon(c(0.1,89.8,89.8,0.1),c(confint_lm[3,2],confint_lm[3,2]+89.8*confint_lm["Age",2],confint_lm[3,1]+89.8*confint_lm["Age",1],confint_lm[3,1]),border = NA,col = "grey90")
points(x=wgs_data_cancer_normal$Age[select],y=wgs_data_cancer_normal$SNV_Burden_adjusted[select],bg=alpha(site_cols[wgs_data_cancer_normal$Sample],1),cex=1.5,lwd=0.6,pch=21)
axis(1)
abline(a=56.336,b=27.830,lwd=2,lty='dashed',col='firebrick')

#Test enrichment for H pylori status in mutation burden model
donor_meta_data=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=1))
rownames(donor_meta_data)=donor_meta_data$Donor

wgs_data$HP_status=donor_meta_data[wgs_data$Donor,"H pylori status"]
wgs_data$HP_status[wgs_data$HP_status=="Previous"]="Positive" #Include "previous" in the known positive group

wgs_data_subset=wgs_data[wgs_data$HP_status%in%c("Positive","Negative")&wgs_data$Feature!="Tumour",]
library(lmerTest)
test1=lmer(SNV_Burden_adjusted~Age+(1|Donor)+HP_status,data=wgs_data_subset)
test2=lmer(SNV_Burden_adjusted~Age+(1|Donor),data=wgs_data_subset)
anova(test1,test2) #p=0.07493

col=rep("steelblue",nrow(mutation_burden_select))
col[mutation_burden_select$HP_status=="Positive"]="firebrick"
plot(mutation_burden_select$SNV_Burden_adjusted,x=mutation_burden_select$Age,bg=col,pch=21, xlab="Age",ylab="SNV burden",ylim=c(0,8000),xlim=c(0,90),cex=1.2,lwd=0.6)





