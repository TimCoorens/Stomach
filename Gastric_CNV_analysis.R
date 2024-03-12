# Script for the analysis of mutation burdens and rates in normal stomach
# Produces elements of Figure 4 and all statistical results presented in the section "Recurrent trisomies in normal gastric glands"
# Tim Coorens
# March 2024

options(stringsAsFactors = F)
library(readxl)
library(lmerTest)
cnvs=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=4))
wgs_data=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=2))
met_patient=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=1))
rownames(met_patient)=met_patient$Donor
nsamples=table(wgs_data$Donor[wgs_data$Feature!="Tumour"])

cnv_type=cnvs$CNV
cnv_type[grepl('gain',cnvs$CNV)]="Trisomy"
cnv_type[grepl('LOH',cnvs$CNV)]="Arm cnn-LOH"

CNV_freq=table(cnv_type,cnvs$Donor)

CNV_freq_full=matrix(0,nrow=nrow(CNV_freq),ncol=length(nsamples))
rownames(CNV_freq_full)=rownames(CNV_freq)
colnames(CNV_freq_full)=names(nsamples)
CNV_freq_full[,colnames(CNV_freq)]=CNV_freq

test1=lm(colSums(CNV_freq_full)/nsamples~met_patient[colnames(CNV_freq_full),"Age"]) #p=0.12
test2=lm(colSums(CNV_freq_full[c("Deletion","Duplication"),])/nsamples~met_patient[colnames(CNV_freq_full),"Age"]) #p-value: 0.008843
test3=lm(CNV_freq_full["Trisomy",]/nsamples~met_patient[colnames(CNV_freq_full),"Age"]) #p=0.4

cnvs$CNV_id=paste0(cnvs$CNV,"_",cnvs$Fragile_site,"_",cnvs$Donor,"_",cnvs$Event)
CNV_uniq=unique(cnvs$CNV_id)
c(sum(cnv_type%in%c("Deletion","Duplication")),)

short=c(sum(grepl("Deletion|Duplication",CNV_uniq)), #total number
        sum(grepl("Deletion_NA",CNV_uniq)),
        sum(grepl("Deletion_FHIT",CNV_uniq)),
        sum(grepl("Deletion_PTPRD",CNV_uniq)),
        sum(grepl("Deletion_IMMP2L",CNV_uniq)),
        sum(grepl("Deletion_MACROD2",CNV_uniq)),
        sum(grepl("Duplication",CNV_uniq)))
LOH=c(sum(grepl("cnnLOH",CNV_uniq)),
      sum(grepl("chr11q_cnnLOH",CNV_uniq)),
      sum(grepl("chr15q_cnnLOH",CNV_uniq)),
      sum(grepl("chr17q_cnnLOH",CNV_uniq)))
tri=c(sum(grepl("whole_gain",CNV_uniq)),
      sum(grepl("chr20_whole_gain",CNV_uniq)),
      sum(grepl("chr13_whole_gain",CNV_uniq)),
      sum(grepl("chr7_whole_gain",CNV_uniq)),
      sum(grepl("chr18_whole_gain",CNV_uniq)),
      sum(grepl("chr10_whole_gain",CNV_uniq)),
      sum(grepl("chr16_whole_gain",CNV_uniq)))

#Fig. 3a
barplot(c(short,0,0,0,
          LOH,0,0,0,
          tri),
        col=c('dodgerblue3',rep('lightskyblue',5),'firebrick',rep('white',3),
              'black',rep("grey70",3),rep('white',3),
              'darkgreen',rep("olivedrab1",6)))

patients=met_patient$Donor[order(met_patient$Sex,met_patient$Age,decreasing = F)]

#Fig. 3b
CNV_freq_full=CNV_freq_full[c('Trisomy','Arm cnn-LOH','Deletion',"Duplication"),patients[patients%in%colnames(CNV_freq_full)]]
grid = expand.grid(x=rownames(CNV_freq_full), y=colnames(CNV_freq_full))
grid$z = c(CNV_freq_full)
#dev.new(width=6.5, height=6)
library(lattice)
color.palette = colorRampPalette(c("white", "darkorange", "darkorchid2","darkorchid4"))
levelplot(z~y*x, grid, col.regions=color.palette, scales = list(tck = c(0,0), y = list(cex=1), x = list(rot=90)), ylab="Genes", xlab="Samples", colorkey=list(space="bottom"), 
          panel=function(...) { arg <- list(...)
          panel.levelplot(...)
          panel.text(arg$x, arg$y, arg$z)})


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


wgs_data_sub=wgs_data[wgs_data$Feature!="Tumour"&wgs_data$Site!="Unknown",]
has_trisomy=wgs_data_sub$Sample%in%cnvs$sample[grepl("whole_gain",cnvs$CNV)]

fisher.test(table(has_trisomy,wgs_data_sub$Site)) #p-value = 0.3188

wgs_data_sub=wgs_data[wgs_data$Feature!="Tumour",]
has_trisomy=wgs_data_sub$Sample%in%cnvs$sample[grepl("whole_gain",cnvs$CNV)]
fisher.test(table(has_trisomy,wgs_data_sub$Cohort)) #p-value = 0.3441
fisher.test(table(has_trisomy,wgs_data_sub$Hypermut)) #p-value = 0.03

drivers=as.data.frame(read_xlsx("~/Desktop/Gastric/R0/Extended_Data_Table_1_5.xlsx",sheet=5))
has_driver=wgs_data_sub$Sample%in%drivers$Sample
fisher.test(table(has_trisomy,has_driver)) #p-value = 0.03
