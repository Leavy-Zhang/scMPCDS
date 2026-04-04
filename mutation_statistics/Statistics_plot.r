library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(pheatmap)
library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(corrplot)

t293Ddat <- read.csv('path/to/293t_vcf_table.csv',row.names=1, header=TRUE)
helaDdat <- read.csv('path/to/hela_vcf_table.csv',row.names=1, header=TRUE)
helabulkDat <- read.csv('path/to/hela_bulk_vcf_table.csv',row.names=1, header=TRUE)
t293bulkDat <- read.csv('path/to/293t_bulk_vcf_table.csv',row.names=1, header=TRUE)
t293editDat <- read.csv('path/to/293t_edit_vcf_table.csv',row.names=1, header=TRUE)
kidneyDat <- read.csv('path/to/kidney_vcf_table.csv',row.names=1, header=TRUE)


schelacountRegion <- data.frame(celltype = "hela",helaDat %>%
  group_by(Effect, Phenotype, cellID) %>%
  dplyr::summarize(n = n()))
sct293countRegion <- data.frame(celltype = "293T",t293Dat %>%
  group_by(Effect, Phenotype, cellID) %>%
    dplyr::summarize(n = n()))
sccountRegion <- rbind(schelacountRegion, sct293countRegion)

pdf('hela.snp.region.count.pdf', width=5, height =  5)

ggplot(subset(sccountRegion, subset = celltype == "hela"), aes(x = Effect,y = n, color = celltype))+geom_boxplot(fill = "white") + ylab("SNP number")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.snp.region.Freq.pdf', width=5, height =  5)
ggplot(subset(scFREQRegion, subset = celltype == "hela"), aes(x = Effect,y = Mean, color = celltype))+geom_boxplot(fill = "white") + ylab("Frequency (%)")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('sc.snp.region.count.pdf', width=5, height =  5)

ggplot(sccountRegion, aes(x = Effect,y = n, color = celltype))+geom_boxplot(fill = "white") + ylab("SNP number")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('sc.snp.region.Freq.pdf', width=5, height =  5)
ggplot(scFREQRegion, aes(x = Effect,y = Mean, color = celltype))+geom_boxplot(fill = "white") + ylab("Frequency (%)")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('sc.snp.region.count.with_phenotype.pdf',width=14, height =  5)
ggplot(sccountRegion, aes(x = Phenotype,y = n, color = celltype))+geom_boxplot(fill = "white") + facet_grid(~Effect) + ylab("SNP number")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('sc.snp.region.Freq.with_phenotype.pdf',width=14, height =  5)
ggplot(scFREQRegion, aes(x = Phenotype,y = Mean, color = celltype))+geom_boxplot(fill = "white") + facet_grid(~Effect) + ylab("Frequency (%)")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('sc.snp.region.count.with_phenotype.pdf',width=14, height =  5)
ggplot(subset(sccountRegion, subset = celltype == "hela"), aes(x = Phenotype,y = n, color = celltype))+geom_boxplot(fill = "white") + facet_grid(~Effect) + ylab("SNP number")  + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


helaoverlapbulk <- lapply(seq_along(unique(helaDat$Cell)), function(m){
    tmp <- helaDat[helaDat$Cell ==unique(helaDat$Cell)[m], ]
    return(data.frame(celltype = "hela", Cell = unique(helaDat$Cell)[m], percent = sum(tmp$POS %in% helabulkDat$POS)/dim(tmp)[1]*100))}) %>% Reduce('rbind',.)

t293overlapbulk <- lapply(seq_along(unique(t293Dat$Cell)), function(m){
    tmp <- t293Dat[t293Dat$Cell ==unique(t293Dat$Cell)[m], ]
    return(data.frame(celltype = "293T", Cell = unique(t293Dat$Cell)[m], percent = sum(tmp$POS %in% t293bulkDat$POS)/dim(tmp)[1]*100))}) %>% Reduce('rbind',.)
scoverlapbulk <- rbind(helaoverlapbulk,t293overlapbulk)
pdf('sc.in.bulk.pdf')
ggplot(scoverlapbulk,aes(x=celltype,y=percent,color=celltype)) + geom_violin(fill="white")+geom_point()+theme_classic()
dev.off()

helaoverlapbulk2 <- lapply(seq_along(unique(helaDat$Cell)), function(m){
    tmp <- helaDat[helaDat$Cell ==unique(helaDat$Cell)[m], ]
    return(data.frame(celltype = "hela", Cell = unique(helaDat$Cell)[m], percent = sum(tmp$POS %in% helabulkDat$POS)/dim(t293bulkDat)[1]*100))}) %>% Reduce('rbind',.)

t293overlapbulk2 <- lapply(seq_along(unique(t293Dat$Cell)), function(m){
    tmp <- t293Dat[t293Dat$Cell ==unique(t293Dat$Cell)[m], ]
    return(data.frame(celltype = "293T", Cell = unique(t293Dat$Cell)[m], percent = sum(tmp$POS %in% t293bulkDat$POS)/dim(t293bulkDat)[1]*100))}) %>% Reduce('rbind',.)
bulkoverlapsc <- rbind(helaoverlapbulk2, t293overlapbulk2)

pdf('bulk.in.sc.pdf')

ggplot(subset(data.frame(celltype = c("hela","293T"), percent = c(sum(t293bulkDat$POS %in% t293Dat$POS)/dim(t293bulkDat)[1]*100, sum(helabulkDat$POS %in% helaDat$POS)/dim(helabulkDat)[1]*100)),subset = celltype == "hela")) +
       geom_bar(stat = "identity",aes(x = celltype,y=percent, fill = celltype)) + theme_classic()

dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Cell,Transcript_BioType) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Cell,Transcript_BioType) %>% summarise(number = n())))

pdf('sc.gene_type.mut_number.pdf')

ggplot(subset(dat, subset =  Transcript_BioType != ""),aes(x=Transcript_BioType, color = celltype, y= number))+ 
geom_boxplot(aes(x=Transcript_BioType, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()
dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())))

pdf('sc.protein-coding_gene.mut_number.pdf',width = 14,height = 5)

ggplot(dat,aes(x=Gene_Name, color = celltype, y= number))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())))

pdf('sc.tRNA_gene.mut_number.pdf',width = 14,height = 5)
ggplot(dat,aes(x=Gene_Name, color = celltype, y= number))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())))
pdf('sc.rRNA_gene.mut_number.pdf',width = 5,height = 5)

ggplot(dat,aes(x=Gene_Name, color = celltype, y= number))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Cell,Transcript_BioType) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Cell,Transcript_BioType) %>% summarise(number = n())))


pdf('hela.sc.gene_type.mut_number.pdf')

ggplot(subset(dat, subset =  (Transcript_BioType != "") & (celltype =="Hela")),aes(x=Transcript_BioType, color = celltype, y= number))+ 
geom_boxplot(aes(x=Transcript_BioType, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()
dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())))

pdf('hela.sc.protein-coding_gene.mut_number.pdf',width = 7,height = 5)

ggplot(subset(dat,subset = (celltype =="Hela")),aes(x=Gene_Name, color = celltype, y= number))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())))


pdf('sc.tRNA_gene.mut_number.pdf',width = 7,height = 5)
ggplot(subset(dat,subset = (celltype =="Hela")),aes(x=Gene_Name, color = celltype, y= number))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Cell,Gene_Name) %>% summarise(number = n())))

pdf('hela.sc.rRNA_gene.mut_number.pdf',width = 5,height = 5)

ggplot(subset(dat,subset = (celltype =="Hela")),aes(x=Gene_Name, color = celltype, y= number))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= number),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Cell,Transcript_BioType) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Cell,Transcript_BioType) %>% summarise(Freq = mean(FREQ))))

pdf('hela.sc.gene_type.mut_freq.pdf')

ggplot(subset(dat, subset =  (Transcript_BioType != "") & (celltype =="Hela")),aes(x=Transcript_BioType, color = celltype, y= Freq))+ 
geom_boxplot(aes(x=Transcript_BioType, color = celltype, y= Freq),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Cell,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Cell,Gene_Name) %>% summarise(Freq = mean(FREQ))))

pdf('hela.sc.rRNA_gene.mut_freq.pdf',width = 7,height = 5)

ggplot(subset(dat,subset = (celltype =="Hela")),aes(x=Gene_Name, color = celltype, y= Freq))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= Freq),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Cell,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Cell,Gene_Name) %>% summarise(Freq = mean(FREQ))))

pdf('hela.sc.tRNA_gene.mut_freq.pdf',width = 7,height = 5)

ggplot(subset(dat,subset = (celltype =="Hela")),aes(x=Gene_Name, color = celltype,   y= Freq))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= Freq),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Cell,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Cell,Gene_Name) %>% summarise(Freq = mean(FREQ))))

pdf('hela.sc.protein-coding_gene.mut_freq.pdf',width = 7,height = 5)

ggplot(subset(dat,subset = (celltype =="Hela")),aes(x=Gene_Name, color = celltype,   y= Freq))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= Freq),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(cellID,Transcript_BioType) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(cellID,Transcript_BioType) %>% summarise(number = n())))

tail(subset(helaDat,subset = Transcript_BioType == "Mt_tRNA"))

dat <- rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(number = n())))


dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)



dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('mo_sc.gene_type.mut_number.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& grepl('mo',Phenotype)),aes(x=Transcript_BioType, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means( method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()

pdf('ros_sc.gene_type.mut_number.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& !grepl('mo',Phenotype)),aes(x=Transcript_BioType, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means( method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()



dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)



pdf('mo_sc.gene_type.mut_freq.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& grepl('mo',Phenotype)),aes(x=Transcript_BioType, color = Phenotype, y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()

pdf('ros_sc.gene_type.mut_freq.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& !grepl('mo',Phenotype)),aes(x=Transcript_BioType, color = Phenotype, y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.mo_sc.gene_type.mut_number.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& grepl('mo',Phenotype)& (celltype =="Hela")),aes(x=Transcript_BioType, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means( method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()

pdf('hela.ros_sc.gene_type.mut_number.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& !grepl('mo',Phenotype)& (celltype =="Hela")),aes(x=Transcript_BioType, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means( method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()



dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = !grepl('-',t293Dat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = !grepl('-',helaDat$Transcript_BioType)) %>% group_by(Phenotype,cellID,Transcript_BioType) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)



pdf('hela.mo_sc.gene_type.mut_freq.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& grepl('mo',Phenotype)& (celltype =="Hela")),aes(x=Transcript_BioType, color = Phenotype, y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()

pdf('hela.ros_sc.gene_type.mut_freq.pdf',height = 5)

ggplot(subset(dat, subset =  (Transcript_BioType != "")& !grepl('mo',Phenotype)& (celltype =="Hela")),aes(x=Transcript_BioType, color = Phenotype, y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('ros_sc.protein-coding_gene.mut_number.pdf',width = 14,height = 5)


#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('mo_sc.protein-coding_gene.mut_number.pdf',width = 14,height = 5)


#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)


pdf('ros_sc.tRNA_gene.mut_number.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('mo_sc.tRNA_gene.mut_number.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('ros_sc.rRNA_gene.mut_number.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('mo_sc.rRNA_gene.mut_number.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()



dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.ros_sc.protein-coding_gene.mut_number.pdf',width =7,height = 5)


#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('hela.mo_sc.protein-coding_gene.mut_number.pdf',width = 7,height = 5)


#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
#theme_classic()
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)


pdf('hela.ros_sc.tRNA_gene.mut_number.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('hela.mo_sc.tRNA_gene.mut_number.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(number = n())))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.ros_sc.rRNA_gene.mut_number.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('hela.mo_sc.rRNA_gene.mut_number.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype, y= number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= number),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.mo_sc.protein-coding_gene.mut_freq.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.protein-coding_gene.mut_freq.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.mo_sc.rRNA_gene.mut_freq.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.rRNA_gene.mut_freq.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.mo_sc.tRNA_gene.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.tRNA_gene.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)


pdf('hela.mo_sc.protein-coding_gene.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.protein-coding.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()




dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('mo_sc.rRNA_gene.mut_freq.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('ros_sc.rRNA_gene.mut_freq.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('mo_sc.tRNA_gene.mut_freq.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('ros_sc.tRNA_gene.mut_freq.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)


pdf('mo_sc.protein-coding_gene.mut_freq.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('ros_sc.protein-coding.mut_freq.pdf',width = 14,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()




dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_rRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.mo_sc.rRNA_gene.mut_freq.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.rRNA_gene.mut_freq.pdf',width = 5,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'Mt_tRNA') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)

pdf('hela.mo_sc.tRNA_gene.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.tRNA_gene.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


dat = rbind(cbind(celltype = "293T",subset(t293Dat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))),
            cbind(celltype = "Hela",subset(helaDat,subset = Transcript_BioType == 'protein_coding') %>% group_by(Phenotype,cellID,Gene_Name) %>% summarise(Freq = mean(FREQ))))
dat$Phenotype <- gsub('(293t-wt-)|(hela)','',dat$Phenotype)
dat$Phenotype <- gsub('fu','--',dat$Phenotype)


pdf('hela.mo_sc.protein-coding_gene.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
pdf('hela.ros_sc.protein-coding.mut_freq.pdf',width = 7,height = 5)
ggplot(subset(dat, subset = !grepl('mo', Phenotype)& (celltype =="Hela")),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()



ggplot(subset(dat, subset = !grepl('mo', Phenotype)),aes(x=Gene_Name, color = Phenotype,  y= Freq))+ 
geom_boxplot(aes(x=Gene_Name, color = celltype, y= Freq),fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
 scale_color_viridis_d()+
facet_grid(~celltype)+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

t293$mut <- paste(helaDat$REF,helaDat$ALT,sep="->")

t293Dat$mut <- paste(t293Dat$REF,t293Dat$ALT,sep="->")


dat = rbind(cbind(celltype = "293T", t293Dat %>% group_by(Cell,mut) %>% summarise(number = n())))
#            cbind(celltype = "Hela", helaDat %>% group_by(Cell,mut) %>% summarise(number = n())))
pdf('293T.sc.nuc_mut_number.pdf',width = 14,height = 5)

ggplot(dat,aes(x=mut,   y=number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
# scale_color_gradient2()+
#facet_grid(~mut)+
#    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


helaDat$mut <- paste(helaDat$REF,helaDat$ALT,sep="->")
dat = rbind(#cbind(celltype = "293T", helaDat %>% group_by(Cell,mut) %>% summarise(number = n())))
            cbind(celltype = "Hela", helaDat %>% group_by(REF, ALT, Cell,mut) %>% summarise(number = n())))
pdf('Hela.sc.nuc_mut_number.pdf',width = 14,height = 5)

ggplot(dat,aes(x=mut,   y=number))+ 
geom_boxplot(fill="white", outlier.shape = NA)+
 geom_point(position=position_jitterdodge(jitter.width=0.2)) +
# scale_color_gradient2()+
#facet_grid(~mut)+
#    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +

#geom_dotplot(aes(x=Transcript_BioType, color = celltype, y= freq),shape=1,size=1,stroke=1.5)+
theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

helaDat$Mutation <- paste(helaDat$REF,helaDat$ALT,sep="->")
dat <- data.frame(helaDat %>% group_by(REF, ALT, Cell,Mutation) %>% summarise(number = n()))



print(length(unique(dat$Cell)))
helaBaseDat <- lapply(seq_along(unique(dat$Cell)), function(m){
    tmp <- subset(dat, subset = Cell %in% unique(dat$Cell)[m])
#    print(dim(tmp))
    tmp$ratio <- tmp$number/sum(tmp$number)*100
    return(tmp)
}) %>% Reduce('rbind',.)
ggplot(subset(helaBaseDat,subset = (Mutation %in% c("A->G", "G->A", "C->T","T->C")) & (REF != "N")), aes(x= Mutation, y= ratio, fill = "#23FF00"))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")


#helascMutCount <- helaDat %>% group_by(Cell, REF, Gene_Name, Transcript_BioType) %>% summarize(count = n())
#
#helascMutRatio <- lapply(seq_along(unique(helascMutCount$Cell)), function(n){
#    helascMutCountThisCell <- subset(helascMutCount, subset = Cell == unique(helascMutCount$Cell)[n])
#    helascMutRatioThisCell <- lapply(seq_along(unique(helascMutCountThisCell$Transcript_BioType)), function(m){
#        tmpDatCount <- subset(helascMutCountThisCell, subset = Transcript_BioType == unique(helascMutCountThisCell$Transcript_BioType)[m])
#        tmpDatRatio <- lapply(seq_along(unique(tmpDatCount$Gene_Name)), function(i){
#            tmpDat <- subset(tmpDatCount, subset = Gene_Name == unique(tmpDatCount$Gene_Name)[i])
#            REFbase <- c('T','C','G','A')[!(c('T','C','G','A') %in% tmpDat$REF)]
#            if(length(REFbase) > 0){
#
#                tmpDatAppend <- data.frame(Cell= rep(tmpDat$Cell[1], length(REFbase)),
#                                       REF = REFbase,
#                                       Gene_Name = rep(tmpDat$Gene_Name[1], length(REFbase)),
#                                       Transcript_BioType = rep(tmpDat$Transcript_BioType[1], length(REFbase)),
#                                       count = rep(0, length(REFbase))
#                                      )
#                tmpDat <- rbind(tmpDat, tmpDatAppend)
#            }
#            tmpDat$ratio <- tmpDat$count/sum(helascMutCountThisCell$count)*100
#            tmpDat <- rbind(tmpDat, data.frame(
#                Cell= tmpDat$'Cell'[1],
#                REF = "WM",
#                Gene_Name = tmpDat$Gene_Name[1],
#                Transcript_BioType = tmpDat$Transcript_BioType[1],
#                count =sum(tmpDat$count),
#                ratio = sum(tmpDat$count)/sum(helascMutCountThisCell$count)*100
#                )
#            )
#            return(tmpDat)
#            
#        }) %>% Reduce('rbind', .)
#    }) %>% Reduce('rbind', .)
#}) %>% Reduce('rbind', .)
helaBaseDat <- subset(helaBaseDat, subset = !grepl('TRUE',Mutation))
p1 <- ggplot(subset(helascMutRatio,subset = (Transcript_BioType == "protein_coding") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p1.1 <- ggplot(subset(helascMutRatio,subset = (Transcript_BioType == "protein_coding") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") 
p2<- ggplot(subset(helascMutRatio,subset = (Transcript_BioType == "Mt_tRNA") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+ylab('')+ ylim(0,4)+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p3 <- ggplot(subset(helascMutRatio,subset = (Transcript_BioType == "Mt_rRNA") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+ylab('')+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p4 <- ggplot(subset(helaBaseDat,subset = (Mutation %in% c("A->G", "G->A", "C->T","T->C")) & (REF != "N")), aes(x= Mutation, y= ratio, fill = "#23FF00"))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p5 <- ggplot(subset(helaBaseDat,subset = !(Mutation %in% c("A->G", "G->A", "C->T","T->C")) & (REF != "N")), aes(x= Mutation, y= ratio, fill = "#23FF00"))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")

grid.arrange(grobs = list(p1.1,p2,p3,p4,p5), layout_matrix = rbind(c(1, 1, 1, 2),
                        c(3, 4, 5)))


pdf('hela.nucleotide.mutation.ratio.pdf',height=5)
grid.arrange(grobs = list(p4,p5), layout_matrix = matrix(c(1,1,2, 2,2), nrow = 1))
dev.off()


spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))

#myQuantile <- quantile((subset(helaDat, subset = !is.na(mo)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$mo,probs = seq(0,1,by=0.2))
myQuantile <- quantile((subset(helaDat, subset = !is.na(mo)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$mo,probs = round(seq(0,1,length.out = 4),2))

helaDatCatPhenoQuantile <- lapply(seq(1,length(myQuantile)-1), function(m){
    tmp <- subset(helaDat, subset = !is.na(mo) & (mo >= myQuantile[m]) & (mo < myQuantile[m+1]))
    tmpDat <- tmp %>% group_by(Cell, Transcript_BioType) %>% summarise(Number = n(),Freq = median(FREQ))
    tmpDat$Cat <- paste0(names(myQuantile)[m],"-",names(myQuantile)[m+1])
    return(tmpDat)
}) %>% Reduce('rbind',.)
helaDatCatPhenoQuantile <- helaDatCatPhenoQuantile %>% group_by(Transcript_BioType, Cat) %>% summarise(NUMBER = median(Number),FREQ = median(Freq))
p1 <- ggplot(helaDatCatPhenoQuantile,aes(x=Transcript_BioType, y = Cat)) + spot.theme+
#    geom_point(colour = "black", fill = "white", aes(size = 100)) +
#    geom_point(colour = "white", aes(size = 90)) +
    geom_point(aes(size = NUMBER, colour = FREQ))  + ylab("quantile of mo") + scale_color_viridis_c()

#p1<- update_labels(p1, list(size="FREQ"))


myQuantile <- quantile((subset(helaDat, subset = !is.na(ros)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$ros,probs = round(seq(0,1,length.out = 4),2))

helaDatCatPhenoQuantile <- lapply(seq(1,length(myQuantile)-1), function(m){
    tmp <- subset(helaDat, subset = !is.na(ros) & (ros >= myQuantile[m]) & (ros < myQuantile[m+1]))
    tmpDat <- tmp %>% group_by(Cell, Transcript_BioType) %>% summarise(Number = n(),Freq = median(FREQ))
    tmpDat$Cat <- paste0(names(myQuantile)[m],"-",names(myQuantile)[m+1])
    return(tmpDat)
}) %>% Reduce('rbind',.)
helaDatCatPhenoQuantile <- helaDatCatPhenoQuantile %>% group_by(Transcript_BioType, Cat) %>% summarise(NUMBER = median(Number),FREQ = median(Freq))

p2 <- ggplot(helaDatCatPhenoQuantile,aes(x=Transcript_BioType, y = Cat)) + spot.theme+
#    geom_point(colour = "black", fill = "white", aes(size = 100)) +
#    geom_point(colour = "white", aes(size = 90)) +
    geom_point(aes(size = NUMBER, colour = FREQ))  + ylab("quantile of ros") + scale_color_viridis_c()
#    geom_point(aes(colour = NUMBER, size = FREQ*0.9))  + ylab("quantile of ros") + scale_color_viridis_c()
#p2 <- update_labels(p2, list(size="FREQ"))

pdf('hela.wt.region.phenotype.quantile.pdf', width=15,height=7)
do.call("grid.arrange", c(list(p1,p2), ncol = 2))
dev.off()

thisDat = subset(helaDat,subset=FREQ >=1)
myquantile = round(seq(0,100,length.out=4),2)
T293Quantile <- data.frame(t(lapply(seq_along(unique(thisDat$Cell)), function(n){
    tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
    cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
        ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
        return(ctmp)
    }) %>% unlist
    cellQuantile <- c(tmpDat$ros[1], tmpDat$mo[1],cellQuantile)
    cellQuantile <- data.frame(cellQuantile)
    return(cellQuantile)
}) %>% Reduce('cbind',.)))
colnames(T293Quantile) <- c('ros','mo', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
rownames(T293Quantile) <- unique(thisDat$Cell)
head(T293Quantile)

T293Quantile <- subset(T293Quantile, subset = !is.na(T293Quantile$ros) & !is.na(T293Quantile$mo))
pdf('hela.wt.total.phenotype.correlation.pdf', width=5,height=3)
corrplot::corrplot(cor(T293Quantile[,1:2], T293Quantile[,3:dim(T293Quantile)[2]]),col = COL2("RdBu",100)[100:1])
dev.off()


cellNumber <- c(1,5,10,50,100,200,350)
cellMeanCount2 <- lapply(cellNumber,function(m){
    cellSNPCount <- subset(t293Dat,subset=FREQ >=1) %>% group_by(Cell) %>% summarise(n = n())
    thisCellMeanCount <- lapply(seq_len(300), function(N){
        tmp <- cellSNPCount$n[sample(seq(1,dim(cellSNPCount)[1]),size = m)]
        return(c(m, sum(tmp)/m))
    }) %>% Reduce('rbind',.)
        }) %>% Reduce('rbind',.)
colnames(cellMeanCount2) <- c("cell_number","mean_count")
cellMeanCount2 <- data.frame(cellMeanCount2)
cellMeanCount2$cell_number <- factor(cellMeanCount2$cell_number)
mycomp <- list(c("1","5"),
c("1","10"),
c("1","50"),
c("1","100"),
c("1","200"),
c("1","350"))
pdf('293T.random.mean_count.pdf',height = 5)
ggplot(cellMeanCount2,aes(x=cell_number,fill = cell_number,color = cell_number , y=mean_count))+ 
ylab("Mutation sites per each cell") + 
xlab("Cell number") + geom_boxplot(fill="white",outlier.shape = NA) + geom_point(position=position_jitterdodge(jitter.width=1.5))  + 
stat_compare_means(method ="wilcox.test",label = paste0("p = ", after_stat("p.format")),paired = TRUE,comparisons = mycomp) + theme_classic()
dev.off()

thisDat = subset(t293Dat,subset=FREQ >=1)
myquantile = seq(0,100,10)
T293mutQuantile <- lapply(seq_along(unique(thisDat$Cell)), function(n){
    tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
    T293mutratio <- lapply(seq_len(length(myquantile)-1),function(m){
        ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))/dim(tmpDat)[2]
        return(data.frame(VAF = paste0("(",myquantile[m],"-",myquantile[(m+1)],")"), Percent = ctmp, Cell = unique(thisDat$Cell)[n]))
    }) %>% Reduce('rbind',.)
    return(T293mutratio)
}) %>% Reduce('rbind',.)

pdf("293T.sc.percentile.VAF.boxplot.pdf",width = 5,height = 5)
ggplot(T293mutQuantile,aes(x=VAF,y=Percent,color=VAF)) + 
geom_boxplot(fill = "white", outlier.shape = NA) + 
#geom_point(position=position_jitterdodge(jitter.width=0.8))  +
ylab("Percentage of mutation sites") + xlab("VAF Quantiles")+theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
dev.off()



head(T293mutQuantile)

T293SNPinCell <- data.frame(table(subset(thisDat %>% group_by(Cell) %>% summarise(Number=n()), select = Number)))
#T293SNPinCell$Freq <- T293SNPinCell$Freq/sum(T293SNPinCell$Freq)*100

T293SNPinCell <- thisDat %>% group_by(Cell) %>% summarise(Number=n())
head(T293SNPinCell)

pdf("293T.sc.SNP_number.in.Cell.barplot.pdf",width = 5,height = 5)

ggplot(T293SNPinCell,aes(x=Number)) + 
geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "lightblue") +
ylab("Number of cells")+
xlab("Number of mutations in a cell")+
  geom_density(alpha=0.3,size=1.5, fill =  "#21918c") + theme_classic()
dev.off()

t293bulkDat$inSc <- t293bulkDat$POS %in% t293Dat$POS
t293bulkDat$POS <- as.numeric(t293bulkDat$POS)
t293bulkVAF <- subset(t293bulkDat,subset = FREQ >=1) %>% group_by(POS) %>% summarise(VAF = mean(FREQ))
t293bulkVAF$inSc <- lapply(seq_along(t293bulkVAF$POS), function(m){return(t293bulkVAF$POS[m] %in% t293Dat$POS)}) %>% unlist

pdf("293T.bulk.dotplot.pdf",width = 9,height = 3)
ggplot(t293bulkVAF, aes(x = POS, color = inSc, y = VAF)) + 
geom_point()+
scale_color_manual(values=c('red',"black"))+ 
scale_x_continuous(name="Chromosome position", limits=c(1, 17000))+
theme_classic()
dev.off()

head(t293Dat)

thisDat = subset(t293Dat,subset=FREQ >=1)
allPos <- unique(c(t293bulkVAF$POS, thisDat$POS))
t293DatVAFMatrix <- lapply(seq_along(unique(thisDat$Cell)), function(m){
    cDat <- subset(thisDat, subset = Cell == unique(thisDat$Cell)[m])
    cVAF <- lapply(seq_along(allPos), function(n){
        if (allPos[n] %in% cDat$POS)
            return(cDat$FREQ[cDat$POS == allPos[n]])
        else
            return(0)
    }) %>% unlist
    cVAF <- data.frame(cVAF)
    colnames(cVAF) <- unique(thisDat$Cell)[m]
    rownames(cVAF) <- allPos
    return(cVAF)
}) %>% Reduce('cbind',.)
t293DatVAFMatrix<- t(t293DatVAFMatrix[order(as.numeric(rownames(t293DatVAFMatrix))),])

table(unique(colnames(t293DatVAFMatrix)) %in% t293bulkVAF$POS)

table(unique(t293Dat$POS) %in% t293bulkDat$POS)/1221

dev.off()

for(i in seq(1,20))dev.off()

ggplot(t293bulkVAF, aes(x = POS, color = inSc, y = VAF)) + 
geom_point()+
scale_color_manual(values=c('red',"black"))+ 
scale_x_continuous(name="Chromosome position", limits=c(1, 17000))+
theme_classic() + theme(legend.position ="none")

t293DatVAFMatrixMelt$value <- as.numeric(t293DatVAFMatrixMelt$value)
ggplot(t293DatVAFMatrixMelt,aes(x=Var2,y=Var1,fill=value)) + geom_tile() + 
scale_fill_gradientn(colors=viridis(101,begin = 0,end = 1,direction = 1),breaks=seq(0,100,25))+
#scale_fill_distiller(name = "Legend title", palette = "Reds", direction = 1, na.value = "transparent") + 
theme(legend.position = "bottom", legend.direction = "horizontal",)

t293sccountRegion <- subset(t293Dat, subset = (FREQ >=1)) %>% group_by(Cell, Transcript_BioType, Functional_Class) %>% summarise(count = n())
t293sccountRegion <- data.frame(complete(data.frame(t293sccountRegion), Cell, Transcript_BioType, Functional_Class))
t293sccountRegion$count <- lapply(seq_along(t293sccountRegion$count), function(m){
    if(is.na(t293sccountRegion$count[m])){
        return(0)
    }else{
        return(t293sccountRegion$count[m])
    }
    }) %>% unlist
t293bulkcountRegion <- subset(t293bulkDat, subset = (FREQ >=1)) %>% group_by(cellID, Transcript_BioType, Functional_Class) %>% summarise(count = n())
colnames(t293bulkcountRegion) <- colnames(t293sccountRegion)
t293bulkcountRegion <- data.frame(complete(data.frame(t293bulkcountRegion), Cell, Transcript_BioType, Functional_Class))
t293bulkcountRegion$count <- lapply(seq_along(t293bulkcountRegion$count), function(m){
    if(is.na(t293bulkcountRegion$count[m])){
        return(0)
    }else{
        return(t293bulkcountRegion$count[m])
    }
    }) %>% unlist
t293bulkcountRegion$libraryType = "bulk"
t293sccountRegion$libraryType = "sc"
t293sccountRegion$Functional_Class <- factor(t293sccountRegion$Functional_Class, levels = c('MISSENSE','NONSENSE','SILENT','Other'))
t293bulkcountRegion$Functional_Class <- factor(t293bulkcountRegion$Functional_Class, levels = c('MISSENSE','NONSENSE','SILENT','Other'))
pdf('293T.sc-bulk.SNP_count.in_regions.pdf',width = 10,height=5)
ggplot(rbind(t293bulkcountRegion,t293sccountRegion),  aes(Functional_Class,  count)) +
geom_boxplot(aes(fill = libraryType, drop = TRUE), na.rm = FALSE,position = position_dodge(preserve = 'single')) + 
facet_grid(~Transcript_BioType)+ stat_compare_means(method = "wilcox.test", label="p.format") +
theme_classic()+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dev.off()

t293sccountRegion

t293countRegion <- rbind(t293bulkcountRegion,t293sccountRegion)

ggplot(subset(t293countRegion, subset = !(Transcript_BioType %in% c("Dloop",'Mt_rRNA',"Mt_tRNA") & Functional_Class %in% c("MISSENSE","NONSENSE","SILENT"))),  aes(Functional_Class,  count)) +
geom_boxplot(aes(fill = libraryType, drop = TRUE), na.rm = FALSE,position = position_dodge(preserve = 'single')) + 
facet_grid(~Transcript_BioType)+ stat_compare_means(method = "wilcox.test", label="p.format") +
theme_classic()+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

t293scFREQRegion <- subset(t293Dat, subset = (FREQ >=1)) %>% group_by(Cell, Transcript_BioType, Functional_Class) %>% summarise(FREQ = mean(FREQ))
#t293scFREQRegion <- data.frame(complete(data.frame(t293scFREQRegion), Cell, Transcript_BioType, Functional_Class))
t293bulkFREQRegion <- subset(t293bulkDat, subset = (FREQ >=1)) %>% group_by(cellID, Transcript_BioType, Functional_Class) %>% summarise(FREQ = mean(FREQ))
colnames(t293bulkFREQRegion) <- colnames(t293scFREQRegion)
t293scFREQRegion$FREQ <- lapply(seq_along(t293scFREQRegion$FREQ), function(m){
    if(is.na(t293scFREQRegion$FREQ[m])){
        return(0)
    }else{
        return(t293scFREQRegion$FREQ[m])
    }
    }) %>% unlist
t293bulkFREQRegion$FREQ <- lapply(seq_along(t293bulkFREQRegion$FREQ), function(m){
    if(is.na(t293bulkFREQRegion$FREQ[m])){
        return(0)
    }else{
        return(t293bulkFREQRegion$FREQ[m])
    }
    }) %>% unlist
#t293bulkFREQRegion <- data.frame(complete(data.frame(t293bulkFREQRegion), Cell, Transcript_BioType, Functional_Class))
t293bulkFREQRegion$libraryType = "bulk"
t293scFREQRegion$libraryType = "sc"
t293scFREQRegion$Functional_Class <- factor(t293scFREQRegion$Functional_Class, levels = c('MISSENSE','NONSENSE','SILENT','Other'))
t293bulkFREQRegion$Functional_Class <- factor(t293bulkFREQRegion$Functional_Class, levels = c('MISSENSE','NONSENSE','SILENT','Other'))
pdf('293T.sc-bulk.SNP_FREQ.in_regions.pdf',width = 10,height=5)
ggplot(rbind(t293bulkFREQRegion,t293scFREQRegion),  aes(Functional_Class,  FREQ)) +
geom_boxplot(aes(fill = libraryType, drop = TRUE), na.rm = FALSE,position = position_dodge(preserve = 'single')) + 
facet_grid(~Transcript_BioType)+ stat_compare_means(method = "wilcox.test", label="p.format") +
theme_classic()+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

dev.off()

wilcox.test(subset(t293bulkcountRegion,subset = (Functional_Class == "Other") & (Transcript_BioType ==   "Mt_tRNA"))$count,
            subset(t293sccountRegion,subset = (Functional_Class ==   "Other")   & (Transcript_BioType == "Mt_tRNA"))$count)

helaBaseDat <- helaDat %>% group_by(POS, REF, Mutation) %>% summarise(n=n())

502/388



t293scMutCount <- t293Dat %>% group_by(Cell, REF, Gene_Name, Transcript_BioType) %>% summarize(count = n())

t293scMutRatio <- lapply(seq_along(unique(t293scMutCount$Cell)), function(n){
    t293scMutCountThisCell <- subset(t293scMutCount, subset = Cell == unique(t293scMutCount$Cell)[n])
    t293scMutRatioThisCell <- lapply(seq_along(unique(t293scMutCountThisCell$Transcript_BioType)), function(m){
        tmpDatCount <- subset(t293scMutCountThisCell, subset = Transcript_BioType == unique(t293scMutCountThisCell$Transcript_BioType)[m])
        tmpDatRatio <- lapply(seq_along(unique(tmpDatCount$Gene_Name)), function(i){
            tmpDat <- subset(tmpDatCount, subset = Gene_Name == unique(tmpDatCount$Gene_Name)[i])
            REFbase <- c('T','C','G','A')[!(c('T','C','G','A') %in% tmpDat$REF)]
            if(length(REFbase) > 0){

                tmpDatAppend <- data.frame(Cell= rep(tmpDat$Cell[1], length(REFbase)),
                                       REF = REFbase,
                                       Gene_Name = rep(tmpDat$Gene_Name[1], length(REFbase)),
                                       Transcript_BioType = rep(tmpDat$Transcript_BioType[1], length(REFbase)),
                                       count = rep(0, length(REFbase))
                                      )
                tmpDat <- rbind(tmpDat, tmpDatAppend)
            }
            tmpDat$ratio <- tmpDat$count/sum(t293scMutCountThisCell$count)*100
            tmpDat <- rbind(tmpDat, data.frame(
                Cell= tmpDat$'Cell'[1],
                REF = "WM",
                Gene_Name = tmpDat$Gene_Name[1],
                Transcript_BioType = tmpDat$Transcript_BioType[1],
                count =sum(tmpDat$count),
                ratio = sum(tmpDat$count)/sum(t293scMutCountThisCell$count)*100
                )
            )
            return(tmpDat)
            
        }) %>% Reduce('rbind', .)
    }) %>% Reduce('rbind', .)
}) %>% Reduce('rbind', .)



p1 <- ggplot(subset(t293scMutRatio,subset = (Transcript_BioType == "protein_coding") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p1.1 <- ggplot(subset(t293scMutRatio,subset = (Transcript_BioType == "protein_coding") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") 
p2<- ggplot(subset(t293scMutRatio,subset = (Transcript_BioType == "Mt_tRNA") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+ylab('')+ ylim(0,4)+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p3 <- ggplot(subset(t293scMutRatio,subset = (Transcript_BioType == "Mt_rRNA") & (REF != "N")), aes(x= Gene_Name, y= ratio, color = REF))+ylab('')+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p4 <- ggplot(subset(t293BaseDat,subset = (Mutation %in% c("A->G", "G->A", "C->T","T->C")) & (REF != "N")), aes(x= Mutation, y= ratio, fill = "#23FF00"))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p5 <- ggplot(subset(t293BaseDat,subset = !(Mutation %in% c("A->G", "G->A", "C->T","T->C")) & (REF != "N")), aes(x= Mutation, y= ratio, fill = "#23FF00"))+
geom_boxplot() +theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")

grid.arrange(grobs = list(p1.1,p2,p3,p4,p5), layout_matrix = rbind(c(1, 1, 1, 2),
                        c(3, 4, 5)))


pdf('293T.sc.genes.nuc.mutRatio.pdf',width = 20,height = 7)
grid.arrange(grobs = list(p1,p2,p3,p4,p5), widths = c(1,1.2,0.3), layout_matrix = rbind(c(1, 2,3),
                        c(4,5, 5)))
dev.off()

grid.arrange(grobs = list(p1,p2,p3, p4, p5), 
             widths = c(10,4,4,6),
             layout_matrix = rbind(c(1, 1, 1,1),
                        c(2, 3, 4, 5)))


t293DatFREQDist <- lapply(seq_along(unique(t293Dat$Cell)), function(m){
    thisDat <- subset(t293Dat, subset = Cell == unique(t293Dat$Cell)[m])
    return(lapply(seq(1,19), function(n){
        tmpDat <- subset(thisDat, subset = (FREQ >= n) & (FREQ <= (n+1))) %>% summarise(count = n())
        tmpDat$Start <- n
        tmpDat$Cell <- unique(t293Dat$Cell)[m]
        return(c(Cell = unique(t293Dat$Cell)[m], Start = n , count = tmpDat$count))
    }) %>% Reduce('rbind',.))  
})  %>% Reduce('rbind',.)

tpos <- colnames(t293DatVAFMatrix)

tVAF <- lapply(seq_along(tpos), function(m){
    tmp <- subset(t293bulkVAF, subset = POS == tpos[m])
    if(dim(tmp)[1] == 1){
        return(tmp$VAF)
    }else{
        return(NA)
    }
}) %>% unlist

tOverlap <- lapply(seq_along(tpos), function(m){
    tmp <- subset(t293bulkVAF, subset = POS == tpos[m])
    if(dim(tmp)[1] == 1){
        return(tmp$inSc)
    }else{
        return(FALSE)
    }
}) %>% unlist

table(colnames(t293DatVAFMatrix) %in% t293bulkDat$POS)

# Actually it can be made easily by `anno_points()`.
library(circlize)
library(ComplexHeatmap)
col_runif = colorRamp2(c(0,60, 80, 100), c("white","purple","purple", "purple"))
ha = HeatmapAnnotation(foo = anno_empty(border = TRUE, height = unit(3, "cm")))
sc = HeatmapAnnotation(ifOverlap = ifelse(colnames(t293DatVAFMatrix) %in% t293bulkDat$POS,"YES","NO"), 
                       col = list(ifOverlap = c( "YES"="red", "NO"="black")))
ht = Heatmap(t293DatVAFMatrix, name = "VAF",show_row_dend = FALSE, col = col_runif,
             border_gp = gpar(col = "black"), cluster_columns = FALSE, 
             show_row_names = FALSE, show_column_names = FALSE, top_annotation = ha,
            bottom_annotation = sc)
pdf("/293T.bulk.dotplot-heatmap.pdf")
ht = draw(ht)
co = column_order(ht)
value = runif(10)
decorate_annotation("foo", {
    # value on x-axis is always 1:ncol(mat)
    #x = 1:10
    # while values on y-axis is the value after column reordering
    #value = value[co]
#    pushViewport(viewport(xscale = range(as.numeric(colnames(t293DatVAFMatrix))), yscale = c(0, 110)))
    pushViewport(viewport(xscale = c(1,1232), yscale = c(-10, 110)))
    #grid.lines(c(0.5, 10.5), c(0.5, 0.5), gp = gpar(lty = 2),
    #    default.units = "native")
    grid.points(seq(1,1232), tVAF, pch = 16, size = unit(2, "mm"),
        gp = gpar(col = ifelse(tOverlap, "black", "red")), 
                default.units = "native")
    grid.yaxis(at = c(0, 50, 100))
    popViewport()
})
dev.off()

t293DatExample <- subset(t293Dat,POS %in% c(15739,15740))

t293editDatExample <- subset(t293editDat,POS %in% c(15739,15740))
t293DatExample$Genotype <- "WT"
t293editDatExample$Genotype <- "edit"
t293Example <- rbind(t293DatExample[,c(1:18,21)],t293editDatExample)
t293Example$Genotype <- factor(t293Example$Genotype,  levels = c('WT',"edit"))
ggplot(subset(t293Example, subset = !grepl('edit1',Cell)),aes(x = Genotype, y = FREQ,color = Genotype)) + facet_grid("~POS") +
geom_beeswarm(cex = 2.5, corral = "wrap")+ stat_compare_means(method = "wilcox.test", label="p.format", label.x = 1.3, label.y = 40) + scale_color_manual(values = RColorBrewer::brewer.pal(n = 2,name = "Dark2"))+ theme_classic()+theme(legend.position = "none")


t293DatExample <- subset(t293Dat,POS %in% c(15739,15740))
t293editDatExample <- subset(t293editDat,POS %in% c(15739,15740))
t293DatExample$Genotype <- "WT"
t293editDatExample$Genotype <- "edit"
t293Example <- rbind(t293DatExample[,c(1:18,21)],t293editDatExample)
t293Example$Genotype <- factor(t293Example$Genotype,  levels = c('WT',"edit"))
ggplot(subset(t293Example, subset = !grepl('edit1',Cell)),aes(x = Genotype, y = FREQ,color = Genotype)) + facet_grid("~POS") +
geom_beeswarm(cex = 2.5, corral = "wrap")+ stat_compare_means(method = "wilcox.test", label="p.format", label.x = 1.3, label.y = 40) + scale_color_manual(values = RColorBrewer::brewer.pal(n = 2,name = "Dark2"))+ theme_classic()+theme(legend.position = "none")

pdf('293T.WT-edit.Example.pdf',height=5)
p1 <- ggplot(subset(t293Example, subset = !grepl('edit1',Cell)),aes(x = Genotype, y = FREQ,color = Genotype)) + 
        facet_grid("~POS") +
        geom_beeswarm(cex = 2.5, corral = "wrap")+ 
        stat_compare_means(method = "wilcox.test", label="p.format", label.x = 1.3, label.y = 40) + 
        scale_color_manual(values = RColorBrewer::brewer.pal(n = 2,name = "Dark2"))+ 
        theme_classic()+theme(legend.position = "none")

t293Stat1 <- t293Example %>% group_by(Genotype, POS) %>% summarise(percent = n())
t293Stat1$percent <- c(t293Stat1$percent[1:2]/length(unique(t293Dat$Cell))*100, t293Stat1$percent[3:4]/length(unique(t293editDat2$Cell))*100)
t293Stat2<- t293Stat1
t293Stat2$percent <- 100 - t293Stat1$percent
t293Stat2$ifDetected <- "NO"
t293Stat1$ifDetected <- "YES"
t293Stat <- rbind(t293Stat1,t293Stat2)
t293Stat$ifDetected <- factor(t293Stat$ifDetected, levels = c("NO","YES"))
p2 <- ggplot(t293Stat, aes(x = Genotype, y = percent, fill = forcats::fct_rev(ifDetected))) + 
geom_bar(stat = "identity",position = "fill")+ 
scale_fill_manual(values = c("darkgreen","gray"))+ facet_grid("~POS") + theme_classic() + ylab("Percentage of mutation detected cell")+guides(fill=guide_legend(title="ifDetected"))
do.call(grid.arrange, c(list(p1,p2), ncol = 2))
dev.off()

t293bulkeditDat$FREQ <- as.numeric(t293bulkeditDat$FREQ)
t293bulkDatExample <- subset(t293bulkDat,POS %in% c("15739","15740"))
t293bulkeditDatExample <- subset(t293bulkeditDat,POS %in% c("15739","15740"))

t293bulkDatExample$Genotype <- "WT"
t293bulkeditDatExample$Genotype <- "edit"
t293bulkeditDat
t293bulkExample <- rbind(t293bulkDatExample,t293bulkeditDatExample)
#t293bulkExample$Genotype <- factor(t293bulkExample$Genotype,  levels = c('WT',"edit"))

pdf('293T.bulk.edit.Example.pdf',width = 4, height=5)
ggplot(t293bulkeditDatExample,aes(x = POS, y = FREQ, fill = "skyblue")) + 
        geom_boxplot()+ theme_classic()
dev.off()

t293bulkeditDatExample

t293Stat <- rbind(t293Stat1,t293Stat2)


for(i in seq(1,20))dev.off()

t293DatSNPdetectedCell <- data.frame(sort(table((t293Dat %>% group_by(POS) %>% summarise(detected_cell = n()))$detected_cell)))#[174:155,]
t293DatSNPdetectedCell<-t293DatSNPdetectedCell[order(t293DatSNPdetectedCell$Var1,decreasing = TRUE),][1:20,]
t293DatSNPdetectedCell$Var1 <- factor(t293DatSNPdetectedCell$Var1, levels = c(1, 2, 3, 4, 5, 6, 8, 7, 11, 10, 9, 373, 20, 376, 18, 14, 13, 36, 32, 26))
#t293DatSNPdetectedCell$Var1 <- as.numeric(t293DatSNPdetectedCell$Var1)
pdf('293T.WT-edit.Example.pdf',width=5, height=5)
ggplot(t293DatSNPdetectedCell, aes(x=Var1,y=Freq)) + 
geom_bar(stat = "identity")+ ylab("Number of mutation sites") + xlab("Cell number of detected mutation sites") +theme_classic()
dev.off()

head(helaDat$Mutation)

head(t293bulkeditDat)

t293bulkCircos <- t293bulkDat %>% group_by(POS) %>% summarise(FREQ = mean(FREQ))
t293Circos <- t293Dat %>% group_by(POS) %>% summarise(FREQ = mean(FREQ))

t293bulkCircos$FREQ <- -t293bulkCircos$FREQ

write.csv(rbind(t293bulkCircos,t293Circos),"293t.freq.for.circos.csv")



head(t293editDat)


#ggscatter(t293Dat %>% group_by(Cell, ros, mo) %>% summarise(Number = n()),x="ros",y = "Number", add = "reg.line") +
#  stat_cor(label.x = 3, label.y = 36) +
#  stat_regline_equation(label.x = 3, label.y = 32)


ggplot(data = t293Dat %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)), aes(x = ros, y = Freq)) +
        geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        geom_point()+ ylim(0,50)+
        stat_cor(label.y = 10)+ #this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 5) + theme_classic()#this means at 30th unit regresion line equation will be shown


thisDat = subset(t293Dat,subset=FREQ >=1)
myquantile = seq(0,100,10)
T293mutQuantile <- lapply(seq_along(unique(thisDat$Cell)), function(n){
    tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
    T293mutratio <- lapply(seq_len(length(myquantile)-1),function(m){
        ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))/dim(tmpDat)[2]
        return(data.frame(VAF = paste0("(",myquantile[m],"-",myquantile[(m+1)],")"), Percent = ctmp, Cell = unique(thisDat$Cell)[n]))
    }) %>% Reduce('rbind',.)
    return(T293mutratio)
}) %>% Reduce('rbind',.)

pdf("293T.sc.percentile.VAF.boxplot.pdf",width = 5,height = 5)
ggplot(T293mutQuantile,aes(x=VAF,y=Percent,color=VAF)) + 
geom_boxplot(fill = "white", outlier.shape = NA) + 
#geom_point(position=position_jitterdodge(jitter.width=0.8))  +
ylab("Percentage of mutation sites") + xlab("VAF Quantiles")+theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
dev.off()

myQuantile <- quantile((subset(t293Dat, subset = !is.na(ros)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$ros,probs = seq(0,1,by=0.2))

t293DatPhenoQuantile1 <- lapply(seq(1,length(myQuantile)-1), function(m){
    tmp <- subset(t293Dat, subset = !is.na(ros) & (ros >= myQuantile[m]) & (ros < myQuantile[m+1]))
    tmpDat <- tmp %>% group_by(Cell) %>% summarise(Number = n())
    tmpDat$Cat <- paste0(names(myQuantile)[m],"-",names(myQuantile)[m+1])
    return(tmpDat)
}) %>% Reduce('rbind',.)
p1 <- ggplot(t293DatPhenoQuantile1, aes(x = Cat, y = Number,color = Cat)) + 
        geom_boxplot(fill="white") + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

myQuantile <- quantile((subset(t293Dat, subset = !is.na(ros)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$ros,probs = seq(0,1,by=0.2))

t293DatPhenoQuantile2 <- lapply(seq(1,length(myQuantile)-1), function(m){
    tmp <- subset(t293Dat, subset = !is.na(ros) & (ros >= myQuantile[m]) & (ros < myQuantile[m+1]))
    tmpDat <- tmp %>% group_by(Cell) %>% summarise(Freq = mean(FREQ))
    tmpDat$Cat <- paste0(names(myQuantile)[m],"-",names(myQuantile)[m+1])
    return(tmpDat)
}) %>% Reduce('rbind',.)
p2 <- ggplot(t293DatPhenoQuantile2, aes(x = Cat, y = Freq,color = Cat)) + 
        geom_boxplot(fill="white") + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pdf('293T.wt.ros.quantile.pdf', width=7,height=5)
do.call("grid.arrange", c(list(p1,p2), ncol = 2))
dev.off()

helaDat$Phenotype <- gsub('hela','',helaDat$Phenotype)
unique(helaDat$Phenotype)


spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))

#myQuantile <- quantile((subset(t293Dat, subset = !is.na(mo)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$mo,probs = seq(0,1,by=0.2))
myQuantile <- quantile((subset(t293Dat, subset = !is.na(mo)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$mo,probs = round(seq(0,1,length.out = 4),2))

t293DatCatPhenoQuantile <- lapply(seq(1,length(myQuantile)-1), function(m){
    tmp <- subset(t293Dat, subset = !is.na(mo) & (mo >= myQuantile[m]) & (mo < myQuantile[m+1]))
    tmpDat <- tmp %>% group_by(Cell, Transcript_BioType) %>% summarise(Number = n(),Freq = median(FREQ))
    tmpDat$Cat <- paste0(names(myQuantile)[m],"-",names(myQuantile)[m+1])
    return(tmpDat)
}) %>% Reduce('rbind',.)
t293DatCatPhenoQuantile <- t293DatCatPhenoQuantile %>% group_by(Transcript_BioType, Cat) %>% summarise(NUMBER = median(Number),FREQ = median(Freq))
p1 <- ggplot(t293DatCatPhenoQuantile,aes(x=Transcript_BioType, y = Cat)) + spot.theme+
#    geom_point(colour = "black", fill = "white", aes(size = 100)) +
#    geom_point(colour = "white", aes(size = 90)) +
    geom_point(aes(size = NUMBER, colour = FREQ))  + ylab("quantile of mo") + scale_color_viridis_c()

#p1<- update_labels(p1, list(size="FREQ"))


myQuantile <- quantile((subset(t293Dat, subset = !is.na(ros)) %>% group_by(Cell, ros, mo) %>% summarise(Freq = mean(FREQ)))$ros,probs = round(seq(0,1,length.out = 4),2))

t293DatCatPhenoQuantile <- lapply(seq(1,length(myQuantile)-1), function(m){
    tmp <- subset(t293Dat, subset = !is.na(ros) & (ros >= myQuantile[m]) & (ros < myQuantile[m+1]))
    tmpDat <- tmp %>% group_by(Cell, Transcript_BioType) %>% summarise(Number = n(),Freq = median(FREQ))
    tmpDat$Cat <- paste0(names(myQuantile)[m],"-",names(myQuantile)[m+1])
    return(tmpDat)
}) %>% Reduce('rbind',.)
t293DatCatPhenoQuantile <- t293DatCatPhenoQuantile %>% group_by(Transcript_BioType, Cat) %>% summarise(NUMBER = median(Number),FREQ = median(Freq))

p2 <- ggplot(t293DatCatPhenoQuantile,aes(x=Transcript_BioType, y = Cat)) + spot.theme+
#    geom_point(colour = "black", fill = "white", aes(size = 100)) +
#    geom_point(colour = "white", aes(size = 90)) +
    geom_point(aes(size = NUMBER, colour = FREQ))  + ylab("quantile of ros") + scale_color_viridis_c()
#    geom_point(aes(colour = NUMBER, size = FREQ*0.9))  + ylab("quantile of ros") + scale_color_viridis_c()
#p2 <- update_labels(p2, list(size="FREQ"))

pdf('293T.wt.region.phenotype.quantile.pdf', width=15,height=7)
do.call("grid.arrange", c(list(p1,p2), ncol = 2))
dev.off()

thisDat = subset(t293Dat,subset=FREQ >=1)
myquantile = round(seq(0,100,length.out=4),2)
T293Quantile <- data.frame(t(lapply(seq_along(unique(thisDat$Cell)), function(n){
    tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
    cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
        ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
        return(ctmp)
    }) %>% unlist
    cellQuantile <- c(tmpDat$ros[1], tmpDat$mo[1],cellQuantile)
    cellQuantile <- data.frame(cellQuantile)
    return(cellQuantile)
}) %>% Reduce('cbind',.)))
colnames(T293Quantile) <- c('ros','mo', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
rownames(T293Quantile) <- unique(thisDat$Cell)
head(T293Quantile)

T293Quantile <- subset(T293Quantile, subset = !is.na(T293Quantile$ros) & !is.na(T293Quantile$mo))
pdf('293T.wt.total.phenotype.correlation.pdf', width=5,height=3)
corrplot::corrplot(cor(T293Quantile[,1:2], T293Quantile[,3:dim(T293Quantile)[2]]),col = COL2("RdBu",100)[100:1])
dev.off()


myquantile = round(seq(0,100,length.out = 4),2)
T293RegionQuantile <- lapply(seq_along(unique(t293Dat$Transcript_BioType)), function(j){
    thisDat = subset(t293Dat,subset=Transcript_BioType == unique(t293Dat$Transcript_BioType)[j])
    print(dim(thisDat))
    T293Quantile <- data.frame(t(lapply(seq_along(unique(thisDat$Cell)), function(n){
        tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
        cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
            ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
            return(ctmp)
        }) %>% unlist
        cellQuantile <- c(tmpDat$ros[1], tmpDat$mo[1],cellQuantile)
        cellQuantile <- data.frame(cellQuantile)
        return(cellQuantile)
    }) %>% Reduce('cbind',.)))
    colnames(T293Quantile) <- c('ros','mo', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
    rownames(T293Quantile) <- unique(thisDat$Cell)
    T293Quantile <- subset(T293Quantile, subset = !is.na(T293Quantile$ros) & !is.na(T293Quantile$mo))
    print(dim(T293Quantile))
    return(T293Quantile)
})


round(seq(0,100,length.out = 4),2)

lapply(seq_along(unique(t293Dat$Transcript_BioType)),function(j){
    pdf(paste0('293T.wt.',unique(t293Dat$Transcript_BioType)[j],'.phenotype.correlation.pdf'), width=5,height=3)
        corrplot::corrplot(cor(T293RegionQuantile[[j]][,1:2], T293RegionQuantile[[j]][,3:dim(T293RegionQuantile[[j]])[2]]),mar = c(0,3,0,3),tl.col = 'black',col = COL2("RdBu",100)[100:1])
        title(ylab = unique(t293Dat$Transcript_BioType)[j])
    dev.off()
})

cor.test(T293Quantile[,2],T293Quantile[,8],method = "spearman")

head(T293Quantile)

tmpDat  =  data.frame(T293Quantile[,c(2,dim(T293Quantile)[2])])
colnames(tmpDat) <- c("Phenotype", "Number")
pdf("293T.wt.top30.mo.phenotype.correlation.lm.pdf", width = 5, height = 5)
ggplot(data = tmpDat , aes(x = Number, y = Phenotype)) +
        geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        geom_point()+
        stat_cor(label.y = 4800)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 5000) + theme_classic()#this means at 30th unit regresion line equation will be shown
dev.off()
tmpDat  =  data.frame(T293Quantile[,c(1,dim(T293Quantile)[2])])
colnames(tmpDat) <- c("Phenotype", "Number")
pdf("293T.wt.top30.ros.phenotype.correlation.lm.pdf", width = 5, height = 5)
ggplot(data = tmpDat , aes(x = Number, y = Phenotype)) +
        geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        geom_point()+
        stat_cor(label.y = 950)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 1000) + theme_classic()#this means at 30th unit regresion line equation will be shown
dev.off()

lapply(seq_along(unique(t293Dat$Transcript_BioType)),function(j){
    tmpDat  =  data.frame(T293RegionQuantile[[j]][,c(2,dim(T293RegionQuantile[[j]])[2])])
    colnames(tmpDat) <- c("Phenotype", "Number")
    pdf(paste0('293T.wt.top30.ros.',unique(t293Dat$Transcript_BioType)[j],'.phenotype.correlation.lm.pdf'), width=5,height=5)
        ggplot(data = tmpDat , aes(x = Number, y = Phenotype)) +
        geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
        geom_point()+
        stat_cor(label.y = 4800)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 5000) + theme_classic()#this means at 30th unit regresion line equation will be shown
    dev.off()
})

pdf('293T.wt.ros-mo.correlation.lm.pdf',height = 5, width = 5.5)
thisDat <- lapply(seq_along(unique(t293Dat$Cell)), function(m){
    tmpDat <- subset(t293Dat, subset = Cell == unique(t293Dat$Cell)[m])
    tmpDat <- cbind(Number = dim(tmpDat)[1], FREQ = mean(tmpDat$FREQ),head(subset(tmpDat, select =  c(ros,mo)),1))
    return(tmpDat)
}) %>% Reduce('rbind',.)
head(thisDat$ros)
ggplot(thisDat, aes(x = ros, y=mo)) +
geom_point(aes(size = Number,color = FREQ, alpha =0.5)) + 
geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
scale_color_viridis_c(begin = 0.2,end = 1)+
scale_radius()+
stat_cor(label.y = 4800)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 5000) +
theme_classic()#this means at 30th unit regresion line equation will be shown
dev.off()

pdf('hela.wt.ros-mo.correlation.lm.pdf',height = 5, width = 5.5)
thisDat <- lapply(seq_along(unique(helaDat$Cell)), function(m){
    tmpDat <- subset(helaDat, subset = Cell == unique(helaDat$Cell)[m])
    tmpDat <- cbind(Number = dim(tmpDat)[1], FREQ = mean(tmpDat$FREQ),head(subset(tmpDat, select =  c(ros,mo)),1))
    return(tmpDat)
}) %>% Reduce('rbind',.)
head(thisDat$ros)
ggplot(thisDat, aes(x = ros, y=mo)) +
geom_point(aes(size = Number,color = FREQ, alpha =0.5)) + 
geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
scale_color_viridis_c(begin = 0.2,end = 1)+
scale_radius()+
stat_cor(label.y = 4800)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 5000) +
theme_classic()#this means at 30th unit regresion line equation will be shown
dev.off()

myquantile = round(seq(0,100,length.out = 4),2)
T293RegionQuantile <- lapply(c(2,3), function(j){
    thisDat = subset(t293Dat,subset=Transcript_BioType == unique(t293Dat$Transcript_BioType)[j])
    T293Quantile <- data.frame(t(lapply(seq_along(unique(thisDat$Cell)), function(n){
        tmpDat1 <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
        genes <- unique(tmpDat1$Gene_Name)
        geneQuantile <- lapply(seq_along(genes), function(i){
            tmpDat <- subset(tmpDat1, subset = Gene_Name == genes[i])
            cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
                ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
                return(ctmp)
            }) %>% unlist
            cellQuantile <- c(tmpDat$Cell[1], genes[i], tmpDat$ros[1], tmpDat$mo[1],cellQuantile)
            cellQuantile <- data.frame(cellQuantile)
            return(cellQuantile)
        }) %>% Reduce('cbind',.)
        return(geneQuantile)
    }) %>% Reduce('cbind',.)))
    colnames(T293Quantile) <- c("Cell","Gene_Name",'ros','mo', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
    T293Quantile <- subset(T293Quantile, subset = !is.na(T293Quantile$ros) & !is.na(T293Quantile$mo))
    T293Quantile[,3] <- as.numeric(T293Quantile[,3])
    T293Quantile[,4] <- as.numeric(T293Quantile[,4])
    T293Quantile[,5] <- as.numeric(T293Quantile[,5])
    T293Quantile[,6] <- as.numeric(T293Quantile[,6])
    T293Quantile[,7] <- as.numeric(T293Quantile[,7])
    return(T293Quantile)
})


head(T293RegionQuantile[[1]])

geneClass <- c('rRNA',"PC")
lapply(seq_along(geneClass), function(m){
    tmpDat1  =  T293RegionQuantile[[m]]
    genes <- unique(tmpDat1$Gene_Name)
    lapply(seq_along(genes), function(n){
        tmpDat = subset(tmpDat1, subset = Gene_Name == genes[n])[,c(3,7)]
        colnames(tmpDat) <- c("Phenotype", "Number")
        pdf(paste0('293T.wt.top30.ros.',geneClass[m],'.', genes[n],'.phenotype.correlation.lm.pdf'), width=5,height=5)
            print(ggplot(data = tmpDat , aes(x = Number, y = Phenotype)) +
            geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
            geom_point()+
            stat_cor(label.y = 4800)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
            stat_regline_equation(label.y = 5000) + theme_classic()
                  )#this means at 30th unit regresion line equation will be shown
        dev.off()
        
    })
                           })


#thisDat = subset(t293Dat,subset=FREQ >=1)
myquantile = round(seq(0,100,length.out = 4),2)
T293GroupedQuantile <- lapply(seq_along(unique(t293Dat$Phenotype)), function(i){
    thisDat <- subset(t293Dat, subset = Phenotype == unique(t293Dat$Phenotype)[i])
    T293Quantile <- data.frame(t(lapply(seq_along(unique(thisDat$Cell)), function(n){
        tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
        cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
            ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
            return(ctmp)
        }) %>% unlist
        cellQuantile <- c(tmpDat$Phenotype[1],cellQuantile)
        cellQuantile <- data.frame(cellQuantile)
        return(cellQuantile)
    }) %>% Reduce('cbind',.)))
#    print(T293Quantile)
    colnames(T293Quantile) <- c('Phenotype', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
    rownames(T293Quantile) <- unique(thisDat$Cell)
#    T293Quantile <- subset(T293Quantile, subset = !is.na(T293Quantile$ros) & !is.na(T293Quantile$mo))
    return(T293Quantile)
    }) %>% Reduce('rbind',.)


head(T293GroupedQuantile)

pList1 <- lapply(seq(2,4), function(m){
    tmpDat <- T293GroupedQuantile[,c(1,m)]
    colnames(tmpDat) <- c("Phenotype","Number")
    tmpDat$Number <- as.numeric(tmpDat$Number)
    tmpDat$Phenotype <- gsub('fu',"--", tmpDat$Phenotype)
    tmpDat <- subset(tmpDat, subset = grepl('mo',Phenotype))
    ggplot(tmpDat, aes(x = Phenotype, y = Number,color = Phenotype)) + ggtitle(colnames(T293GroupedQuantile)[m])+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
    geom_boxplot() + scale_color_viridis_d()+ theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

})
pList2 <- lapply(seq(2,4), function(m){
    tmpDat <- T293GroupedQuantile[,c(1,m)]
    colnames(tmpDat) <- c("Phenotype","Number")
    tmpDat$Number <- as.numeric(tmpDat$Number)
    tmpDat$Phenotype <- gsub('fu',"--", tmpDat$Phenotype)
    tmpDat <- subset(tmpDat, subset = !grepl('mo',Phenotype))
    ggplot(tmpDat, aes(x = Phenotype, y = Number,color = Phenotype)) + ggtitle(colnames(T293GroupedQuantile)[m])+
    stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
    geom_boxplot() + scale_color_viridis_d()+ theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

})
pdf('293T.WT.mo_grouped_quantile.number.pdf',width = 20,height=5)
do.call("grid.arrange", c(pList1, nrow = 1))
dev.off()
pdf('293T.WT.ros_grouped_quantile.number.pdf',width = 20,height=5)
do.call("grid.arrange", c(pList2, nrow = 1))
dev.off()

#thisDat = subset(t293Dat,subset=FREQ >=1)
myquantile = round(seq(0,100,length.out = 4),2)
T293GroupedRegionQuantile <- lapply(seq(1,4), function(j){
    that293Dat <- subset(t293Dat, subset = Transcript_BioType == unique(t293Dat$Transcript_BioType)[j])
    regionQuantile <- lapply(seq_along(unique(t293Dat$Phenotype)), function(i){
    thisDat <- subset(that293Dat, subset = Phenotype == unique(t293Dat$Phenotype)[i])
    T293Quantile <- data.frame(t(lapply(seq_along(unique(thisDat$Cell)), function(n){
        tmpDat <- subset(thisDat, subset = Cell %in% unique(thisDat$Cell)[n])
        cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
            ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
            return(ctmp)
        }) %>% unlist
        cellQuantile <- c(unique(t293Dat$Transcript_BioType)[j], tmpDat$Phenotype[1],cellQuantile)
        cellQuantile <- data.frame(cellQuantile)
        return(cellQuantile)
    }) %>% Reduce('cbind',.)))
#    print(T293Quantile)
    colnames(T293Quantile) <- c("Transcript_Biotype",'Phenotype', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
    rownames(T293Quantile) <- unique(thisDat$Cell)
#    T293Quantile <- subset(T293Quantile, subset = !is.na(T293Quantile$ros) & !is.na(T293Quantile$mo))
    return(T293Quantile)
    }) %>% Reduce('rbind',.)
    return(regionQuantile)
}) %>% Reduce('rbind',.)


head(T293GroupedRegionQuantile)

lapply(seq(1,4), function(j){
    thisRegionQuantile <- subset(T293GroupedRegionQuantile, subset = Transcript_Biotype == unique(T293GroupedRegionQuantile$Transcript_Biotype)[j])
    pList1 <- lapply(seq(3,5), function(m){
        tmpDat <- thisRegionQuantile[,c(2,m)]
        colnames(tmpDat) <- c("Phenotype","Number")
        tmpDat$Number <- as.numeric(tmpDat$Number)
        tmpDat$Phenotype <- gsub('fu',"--", tmpDat$Phenotype)
        tmpDat <- subset(tmpDat, subset = grepl('mo',Phenotype))
        ggplot(tmpDat, aes(x = Phenotype, y = Number,color = Phenotype)) + ggtitle(colnames(thisRegionQuantile)[m])+
        stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
        geom_boxplot() + scale_color_viridis_d()+ theme_classic() +
        theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    })
    pList2 <- lapply(seq(3,5), function(m){
        tmpDat <- thisRegionQuantile[,c(2,m)]
        colnames(tmpDat) <- c("Phenotype","Number")
        tmpDat$Number <- as.numeric(tmpDat$Number)
        tmpDat$Phenotype <- gsub('fu',"--", tmpDat$Phenotype)
        tmpDat <- subset(tmpDat, subset = !grepl('mo',Phenotype))
        ggplot(tmpDat, aes(x = Phenotype, y = Number,color = Phenotype)) + ggtitle(colnames(thisRegionQuantile)[m])+
        stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
        geom_boxplot() + scale_color_viridis_d()+ theme_classic() +
        theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    })
    pdf(paste0('293T.WT.',unique(T293GroupedRegionQuantile$Transcript_Biotype)[j],'.mo_grouped_quantile.number.pdf'),width = 20,height=5)
    do.call("grid.arrange", c(pList1, nrow = 1))
    dev.off()
    pdf(paste0('293T.WT.',unique(T293GroupedRegionQuantile$Transcript_Biotype)[j],'.ros_grouped_quantile.number.pdf'),width = 20,height=5)
    do.call("grid.arrange", c(pList2, nrow = 1))
    dev.off()
})

thisDat <- subset(t293Dat, subset = FREQ >= myquantile[3])
t293DatTop30Cor <- lapply(seq_along(unique(thisDat$POS)),function(m){
    tmpDat <- subset(thisDat, subset = POS == unique(thisDat$POS)[m])
    if(dim(tmpDat)[1] > 2){
        print(dim(tmpDat)[1])
        p1 <- cor.test(tmpDat$FREQ,tmpDat$ros)
        p2 <- cor.test(tmpDat$FREQ,tmpDat$mo)
        return(cbind(tmpDat[1,], ros_pv = p1$p.value, ros_cor = p1$estimate, mo_pv = p2$p.value, mo_cor = p2$estimate))
    }else{
        return(NULL)
    }
}) %>% Reduce('rbind',.)

subset(t293DatTop30Cor, subset = tmp & (ros_cor >= 0) & (Gene_Name %in% c("MT-CYB",
"MT-ND5",
"MT-RNR2",
"MT-ND4",
"MT-RNR1",
"MT-CO1",
"MT-CO3",
"MT-ND2",
"MT-TL2")), select = c(POS,REF,ALT,Functional_Class, Effect,Gene_Name, Codon_Change,Amino_Acid_change))

###### off target
setwd('../off-target-vcf/')


thisDir = "293t//noPheno/"
dirlist = grep('vcf$',dir(thisDir),value=TRUE)
t293offTargetDat <- data.frame(lapply(seq_along(dirlist), function(n){
#t293Dat <- data.frame(lapply(c(4), function(n){
                flist <- grep('anno.vcf',dir(file.path(thisDir, dirlist[n])), value = TRUE)
                dat <- lapply(seq_along(flist), function(m){
#                dat <- lapply(seq_along(flist), function(m){
#                    print(flist[m])
                    x<- read.delim(file.path(thisDir,dirlist[n],flist[m]),sep="\t",header=TRUE)
                    if(dim(x)[1] != 0){
                        tmpdat <- cbind(
                           rep(gsub("_?vcf","",dirlist[n]),length(x$POS)),
                           rep(gsub("\\.anno.*","",flist[m]),length(x$POS)),
                           x$POS,x$REF,x$ALT,
                           lapply(seq_along(x$INFO), function(i){
                               tmp <- unlist(strsplit(x$INFO[i],split = "="))
                               tmp <- tmp[length(tmp)]
                               tmp <- gsub(')','',tmp)
                               tmp <- unlist(strsplit(tmp,split = "\\("))
                               if(tmp[1] == "INTERGENIC"){
                                   return(t(c(tmp[1],rep('-',10))))
                                }else{
                                   if(length(unlist(strsplit(tmp[2], split = '\\|'))) <10){
                                       return(t(c(tmp[1], unlist(strsplit(tmp[2], split = '\\|')),rep('-',10-length(unlist(strsplit(tmp[2], split = '\\|')))))))
                                    }else{
                                   return(t(c(tmp[1], unlist(strsplit(tmp[2], split = '\\|')))))
                                    }
                                }
                           }) %>% Reduce('rbind',.),
                           lapply(seq_along(x$INFO), function(i){
                               return(as.numeric(gsub(pattern = "%",replacement = "",unlist(strsplit(x[i,dim(x)[2]][1], split = ":"))[7])))
                           }) %>% unlist
                          )
#                        print(paste0(flist[m], '->',dim(tmpdat)[2]))
                        return(tmpdat)
                    }}) %>% Reduce('rbind',.)
                return(dat)
#    }))
            }) %>% Reduce('rbind', .))
colnames(t293offTargetDat) <- c("Phenotype",'cellID',"POS", "REF","ALT","Effect", "Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_change", "Gene_Name", "Transcript_BioType", "Gene_Coding", "Transcript_ID", "Exon", "INFO","FREQ")

dim(t293offTargetDat)

length(unique(subset(t293editDat2, subset = POS %in% t293Dat$POS)$POS))

length(unique(t293editDat2$POS))

POS

length(unique(subset(t293editDat2, subset = !(POS %in% t293Dat$POS))$POS))

hist(abs(as.integer(subset(t293editDat2, subset = !(POS %in% t293Dat$POS), select = POS)$POS) - 14867))

pdf('293T.edit1.distance.SNP.histogram.pdf')
ggplot(data.frame(distance = abs(as.integer(subset(t293editDat1, subset = !(POS %in% t293Dat$POS), select = POS)$POS) - 14867)),aes(x=distance)) + 
geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "lightblue") +
ylab("Number of sites")+
xlab("distance to the editted site")+
  geom_density(alpha=0.3,size=1.5, fill =  "#21918c") + theme_classic()
dev.off()

head(t293editDat1)

pdf('293T.edit1.distance.SNP.FREQ.histogram.pdf')
ggplot(data.frame(distance = abs(as.integer(subset(t293editDat1, subset = !(POS %in% t293Dat$POS), select = POS)$POS) - 14867),
                  FREQ=as.numeric(subset(t293editDat1, subset = !(POS %in% t293Dat$POS))$FREQ)),
                  aes(x=distance,y=FREQ)) + 
geom_point() +
ylab("FREQ")+
xlab("distance to the editted site")+ theme_classic()
dev.off()

pdf('293T.edit2.distance.SNP.histogram.pdf')
ggplot(data.frame(distance = abs(as.integer(subset(t293editDat2, subset = !(POS %in% t293Dat$POS), select = POS)$POS) - 15739)),aes(x=distance)) + 
geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "lightblue") +
ylab("Number of sites")+
xlab("distance to the editted site")+
  geom_density(alpha=0.3,size=1.5, fill =  "#21918c") + theme_classic()
dev.off()

pdf('293T.edit2.distance.SNP.FREQ.histogram.pdf')
ggplot(data.frame(distance = abs(as.integer(subset(t293editDat2, subset = !(POS %in% t293Dat$POS), select = POS)$POS) - 15739),
                  FREQ=as.numeric(subset(t293editDat2, subset = !(POS %in% t293Dat$POS))$FREQ)),
                  aes(x=distance,y=FREQ)) + 
geom_point() +
ylab("FREQ")+
xlab("distance to the editted site")+ theme_classic()
dev.off()

t293DatExample <- subset(t293Dat,POS %in% c(14867))
t293editDatExample <- subset(t293editDat1,POS %in% c(14867))
t293editDatExample$FREQ <- as.numeric(t293editDatExample$FREQ)
t293DatExample$Genotype <- "WT"
t293editDatExample$Genotype <- "edit"
t293Example <- rbind(t293DatExample[,c(1:18,21)],t293editDatExample)
t293Example$Genotype <- factor(t293Example$Genotype,  levels = c('WT',"edit"))

pdf('293T.WT-edit1.Example.pdf',height=5)



p1 <- ggplot(subset(t293Example, subset = !grepl('edit2',Cell)),aes(x = Genotype, y = FREQ,color = Genotype)) + 
        facet_grid("~POS") +
        geom_beeswarm(cex = 2.5, corral = "wrap")+ 
        stat_compare_means(method = "wilcox.test", label="p.format", label.x = 1.3, label.y = 70) + 
        scale_color_manual(values = RColorBrewer::brewer.pal(n = 2,name = "Dark2"))+ 
        theme_classic()+theme(legend.position = "none")
p1


t293Stat1 <- t293Example %>% group_by(Genotype, POS) %>% summarise(percent = n())
head(t293Stat1)

t293Stat1$percent <- c(t293Stat1$percent[1]/length(unique(t293Dat$Cell))*100, t293Stat1$percent[2]/length(unique(t293editDat1$Cell))*100)
t293Stat2<- t293Stat1
head(t293Stat1)

t293Stat2$percent <- 100 - t293Stat1$percent
t293Stat2$ifDetected <- "NO"
t293Stat1$ifDetected <- "YES"
t293Stat <- rbind(t293Stat1,t293Stat2)
t293Stat$ifDetected <- factor(t293Stat$ifDetected, levels = c("NO","YES"))

length(unique(t293Dat$Cell))

p2 <- ggplot(t293Stat, aes(x = Genotype, y = percent, fill = forcats::fct_rev(ifDetected))) + 
geom_bar(stat = "identity",position = "fill")+ 
scale_fill_manual(values = c("darkgreen","gray"))+ facet_grid("~POS") + theme_classic() + ylab("Percentage of mutation detected cell")+guides(fill=guide_legend(title="ifDetected"))
do.call(grid.arrange, c(list(p1,p2), ncol = 2))


pdf('293T.WT-edit1.Example.pdf',height=5,width=5)
do.call(grid.arrange, c(list(p1,p2), ncol = 2))

    dev.off()

t293editDat2$relDist <- abs(as.integer(t293editDat2$POS) - 15739)
t293editDat2$FREQ <- as.numeric(t293editDat2$FREQ)
#as.integer(subset(t293editDat1, subset = !(POS %in% t293Dat$POS), select = POS)$POS) - 15739)
ggplot(subset(t293editDat1, subset = !(POS %in% t293Dat$POS)), aes(x = relDist, y= FREQ)) + geom_point()+ theme_classic()

length(unique(subset(t293editDat2, subset = (POS %in% t293editDat1$POS) & !(POS %in% t293Dat$POS))$POS)) - 2584

t293editDat1$Mut <- NULL
t293editDat2$Mut <- NULL
t293editDat1$Mutation <- paste(t293editDat1$REF, t293editDat1$ALT, sep= "->")
t293editDat2$Mutation <- paste(t293editDat2$REF, t293editDat2$ALT, sep= "->")


t293scEditMutCount1 <- subset(t293editDat1, subset = !grepl('N', Mutation) & !grepl('TRUE', Mutation)& !(POS %in% t293Dat$POS)) %>% group_by(Cell, Mutation) %>% summarize(count = n())

t293scEditMutCount1 <- lapply(seq_along(unique(t293scEditMutCount1$Cell)), function(n){
    t293scMutCountThisCell <- subset(t293scEditMutCount1, subset = Cell == unique(t293scEditMutCount1$Cell)[n])
    t293scMutCountThisCell$ratio = t293scMutCountThisCell$count/sum(t293scMutCountThisCell$count)*100
    return(t293scMutCountThisCell)
}) %>% Reduce('rbind', .)

t293scEditMutCount2 <- subset(t293editDat2, subset = !grepl('N', Mutation) & !grepl('TRUE', Mutation)& !(POS %in% t293Dat$POS)) %>% group_by(Cell, Mutation) %>% summarize(count = n())

t293scEditMutCount2 <- lapply(seq_along(unique(t293scEditMutCount2$Cell)), function(n){
    t293scMutCountThisCell <- subset(t293scEditMutCount2, subset = Cell == unique(t293scEditMutCount2$Cell)[n])
    t293scMutCountThisCell$ratio = t293scMutCountThisCell$count/sum(t293scMutCountThisCell$count)*100
    return(t293scMutCountThisCell)
}) %>% Reduce('rbind', .)




p1 <- ggplot(subset(t293scEditMutCount1,subset = Mutation %in% c("A->G", "G->A", "C->T","T->C")),
                     aes(x= Mutation, y= ratio, colour = Mutation))+
geom_boxplot() + 
#geom_area(aes(fill = variable, group = variable),alpha = 0.5, position = 'identity') +
theme_classic()  +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p2 <- ggplot(subset(t293scEditMutCount1,subset = !(Mutation %in% c("A->G", "G->A", "C->T","T->C"))),
                     aes(x= Mutation, y= ratio, colour = Mutation))+
geom_boxplot() + 
#geom_area(aes(fill = variable, group = variable),alpha = 0.5, position = 'identity') +
theme_classic()  +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
pdf('293T.edit1.nucleotide.ratio.pdf',height = 5)
grid.arrange(grobs=list(p1,p2),layout_matrix = matrix(c(1,2,2),nrow=1),nrow = 1)
dev.off()



p1.2 <- ggplot(subset(t293scEditMutCount2,subset = Mutation %in% c("A->G", "G->A", "C->T","T->C")),
                     aes(x= Mutation, y= ratio, colour = Mutation))+
geom_boxplot() + 
#geom_area(aes(fill = variable, group = variable),alpha = 0.5, position = 'identity') +
theme_classic()  +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p2.2 <- ggplot(subset(t293scEditMutCount2,subset = !(Mutation %in% c("A->G", "G->A", "C->T","T->C"))),
                     aes(x= Mutation, y= ratio, colour = Mutation))+
geom_boxplot() +
#geom_area(aes(fill = variable, group = variable),alpha = 0.5, position = 'identity') +
theme_classic()  +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
grid.arrange(grobs=list(p1.2,p2.2),width=c(1,2),nrow = 1)
pdf('293T.edit2.nucleotide.ratio.pdf',height = 5)
grid.arrange(grobs=list(p1,p2),layout_matrix = matrix(c(1,2,2),nrow=1),nrow = 1)
dev.off()

t293editDat2Both <- subset(t293editDat2, subset = (POS %in% t293editDat1$POS) & !(POS %in% t293Dat$POS))

t293editDat2BothCount <- t293editDat2Both %>% group_by(Cell, Mutation, Gene_Name, Transcript_BioType) %>% summarize(count = n())
t293editDat1Both <- subset(t293editDat1, subset = (POS %in% t293editDat2$POS) & !(POS %in% t293Dat$POS))

t293editDat1BothCount <- t293editDat1Both %>% group_by(Cell, Mutation, Gene_Name, Transcript_BioType) %>% summarize(count = n())

#t293scMutRatio <- lapply(seq_along(unique(t293scMutCount$Cell)), function(n){
#    t293scMutCountThisCell <- subset(t293scMutCount, subset = Cell == unique(t293scMutCount$Cell)[n])
#    t293scMutRatioThisCell <- lapply(seq_along(unique(t293scMutCountThisCell$Transcript_BioType)), function(m){
#        tmpDatCount <- subset(t293scMutCountThisCell, subset = Transcript_BioType == unique(t293scMutCountThisCell$Transcript_BioType)[m])
#        tmpDatRatio <- lapply(seq_along(unique(tmpDatCount$Gene_Name)), function(i){
#            tmpDat <- subset(tmpDatCount, subset = Gene_Name == unique(tmpDatCount$Gene_Name)[i])
#            REFbase <- c('T','C','G','A')[!(c('T','C','G','A') %in% tmpDat$REF)]
#            if(length(REFbase) > 0){
#
#                tmpDatAppend <- data.frame(Cell= rep(tmpDat$Cell[1], length(REFbase)),
#                                       REF = REFbase,
#                                       Gene_Name = rep(tmpDat$Gene_Name[1], length(REFbase)),
#                                       Transcript_BioType = rep(tmpDat$Transcript_BioType[1], length(REFbase)),
#                                       count = rep(0, length(REFbase))
#                                      )
#                tmpDat <- rbind(tmpDat, tmpDatAppend)
#            }
#            tmpDat$ratio <- tmpDat$count/sum(t293scMutCountThisCell$count)*100
#            tmpDat <- rbind(tmpDat, data.frame(
#                Cell= tmpDat$'Cell'[1],
#                REF = "WM",
#                Gene_Name = tmpDat$Gene_Name[1],
#                Transcript_BioType = tmpDat$Transcript_BioType[1],
#                count =sum(tmpDat$count),
#                ratio = sum(tmpDat$count)/sum(t293scMutCountThisCell$count)*100
#                )
#            )
#            return(tmpDat)
#            
#        }) %>% Reduce('rbind', .)
#    }) %>% Reduce('rbind', .)
#}) %>% Reduce('rbind', .)


colorList = c('C->T' = '#A6CEE3',
 'G->A'= '#1F78B4',
'C->A' = '#B2DF8A',
'G->T' = '#33A02C',
'C->G' = '#FB9A99',
'A->C' = '#E31A1C',
'T->C' = '#FDBF6F',
'A->G' = '#FF7F00',
'T->A' = '#CAB2D6',
 'T->G'= '#6A3D9A'
)
p1 <- ggplot(subset(t293editDat2BothCount,subset = Transcript_BioType == "protein_coding"), aes(x = Gene_Name, colour = Mutation,y = count))+ 
geom_boxplot() + theme_classic() + scale_colour_manual(values = colorList) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2 <- ggplot(subset(t293editDat2BothCount,subset = Transcript_BioType == "Mt_tRNA"), aes(x = Gene_Name, colour = Mutation,y = count))+ 
geom_boxplot() + theme_classic() +  scale_colour_manual(values = colorList) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p3 <- ggplot(subset(t293editDat2BothCount,subset = Transcript_BioType == "Mt_rRNA"), aes(x = Gene_Name, colour = Mutation,y = count))+ 
geom_boxplot() + theme_classic() + scale_colour_manual(values = colorList) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
pdf('293T.edit2-offtarget-in-both.boxplot.pdf',width = 10,height = 5)
grid.arrange(grobs= list(p1, p2, p3), nrow = 1, widths =c(12,5,2) )
dev.off()

p1 <- ggplot(subset(t293editDat1BothCount,subset = Transcript_BioType == "protein_coding"), aes(x = Gene_Name, colour = Mutation,y = count))+ 
geom_boxplot() + theme_classic() + scale_colour_manual(values = colorList) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2 <- ggplot(subset(t293editDat1BothCount,subset = Transcript_BioType == "Mt_tRNA"), aes(x = Gene_Name, colour = Mutation,y = count))+ 
geom_boxplot() + theme_classic() +  scale_colour_manual(values = colorList) +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
p3 <- ggplot(subset(t293editDat1BothCount,subset = Transcript_BioType == "Mt_rRNA"), aes(x = Gene_Name, colour = Mutation,y = count))+ 
geom_boxplot() + theme_classic() + scale_colour_manual(values = colorList) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
pdf('293T.edit1-offtarget-in-both.boxplot.pdf',width = 10,height = 5)
grid.arrange(grobs= list(p1, p2, p3), nrow = 1, widths =c(12,5,2) )
dev.off()

table(unique(t293editDat1$POS) %in% unique(t293bulkeditDat$POS))

table(unique(t293editDat2$POS) %in% unique(t293bulkeditDat$POS))

table(unique(t293bulkeditDat$POS) %in% unique(t293editDat1$POS))

table(unique(t293bulkeditDat$POS) %in% unique(t293editDat2$POS))

table(unique(t293editDat2$POS)[unique(t293editDat2$POS) %in% unique(t293bulkeditDat$POS)] %in% unique(t293editDat1$POS))

library(ggVennDiagram)
ggVennDiagram(list(edit1=unique(t293editDat1$POS), edit2= unique(t293editDat2$POS), bulk =unique(t293bulkeditDat$POS)))

head(t293editDat)

editPheno <- read.table('../../293t.edit.phenotype.txt', header=TRUE, row.names = 1,sep="\t")
editPheno2 <- lapply(seq_along(t293editDat$Cell), function(m){
    tmp <- editPheno[rownames(editPheno) == t293editDat$Cell[m],]
    if(dim(tmp)[1] == 0){
        return(c(NA,NA))
    }else{
        return(tmp)
    }
}) %>% Reduce('rbind',.)
t293editDat <- cbind(t293editDat, editPheno2)

pos1 <- unique(t293editDat$POS[!(t293editDat$POS %in% t293Dat$POS) & grepl('edit1', t293editDat$Cell)])
pos2 <- unique(t293editDat$POS[!(t293editDat$POS %in% t293Dat$POS) & grepl('edit2', t293editDat$Cell)])

pos <- pos1[pos1 %in% pos2]
length(pos)

cor.test(subset(t293editDat , subset = (POS %in% pos) & grepl('edit1', t293editDat$Cell))$FREQ, 
         subset(t293editDat , subset = (POS %in% pos) & grepl('edit1', t293editDat$Cell))$mo)
cor.test(subset(t293editDat , subset = (POS %in% pos) & grepl('edit1', t293editDat$Cell))$FREQ, 
         subset(t293editDat , subset = (POS %in% pos) & grepl('edit1', t293editDat$Cell))$ros)
cor.test(subset(t293editDat , subset = (POS %in% pos) & grepl('edit2', t293editDat$Cell))$FREQ, 
         subset(t293editDat , subset = (POS %in% pos) & grepl('edit2', t293editDat$Cell))$mo)
cor.test(subset(t293editDat , subset = (POS %in% pos) & grepl('edit2', t293editDat$Cell))$FREQ, 
         subset(t293editDat , subset = (POS %in% pos) & grepl('edit2', t293editDat$Cell))$ros)

thisstat <-lapply(seq(0,60,by = 10),function(n){
    step <- 10
    low <- n+0
    high <- low+10
    print(c(low,high))
    corr <- c(cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$mo)$estimate,
    cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$ros)$estimate,
    cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$mo)$estimate,
    cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$ros)$estimate
    )
    pv <- c(cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$mo)$p.value,
    cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit1', t293editDat$Cell))$ros)$p.value,
    cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$mo)$p.value,
    cor.test(subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$FREQ, 
             subset(t293editDat , subset = (POS %in% pos) & (FREQ >= low) & (FREQ < high) & grepl('edit2', t293editDat$Cell))$ros)$p.value
    )
    return(data.frame(group = paste0(low,'-',high), edit = c(1,1,2,2),pheno = c('mo','ros'), corr,pv))
}) %>% Reduce('rbind',.)


#colnames(thisstat) <- c('FREQ','edit','pheno','corr','pv')
#thisstat$group <- paste(paste0('edit-',thisstat$edit), thisstat$pheno,sep=":")
#thisstat$edit <- NULL
#thisstat$pheno <- NULL
thisstat$sig <- lapply(seq_along(thisstat$pv), function(m){
    if(thisstat$pv[m] <= 0.001){
        return('***')
    }else if(thisstat$pv[m] <= 0.01){
        return('**')
    }else if(thisstat$pv[m] <= 0.05){
        return('*')
    }else {
        return('ns')
    }
}) %>% unlist

library(reshape2)
cordat <- reshape2::dcast(thisstat[,c(1,2,4)],formula = FREQ~group, value.var = "corr")
sigdat <- reshape2::dcast(thisstat[,c(1,4,5)],formula = FREQ~group, value.var = "sig")


rownames(cordat) <- cordat$FREQ
rownames(sigdat) <- sigdat$FREQ
cordat$FREQ <- NULL
sigdat$FREQ <- NULL
pheatmap(cordat)

pdf("edit.cor.pheno.pdf")
pheatmap(cordat, display_numbers = sigdat,cluster_cols = F,cluster_rows = F,border_color = "white",fontsize_number = 20,breaks = seq(-1,1,length.out = 100))
dev.off()

write.table(cordat, "~/jjluo/edit.cor.txt",row.names = T,col.names = T,sep="\t")
write.table(sigdat, "~/jjluo/edit.cor-sig.txt",row.names = T,col.names = T,sep="\t")

tmp1 <- unique(subset(t293editDat,select = c(Phenotype, cellID, ros, mo)))
tmp2 <- unique(subset(t293Dat,select = c(Phenotype, cellID, ros, mo)))

tmpDat <- rbind(tmp1, tmp2)

tmpDat$Phenotype <- gsub('fu','--',tmpDat$Phenotype)
tmpDat$Phenotype <- gsub('zheng','\\+\\+',tmpDat$Phenotype)

tmpDat <- tmpDat[!is.na(tmpDat$ros),]

pdf('293T.edit-bulk-sc.ros.boxplot.pdf',width = 5, height=5)
ggplot(tmpDat, aes(x=Phenotype,y=ros, color = Phenotype)) + geom_boxplot(fill="white") + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf('293T.edit-bulk-sc.mo.boxplot.pdf',width = 5, height=5)
ggplot(tmpDat, aes(x=Phenotype,y=mo, color = Phenotype)) + geom_boxplot(fill="white") + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

head(helaDat)

helaDat <- subset(helaDat, subset = Cell %in% editPheno$FORMAT)

editPheno <- read.table('../../hela.phenotype.txt', header=TRUE, row.names = 1,sep="\t")
editPheno2 <- lapply(seq_along(helaDat$Cell), function(m){
    tmp <- editPheno[rownames(editPheno) == helaDat$Cell[m],]
    if(dim(tmp)[1] == 0){
        return(c(NA,NA))
    }else{
        return(tmp)
    }
}) %>% Reduce('rbind',.)
helaDat <- cbind(helaDat, editPheno2)

editPheno <- read.table('../../hela.phenotype.txt', header=TRUE, sep="\t")


tmpDat <- rbind(helaDat,t293Dat)

tmpDat$cellType <- gsub('hela.*', 'Hela', tmpDat$Phenotype)
tmpDat$cellType <- gsub('293t.*', '293T', tmpDat$cellType)
tmpDat$Phenotype <- gsub('hela','',tmpDat$Phenotype)
tmpDat$Phenotype <- gsub('293t-wt-','',tmpDat$Phenotype)
tmpDat$Phenotype <- gsub('fu','--',tmpDat$Phenotype)
tmpDat <- subset(tmpDat, subset = !is.na(tmpDat$ros))

p1 <- ggplot(tmpDat %>% group_by(Cell, cellType, Phenotype) %>% summarise(number = n(), Freq = median(FREQ)), aes(x = cellType, y = number,  colour = Freq)) + 
geom_violin(fill = NA, outlier.shape = NA) +
geom_beeswarm(cex = 2.5, corral = "wrap") + scale_colour_viridis_c(begin = 0.4, end = 1)+stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1.5 ) + theme_classic() 

p2 <- ggplot(tmpDat %>% group_by(Cell, Phenotype, cellType, Transcript_BioType) %>% 
           summarise(number = n(), Freq = median(FREQ)),
       aes(x = cellType, y = number, colour = Freq)) +
       geom_violin(fill = NA, outlier.shape = NA) +
       geom_beeswarm(cex = 2.5, corral = "wrap") +
       scale_colour_viridis_c(begin = 0.3, end = 1) +
       stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1.5 ) +
       facet_grid('~Transcript_BioType') +
    theme_classic() 

p3.0 <- ggplot(subset(tmpDat %>% group_by(Cell, Phenotype, cellType, Gene_Name, Transcript_BioType) %>% 
           summarise(number = n(), Freq = median(FREQ)), subset = Transcript_BioType == "protein_coding"),
       aes(x = cellType, colour = number, y = Freq))+
       geom_boxplot(fill = NA, outlier.shape = NA) +
       geom_beeswarm(cex = 2.5, corral = "wrap") +
       scale_colour_viridis_c(begin = 0.3, end = 1) +
       stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
       facet_grid('~Gene_Name') +
    theme_classic()
p3.1 <- ggplot(subset(tmpDat %>% group_by(Cell, Phenotype, cellType, Gene_Name, Transcript_BioType) %>% 
           summarise(number = n(), Freq = median(FREQ)), subset = Transcript_BioType == "protein_coding"),
       aes(x = cellType, colour = number, y = Freq))+
       geom_boxplot(fill = NA, outlier.shape = NA) +
       geom_beeswarm(cex = 2.5, corral = "wrap") +
       scale_colour_viridis_c(begin = 0.3, end = 1) +
       stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
       facet_grid('~Gene_Name') +
    theme_classic() + theme(legend.position = "none")
p3.2 <- ggplot(subset(tmpDat %>% group_by(Cell, Phenotype, cellType, Gene_Name, Transcript_BioType) %>% 
           summarise(number = n(), Freq = median(FREQ)), subset = Transcript_BioType == "Mt_rRNA"),
       aes(x = cellType, colour = number, y = Freq))+
       geom_boxplot(fill = NA, outlier.shape = NA) +
       geom_beeswarm(cex = 2.5, corral = "wrap") +
       scale_colour_viridis_c(begin = 0.3, end = 1) +
       stat_compare_means(method ="wilcox.test", label="p.format", label.x = 1 ) +
       facet_grid('~Gene_Name') +
    theme_classic() 

pdf('Hela-293T.global.pdf',width = 12, height = 5)
grid.arrange(grobs = list(p1, p2),widths=c(1,2))
dev.off()

p3.2

pdf('Hela-293T.category0.pdf',width = 20, height = 5)
grid.arrange(grobs = list(p3.0, p3.2),widths=c(4, 1))
dev.off()

unlist(tmpDat[1,])

tmpDat2 <- subset(tmpDat, subset = Gene_Name %in% c('MT-RNR2', 'MT-ATP6', 'MT-CO1','MT-CYB', 'MT-ND1', 'MT-ND3', 'MT-ND4', 'MT-ND4L','MT-RNR2', 'MT-ND5', 'MT-ND6'))

length(pos)

tmpDat3 <- lapply(seq_along(unique(tmpDat$POS)),function(m){
                    t1 <- subset(tmpDat, subset = (POS == unique(tmpDat$POS)[m]) & (cellType == "Hela"))
                    t2 <- subset(tmpDat, subset = (POS == unique(tmpDat$POS)[m]) & (cellType == "293T"))
                    tt <- data.frame(POS = unique(tmpDat$POS)[m],
                                     REF = t1$REF[1],
                                     ALT = t1$ALT[1],
                                     Effect = t1$Effect[1],
                                     Effect_Impact	= t1$Effect_Impact[1],
                                     Functional_Class = t1$Functional_Class[1],
                                     Codon_Change = t1$Codon_Change[1],
                                     Amino_Acid_change = t1$Amino_Acid_change[1],
                                     Hela_number = length(t1$cellType), Hela_median_Freq = median(t1$FREQ), Hela_mean_Freq = mean(t1$FREQ),
                                     T293_number = length(t2$cellType), T293_median_Freq = median(t2$FREQ), T293_mean_Freq = mean(t2$FREQ),
                                     Gene_Name=t1$Gene_Name[1], 
                                     Transcript_BioType=t1$Transcript_BioType[1]
                                    )
                    if(dim(t1)[1] > 2){
                        tt$ros_corr = cor.test(as.numeric(t1$FREQ), as.numeric(t1$ros))$estimate
                        tt$ros_pv = cor.test(as.numeric(t1$FREQ), as.numeric(t1$ros))$p.value
                        tt$mo_corr = cor.test(as.numeric(t1$FREQ), as.numeric(t1$mo))$estimate
                        tt$mo_pv = cor.test(as.numeric(t1$FREQ), as.numeric(t1$mo))$p.value
                    }else{
                        tt$ros_corr = NA
                        tt$ros_pv = NA
                        tt$mo_corr = NA
                        tt$mo_pv = NA
                    }
                    if(min(dim(t1)[1],dim(t2)[1]) == 0)
                        p1 <- NA
                    else
                        p1 <- wilcox.test(t1$FREQ,t2$FREQ)$p.value
                    tt$pv = p1
                    return(tt)
                    
}) %>% Reduce('rbind',.)

write.csv(subset(tmpDat, subset= Transcript_BioType == "protein_coding"),'../../Hela.significant.SNP.csv')

tmpDat3$fdr <- p.adjust(tmpDat3$pv, method = "fdr")

write.csv(tmpDat3,'../../Hela.SNP.csv')



#kidneyDat$tmp <- paste0(kidneyDat$Tissue, ":", kidneyDat$MutType)
if(!dir.exists('figures'))
    dir.create("figures")
pdf(file.path('figures',"overall.with_MO.freq2.pdf"))
my_comparisons <-list(c( "MO--","MO++"),
c("RCC:EGFR−","RCC:EGFR+"),
c( "NC:EGFR−","RCC:EGFR+"))

p1 <- ggplot(subset(kidneyDat %>% group_by(tmp,Functional_Class,Phenotype, cellID) %>% summarise(number = n()), 
                    subset= Functional_Class != "Other"), 
             aes(x=tmp, y=number)) +
    geom_violin(aes(fill = tmp),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p2 <- ggplot(subset(kidneyDat %>% group_by(tmp,Functional_Class,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
                    subset= Functional_Class != "Other"), aes(x = tmp, y = Freq)) + 
    geom_violin(aes(fill = tmp),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
grid.arrange(grobs = list(p1,p2), ncols =1)
dev.off()

write.table(kidneyDat,'kidney-table.txt',col.names = T, sep="\t",quote=F)

pdf(file.path('figures',"rRNA.freq.pdf"),width = 7)
my_comparisons <-list(c( "NC:EGFR−","RCC:EGFR−"),
c("RCC:EGFR−","RCC:EGFR+"),
c( "NC:EGFR−","RCC:EGFR+"))

p1 <- ggplot(subset(kidneyDat, subset = Transcript_BioType == "Mt_rRNA") %>% group_by(tmp,Gene_Name,cellID) %>% summarise(number = n()), 
             aes(x=tmp, y=number)) +
    geom_violin(aes(fill = tmp),drop=FALSE) + 
    facet_wrap(facets = "~Gene_Name") + stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p2 <- ggplot(subset(kidneyDat, subset = Transcript_BioType == "Mt_rRNA") %>% group_by(tmp,Gene_Name,cellID) %>% summarise(Freq = mean(FREQ)),
             aes(x = tmp, y = Freq)) + 
    geom_violin(aes(fill = tmp),drop=FALSE) + 
    facet_wrap(facets = "~Gene_Name") + stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
grid.arrange(grobs = list(p1,p2), ncols =1)
dev.off()

pdf(file.path('figures',"overall.with_phenotype.freq2.pdf"),width = 12)

p1<- ggplot(subset(subset(kidneyDat, subset = tmp == 'NC:EGFR-') %>% group_by(Functional_Class,Phenotype, cellID) %>% summarise(number = n()), 
                    subset= Functional_Class != "Other"), 
             aes(x=Phenotype, y=number)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2 <- ggplot(subset(subset(kidneyDat, subset = tmp == 'NC:EGFR-') %>% group_by(Functional_Class,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
                    subset= Functional_Class != "Other"), 
             aes(x=Phenotype, y=Freq)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p3<- ggplot(subset(subset(kidneyDat, subset = tmp == 'RCC:EGFR-') %>% group_by(Functional_Class,Phenotype, cellID) %>% summarise(number = n()), 
                    subset= Functional_Class != "Other"), 
             aes(x=Phenotype, y=number)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p4 <- ggplot(subset(subset(kidneyDat, subset = tmp == 'RCC:EGFR-') %>% group_by(Functional_Class,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
                    subset= Functional_Class != "Other"), 
             aes(x=Phenotype, y=Freq)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() +    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


p5<- ggplot(subset(subset(kidneyDat, subset = tmp == 'RCC:EGFR+') %>% group_by(Functional_Class,Phenotype, cellID) %>% summarise(number = n()), 
                    subset= Functional_Class != "Other"), 
             aes(x=Phenotype, y=number)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p6 <- ggplot(subset(subset(kidneyDat, subset = tmp == 'RCC:EGFR+') %>% group_by(Functional_Class,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
                    subset= Functional_Class != "Other"), 
             aes(x=Phenotype, y=Freq)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~Functional_Class") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
grid.arrange(grobs = list(p1,p2,p3,p4,p5,p6), ncol =2)
dev.off()


pdf(file.path('figures',"ND4L.with_phenotype.freq2.pdf"),width = 12)

p1<- ggplot(subset(kidneyDat, subset = (Transcript_BioType == "protein_coding")) %>% group_by(tmp,Phenotype, cellID) %>% summarise(number = n()), 
             aes(x=Phenotype, y=number)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~tmp") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p2 <- ggplot(subset(kidneyDat, subset =  (Transcript_BioType == "protein_coding")) %>% group_by(tmp,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
             aes(x=Phenotype, y=Freq)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~tmp") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p3<- ggplot(subset(kidneyDat, subset =  (Transcript_BioType == "protein_coding")) %>% group_by(tmp,Phenotype, cellID) %>% summarise(number = n()), 
             aes(x=Phenotype, y=number)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~tmp") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p4 <- ggplot(subset(kidneyDat, subset =  (Transcript_BioType == "protein_coding")) %>% group_by(tmp,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
             aes(x=Phenotype, y=Freq)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~tmp") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() +    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


p5<- ggplot(subset(kidneyDat, subset =  (Transcript_BioType == "protein_coding")) %>% group_by(tmp,Phenotype, cellID) %>% summarise(number = n()), 
             aes(x=Phenotype, y=number)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~tmp") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p6 <- ggplot(subset(kidneyDat, subset =  (Transcript_BioType == "protein_coding")) %>% group_by(tmp,Phenotype, cellID) %>% summarise(Freq = mean(FREQ)), 
             aes(x=Phenotype, y=Freq)) +
    geom_violin(aes(fill = Phenotype),drop=FALSE) + 
    facet_wrap(facets = "~tmp") + 
    #stat_compare_means(comparisons =  my_comparisons, method = "wilcox.test", label="p.format")+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

grid.arrange(grobs = list(p1,p2,p3,p4,p5,p6), ncol =2)
dev.off()

tmp <- subset(kidneyDat, subset = Gene_Name == "MT-ND4L") %>% group_by(Tissue, MutType,Phenotype,cellID) %>% summarise(number = n())
tmp$group <- paste0(tmp$MutType,":",tmp$Phenotype)
pdf("kidney.ND4L.detailed.boxplot.pdf",height=5)
ggplot(tmp, aes(x=group, y=number)) + geom_boxplot(aes(fill = group)) + facet_wrap('~Tissue') + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

tmp <- subset(kidneyDat, subset = Gene_Name == "MT-ND4L") %>% group_by(Tissue, POS,Functional_Class,Amino_Acid_change, MutType,Phenotype) %>% summarise(number = n())


tmp2 <- tmp %>% group_by(Tissue, Functional_Class,		MutType,	Phenotype,POS) %>% summarise(number = sum(number)) %>% top_n(wt = number,n = 3)

tmp2 <- tmp2[order(tmp2$number,decreasing = TRUE),]

tmp <- subset(kidneyDat, subset = POS %in% c(10750,10708,10683)) %>% group_by(Tissue, POS, MutType,Phenotype ,cellID) %>% summarise(Freq = mean(FREQ))

tmp$class <- paste(tmp$Tissue,tmp$MutType,tmp$Phenotype, sep=":")

ggplot(tmp, aes(x=class,y= Freq)) + geom_boxplot() + facet_wrap('~POS') + theme_classic() +theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tmp <- subset(kidneyDat, subset = Gene_Name == "MT-TE") %>% group_by(Tissue, POS, MutType,Phenotype,cellID) %>% summarise(number = n())


tmp2<- tmp %>% group_by(POS,MutType,Phenotype) %>% summarise(number = sum(number))

tmp2<- tmp2[order(tmp2$number,decreasing = T),]

head(kidneyDat[grep("NO4RCCE-ROSfu-",kidneyDat$cellID),])

gsub('vcf.*-','_',head(kidneyDat$cellID))

kidneyDat$cellID <- gsub('L','',kidneyDat$cellID)

kidney_pheno <- read.table('phenotype.txt', row.names = 1,sep="\t", header=TRUE)

rownames(kidney_pheno)[!(rownames(kidney_pheno) %in% kidneyDat$cellID)]

kidneyDat <- subset(kidneyDat, subset =cellID %in% rownames(kidney_pheno))

kidney_pheno <- lapply(seq_along(kidneyDat$cellID), function(m){return(subset(kidney_pheno, subset = rownames(kidney_pheno) == kidneyDat$cellID[m]))}) %>% Reduce('rbind',.)

kidneyDat <- cbind(kidneyDat, kidney_pheno)

head(kidneyDat %>% group_by(cellID, mo, ros) %>%  top_n(n = 1))

#thisDat = subset(t293Dat,subset=FREQ >=1)
myquantile = round(seq(0,100,length.out=4),2)
cells <- unique(kidneyDat$cellID)
kidneyQuantile <- data.frame(t(lapply(seq_along(cells), function(n){
    tmpDat <- subset(kidneyDat, subset = cellID %in% cells[n])
    cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
        ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
        return(ctmp)
    }) %>% unlist
    cellQuantile <- c(tmpDat$ros[1], tmpDat$mo[1],cellQuantile)
    cellQuantile <- data.frame(cellQuantile)
    return(cellQuantile)
}) %>% Reduce('cbind',.)))
colnames(kidneyQuantile) <- c('ros','mo', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
rownames(kidneyQuantile) <- unique(thisDat$Cell)
#head(T293Quantile)


head(kidneyDat)

pdf('../kidney/kidney.wt.ros-mo.correlation.lm.pdf',height = 5, width = 5.5)
thisDat <- lapply(seq_along(unique(kidneyDat$cellID)), function(m){
    tmpDat <- subset(kidneyDat, subset = cellID == unique(kidneyDat$cellID)[m])
    tmpDat <- cbind(Number = dim(tmpDat)[1], FREQ = mean(tmpDat$FREQ),head(subset(tmpDat, select =  c(ros,mo)),1))
    return(tmpDat)
}) %>% Reduce('rbind',.)
head(thisDat$ros)
ggplot(thisDat, aes(x = ros, y=mo)) +
geom_point(aes(size = Number,color = FREQ, alpha =0.5)) + 
geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x) +
scale_color_viridis_c(begin = 0.2,end = 1)+
scale_radius()+
stat_cor(label.y = 4800)+ ylab("mo")+#this means at 35th unit in the y axis, the r squared and p value will be shown
        stat_regline_equation(label.y = 5000) +
theme_classic()#this means at 30th unit regresion line equation will be shown
dev.off()

myquantile = round(seq(0,100,length.out=4),2)
condition <- unique(subset(kidneyDat, select = c(Tissue, MutType)))
kidneyQuantileList <- lapply(seq_len(dim(condition)[1]), function(k){
    thisDat <- subset(kidneyDat, subset = (Tissue == condition[k,1]) & (MutType == condition[k,2]))
    cells <- unique(thisDat$cellID)
    kidneyQuantile <- data.frame(t(lapply(seq_along(cells), function(n){
        tmpDat <- subset(thisDat, subset = cellID %in% cells[n])
        cellQuantile <- lapply(seq_len(length(myquantile)-1),function(m){
            ctmp = length(which((tmpDat$FREQ > myquantile[m]) & (tmpDat$FREQ <= myquantile[(m+1)])))
            return(ctmp)
        }) %>% unlist
        cellQuantile <- c(tmpDat$ros[1], tmpDat$mo[1],cellQuantile)
        cellQuantile <- data.frame(cellQuantile)
        return(cellQuantile)
    }) %>% Reduce('cbind',.)))
    colnames(kidneyQuantile) <- c('ros','mo', lapply(seq_len(length(myquantile)-1),function(m){return(paste0(myquantile[m],"-",myquantile[m+1],"%"))}) %>% unlist)
    rownames(kidneyQuantile) <- cells
    kidneyQuantile <- subset(kidneyQuantile, subset = !is.na(kidneyQuantile$ros) & !is.na(kidneyQuantile$mo))
    
    return(kidneyQuantile)
})
#head(T293Quantile)


pdf('kidney.tissue.phenotype.correlation.pdf', width=7,height=3)
p <- corrplot::corrplot(cor(kidneyQuantile[,1:2], kidneyQuantile[,3:dim(kidneyQuantile)[2]]),col = COL2("RdBu",100)[100:1])
grid.arrange(grobs = pList, ncol=3)
dev.off()


pval <- psych::corr.test(kidneyQuantileList[[1]][,1:2],kidneyQuantileList[[1]][3:5], method = "pearson")$p
pdf('kidney.NC.EGFR-.quantile.phenotype-correlation.quantiled.pdf',width = 5, height  = 5)
corrplot::corrplot(cor(kidneyQuantileList[[1]][,1:2], kidneyQuantileList[[1]][,3:dim(kidneyQuantileList[[1]])[2]]),col = COL2("RdBu",100)[100:1])
dev.off()
#                   p.mat=pval, insig="p-value",  tl.pos="n", sig.level=0)
pval <- psych::corr.test(kidneyQuantileList[[2]][,1:2],kidneyQuantileList[[2]][3:5], method = "pearson")$p
pdf('kidney.RCC.EGFR-.quantile.phenotype-correlation.quantiled.pdf',width = 5, height  = 5)
corrplot::corrplot(cor(kidneyQuantileList[[2]][,1:2], kidneyQuantileList[[2]][,3:dim(kidneyQuantileList[[2]])[2]]))
dev.off()
                   #,col = COL2("RdBu",100)[100:1],p.mat=pval, insig="p-value",   tl.pos="n", sig.level=0)
pval <- psych::corr.test(kidneyQuantileList[[3]][,1:2],kidneyQuantileList[[3]][3:5], method = "pearson")$p
pdf('kidney.RCC.EGFR+.quantile.phenotype-correlation.quantiled.pdf',width = 5, height  = 5)

corrplot::corrplot(cor(kidneyQuantileList[[3]][,1:2], kidneyQuantileList[[3]][,3:dim(kidneyQuantileList[[3]])[2]]))
dev.off()#,col = COL2("RdBu",100)[100:1],p.mat=pval, insig="p-value", 
                                               #tl.pos="n", sig.level=0)


## SNP stat
tmp <- (subset(kidneyDat, subset = (FREQ >= 66.67)) %>% group_by(Tissue, MutType, cellID, POS) %>% summarise(Number = n())) %>%
group_by(Tissue, MutType,POS) %>% summarise(Number = sum(Number))

# core snp stat
tmp <- subset(tmp, subset = Number >= 50)
       #Number >= 100)

data.frame(table(unique(subset(kidneyDat, subset = POS %in% tmp$POS,select =c(POS, Gene_Name)))$Gene_Name))

library(ggplot2)
library(plotly)
library(dplyr)
corsnpInfo <- subset(kidneyDat, subset = POS %in% tmp$POS)
corsnpInfo$Group <- paste(corsnpInfo$Tissue, corsnpInfo$MutType, sep = ":")
corsnpInfoStat1 <- unique(subset(corsnpInfo, select = c(POS,Transcript_BioType))) %>% group_by(Transcript_BioType) %>% summarise(number = n())
corsnpInfoStat2 <- unique(subset(corsnpInfo, select = c(POS,Group,Transcript_BioType))) %>% group_by(Group,Transcript_BioType) %>% summarise(number = n())

pie(table(unique(subset(corsnpInfo, select = c(POS,Transcript_BioType)))$Transcript_BioType),labels = TRUE)

plot_ly(corsnpInfoStat) %>%
  add_pie(labels = ~`group`, values = ~`value1`, 
          type = 'pie', hole = 0.7, sort = F,
          marker = list(line = list(width = 2))) %>%
  add_pie(sample_data, labels = ~`group`, values = ~`value2`, 
          domain = list(
            x = c(0.15, 0.85),
            y = c(0.15, 0.85)),
          sort = F)
  
#ggplot(corsnpInfoStat, aes(x = factor(Transcript_BioType), y = number, fill = factor(Group))) +
#          geom_col() +
#         scale_x_discrete()+#limits = c("Group1", "Group2")) +
#          coord_polar("y")

corsnpInfoStat1

plot2 <- plot_ly(corsnpInfoStat1) %>%
  add_pie(corsnpInfoStat,labels = ~Transcript_BioType, values = ~number,hole=0.7,name="Category",sort = F,type='pie',)
plot2.update_traces(textinfo='value')
#plot1 <- plot_ly(corsnpInfoStat) %>%
#  add_pie(corsnpInfoStat,labels = ~Group,values = ~number,hole=0.7,name = "Sub-Category",domain = list(x = c(0.15,0.85),y=c(0.15,0.85),sort=F))
#
#combined_plot <- subplot(plot2, plot1, nrows = 1, margin = 0.03) %>%
#  layout(layout)
#
#combined_plot

ggplot(data.frame(Number = quantile(tmp$Number,probs = seq(0,1,by=0.1)), percentile =seq(0,1,by=0.1)*100), aes(x= percentile,y = Number, color= Number)) +
geom_point(size=6, alpha=0.3) + geom_line(size=1.5) +theme_classic()

write.table(kidneyDat, 'kidney_table.txt',col.names = TRUE,)

#
#gg(p1)
#ggsave(filename = "kidney.NC.EGFR-.quantile.phenotype-correlation.quantiled.pdf", plot = replayPlot(p1))#print(p2)
#ggsave(filename = "kidney.RCC.EGFR-.quantile.phenotype-correlation.quantiled.pdf", plot = replayPlot(p2))
#ggsave(filename = "kidney.RCC.EGFR+.quantile.phenotype-correlation.quantiled.pdf", plot = replayPlot(p3))
#print(p3)
#dev.off()
#pdf('kidney.RCC.EGFR-.quantile.phenotype-correlation.quantiled.pdf',width = 5, height  = 5)
##print(p1)
#print(p2)
##print(p3)
#
#dev.off()
pdf('kidney.RCC.EGFR+.quantile.phenotype-correlation.quantiled.pdf',width = 5, height  = 5)
##print(p1)
##print(p2)
p3
#
dev.off()

tmp2 <- subset(kidneyDat, subset =(FREQ >= 33.33) & (FREQ <66.67))

tmp2.1 <-  tmp2 %>% group_by(Tissue, MutType, cellID, POS) %>% summarise(Number = n()) %>%
group_by(Tissue, MutType,POS) %>% summarise(Number = sum(Number))
ggplot(data.frame(Number = quantile(tmp2.1$Number,probs = seq(0,1,by=0.1)), percentile =seq(0,1,by=0.1)*100), aes(x= percentile,y = Number, color= Number)) +
geom_point(size=6, alpha=0.3) + geom_line(size=1.5) +theme_classic()

tmp2 <- subset(tmp2, subset = POS %in% tmp2.1$POS[tmp2.1$Number >= 50])

data.frame(table(unique(subset(kidneyDat,subset = POS %in% tmp2.1$POS,select=c(POS,Transcript_BioType)))$Transcript_BioType))

pos <- c(
10683
)
tmpDat<- subset(kidneyDat, subset = POS %in% pos) %>% group_by(POS,Tissue, MutType,cellID, FREQ) %>% summarise(n = n())
unique(subset(kidneyDat, subset = (POS %in% pos) & (Functional_Class !="Other"), select = c(POS,Gene_Name, Functional_Class, Amino_Acid_change)))
tmpDat$group <- paste(tmpDat$Tissue, tmpDat$MutType, sep=":")
pdf('tmp.kidney_core.10683.snp.stat.boxplot.pdf',height= 12,width=14)
ggplot(tmpDat, aes(x=group,y= FREQ)) +geom_boxplot() + 
stat_compare_means(comparisons = list(c("NC:EGFR-","RCC:EGFR-"), c("NC:EGFR-","RCC:EGFR+"),c("RCC:EGFR+","RCC:EGFR-")),
                   method = "wilcox.test", 
                   , label="p.format"
                  )+ facet_wrap('~POS',ncol = 4) + theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

kidneyDat$group <- paste(kidneyDat$Tissue,kidneyDat$MutType,sep=":")
head(kidneyDat)

pos <- unique(c(
12308,
12372,
12801,
13617,
16399,
195,
9477,
    872,11172,8563,
    10750,
    14793,
    12335,
    12329,
    16266,
    16290,
    16319,
    16231,
    16235,
    16189,
    3571,
    2699,
    1736,
    539))
length(pos)

pos <- unique(c(
12308,
12372,
12801,
13617,
16399,
195,
9477,
    8772,11172,8563,
    10750,
    14793,
    12335,
    12329,
    16266,
    16290,
    16319,
    16231,
    16235,
    16189,
    3571,
    2699,
    1736,
    539))
tmpDat<- subset(kidneyDat, subset = POS %in% pos) %>% group_by(POS,Tissue, MutType,cellID, FREQ) %>% summarise(n = n())
tmpDat$group <- paste(tmpDat$Tissue, tmpDat$MutType, sep=":")
pdf('kindey.freq.pos.3.pdf',height=13,width=13)
ggplot(tmpDat, aes(x=group,y= FREQ,color = group)) +geom_boxplot() + 
stat_compare_means(comparisons = list(c("NC:EGFR-","RCC:EGFR-"), c("NC:EGFR-","RCC:EGFR+"),c("RCC:EGFR+","RCC:EGFR-")),
                   method = "wilcox.test", 
                   , label="p.format"
                  )+ facet_wrap("~POS",nrow = 4)+ theme_classic()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

write.csv(subset(kidneyDat, subset = POS %in% c(12308,
12372,2523,
16189,
53,
16234
,16243,8772,1736,16126,235,11172,8563,10750)) %>% group_by(POS, REF, ALT, Gene_Name, Functional_Class, Amino_Acid_change) %>% summarise(Number = n(), FREQ=mean(FREQ)), 'mut_site.info2.csv')

pos <- c(8772,1736,16126,235,11172,8563,10750)
corDat <- lapply(seq_along(pos), function(m){
    tmpDat <- subset(kidneyDat, subset = POS == pos[m])
    cor1 <- cor.test(tmpDat$FREQ, tmpDat$ros, method="spearman")
    cor2 <- cor.test(tmpDat$FREQ, tmpDat$mo , method="spearman")
    return(c(cor1$estimate,cor2$estimate))
    }) %>% Reduce('rbind',.)
pvDat <- lapply(seq_along(pos), function(m){
    tmpDat <- subset(kidneyDat, subset = POS == pos[m])
    cor1 <- cor.test(tmpDat$FREQ, tmpDat$ros, method="spearman")
    cor2 <- cor.test(tmpDat$FREQ, tmpDat$mo , method="spearman")
    if(cor1$p.value <= 0.001){
        pv1 <- "***"
    }else if(cor1$p.value <= 0.01){
        pv1 <- "**"
    }else if(cor1$p.value <= 0.05){
        pv1 <- "*"
    }else{
        pv1 <- "ns"
    }
    if(cor2$p.value <= 0.001){
        pv2 <- "***"
    }else if(cor2$p.value <= 0.01){
        pv2 <- "**"
    }else if(cor2$p.value <= 0.05){
        pv2 <- "*"
    }else{
        pv2 <- "ns"
    }
    return(c(pv1, pv2))
    }) %>% Reduce('rbind',.)

rownames(corDat) <- pos
colnames(corDat) <- c("ros","mo")
rownames(pvDat) <- pos
colnames(pvDat) <- c("ros","mo")

pdf('kidney.selected.site.cor.pheatmap.pdf')
pheatmap(corDat,display_numbers =  pvDat,
               fontsize_number = 20,cellwidth = 35,
               breaks=seq(-0.5,0.5, length.out = 101),
               cellheight = 35, border_color = "white")

dev.off()
