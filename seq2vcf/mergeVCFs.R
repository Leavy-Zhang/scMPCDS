library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(pheatmap)
library(tidyverse)
library(ggbeeswarm)
library(cowplot)
library(corrplot)

setwd('/path/to/project')


thisDir = "vcf_dir"
dirlist = grep('anno.vcf',dir(thisDir),value=TRUE)
dat <- data.frame(lapply(seq_along(dirlist), function(n){
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
colnames(dat) <- c("Phenotype",'cellID',"POS", "REF","ALT","Effect", "Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_change", "Gene_Name", "Transcript_BioType", "Gene_Coding", "Transcript_ID", "Exon", "INFO","FREQ")


dat$Cell <- paste(dat$Phenotype,dat$cellID, sep=":")

fixEff <- function(dat){
    dat$Transcript_BioType <- lapply(seq_along(dat$Transcript_BioType),function(m){
        if((as.numeric(dat$POS[m]) <= 576)){
            return("Dloop")
        }else if((as.numeric(dat$POS[m]) >= 16024)){
            return("Dloop")
        }else{
            return(dat$Transcript_BioType[m])
        }
    }) %>% unlist
    dat$Transcript_BioType <- lapply(seq_along(dat$Transcript_BioType),function(m){
        if(dat$Transcript_BioType[m] == "-"){
            return("Intergenic")
        }else{
            return(dat$Transcript_BioType[m])
        }
    }) %>% unlist
    dat$Functional_Class <- lapply(seq_along(dat$Transcript_BioType),function(m){
        if(dat$Functional_Class[m] == "-"){
            return("Other")
        }else if(dat$Functional_Class[m] == ""){
            return("Other")
        }else{
            return(dat$Functional_Class[m])
        }
    }) %>% unlist
    return(dat)
}


dat <- fixEff(dat)
dat$FREQ <- as.numeric(dat$FREQ)
pheno <- read.table('phenotype.txt',sep="\t",header=TRUE,row.names = 1)
pheno2 <- lapply(seq_along(dat$Cell), function(m){
    tmp <- pheno[rownames(pheno) == dat$Cell[m],]
    if(dim(tmp)[1] == 0){
        return(c(NA,NA))
    }else{
        return(tmp)
    }
}) %>% Reduce('rbind',.)
dat <- cbind(dat, pheno2)
write.csv(dat, "path/to/vcf_table.csv")
