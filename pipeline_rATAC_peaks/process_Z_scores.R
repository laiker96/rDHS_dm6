library(tidyverse)

histone_left <- read.delim("dataset_left.bed", header=F)
histone_right <- read.delim("dataset_right.bed", header=F)
histone_center <- read.delim("dataset_central.bed", header=F)
histone <- data.frame(V1 = histone_left$V1, 
                      V2 = histone_left$V2, 
                      V3 = histone_left$V3, 
                      V4 = histone_left$V4, 
                      V5 = pmax(histone_left$V5, histone_right$V5, histone_center$V5), 
                      V6 = pmax(histone_left$V6, histone_right$V6, histone_center$V6),
                      V7 = pmax(histone_left$V7, histone_right$V7, histone_center$V7),
                      V8 = pmax(histone_left$V8, histone_right$V8, histone_center$V8),
                      V9 = pmax(histone_left$V9, histone_right$V9, histone_center$V9),
                      V10 = pmax(histone_left$V10, histone_right$V10, histone_center$V10),
                      V11 = pmax(histone_left$V11, histone_right$V11, histone_center$V11),
                      V12 = pmax(histone_left$V12, histone_right$V12, histone_center$V12),
                      V13 = pmax(histone_left$V13, histone_right$V13, histone_center$V13))

atac <- read.delim("reference.bed", header=F)

histone_ids <- read.delim("metadata/histone_data_ids.txt", header=F)
atac_ids <- read.csv("metadata/contexts.csv", header=T)

histone <- histone[,-c(1:4)]
names(histone) <- histone_ids$V1

names(atac) <- c("chr", "start", "end", "id", atac_ids$context)
histone[histone == 0] <- 1e-10
histone <- log10(histone)

calc_z_scores <- function(column) {
  
  return((column - mean(column)) / sd(column))
  
}

histone <- as.data.frame(apply(histone, 2, calc_z_scores))
str(histone)
str(atac)
summary(histone)

data.frame(chr = atac$chr, 
           star = atac$start, 
           end = atac$end, 
           id = paste("cCRE_", seq(1, dim(atac)[1]), sep = ""), 
           AB = ifelse(test = histone$AB_H3K27ac >= 0.4 & atac$AB == 1, yes = 1, no = 0), 
           LB = ifelse(test = histone$LB_H3K27ac >= 0.4 & atac$LB == 1, yes = 1, no = 0), 
           E5 = ifelse(test = histone$E5_H3K27ac >= 0.4 & atac$E5 == 1, yes = 1, no = 0), 
           E11 = ifelse(test = histone$E11_H3K27ac >= 0.4 & atac$E11 == 1, yes = 1, no = 0), 
           E13 = ifelse(test = histone$E13_H3K27ac >= 0.4 & atac$E13 == 1, yes = 1, no = 0), 
           EAD = ifelse(test = histone$EAD_H3K27ac >= 0.4 & atac$EAD == 1, yes = 1, no = 0), 
           WID = ifelse(test = histone$WID_H3K27ac >= 0.4 & atac$WID == 1, yes = 1, no = 0)) -> cCREs

cCREs %>% filter(rowSums(cCREs[,-c(1:4)]) > 0) -> cCREs
colSums(cCREs[,-c(1:4)])

probs <- 0.7
data.frame(chr = atac$chr, 
           star = atac$start, 
           end = atac$end, 
           id = paste("cCRE_", seq(1, dim(atac)[1]), sep = ""), 
           AB = ifelse(test = (histone$AB_H3K27ac >= tapply(histone$AB_H3K27ac, atac$AB, quantile, probs = probs)[1]) 
                       & atac$AB == 1, yes = 1, no = 0), 
           E11 = ifelse(test = (histone$E11_H3K27ac >= tapply(histone$E11_H3K27ac, atac$E11, quantile, probs = probs)[1]) 
                        & atac$E11 == 1, yes = 1, no = 0), 
           E13 = ifelse(test = (histone$E13_H3K27ac >= tapply(histone$E13_H3K27ac, atac$E13, quantile, probs = probs)[1]) 
                        & atac$E13 == 1, yes = 1, no = 0), 
           E5 = ifelse(test = (histone$E5_H3K27ac >= tapply(histone$E5_H3K27ac, atac$E5, quantile, probs = probs)[1]) 
                       & atac$E5 == 1, yes = 1, no = 0),
           EAD = ifelse(test = (histone$EAD_H3K27ac >= tapply(histone$EAD_H3K27ac, atac$EAD, quantile, probs = probs)[1]) 
                        & atac$EAD == 1, yes = 1, no = 0), 
           HID = ifelse(test = (histone$HID_H3K27ac >= tapply(histone$HID_H3K27ac, atac$HID, quantile, probs = probs)[1]) 
                        & atac$HID == 1, yes = 1, no = 0), 
           LB = ifelse(test = (histone$LB_H3K27ac >= tapply(histone$LB_H3K27ac, atac$LB, quantile, probs = probs)[1]) 
                       & atac$LB == 1, yes = 1, no = 0),
           O = ifelse(test = (histone$O_H3K27ac >= tapply(histone$O_H3K27ac, atac$O, quantile, probs = probs)[1]) 
                      & atac$O == 1, yes = 1, no = 0),
           WID = ifelse(test = (histone$WID_H3K27ac >= tapply(histone$WID_H3K27ac, atac$WID, quantile, probs = probs)[1]) 
                        & atac$WID == 1, yes = 1, no = 0)) -> cCREs

cCREs %>% filter(rowSums(cCREs[,-c(1:4)]) > 0) -> cCREs
colSums(cCREs[,-c(1:4)])
write.table(cCREs, "cCREs.bed", quote=F, sep="\t", row.names = F, col.names=F)
