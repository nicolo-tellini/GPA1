# Mon Dec  2 15:55:49 2024 

# Title:
# Author: Nicolò T.
# Status: Draft

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)
gc()
gcinfo(FALSE)
options(scipen=999)

# Variables ----

baseDir <- '/home/ntellini/proj/collaborations/chiara-GPA1'
# allDir <- list.dirs()
# allFiles <- list.files()
setwd(baseDir)

# Libraries ----

library(seqinr)

# body ----


# these code below does the reverse complement 
vcf <- fread('grep -v "##" GPA1.codonAA469.vcf',data.table = F)

vcf <- vcf[-2,]  

vcf <- vcf[nrow(vcf):1,]

vcf[,4] <- chartr("ACGT","TGCA",vcf[,4])

vcf[,5] <- chartr("ACGT","TGCA",vcf[,5])

write.table(vcf,file = "GPA1.codonAA469.revcompl.vcf",append = F,quote = F,sep = "\t",col.names = T,row.names = F)

# these code below take the fasta codon and count them
AA82 <- read.fasta("GPA1.codonAA82.revcompl.min0.fasta",seqtype = "DNA",as.string = T,forceDNAtolower = F,set.attributes = F)

AA469 <- read.fasta("GPA1.codonAA469.revcompl.min0.fasta",seqtype = "DNA",as.string = T,forceDNAtolower = F,set.attributes = F)

AA82 <- AA82[grep("scer-",names(AA82),invert = T)]

AA469 <- AA469[grep("scer-",names(AA469),invert = T)]

result <- mapply(function(x, y) paste(x, y, sep = "_"), AA82, AA469, SIMPLIFY = T)  

df <- as.data.frame(grep("N",invert = T,table(result),value = T))
df$couples <- rownames(df)

write.table(x = df,file = "codonaa82_codonaa469.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

AA82_YGG <- as.data.frame(t(as.data.frame(AA82[AA82[] == "YGG"])))
AA82_YGG$name <- rownames(AA82_YGG)

AA469_AKT <- as.data.frame(t(as.data.frame(AA469[AA469[] == "AKT"])))
AA469_AKT$name <- rownames(AA469_AKT)

ven_d <- gplots::venn(data = list(
  "AA82_YGG" = AA82_YGG$name,
  "AA469_AKT" = AA469_AKT$name
))

attr_vend <- attributes(ven_d)

list_strains <- as.data.frame(names(AA469))

write.table(x = list_strains,file = "list_scer_strains_wo_sgrp5.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

### part II ----

# Mon Dec  2 15:55:49 2024 

# Title:
# Author: Nicolò T.
# Status: Draft

# Comments:

# Options ----

rm(list = ls())
options(warn = 1)
options(stringsAsFactors = F)
gc()
gcinfo(FALSE)
options(scipen=999)

# Variables ----

baseDir <- '/home/ntellini/proj/collaborations/chiara-GPA1'
# allDir <- list.dirs()
# allFiles <- list.files()
setwd(baseDir)

# Libraries ----

library(seqinr)

# body ----


# these code below does the reverse complement 
vcf <- fread('grep -v "##" GPA1.codonAA469.vcf',data.table = F)

vcf <- vcf[-2,]  

vcf <- vcf[nrow(vcf):1,]

vcf[,4] <- chartr("ACGT","TGCA",vcf[,4])

vcf[,5] <- chartr("ACGT","TGCA",vcf[,5])

write.table(vcf,file = "GPA1.codonAA469.revcompl.vcf",append = F,quote = F,sep = "\t",col.names = T,row.names = F)

# these code below take the fasta codon and count them
AA82 <- read.fasta("GPA1.codonAA82.revcompl.min0.fasta",seqtype = "DNA",as.string = T,forceDNAtolower = F,set.attributes = F)

AA469 <- read.fasta("GPA1.codonAA469.revcompl.min0.fasta",seqtype = "DNA",as.string = T,forceDNAtolower = F,set.attributes = F)

AA82 <- AA82[grep("scer-",names(AA82),invert = T)]

AA469 <- AA469[grep("scer-",names(AA469),invert = T)]


aa82df <- as.data.frame(do.call(rbind,AA82))
AA469df <-  as.data.frame(do.call(rbind,AA469))

aa82df$names <- rownames(aa82df)
AA469df$names <- rownames(AA469df)

aa_all <- base::merge(aa82df,AA469df,by="names",all.x = T)

aa_all[aa_all[,2] == "TGG","aa82"] <- "W"
aa_all[aa_all[,2] == "CGG","aa82"] <- "R"
aa_all[aa_all[,2] == "YGG","aa82"] <- "W-R"
aa_all[aa_all[,2] == "NGG","aa82"] <- "-"
aa_all[aa_all[,2] == "TGN","aa82"] <- "-"

aa_all[aa_all[,3] == "AGT","aa469"] <- "S"
aa_all[aa_all[,3] == "ATT","aa469"] <- "I"
aa_all[aa_all[,3] == "AKT","aa469"] <- "S-I"
aa_all[aa_all[,3] == "ANT","aa469"] <- "-"
aa_all[aa_all[,3] == "NNN","aa469"] <- "-"

aa_all$combination <- paste0(aa_all$aa82,"_",aa_all$aa469)

df_clades_1011 <- fread("1011-clades",data.table = F)
df_clades_pub <- fread("pub-clades",data.table = F)

colnames(df_clades_1011) <- c("str","cld")
colnames(df_clades_pub) <- c("str","cld")

clades <- rbind(df_clades_1011,df_clades_pub)

aa_all_clade <- base::merge(aa_all,clades,by.x = "names",by.y ="str",all.x = T)
aa_all_clade[is.na(aa_all_clade$cld),"cld"] <- ""

df_final <- as.data.frame(table(aa_all_clade$combination,aa_all_clade$cld))        
df_final$Var2 <- as.character(df_final$Var2)

for (i in 1:length(unique(df_final$Var2))) {
  df_final[df_final[,"Var2"] == unique(df_final$Var2)[i],"fraction"] <- round(df_final[df_final[,"Var2"] == unique(df_final$Var2)[i],"Freq"] / sum(df_final[df_final[,"Var2"] == unique(df_final$Var2)[i],"Freq"] ),digits = 3)
}
colnames(df_final) <- c("combination aa82_aa469","clade/grouping","counts","frequencies")
colnames(aa_all_clade) <- c("strain name/sequencing code","codon aa82","codon aa469","aa82","aa469","combination aa82_aa469","clade")
write.table(df_final,file = "summary_table_perclade_GPA1_aa82aa469.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
write.table(aa_all_clade,file = "summary_table_perstrain_GPA1_aa82aa469.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
