###############       Libraries       ###############

library(tidyverse)
library(ggpubr)
library(viridis)
library(reshape)
library(ComplexHeatmap)
library(circlize)

#####################################################
###############      data loading     ###############
#####################################################
# TSSs class
TSS_Class <- as_tibble(read.csv("../1_Files/15.stable.variable.tsv", sep = '\t'))
TSS_Class$pnl <- gsub("stable \\(566\\)", "stable", TSS_Class$pnl)
TSS_Class$pnl <- gsub("variable \\(623\\)", "variable", TSS_Class$pnl)
TSS_Class$gid <- as.character(TSS_Class$gid)

# Genes-Transcripts
GenTran <- as_tibble(read.table("Zm_Gene_Transcript.txt", h=F))
colnames(GenTran) <- c("gid", "type", "t")
GenTran[,1:3] <- apply(GenTran[,1:3], 2, as.character)

# 
TranscriptByGeneFreq <- as_tibble(read.table("Zm_Freq.Transcript.by.Gene.txt", h=F))
TranscriptByGeneFreq$V1 <- as.character(TranscriptByGeneFreq$V1)

TranscriptByGeneFreq$V2[(TranscriptByGeneFreq$V2 >50)==TRUE] <- 50

TranscriptByGeneFreq_Interest <-  subset(TranscriptByGeneFreq, V1 %in% TSS_Class$gid)
TranscriptByGeneFreq_Interest <- left_join(TranscriptByGeneFreq_Interest, 
                                           TSS_Class[,1:2], by=c("V1"= "gid"))

# Transcript vs TSS distance
Exon2 <- as_tibble(read.table("Zm_Exon2.Start_by_Transcript.txt", h=F))
colnames(Exon2) <- c('t', 'Pos.t', 'strand')
TSS.annotated <- as_tibble(read.table("Zm_Annotated_TSS.txt", h=F))
colnames(TSS.annotated) <- c('gid', 'Pos.tss', 'strand')



#####################################################
###########    Processing data        ###############
#####################################################

TranscriptFreq <- TranscriptByGeneFreq$V2

Exon2  <- left_join(Exon2, GenTran, by='t')
Exon2 <- left_join(Exon2 , TSS.annotated[,1:2], by='gid')
Exon2["Dis"] <- abs(Exon2$Pos.tss - Exon2$Pos.t)

# add classes
Exon2["Class"] <- Exon2$gid %in%  TSS_Class$gid

Exon2$Class <- gsub("FALSE", "Expressed", Exon2$Class)

Exon2$Class[Exon2$gid %in% subset(TSS_Class, pnl=="stable")$gid] <- "stable"
Exon2$Class[Exon2$gid %in% subset(TSS_Class, pnl=="variable")$gid] <- "variable"

table(unique(Exon2[,c(4,8)])$Class)
table(TSS_Class$pnl)

#gghistogram(Exon2, x = "Dis", y="..density..", add = "mean",  color = "Class")
Exon_TSS_dis <- gghistogram(subset(Exon2, Class !='Expressed'), 
            x = "Dis", y="..density..", add = "mean",  color = "Class")

Exon_TSS_dis <- ggpar(Exon_TSS_dis, xlab = "Dis. TSS - Start.Exon2 (bps)")

Exon_TSS_dis_violin <- ggplot(subset(Exon2, Class !='Expressed'),
                              aes(x=Class, y=Dis, fill=Class)) +
  geom_violin() + stat_compare_means() +
  theme_pubr()# +

Exon_TSS_dis_violin <- ggpar(Exon_TSS_dis_violin, ylab = "Dis. TSS - Start.Exon2 (bps)", legend = "none")
ggarrange(Exon_TSS_dis, Exon_TSS_dis_violin,  ncol = 2)

####################  ####################

N.trans.by.gene.plot <- ggplot(TranscriptByGeneFreq_Interest, aes(x=pnl, y=V2, fill=pnl)) +
  geom_violin() +
  stat_compare_means() +
  theme_pubr()# +
  # stat_compare_means(method = "t.test", label.y = 150)

N.trans.by.gene.plot <- ggpar(N.trans.by.gene.plot,  xlab = "TSS class", 
                             ylab = "# of transcript", legend = "none")
N.trans.by.gene.plot



TSS_Class.matrix.cor <- left_join(TSS_Class, TranscriptFreq, by=c("gid"="V1"))[,-c(1,2)]
TSS_Class.matrix.cor <- TSS_Class.matrix.cor[rowSums(is.na(TSS_Class.matrix.cor)*1) == 0,]

TSS_Class.matrix.cor$V2 <- scale(TSS_Class.matrix.cor$V2)

plot(TSS_Class.matrix.cor$V2 ~ TSS_Class.matrix.cor$SS_Hs_Ms)

TSS_Class.matrix.cor <- cor(as.matrix(TSS_Class.matrix.cor))
TSS_Class.matrix.cor

#
my_palette <- colorRamp2(seq(-1, 1, 0.1), viridis(21, direction = 1, option = "A"))

Plot <- Heatmap(TSS_Class.matrix.cor, 
                #heatmap_width = unit(9, "cm"),
                name = "log2(Exp)", 
                col=my_palette,
                cluster_columns = TRUE, show_column_dend = FALSE,
                show_row_dend = TRUE, column_dend_reorder = TRUE,
                column_names_gp = gpar(fontsize = 10),
                row_names_gp = gpar(fontsize = 10),
                show_heatmap_legend = T)  
Plot


