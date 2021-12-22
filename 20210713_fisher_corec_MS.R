#Formatting data for testing 
library(dplyr)
library(tidyverse)
#read in
d <- read.csv("~/James_SupplementalFile_1.csv")

#clean
d$wellID <- seq_len(nrow(d)) #add wellID for mapping in downstream visualization
d$wellcount.sum <- (d$wellcount + d$wellcount.pheno) #sum read count for QC
d <- subset(d, d$wellcount.sum >5000) #Remove wells that don't pass threshold read count (5000 reads)

#Clean and format data table 
d[d==""] <- NA
temp <- d[,(16:38)]
temp$BCL6 <- sub(".*? (.+)", "\\1", temp$BCL6)
temp$GATA3 <- sub(".*? (.+)", "\\1", temp$GATA3)
temp$GZMB <- sub(".*? (.+)", "\\1", temp$GZMB)
temp$IFNG <- sub(".*? (.+)", "\\1", temp$IFNG)
temp$IL10 <- sub(".*? (.+)", "\\1", temp$IL10)
temp$IL12A <- sub(".*? (.+)", "\\1", temp$IL12A)
temp$IL13 <- sub(".*? (.+)", "\\1", temp$IL13)
temp$IL2 <- sub(".*? (.+)", "\\1", temp$IL2)
temp$PRF1 <- sub(".*? (.+)", "\\1", temp$PRF1)
temp$RORC <- sub(".*? (.+)", "\\1", temp$RORC)
temp$RUNX1 <- sub(".*? (.+)", "\\1", temp$RUNX1)
temp$RUNX3 <- sub(".*? (.+)", "\\1", temp$RUNX3)
temp$TBET <- sub(".*? (.+)", "\\1", temp$TBET)
temp$TGFB1 <- sub(".*? (.+)", "\\1", temp$TGFB1)
temp$TNF <- sub(".*? (.+)", "\\1", temp$TNF)
temp$CD4 <- sub(".*? (.+)", "\\1", temp$CD4)
temp$CD8 <- sub(".*? (.+)", "\\1", temp$CD8)
temp$EOMES <- sub(".*? (.+)", "\\1", temp$EOMES)
temp$CTLA4 <- sub(".*? (.+)", "\\1", temp$CTLA4)
temp$MKI67 <- sub(".*? (.+)", "\\1", temp$MKI67)
temp$PDCD1 <- sub(".*? (.+)", "\\1", temp$PDCD1)
temp$FOXP3 <- sub(".*? (.+)", "\\1", temp$FOXP3)


temp$FOXP3 <- as.numeric(temp$FOXP3)
temp$BCL6 <- as.numeric(temp$BCL6)
temp$GATA3 <- as.numeric(temp$GATA3)
temp$GZMB <- as.numeric(temp$GZMB)
temp$IFNG <- as.numeric(temp$IFNG)
temp$IL10 <- as.numeric(temp$IL10)
temp$IL12A <- as.numeric(temp$IL12A)
temp$IL13 <- as.numeric(temp$IL13)
temp$IL2 <- as.numeric(temp$IL2)
temp$PRF1 <- as.numeric(temp$PRF1)
temp$RORC <- as.numeric(temp$RORC)
temp$RUNX1 <- as.numeric(temp$RUNX1)
temp$RUNX3 <- as.numeric(temp$RUNX3)
temp$TBET <- as.numeric(temp$TBET)
temp$TGFB1 <- as.numeric(temp$TGFB1)
temp$TNF <- as.numeric(temp$TNF)
temp$CD8 <- as.numeric(temp$CD8)
temp$CD4 <- as.numeric(temp$CD4)
temp$EOMES <- as.numeric(temp$EOMES)
temp$CTLA4 <- as.numeric(temp$CTLA4)
temp$MKI67 <- as.numeric(temp$MKI67)
temp$PDCD1 <- as.numeric(temp$PDCD1)

#set all genes with >5 reads to 1, all else, set to 0 
temp[temp >= 5] <- 1
temp[is.na(temp)] <- 0

#recombine with metadata 
subs <- d[,-(16:38)]
d <- cbind(subs, temp)

head(d)

#subset to only include cells of interest where CD4 or CD8 is annotated
d <-subset(d, CD4==1 | CD8 == 1)


# where both CD4 and CD8 are called, set to CD4 (biological interpretation)
d$Coreceptor <- d$CD4
d$Coreceptor[d$CD8 == 1] <- 'CD8'
d$Coreceptor[d$CD4 == 1] <- 'CD4'

#Running Tests and correcting for multiple comparisons 
# move to matrix
cytokine_cols <- c("PTID", "Antigen", "BCL6" ,"CD4" ,"CD8", "CTLA4" , "EOMES" ,"FOXP3", "GATA3",            
                   "GZMB", "IFNG","IL10","IL12A","IL13","IL17A","IL2","MKI67","PDCD1","PRF1","RORC",
                   "RUNX1","RUNX3","TBET","TGFB1", "TNF") 



d$AgGroup <- d$CD4
d$AgGroup[d$Antigen == "SGL"] <- "1"
d$AgGroup[d$Antigen == "GMM"] <- '1'
d$AgGroup[d$Antigen == "Tet (-)"] <- '0'
d$AgGroup[d$Antigen == "Cell Line"] <- '2'

d4tsne_withmeta <- d %>% 
  dplyr::select(c("Antigen", "PTID", "wellID", "AgGroup", "Coreceptor", cytokine_cols)) 
d4tsne_numericonly <- d4tsne_withmeta %>% 
  dplyr::select(-Antigen, -PTID, -AgGroup, -wellID, -CD4, -CD8)


dim(d4tsne_withmeta)
dim(d4tsne_numericonly)

head(d)

#fisher tests
pval <- d4tsne_numericonly %>% 
  gather(gene, expression, -c(Coreceptor)) %>%
  group_by(Coreceptor, gene) %>%
  summarise(NE = n() - sum(expression),  E = n() - NE) %>%
  ungroup() %>%
  select(-Coreceptor) %>%
  nest(data = c(NE, E)) %>%
  mutate(p_value = data %>%
           map(~fisher.test(.)) %>%
           map_dbl(pluck, "p.value")) %>%
  select(-data)

pval$p_adj <- p.adjust(pval$p_value, method = "BH")
print(pval)

#Plots for manuscript 
  
  d4tsne_numericonly <- d4tsne_withmeta %>% 
    dplyr::select(-Antigen, -PTID, -AgGroup, -wellID, -Coreceptor)
  
  d4tsne_numericonly <- d4tsne_numericonly[c("CD4", "CD8", "BCL6", "CTLA4" , "EOMES" ,"FOXP3", "GATA3",  "TBET",  "RORC",
                                             "RUNX1","RUNX3","GZMB", "IFNG","IL10","IL12A","IL13","IL17A","IL2","MKI67",
                                             "PDCD1","PRF1","TGFB1", "TNF")]

  
# distance & hierarchical clustering
mtx <- as.matrix(d4tsne_numericonly)
dist <- dist(d4tsne_numericonly, method = "binary")
cairo_pdf('~/Desktop/20210802_Figure4A.pdf')
heatmap(mtx, distfun = function (x) dist(x, method = "binary"),  Colv = NA, scale = "none")
dev.off()

cairo_pdf('~/Desktop/20210927_Figure4A_xy_clust.pdf')
heatmap(mtx, distfun = function (x) dist(x, method = "binary"), scale = "none")
dev.off()

mtxCD4 <- mtx[mtx[, "CD4"] == 1,]
cairo_pdf('~/Desktop/20210802_Figure4_CD4only.pdf')
heatmap(mtxCD4, distfun = function (x) dist(x, method = "binary"), Colv = NA, scale = "none")
dev.off()

mtxCD8 <- mtx[mtx[, "CD8"] == 1,]
cairo_pdf('~/Desktop/20210802_Figure4_CD8only.pdf')
heatmap(mtxCD8, distfun = function (x) dist(x, method = "binary"), Colv = NA, scale = "none")
dev.off()

cairo_pdf('~/Desktop/20210927_Figure4A_xy_clust.pdf')
heatmap(mtx, distfun = function (x) dist(x, method = "binary"), scale = "none")
dev.off()

txnfactor <- d4tsne_numericonly[c("BCL6", "EOMES", "CTLA4","FOXP3", "GATA3",  "TBET",  "RORC",
                                           "RUNX1","RUNX3")]
txnmtx <- as.matrix(txnfactor)
cairo_pdf('~/Desktop/20210802_Figure4B.pdf')
heatmap(txnmtx, distfun = function (x) dist(x, method = "binary"),  Colv = NA, scale = "none")
dev.off()
  
#both saved as cairo.pdf using R interface
  
temp <- d4tsne_numericonly %>% 
  gather(gene, expression, -c(Coreceptor)) %>%
  group_by(Coreceptor, gene) %>%
  summarise(NE = n() - sum(expression),  E = n() - NE) %>%
  ungroup() %>%
  select(-Coreceptor) %>%
  nest(data = c(NE, E)) %>%
  filter(gene%in% c("BCL6", "EOMES", "PRF1","TBET"))
print(temp$data)

#[[1]]
#BCL6
#NE     E
#<dbl> <dbl>
 # 1    37    14  
#   2    51     2

14/(37+14) #27.4 CD4 
2/(51+2) #3.8% CD8

BCL6 <- data.frame(per=c(27.4,  3.8),
                 group=c("CD4", "CD8"))

svg("~/Desktop/20210802_Figure4F.svg")
ggplot(BCL6, aes(group, per, fill = group)) + geom_bar(stat = "identity") +  scale_color_manual(c("black")) +
  theme_classic() + scale_fill_manual(values=c('black','white')) + labs(x = "Coreceptor", y = "% of Cells")
dev.off()

#[[2]]
#EOMES
#NE     E
#<dbl> <dbl>
#  1    38    13
#2    27    26

13/(38+13) #25.5% CD4
26/(27+26) #49.0% CD8

EOMES <- data.frame(per=c(25.5,  49.0),
                   group=c("CD4", "CD8"))

svg("~/Desktop/20210802_Figure4E.svg")
ggplot(EOMES, aes(group, per, fill = group)) + geom_bar(stat = "identity") +  scale_color_manual(c("black")) +
  theme_classic() + scale_fill_manual(values=c('black','white')) + labs(x = "Coreceptor", y = "% of Cells")
dev.off()

#[[3]]
#PRF1 
#NE     E
#<dbl> <dbl>
#  1    36    15
#2    23    30

15/(36+15) #29.4% CD4
30/(23+30) #56.5% CD8

PRF1 <- data.frame(per=c(29.4,  56.5),
                    group=c("CD4", "CD8"))

svg("~/Desktop/20210802_Figure4D.svg")
ggplot(PRF1, aes(group, per, fill = group)) + geom_bar(stat = "identity") +  scale_color_manual(c("black")) +
  theme_classic() + scale_fill_manual(values=c('black','white')) + labs(x = "Coreceptor", y = "% of Cells")
dev.off()

#[[4]]
#TBET
#NE     E
#<dbl> <dbl>
#  1    29    22
#2    14    39

22/(29+22) #43.1% CD4
39/(14+39) #73.6% CD8

TBET <- data.frame(per=c(43.1,  73.6),
                   group=c("CD4", "CD8"))

svg("~/Desktop/20210802_Figure4C.svg")
ggplot(TBET, aes(group, per, fill = group)) + geom_bar(stat = "identity") +  scale_color_manual(c("black")) +
  theme_classic() + scale_fill_manual(values=c('black','white')) + labs(x = "Coreceptor", y = "% of Cells")
dev.off()

