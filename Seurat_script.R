# Single-cell Transcriptional Profile Data Analysis
# Author: Charlotte James
# Date: 4/21/2020

library(here)
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)

##########################################
# First analyze just the CD4 or CD8 cells
##########################################

# Read in data
d <- read.csv(here::here("singlecell_data.csv"), header=T, na.strings=c("NA", "null", "<NA>", ""))
dput(colnames(d)) 

d$wellID <- seq_len(nrow(d)) #add wellID for mapping in downstream visualization
d <- subset(d, d$wellcount >5000) #Remove wells that don't pass threshold read count (5000 reads)

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

#set all genes with >= 5 reads to 1, all else, set to 0 
temp <- temp %>% mutate_all( ~ ifelse( !is.na(.) & . >= 5, 1, 0))

#recombine with metadata 
subs <- d[,-(16:38)]
d <- cbind(subs, temp)

# Subset data frame to inclde only cells of interest 

CD4 <- subset(d, CD4 == 1)
CD8 <- subset(d, CD8 == 1)

alphapheno <- rbind(CD4, CD8)


#Add columns for TCR V genes of interest to include in analysis 
table(alphapheno$Va)
#Included due to high prevalence 
alphapheno$TRAV12 <- ifelse(grepl("TRAV12",alphapheno$Va),1,0)
alphapheno$TRAV1.2 <- ifelse(grepl("TRAV1-2",alphapheno$Va),1,0)
alphapheno$TRAV8 <- ifelse(grepl("TRAV8",alphapheno$Va),1,0)
#Included due to biological interest (LDN5-like and A01/A05-like)
alphapheno$TRAV17 <- ifelse(grepl("TRAV17", alphapheno$Va),1,0)
alphapheno$TRAV13 <- ifelse(grepl("TRAV13", alphapheno$Va),1,0)

alphapheno$Coreceptor <- alphapheno$CD4
alphapheno$Coreceptor[alphapheno$CD8 == 1] <- 'CD8'
alphapheno$Coreceptor[alphapheno$CD4 == 1] <- 'CD4'
alphapheno$Coreceptor[alphapheno$Coreceptor == 0] <- 'Null'

alphapheno$AgGroup <- alphapheno$CD4
alphapheno$AgGroup[alphapheno$Antigen == "SGL"] <- "1"
alphapheno$AgGroup[alphapheno$Antigen == "GMM"] <- '1'
alphapheno$AgGroup[alphapheno$Antigen == "Tet (-)"] <- '0'
alphapheno$AgGroup[alphapheno$Antigen == "Cell Line"] <- '2'

cols_to_analyze <- c("BCL6" ,"CD4" ,"CD8", "CTLA4", "EOMES" ,"FOXP3", "GATA3",            
                   "GZMB", "IFNG","IL10","IL12A","IL13","IL17A","IL2","MKI67","PDCD1","PRF1","RORC",
                   "RUNX1","RUNX3","TBET","TGFB1", "TNF", "TRAV1.2", "TRAV13", "TRAV17", "TRAV8", "TRAV12")             

################

# Prepare data for Seurat
d4tsne_withmeta <- alphapheno %>% 
  dplyr::select("Antigen", "PTID", "wellID", "AgGroup", "Coreceptor", all_of(cols_to_analyze)) 
d4tsne_numericonly <- d4tsne_withmeta %>% 
  dplyr::select(-Antigen, -PTID, -wellID, -AgGroup, -Coreceptor) %>% 
  data.matrix()

dim(d4tsne_withmeta)
dim(d4tsne_numericonly)

#Seurat Pipeline to determine which features are defining clusters 
d4_seurat <- t(d4tsne_numericonly) # data has features as rows, not columns
seuratOb <- CreateSeuratObject(counts = d4_seurat,
                               project = "GMM_SGLseurat", min.cells = 3, min.features = 0)
# min.cells = 3 --> includes features detected in at least 3 cells

seuratOb[["Coreceptor"]] <- d4tsne_withmeta$Coreceptor
seuratOb[["Antigen"]] <- d4tsne_withmeta$Antigen
seuratOb[["AgGroup"]] <- d4tsne_withmeta$AgGroup

# We run FindVariableFeatures as part of the pipeline, but all the features will be selected since there are so few
seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(seuratOb), points = head(VariableFeatures(seuratOb), 10), repel = TRUE)

VariableFeatures(object = seuratOb) 

# This is boolean data, so scaling is mostly a formality
seuratOb <- ScaleData(seuratOb, features = NULL)
seuratOb <- RunPCA(seuratOb, features = VariableFeatures(object = seuratOb))

# Cluster the cells
seuratOb <- FindNeighbors(seuratOb, dims = 1:5) # Create Shared Nearest Neighbor Graph
# Next run Louvain algorithm to group cells together
seuratOb <- FindClusters(seuratOb, resolution = 1, random.seed = 20200420) # lower res -> fewer clusters

seuratOb <- RunUMAP(seuratOb, dims = 1:3, seed.use = 20200420)
DimPlot(seuratOb, reduction = "umap")

seuratOb <- RunTSNE(seuratOb, dims = 1:3, seed.use = 20200420, check_duplicates = FALSE, perplexity = 20)
DimPlot(seuratOb, reduction = "tsne") 

# Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 0
cluster0.markers <- FindMarkers(seuratOb, ident.1 = 0, min.pct = 0.1)
cluster0.markers

# find all markers of cluster 1
cluster1.markers <- FindMarkers(seuratOb, ident.1 = 1, min.pct = 0.1)
# head(cluster0.markers, n = 10)
cluster1.markers

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seuratOb, ident.1 = 2, min.pct = 0.1)
# head(cluster0.markers, n = 10)
cluster2.markers

#Plots in Manuscript 
svg(here::here("CD4_CD8_Subset_Plots/umap.svg"))
DimPlot(seuratOb, reduction = "umap")
dev.off()

svg(here::here("CD4_CD8_Subset_Plots/CD4CD8.svg"))
FeaturePlot(seuratOb, reduction = "umap", features = c("CD4", "CD8"))
dev.off()

svg(here::here("CD4_CD8_Subset_Plots/txnfactor_2.svg"))
FeaturePlot(seuratOb, reduction = "umap", features = c("TBET", "GATA3", "RORC", "BCL6", "FOXP3", "CTLA4", "RUNX3"), pt.size = 2)
dev.off()

svg(here::here("CD4_CD8_Subset_Plots/ki67.svg"))
FeaturePlot(seuratOb, reduction = "umap", features = c("MKI67"), pt.size = 5)
dev.off()


svg(here::here("CD4_CD8_Subset_Plots/ki67facet.svg"))
FeaturePlot(seuratOb, reduction = "umap", features = c("MKI67"), pt.size = 3, split.by = c("Coreceptor") ) 
dev.off()

svg(here::here("CD4_CD8_Subset_Plots/txnfactor.svg"))
txnfactor <- select(alphapheno, c("TBET", "GATA3", "RORC", "RUNX3", "BCL6", "FOXP3", "CTLA4"))
txnfactor <- as.matrix(txnfactor)
Heatmap(txnfactor)
dev.off()

##########################################
# Seurat run with all cells, not just CD4 or CD8
##########################################

# Read in and format data
d <- read.csv(here::here("singlecell_data.csv"), header=T, na.strings=c("NA", "null", "<NA>", ""))
dput(colnames(d))

#add column with a wellID for plotting later 
d$wellID <- seq_len(nrow(d)) 

#filter wells based on read count (as defined by Han et al., Nature Biotechnology, 2015)
d <- subset(d, d$wellcount >5000)

#format and booleanize read counts for genes
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

temp <- temp %>% mutate_all( ~ ifelse( !is.na(.) & . >= 5, 1, 0))

subs <- d[,-(16:38)]
d <- cbind(subs, temp)

# no subsetting for this run. 
alphapheno <- d

#Add columns for TCR V genes of interest to include in t-SNE 
table(alphapheno$Va)
#Included due to high prevalence 
alphapheno$TRAV12 <- ifelse(grepl("TRAV12",alphapheno$Va),1,0)
alphapheno$TRAV1.2 <- ifelse(grepl("TRAV1-2",alphapheno$Va),1,0)
alphapheno$TRAV8 <- ifelse(grepl("TRAV8",alphapheno$Va),1,0)
#Included due to biological interest (LDN5-like and A01/A05-like)
alphapheno$TRAV17 <- ifelse(grepl("TRAV17", alphapheno$Va),1,0)
alphapheno$TRAV13 <- ifelse(grepl("TRAV13", alphapheno$Va),1,0)

alphapheno$Coreceptor <- alphapheno$CD4

alphapheno$Coreceptor[alphapheno$CD4 == 1] <- 'CD4'
alphapheno$Coreceptor[alphapheno$CD8 == 1] <- 'CD8'
alphapheno$Coreceptor[alphapheno$Coreceptor == 0] <- 'Null'

alphapheno$AgGroup <- alphapheno$CD4
alphapheno$AgGroup[alphapheno$Antigen == "SGL"] <- "1"
alphapheno$AgGroup[alphapheno$Antigen == "GMM"] <- '1'
alphapheno$AgGroup[alphapheno$Antigen == "Tet (-)"] <- '0'
alphapheno$AgGroup[alphapheno$Antigen == "Cell Line"] <- '2'

cols_to_analyze <- c("BCL6" ,"CD4" ,"CD8", "CTLA4", "EOMES" ,"FOXP3", "GATA3",            
                   "GZMB", "IFNG","IL10","IL12A","IL13","IL17A","IL2","MKI67","PDCD1","PRF1","RORC",
                   "RUNX1","RUNX3","TBET","TGFB1", "TNF", "TRAV1.2", "TRAV13", "TRAV17", "TRAV8", "TRAV12")    

################

# Prepare data for Seurat
d4tsne_withmeta <- alphapheno %>% 
  dplyr::select("Antigen", "PTID", "wellID", "Coreceptor","AgGroup", all_of(cols_to_analyze)) 
d4tsne_numericonly <- d4tsne_withmeta %>% 
  dplyr::select(-Antigen, -PTID, -wellID, -Coreceptor, -AgGroup) %>% 
  data.matrix()

#Check to ensure data was prepared correctly
dim(d4tsne_withmeta)
dim(d4tsne_numericonly)

#Seurat Pipeline to determine which features are defining clusters 

d4_seurat <- t(d4tsne_numericonly) # counts data has features as rows, not columns
seuratOb <- CreateSeuratObject(counts = d4_seurat,
                               project = "GMM_SGLseurat", min.cells = 3, min.features = 0)

d4tsne_withmeta$AgGroup <- as.factor(d4tsne_withmeta$AgGroup)
seuratOb[["Coreceptor"]] <- d4tsne_withmeta$Coreceptor
seuratOb[["Antigen"]] <- d4tsne_withmeta$Antigen
seuratOb[["AgGroup"]] <- d4tsne_withmeta$AgGroup

seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures = 2000)
LabelPoints(plot = VariableFeaturePlot(seuratOb), points = head(VariableFeatures(seuratOb), 10), repel = TRUE)

VariableFeatures(object = seuratOb)

seuratOb <- ScaleData(seuratOb, features = NULL)
seuratOb<- RunPCA(seuratOb, features = VariableFeatures(object = seuratOb))

# Now cluster the cells
seuratOb <- FindNeighbors(seuratOb, dims = 1:5)
seuratOb <- FindClusters(seuratOb, resolution = 1.2, random.seed = 20200420) # lower res -> fewer clusters

#Plot cell clusters
seuratOb <- RunUMAP(seuratOb, dims = 1:5, seed.use = 20200420)

svg(here::here("All_Cells_Plots/20200614_allcells_seurat.svg"))
DimPlot(seuratOb, reduction = "umap", pt.size = 3)
dev.off()

svg(here::here("All_Cells_Plots/facetptid.svg"))
DimPlot(seuratOb, reduction = "umap", pt.size = 3) + facet_grid(. ~ d4tsne_withmeta[match(colnames(seuratOb), alphapheno$wellID), "PTID"]) 
dev.off()

svg(here::here("All_Cells_Plots/facetantigen.svg"))
DimPlot(seuratOb, reduction = "umap", pt.size = 3) + facet_grid(. ~ d4tsne_withmeta[match(colnames(seuratOb), alphapheno$wellID), "Antigen"]) 
dev.off()

svg(here::here("All_Cells_Plots/20200614_co-receptor_allcells.svg"))
FeaturePlot(seuratOb, features = c("CD4", "CD8"), reduction = "umap", pt.size = 3)
dev.off()

svg(here::here("All_Cells_Plots/facetKI67.svg"))
FeaturePlot(seuratOb, feature = c("MKI67"), reduction = "umap", pt.size = 3, split.by = c("AgGroup"), shape.by  = c("Coreceptor"))
#  facet_wrap(d4tsne_withmeta[match(colnames(seuratOb), alphapheno$wellID), "Coreceptor"]~.)
dev.off()
