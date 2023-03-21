library(Seurat)
setwd('/Users/johnny/Documents/singlecell')

data_dir <- '/Users/johnny/Downloads/sc_result/filtered_feature_bc_matrix'
data <- Read10X(data.dir = data_dir)
seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
levels(x = seurat_object)

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
head(seurat_object@meta.data, 5)

VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 5)

seurat_object <- NormalizeData(seurat_object)

seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

DimHeatmap(seurat_object, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(seurat_object, ndims=50)

seurat_object <- FindNeighbors(seurat_object, dims = 1:40)
seurat_object <- FindClusters(seurat_object, resolution = 0.3)

seurat_object <- RunUMAP(seurat_object, dims = 1:40)
DimPlot(seurat_object, reduction = "umap", label = TRUE)

FeaturePlot(seurat_object, features = c("SERPINE1", "MCM5", "MT1E"))




see <- RunUMAP(seurat_object, dims = 1:10)

write.table(as.matrix(GetAssayData(object = seurat_object, slot = "data")), 
            '/Users/johnny/Downloads/r_sc_export_norm.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)




save_seurat(seurat_object, prefix = "unfiltered", proj_dir = getwd())
saveRDS(seurat_object, "sc_raw_count.rds")

#HTO tag process
see <- read.table("raw_counts_no_hto.csv", header=TRUE, row.names = 1, sep=",")

sc_data.umis <- read.table("raw_counts_no_hto_symbol.csv", header=TRUE, row.names = 1, sep=",")
sc_data.htos <- read.table("hto_counts.csv", header=TRUE, row.names = 1, sep=",")

joint.bcs <- intersect(colnames(sc_data.umis), colnames(sc_data.htos))

sc_data.umis <- sc_data.umis[, joint.bcs]
sc_data.htos <- as.matrix(sc_data.htos[, joint.bcs])
rownames(sc_data.htos)

sc_data.hashtag <- CreateSeuratObject(counts = sc_data.umis)

# Normalize RNA data with log normalization
sc_data.hashtag <- NormalizeData(sc_data.hashtag)
# Find and scale variable features
sc_data.hashtag <- FindVariableFeatures(sc_data.hashtag, selection.method = "mean.var.plot")
sc_data.hashtag <- ScaleData(sc_data.hashtag, features = VariableFeatures(sc_data.hashtag))

# Add HTO data as a new assay independent from RNA
sc_data.hashtag[["HTO"]] <- CreateAssayObject(counts = sc_data.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
sc_data.hashtag <- NormalizeData(sc_data.hashtag, assay = "HTO", normalization.method = "CLR")

sc_data.hashtag <- HTODemux(sc_data.hashtag, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(sc_data.hashtag$HTO_classification.global)


# Group cells based on the max HTO signal
Idents(sc_data.hashtag) <- "HTO_maxID"
RidgePlot(sc_data.hashtag, assay = "HTO", features = rownames(sc_data.hashtag[["HTO"]])[1:2], ncol = 2)

final_tags=sc_data.hashtag$HTO_classification


write.table(as.matrix(sc_data.hashtag$HTO_classification), 
            '/Users/johnny/Downloads/tags_final.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)





Idents(sc_data.hashtag) <- "HTO_classification.global"
# Extract the singlets
sc_data.singlet <- subset(sc_data.hashtag, idents = "Singlet")

#sc_data.singlet=subset(x = sc_data.singlet, subset = (HTO_classification == "TH-2" | HTO_classification == "TH-4"))

# Select the top 1000 most variable features
sc_data.singlet <- FindVariableFeatures(sc_data.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
sc_data.singlet <- ScaleData(sc_data.singlet, features = VariableFeatures(sc_data.singlet))

# Run PCA
sc_data.singlet <- RunPCA(sc_data.singlet, features = VariableFeatures(sc_data.singlet))

sc_data.singlet <- FindNeighbors(sc_data.singlet, reduction = "pca", dims = 1:30)#30
sc_data.singlet <- FindClusters(sc_data.singlet, resolution = 0.6) #0.6



#sc_data.singlet <- RunTSNE(sc_data.singlet, reduction = "pca", dims = 1:30)
sc_data.singlet <- RunUMAP(sc_data.singlet, reduction = "pca", dims = 1:30)

# Projecting singlet identities on TSNE visualization
DimPlot(sc_data.singlet, group.by = "HTO_classification", reduction = "umap")

FeaturePlot(sc_data.singlet, features = c('KLF4','MYC', 'EZH2'))
FeaturePlot(sc_data.singlet, features = 'MYC', cols=c("lightgrey", "#ff0000"), min.cutoff = 'q20', max.cutoff = 'q80')
FeaturePlot(sc_data.singlet, features = 'KLF4', cols=c("lightgrey", "#ff0000"), min.cutoff = 'q10', max.cutoff = 'q90')
FeaturePlot(sc_data.singlet, features = 'EZH2', cols=c("lightgrey", "#ff0000"), min.cutoff = 'q10', max.cutoff = 'q90')






#code for DE genes between: 12,13,24

markers1 <- FindMarkers(sc_data.singlet, ident.1 = "TH-1", ident.2 = "TH-2", group.by = 'HTO_classification')
markers2 <- FindMarkers(sc_data.singlet, ident.1 = "TH-1", ident.2 = "TH-3", group.by = 'HTO_classification')
markers3 <- FindMarkers(sc_data.singlet, ident.1 = "TH-2", ident.2 = "TH-4", group.by = 'HTO_classification')
head(x = markers)
write.table(as.matrix(markers3), 
            '/Users/johnny/Downloads/DE_sc_2v4.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)







out=GetAssayData(object = sc_data.singlet, slot = "data")
write.table(as.matrix(out), 
            '/Users/johnny/Downloads/singlet_data_for_stemness_out.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

data_in=read.table("singlet_data_for_stemness_in.csv", header=TRUE, row.names = 1, sep=",")

new_ob <- SetAssayData(object = sc_data.singlet, slot = "scale.data", new.data = as.matrix(data_in))
#new_ob =subset(x = new_ob, subset = (HTO_classification == "TH-2" | HTO_classification == "TH-4"))

FeaturePlot(sc_data.singlet, features = 'LOC102724788', cols=c("lightgrey", "#ff0000"), min.cutoff = 'q10', max.cutoff = 'q90', slot="data")
FeaturePlot(new_ob, features = 'LOC102724788', cols=c("lightgrey", "#ff0000"), min.cutoff = 'q01', max.cutoff = 'q10', slot="scale.data")





#