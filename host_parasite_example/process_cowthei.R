library(Seurat)
library(SeuratDisk)

#Read in Raw HP Data
#Host = Bovine
#Parasite = Theileria 

cowthei_raw = readRDS("cowthei_raw.rds")

#Reformat Gene Names + generate host_genes and parasite_genes

gene_names = rownames(cowthei_raw)

head(gene_names) # parasite genes
#[1] "Theileria-UTR-Tap370b08.q2ca38.01"  "Theileria-UTR-Tap370b08.q2ca38.02c" "Theileria-UTR-Tap370b08.q2ca38.03c"
#[4] "Theileria-UTR-TA16055"  

tail(gene_names) # host genes
#"Cow-----------ARSD" "Cow-----------GYG2" "Cow-----------XG" "Cow-----------CD99"              
#[5] "Cow-----------ENSBTAG00000049411" "Cow-----------ENSBTAG00000044391"

new_gene_names = vector()

host_genes = vector()
parasite_genes = vector()

for (x in gene_names){
  if (startsWith(x,"Cow")){ # if host
    
    # Reformat Gene Name
    genelist = as.list(strsplit(x, "-----------")[[1]])
    new_name = genelist[2]
    new_gene_names = c(new_gene_names,new_name)
    
    # Add to host_genes
    host_genes = c(host_genes,new_name)
    
  }else{ # if parasite
    
    #Reformat Gene name
    genelist = as.list(strsplit(x, "-UTR-")[[1]])
    new_name = genelist[2]
    new_gene_names = c(new_gene_names,new_name)
    
    # Add to parasite_genes
    parasite_genes = c(parasite_genes,new_name)
    
  }
}

new_gene_names = as.character(new_gene_names)

# Create New Seurat Object with updated Gene Names

hp_counts <- GetAssayData(cowthei_raw,assay = "RNA",slot = "counts")
rownames(hp_counts) = new_gene_names

cowthei_new = CreateSeuratObject(counts = hp_counts, meta.data = rawHP@meta.data)

#Add Host-Parasite Metadata

parasite_genes = as.character(parasite_genes)
host_genes = as.character(host_genes)

cowthei_new[["percent_parasite"]] <- PercentageFeatureSet(cowthei_new, features = parasite_genes)
cowthei_new[["parasiteUMI"]] <- PercentageFeatureSet(cowthei_new, features = parasite_genes) * cowthei_new$nCount_RNA/100
cowthei_new[["percent_host"]] <- PercentageFeatureSet(cowthei_new, features = host_genes)
cowthei_new[["hostUMI"]] <- PercentageFeatureSet(cowthei_new, features = host_genes) * cowthei_new$nCount_RNA/100

FeatureScatter(cowthei_new, feature1 = "percent.parasite", feature2 = "percent.host") #perfect correlation ,as expected

#Quality Control

VlnPlot(cowthei_new, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2) 

FeatureScatter(cowthei_new, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cowthei_new <- subset(cowthei_new, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500)

#Basic Processing/Analysis

cowthei_new = NormalizeData(cowthei_new)
cowthei_new = FindVariableFeatures(cowthei_new, selection.method = "vst", nfeatures = 2000)

cowthei_new = ScaleData(cowthei_new, verbose = FALSE)

cowthei_new <- RunPCA(cowthei_new, npcs = 50, verbose = FALSE)

cowthei_new <- RunUMAP(cowthei_new, reduction = "pca", dims = 1:40)
cowthei_new <- FindNeighbors(cowthei_new, reduction = "pca", dims = 1:40)
cowthei_new <- FindClusters(cowthei_new, resolution = 0.5)

DimPlot(cowthei_new, reduction = "umap")

#Reformat Metadata

cowthei_new$orig.ident = NULL
cowthei_new$RNA_snn_res.0.5 = NULL
cowthei_new$seurat_clusters = as.character(cowthei_new$seurat_clusters)

#Save Host-Parasite Gene Lists

write.csv(parasite_genes, "parasite_genes.csv", row.names = FALSE)
write.csv(host_genes, "host_genes.csv", row.names = FALSE)

#Slim Down Object

cowthei_new = DietSeurat(cowthei_new, counts = TRUE, data = TRUE, scale.data = FALSE, dimreducs = c("pca","umap"))

# Convert to annData

SaveH5Seurat(cowthei_new, filename = "cowthei_example.h5Seurat")
Convert("cowthei_example.h5Seurat", dest = "h5ad")



