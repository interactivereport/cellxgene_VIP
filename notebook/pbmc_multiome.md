PBMC10k\_multiome data prepare for cellxgene VIP
================

To generate a multiome instance in cellxgene VIP, three items are required. You can generate them by following this notebook.

+ One h5ad file. Please follows **part 1**.
+ Three .rds file (annotation.rds, links.rds, peaks.rds). Please follows **part 1**.
+ bigwig files (.bw files) and one index file (named bw.cluster). Please follows **part 2**.


## Part 1 
Following R code is based on
<https://satijalab.org/signac/articles/pbmc_multiomic.html>

``` r
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)
```

Download the pbmc multiome data from
<https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0>

or by downloading use curl in shell

```sh
# Input Files

curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_fastqs.tar
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_library.csv

# Output Files
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_cloupe.cloupe
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_analysis.tar.gz
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_web_summary.html
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_summary.csv
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.tar.gz
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_raw_feature_bc_matrix.h5
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_possorted_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_possorted_bam.bam.bai
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_gex_molecule_info.h5
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam.bai
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peaks.bed
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_cut_sites.bigwig
curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_peak_annotation.tsv


```
Requires:

1.  Filtered feature barcode matrix (H5)

2.  ATAC Per fragment information file (TSV.GZ)

3.  ATAC Per fragment information index (TSV.GZ.index)

4.  ATAC Position-sorted alignments (BAM) - for BIGWIG

5.  ATAC Position-sorted alignments (BAM index) - for BIGWIG

``` r
counts <- Read10X_h5("../pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "../pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
```

``` r
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
```

Quality Control

``` r
DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
```

``` r
# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc
```

Peak calling

``` r
# call peaks using MACS
peaks <- CallPeaks(pbmc, macs2.path ='../opt/anaconda3/bin/macs3')

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
```

Gene expression data processing

``` r
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
```

DNA accessibility data processing

``` r
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
```

Annotating cell types PBMC reference downloaded from:
<https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat>

``` r
library(SeuratDisk)
rm(counts)
rm(peaks)
# load PBMC reference
reference <- LoadH5Seurat("../pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                 "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")
```

Joint UMAP visualization

``` r
# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()


# label clusters based on SCT assay 
DefaultAssay(pbmc) <- "SCT"
pbmc <- FindNeighbors(pbmc, dims=1:50)
pbmc <- FindClusters(pbmc)
pbmc@meta.data$SCT_snn_res.0.8 <- as.character(pbmc@meta.data$SCT_snn_res.0.8)
```

UMAP visualization for ATAC assay

``` r
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RUNTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff='q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction ='lsi', dims=2:10, reduction.name ='umap.atac', reduction.key='atacUMAP_')

# label clusters based on ATAC assay 
pbmc <- FindNeighbors(pbmc, reduction ='lsi', dims=2:10)
pbmc <- FindClusters(pbmc, verbose=T, resolution = 0.3)
pbmc@meta.data$SCT_snn_res.0.3 <- as.character(pbmc@meta.data$SCT_snn_res.0.3)
```

Linking peaks to genes

``` r
DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
)
```

Clean the metadata information

``` r
pbmc@meta.data$CellType <- pbmc@meta.data$predicted.id
columns_to_remove <- grep("prediction", names(pbmc@meta.data))
pbmc@meta.data <- pbmc@meta.data[, -columns_to_remove]
```

Save links, peaks and h5ad

``` r
saveRDS(Annotation(pbmc[['peaks']]), file ="annotation.rds")
saveRDS(granges(pbmc[['peaks']], file ='peaks.rds'))
saveRDS(Links(pbmc[["peaks"]], file='links.rds'))

SaveH5Seurat(pbmc, filename='pbmc10k_multiome.h5Seurat')
Convert('pbmc10k_multiome.h5Seurat', dest='h5ad')
```
## Part 2
Generate celltype.txt for bam splitting (In Python)

``` python
import scanpy as sc
import pandas as pd
import numpas as np
data = sc.read_h5ad("../pbmc10k_multiome.h5ad")
data.obs['CellType'].unique()
rename_annotation ={"CD4 Naive": "CD4.Naive", 
                    "CD4 TCM" : "CD4.TCM", 
                    "CD4 CTL" : "CD4.CTL", 
                    "CD4 TEM" : "CD4.TEM", 
                    "CD4 Proliferating" : "CD4.Proliferating",
                    "CD8 Naive" : "CD8.Naive", 
                    "CD8 TEM" : "CD8.TEM", 
                    "CD8 TCM" : "CD8.TCM", 
                    "CD8 Proliferating" : "CD8.Proliferating", 
                    "NK Proliferating" : "NK.Proliferating",
                    "B naive" : "B.naive", 
                    "B intermediate" : "B.intermediate", 
                    "B memory" : "B.memory",
                    "CD14 Mono" : "CD14.Mono", 
                    "CD16 Mono" : "CD16.Mono"
}


data.obs['CellType'] = data.obs['CellType'].map(rename_annotation).astype('category')
df = data.obs['CellType']
df = df.cat.add_categories("Other")
for i, line in enumeratate(df):
  if line is np.nan:
    df[i] = "Other"

df.to_csv('../pbmc10k_celltype.txt', header=None, sep ="\t", mode='a')
```

BIGWIG file generation

Python package requirement

1.  Sinto: <https://timoast.github.io/sinto/index.html>

2.  deepTools: <https://deeptools.readthedocs.io/en/develop/index.html>

``` bash
#! bin/bash
#$ -pe thread 16
#$ -cwd

# split bam file based on pbmc10k_celltype.txt
source ../miniconda3/etc/profile.d/conda.sh
conda active multiome

sinto filterbarcodes pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam -c pbmc10k_celltype.txt -p 16
```

``` bash
#! bin/bash
#$ -cwd

# generate BIGWIG on .bam files
source ../miniconda3/etc/profile.d/conda.sh
conda activate multiome

for val in *.bam; do 
    val2 ="${val%.bam}.bw"
    samtools index $val
    bamCoverage --bam $val -o $val2 --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --ignoreForNormalization chrX
done
```

Generate the bw.cluster file

``` bash
ls *bw > bw.cluster_c1 
ls *bw | sed -e 's/\.bam$//' > bw.cluster_c2
sed -i 's/./ /g' bw.cluster_c2
paste bw.cluster_c1 bw.cluster_c2 > bw.cluster
```
Sample of bw.cluster:
CellType  | |
------------- | -------------
ASDC.bw  | ASDC
B.intermediate.bw | B intermediate
B.memory.bw | B memory
B.naive.bw  | B naive

#### Data files and directory structures for loading into cellxgene VIP
``` bash
cellxgene launch --host 0.0.0.0 --port 8802 --disable-annotations --verbose  ./pbmc10k_multiome.h5ad
```
<img src=https://user-images.githubusercontent.com/15882624/147169763-d6999567-db31-4ea9-950e-16b81252afa5.jpg width=240>
