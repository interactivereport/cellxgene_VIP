# cellxgene VIP unleashes full power of interactive visualization, plotting and analysis of scRNA-seq data in the scale of millions of cells


Since the first single-cell RNA sequencing (scRNA-seq) study was debuted in 2009[1], over 1050 scRNA-seq studies have been published to date. At least 25 studies reported profiles in 200k or more cells[2]. The largest scRNA-seq study reported 2.5 million mouse cells[3]. It is foreseen that there will be a trend of increasing cell size of 500k or more in scRNA-seq studies. The sheer amount of data had brought a challenge in visualizing and interactively exploring such big data set for scientists, even computational savvy specialists.

Cellxgene[4] is a leading open source scRNA-Seq data visualization tool recommended in a recent evaluation[5], which scales well in millions of cells and scores high in user experience by leveraging modern web techniques with interactive features. We tested with a large community of biologists, cellxgene works the best in meeting the need of handling data themselves except lack of some essential plotting seen in scRNA-seq related publications and analysis functions that biologists are used to, hinders its utility and limits scientists from taking advantage of ever accumulating scRNA-seq data in the public to its full potential.

To fill the gap, we developed a plugin of cellxgene named Visualization in Plugin, in short VIP, to address the urgent needs of such essential functions for interactive visual exploration and generation of publication-ready figures. Notably, it enhanced plotting functions significantly to generate violin, stacked violin, stacked bar, heatmap, volcano, embedding, dot, track, density, 2D density, sankey and dual-gene plot in high-resolution by calling server-side scanpy’s[6] plotting functions and general plotting libraries as illustrated in Figure 1 and Supplementary Tutorial.

![cellxgene_VIP](https://interactivereport.github.io/cellxgene_VIP/cellxgene_VIP.png?raw=true "cellxgene_VIP")
**Figure 1 | cellxgene VIP serves as an ecosystem of plugins which provide essential functions for
publication-ready, interactive visualization, as well as CLI for analytics.** (svg files were assembled by
figureComposer[7] with zoomable version available at https://bit.ly/2QqdMg3 that is best viewed by Chrome)
**(a) Multi-tSNE/UMAP plot** visually highlights which cells expressing cell markers on selected embedding (UMAP
based on harmony batch correction in this example). **(b) Dual-gene plot** highlights cells express SYT1 and GAD1
(green SYT1 only, red GAD1 only, yellow co-expression of STY1 and GAD1), expression cutoff 2.2. **(c) Stacked
barplot** demonstrates the fraction of each major cell type across each sample (C are Control and MS are MS
patients). **(d) Trackplot** shows expression of lineage marker genes across individual cells in annotated clusters.
**(e) Violin plot** shows the AQP4 gene expression across cell types. **(f) Sankey diagram** (a.k.a. Riverplot) provides
quick and easy way to explore the inter-dependent relationship of variables in the MS snRNAseq dataset8. **(g)
Density plots** shows expression of marker genes across annotated clusters and split across cell types. **(h)
Stacked violin and Dot plot** are the key visualizations of selected cell markers across cell types. They highlight
their selective expression and validates the scRNAseq approach and visualization method. **(i) Command Line
Interface (CLI)** exposed by mini Jupyter Notebook to provide maximal flexibility of doing various analytics on the
whole or sliced single cell dataset.

1.	Tang, F. et al. mRNA-Seq whole-transcriptome analysis of a single cell. Nature Methods 6, 377-382 (2009).
2.	Svensson, V., da Veiga Beltrame, E. & Pachter, L. A curated database reveals trends in single-cell transcriptomics. bioRxiv, 742304 (2019).
3.	Rodriques, S.G. et al. Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution. Science 363, 1463 (2019).
4.	Chan Zuckerberg Initiative chanzuckerberg/cellxgene: An interactive explorer for single-cell transcriptomics data. Available at https://chanzuckerberg.github.io/cellxgene/ ([Accessed: 5 May 2020]).
5.	Cakir, B. et al. Comparison of visualization tools for single-cell RNAseq data. NAR Genomics and Bioinformatics 2 (2020).
6.	Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19, 15 (2018).
7. Li, K. et al. figureComposer: A web-based interactive multi-panel bio-infographic designing tool. bioRxiv, 2020.2003.2004.976589 (2020).

# Demo site: https://cellxgenevip-ms.bxgenomics.com/

# installation instruction

## 0. Install anaconda python 3.7 and nodejs if not available on server
    - Anaconda3-2020.02-Linux-x86_64.sh
    - node-v12.16.2-linux-x64

## 1. create and enable conda environment
``` bash
conda create -n <env name> python=3.7
conda activate <env name>
```
## 2. Install cellxgene by run the config.sh in the folder
```bash
git clone https://github.com/interactivereport/cellxgene_VIP.git
cd cellxgene_VIP
./config.sh
```
## 3. Run cellxgene by specifiying the single cell h5ad file along with the host and the port, use "ps" to find used ports to spare
```bash
ps -ef | grep cellxgene
cellxgene launch --host <xxx> --port <xxx> --disable-annotations --verbose <h5ad file>
```
## 4. From web browser (Chrome is preferred), access http(s)://host:port
*note: while spinning up the cellxgene from HPC, do **NOT** use qlogin. **ssh directly to the server**.*

# Updating
## run "update.index_template.sh" if jsPanel is modified, seldom.
## run "update.VIPInterface.sh" if interface.html or VIPInterface.py is changed, often.
