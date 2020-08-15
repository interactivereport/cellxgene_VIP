# cellxgene VIP unleashes full power of interactive visualization, plotting and analysis of scRNA-seq data in the scale of millions of cells


Since the first single-cell RNA sequencing (scRNA-seq) study was debuted in 2009[1], over 1050 scRNA-seq studies have been published to date. At least 25 studies reported profiles in 200k or more cells[2]. The largest scRNA-seq study reported 2.5 million mouse cells[3]. It is foreseen that there will be a trend of increasing cell size of 500k or more in scRNA-seq studies. The sheer amount of data had brought a challenge in visualizing and interactively exploring such big data set for scientists, even computational savvy specialists.

   Cellxgene[4] is a leading open source scRNA-Seq data visualization tool recommended in a recent evaluation[5], which can handle millions of cells by leveraging modern web techniques with user-friendly interactive features. However, lack of some essential plotting and analysis functions hinders its utility and limits scientists from taking advantage of ever accumulating scRNA-seq data in the public to its full potential. 
   
   To fill the gap, we developed a plugin of cellxgene named Visualization in Plugin, in short VIP, to address the urgent needs for such essential functions for interactive visual exploration and generation of publication-ready figures. Notably, it greatly extended plotting functions to generate violin, stacked violin, stacked bar, heatmap, volcano, embedding, dot, track, density, 2D density, sankey and dual-gene plot in high-resolution by calling server-side scanpy’s[6] plotting functions and general Python plotting libraries as illustrated in Figure 1 and Supplementary Tutorial.

1. Tang, F. et al. mRNA-Seq whole-transcriptome analysis of a single cell. Nature Methods 6, 377-382 (2009).
2.	Svensson, V., da Veiga Beltrame, E. & Pachter, L. A curated database reveals trends in single-cell transcriptomics. bioRxiv, 742304 (2019).
3.	Rodriques, S.G. et al. Slide-seq: A scalable technology for measuring genome-wide expression at high spatial resolution. Science 363, 1463 (2019).
4.	Chan Zuckerberg Initiative chanzuckerberg/cellxgene: An interactive explorer for single-cell transcriptomics data. Available at https://chanzuckerberg.github.io/cellxgene/ ([Accessed: 5 May 2020]).
5.	Cakir, B. et al. Comparison of visualization tools for single-cell RNAseq data. NAR Genomics and Bioinformatics 2 (2020).
6.	Wolf, F.A., Angerer, P. & Theis, F.J. SCANPY: large-scale single-cell gene expression data analysis. Genome Biology 19, 15 (2018).


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
## 3. Run cellxgene by specifiying the single cell h5ad file along with the host and the port, use ps to find used ports
```bash
ps -ef | grep cellxgene
cellxgene launch --host <xxx> --port <xxx> --disable-annotations --verbose <h5ad file>
```
*note: while spinning up the cellxgene from HPC, do **NOT** use qlogin. **ssh directly to the server**.*

# Updating
## update.index_template.sh if jsPanel is modified, seldom.
## update.VIPInterface.sh if interface.html or VIPInterface.py is changed, often.
