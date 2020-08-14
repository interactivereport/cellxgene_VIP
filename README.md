# cellxgene VIP unleashes full power of interactive visualization, plotting and analysis of scRNA-seq data in the scale of millions of cells


Since the first single-cell RNA sequencing (scRNA-seq) study was debuted in 20091, over 1050 scRNA-seq studies have been published to date. At least 25 studies reported profiles in 200k or more cells2. The largest scRNA-seq study reported 2.5 million mouse cells3. It is foreseen that there will be a trend of increasing cell size of 500k or more in scRNA-seq studies. The sheer amount of data had brought a challenge in visualizing and interactively exploring such big data set for scientists, even computational savvy specialists.

   Cellxgene4 is a leading open source scRNA-Seq data visualization tool recommended in a recent evaluation5, which can handle millions of cells by leveraging modern web techniques with user-friendly interactive features. However, lack of some essential plotting and analysis functions hinders its utility and limits scientists from taking advantage of ever accumulating scRNA-seq data in the public to its full potential. 
   
   To fill the gap, we developed a plugin of cellxgene named Visualization in Plugin, in short VIP, to address the urgent needs for such essential functions for interactive visual exploration and generation of publication-ready figures. Notably, it greatly extended plotting functions to generate violin, stacked violin, stacked bar, heatmap, volcano, embedding, dot, track, density, 2D density, sankey and dual-gene plot in high-resolution by calling server-side scanpy’s6 plotting functions and general Python plotting libraries as illustrated in Figure 1 and Supplementary Tutorial.

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
