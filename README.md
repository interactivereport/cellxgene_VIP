# paraCell: A novel software tool for the interactive analysis and visualization of standard and dual host-parasite scRNA-seq data

To meet the growing demands from scientists to effectively extract deep insights from single cell RNA sequencing, spatial transcriptomics, and emerging multiome datasets, we developed cellxgene VIP (Visualization In Plugin), a frontend interactive visualization plugin of cellxgene framework, which greatly expanded capabilities of the base tool in the following aspects. First, it generates a comprehensive set of over eighteen commonly used quality control and analytical plots in high resolution with highly customizable settings in real time. Second, it provides more advanced analytical functions to gain insights on cellular compositions and deep biology, such as marker gene identification, differential gene expression analysis, and gene set enrichment analysis. Third, it empowers advanced users to perform analysis in a Jupyter Notebook like environment, dubbed Command Line Interface (CLI) by programming in Python and/or R directly without limiting themselves to functional modules available via graphical user interface (GUI). Finally, it pioneers methods to visualize multi-modal data, such as spatial transcriptomics embedding aligned with histological image on one slice or multiple slices in a grid format, and the latest 10x Genomic Multiome dataset where both DNA accessibility and gene expression in the same cells are measured, under the same framework in an integrative way to fully leverage the functionalities mentioned above. Taken together, the open-source tool makes large scale single cell data visualization and analysis more accessible to biologists in a user-friendly manner and fosters computational reproducibility by simplifying data and code reuse through the CLI.  Going forward, it has the potential to become a crowdsourcing ecosystem for the scientific community to contribute even more modules to the Swiss Army knife of single cell data exploration tools.

![paraCell](fig1.png "paraCell")

**Figure 1 | Overview of the four main options in paraCell.** **_Advanced Gene Search_** - a custom gene-set is used to create a cellxgene VIP multi-gene violin plot (_Plasmodium berghei_ dataset). **_Host-Parasite Interactions_** - selection made on the Parasite UMAP is used to update the Host UMAP (Cow-_Theileria_ dataset). **_Database Integration_** - paraCell results table is shown alongside links to relevant external database systems (_P. berghei_ dataset). **_Trajectory Inference_** - [Slingshot](https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html) trajectories drawn over [PHATE](https://phate.readthedocs.io/en/stable/) space of single-cell data, enables [tradeSeq](https://bioconductor.org/packages/release/bioc/html/tradeSeq.html) to plot expression of genes over their “pseudo time” for different conditions - example relates the expression of the gene PAD2 to the progression of both trajectories in the WT and KO condition (_Trypanosoma brucei_ dataset).

## Demo datasets:
- **_Plasmodium berghei_** - https://cellatlas-cxg.mvls.gla.ac.uk/view/Pb.Combined.h5ad/
- **_Trypanosoma brucei_** - https://cellatlas-cxg.mvls.gla.ac.uk/view/Tbrucei.h5ad/

## Online tutorial: https://interactivereport.github.io/cellxgene_VIP/tutorial/docs

## Video tutorial: https://youtu.be/s3Cdz77ioGk

## Pre-print: https://www.biorxiv.org/content/10.1101/2020.08.28.270652v2

# Installation instruction

## 1. Install miniconda if not available on server (https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
``` bash
bash ~/Downloads/Miniconda3-latest-Linux-x86_64.sh
conda install mamba -n base -c conda-forge
```

## 2. Create and enable conda environment
``` bash
git clone https://github.com/sii-cell-atlas/paraCell.git
cd paraCell

source <path to Miniconda3>/etc/profile.d/conda.sh 
conda config --set channel_priority flexible
mamba env create -n <env name, such as: paraCell> -f paraCell_conda_R.yml (local R under conda, no root privilege needed)

For Mac User, conda env create -n <env name, such as: paraCell> -f paraCell.macOS.yml

mamba env update -f r_dependencies.yml --name paraCell

conda activate <env name, such as: paraCell>
or
source activate <env name>
```
## 3. Install cellxgene by running config.sh in "paraCell" directory
```bash
./config.sh
For Mac User, ./config.macOS.sh
```
## 4. Run cellxgene by specifiying a h5ad file storing scRNA-seq data along with a host and a port, use "ps" to find used ports to spare, see https://chanzuckerberg.github.io/cellxgene/posts/launch for details.
```bash
ps -ef | grep cellxgene
Rscript -e 'reticulate::py_config()'
# Run the following command if the output of the above command doesn't point to the Python in your env.
export RETICULATE_PYTHON=`which python`
cellxgene launch --host <xxx> --port <xxx> --disable-annotations --verbose <h5ad file>
```
## 5. From web browser (Chrome is preferred, Version 87.0.4280.88 or 87.0.4280.141 is used), access http(s)://host:port

You should be able to see this in Console of Chrome Developer Tools if everything is right.
![VIP_ready](https://user-images.githubusercontent.com/29576524/92059839-46482d00-ed60-11ea-8890-8e1b513a1656.png)

*note: while spinning up the cellxgene from HPC, do **NOT** use qlogin. **ssh directly to the server**.*

# Updating
```bash
./update.VIPInterface.sh all # if "interface.html" or "VIPInterface.py" is modified or new code needs to go to right location, often.

./update.index_template.sh # if jsPanel is modified, very rare.
```
