# Cellxgene VIP: Visualization in Plugin scales interactive plots of scRNA-Seq data to millions of cells
Development for Cellxgene VIP

Since the first single cell RNA-seq (scRNA-seq) study was published in 2009 (Tang, et al., 2009), there are close to 1000 scRNA-seq studies have been publised to date. At least 25 studies reported 200k or more cells (Svensson, et al., 2019). The largest scRNA-seq study reported 2.5 million mouse cells (Rodriques, et al., 2019). We foresee scRNA-seq size of 500k or more cells would become the new norm. With that, it is getting more and more challenge to visualize and interactively explore such big datasets for both computational and wet-lab bench scientists. 

Cellxgene (Chan Zuckerberg Initiative, [Accessed: 5 May 2020]) is the only web based open source interactive explorer which can handle millions of cells by leveraging modern web de-velopment techniques (Çakır, et al., 2020). However, after beta testing by biologists, we quickly realized that essential plotting /data downloading functions are missing from current Cellxgene release. 

We present a new plugin Cellxgene VIP to address the urgent needs for such additional functions. Cellxgene VIP also provides easy to use server end installation of Cellxgene with customized local modifications. We made Cellxgene VIP open source tool and free of charge.


# installation instruction

## 0. create anaconda python 3.7 and enable conda environment and install nodejs
``` bash
conda create -n <env name> python=3.7
```
## 1. Install cellxgene by run the config.sh in the folder
```bash
./config.sh
```
## 2. Run cellxgene by specifiying the single cell h5ad file along with the host and the port:
```bash
cellxgene launch ...h5ad --host ... -p ... --disable-annotations
```
*note: while spinning up the cellxgene from HPC, do **NOT** use qlogin. **ssh directly to the server**.*
