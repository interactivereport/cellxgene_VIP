# Methods

## Client-side Integration by a jsPanel Window (VIP)

Related lines in config.sh file.




```bash
sed -i "s|<div id=\"root\"></div>|$(sed -e 's/[&\\/]/\\&/g; s/|/\\|/g; s/$/\\/;' -e '$s/\\$//' index_template.insert)\n&|" "cellxgene/client/index_template.html"
```

The content of the index_template.insert file.


```js
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jquery.min.js"></script>
<script src="https://d3js.org/d3.v4.min.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/stackedbar/d3.v3.min.js"></script>
<link href="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/jspanel.css" rel="stylesheet">
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/jspanel.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/modal/jspanel.modal.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/tooltip/jspanel.tooltip.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/hint/jspanel.hint.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/layout/jspanel.layout.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/contextmenu/jspanel.contextmenu.js"></script>
<script src="https://interactivereport.github.io/cellxgene_VIP/static/jspanel/dist/extensions/dock/jspanel.dock.js"></script>
<script>
// execute JavaScript code in panel content
var setInnerHTML = function(elm, html) {
    elm.innerHTML = html;
    Array.from(elm.querySelectorAll('script')).forEach( oldScript => {
        const newScript = document.createElement('script');
        Array.from(oldScript.attributes)
        .forEach( attr => newScript.setAttribute(attr.name, attr.value) );
        newScript.appendChild(document.createTextNode(oldScript.innerHTML));
        oldScript.parentNode.replaceChild(newScript, oldScript);
    });
}
var plotPanel = jsPanel.create({
    panelSize: '190 0',
    position: 'left-top 160 6',
    dragit: { containment: [-10, -2000, -4000, -2000] }, // set dragging range of VIP window
    boxShadow: 1,
    border: "solid #D4DBDE thin",
    contentOverflow: 'scroll scroll', // adding scrolling bars
    headerControls:{
      close: 'remove',
      minimize: 'remove',
      maximize: 'remove'
    },
    headerTitle: function () {return '<strong>Visualization in Plugin</strong>'},
    contentAjax: {
        url: window.location.href.replace(/\\\/+$/,'')+'/static/interface.html',
        done: function (panel) {
               setInnerHTML(panel.content, this.responseText);
        }
    },
    onwindowresize: function(event, panel) {
        var jptop = parseInt(this.currentData.top);
        var jpleft = parseInt(this.currentData.left);
        if (jptop<-10 || window.innerHeight-jptop<10 || window.innerWidth-jpleft<10 || jpleft+parseInt(this.currentData.width)<10) {
            this.reposition("left-top 160 6");
        }
    },
    onunsmallified: function (panel, status) {
        this.reposition('center-top -370 180');
        this.resize({ width: 740, height: function() { return Math.min(480, window.innerHeight*0.6);} });
    },
    onsmallified: function (panel, status) {
        this.reposition('left-top 160 6');
        this.style.width = '190px';
    }
}).smallify();
plotPanel.headerbar.style.background = "#D4DBDE";
</script>
```

All functional VIP HTML and JavaScript code will be in a new file called “interface.html”, which is out of cellxgene code base.

## Server-side Integration

Related in config.sh file.

```bash
echo '
from server.app.VIPInterface import route
@webbp.route("/VIP", methods=["POST"])
def VIP():
    return route(request.data,current_app.app_config)' >> cellxgene/server/app/app.py

cd cellxgene
make pydist
make install-dist
cd ..

## finished setting up ------
./update.VIPInterface.sh all
```

The content of update.VIPInterface.sh file.

```bash
#!/usr/bin/env bash
if [ -n "$1" ]; then
echo "usually update once"
fi

## finished setting up ------
strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb
cp vip.env $strPath/server/app/. 2>/dev/null | true

cp fgsea.R $strPath/server/app/.
mkdir -p $strPath/server/app/gsea
cp gsea/*gmt $strPath/server/app/gsea

if [ -n "$1" ]; then
  cp Density2D.R $strPath/server/app/.
  cp bubbleMap.R $strPath/server/app/.
  cp violin.R $strPath/server/app/.
  cp volcano.R $strPath/server/app/.
  cp browserPlot.R $strPath/server/app/.
  if [ "$(uname -s)" = "Darwin" ]; then
    sed -i .bak "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
    sed -i .bak "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
  else
    sed -i "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
    sed -i "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
  fi

  find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;
fi

echo -e "\nls -l $strweb\n"
ls -l $strweb
```



## Communication between VIP and cellxgene web GUI

Cellxgene client utilizes React Redux that is the official React binding for Redux. It lets your React components read data from a Redux store, and dispatch actions to the store to update data.

So, this following code in config.sh appends "window.store = store;" to the end of client/src/reducers/index.js of cellxgene source code to expose the store to the browser.


```bash
echo -e "\nwindow.store = store;" >> cellxgene/client/src/reducers/index.js
```


By doing this, Redux store holding client data and user selections are visible to VIP to access variables and dispatch actions to control cellxgene user interface. For example,

- Unselect / select a feature. GUI is refreshed automatically after dispatching.



```js
window.store.dispatch({type: "categorical metadata filter deselect", metadataField: "louvain", categoryIndex: 5})
window.store.dispatch({type: "categorical metadata filter select", metadataField: "louvain", categoryIndex: 5})
```

- Get the state of a just finished action and synchronize gene input and cell selections from main window to VIP if corresponding action was performed.



```js
const unsubscribe = window.store.subscribe(() => {
  if (window.store.getState()["@@undoable/filterState"].prevAction) {
    actionType = window.store.getState()["@@undoable/filterState"].prevAction.type;
    if (actionType.includes("user defined gene success") ||
    actionType.includes("store current cell selection as differential set")) {
      sync();
      }
  }
});
```

## Diffxpy Integration

This is the sample pseudocode, please see VIPInterface.py for actual implementation.


```python
import scanpy as sc
import pandas as pd
import diffxpy.api as app
# set 1 of cells as cell1; set 2 of cells as cell2


with app.get_data_adaptor() as data_adaptor:
  X1 = data_adaptor.data.X[cell1]
  X2 = data_adaptor.data.X[cell2]


adata = sc.AnnData(pd.concat([X1,X2]),pd.DataFrame(['grp1' for i in range(X1.shape[0])]+['grp2' for i in range(X2.shape[0])],columns=['comGrp']))
deg = de.test.two_sample(adata,'comGrp').summary()
#deg is a dataframe contains the folloing columns ['gene','log2fc','pval','qval']

```


## Create a h5ad file from Seurat object

First, export the following from Seurat object in R: **expression matrix (assume normalized), metadata and coordinates (pca, tsne, umap) as separate txt files.**

Next in Python, create an AnnData object from 10x (scanpy.read_h5ad function) as a starting point. Then replace the expression matrix, meta data and coordinates as shown in the following Python code block to generate a h5ad file.


```python
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
from numpy import ndarray, unique
from scipy.sparse.csc import csc_matrix

adata= sc.read_h5ad("previous generated .h5ad")

# read clustering res
xpca = pd.read_csv(“./data/harmony_clustered.h5ad.pca_coordinates.txt", sep='\t', encoding='utf-8')
xtsne = pd.read_csv(“./data/harmony_clustered.h5ad.tsne_coordinates.txt", sep='\t', encoding='utf-8')
xumap = pd.read_csv(“./data/harmony_clustered.h5ad.umap_coordinates.txt", sep='\t', encoding='utf-8')
xobs = pd.read_csv(“./data/harmony_clustered.h5ad.meta_data.txt", sep='\t', encoding='utf-8')

xpca.set_index('index', inplace=True)
xtsne.set_index('index', inplace=True)
xumap.set_index('index', inplace=True)
xobs.set_index('index', inplace=True)

adata.obsm['X_pca'] = np.array(xpca.loc[adataRaw.obs.index])

adata.obsm['X_tsne'] = np.array(xtsne.loc[adataRaw.obs.index])
adata.obsm['X_umap'] = np.array(xumap.loc[adataRaw.obs.index])
adata.obs = xobs.loc[adataRaw.obs.index] # this is a pandas dataframe

# read in expression matrix as numpy.ndarray
exp_mat = np.loadtxt(fname =”expression matrix .txt")
adata.X = exp_mat

# convert dense matrix into sparse matrix to save storage space and memory usage
adata.X = csc_matrix(adata.X)_matrix

# add short description and initial graph settings. “|” and “by” are delimiters for VIP to parse the initial settings. Please follow the same rule for your own h5ad files.
adata.obs['>Description'] = ['Human brain snRNAseq 46k cells (MS Nature 2019 Schirmer et al.); data normalized, log transformed and scaled UMI; platform - 10X v2 chemistry | embedding by umap; color by cell_type']*adata.n_obs

# Then last step to save h5ad:
adata.write_h5ad("final output.h5ad")
```

When the h5ad file is uploaded to cellxgeneVIP, AnnData.X matrix is to be used for visualization and DEG analysis. By default, the data (e.g, raw count matrix) is assumed to be unscaled , however, if the data have been scaled or normalized, the user needs to turn off the option ‘Scale data to unit variance for plotting:’ in ‘Global Setting’.


## Prepare Spatial Data for Visualization
### 10X visium data

This is the sample script for the spatial demo data set adapted from [cellxgene_VIP spatial transcriptomics notebook](https://github.com/interactivereport/cellxgene_VIP/blob/master/notebook/spatial_transcriptomics.ipynb). Please go to [cellxgene_VIP spatial transcriptomics notebook](https://github.com/interactivereport/cellxgene_VIP/blob/master/notebook/spatial_transcriptomics.ipynb) to check the intermediate outputs.

#### Download input demo data


```bash
wget -O 10X_ST_demo.tar.gz https://zenodo.org/record/5765589/files/10X_ST_demo.tar.gz?download=1

tar -zxvf 10X_ST_demo.tar.gz
```

#### Using python to prepare 10X visium spatial data


```python
import pandas as pd
import numpy as np
import sys, getopt
import matplotlib.pyplot as plt
import seaborn as sns
import copy
import numpy as np
from PIL import Image, ImageOps
import scanpy as 
import anndata

# change directory to 10X_ST_demo which is just downloaded, if in current directory
os.chdir('./10X_ST_demo')
inputfile='inputfiles.txt'
samples = pd.read_csv(inputfile, header = None)[0]

def read_each(i):
    adata = sc.read_visium(i)
    adata.var_names_make_unique()
    # flip Y axis to show correctly in cellxgene VIP
    adata.obsm['spatial'][:,1] = -adata.obsm['spatial'][:,1]
    return(adata)

adatals = [read_each(i) for i in samples]

sampleIDs = samples.str.extract(r'10X_demo_data_(.*)')
sampleIDs = "V1_"+ sampleIDs

adata_merge = sc.AnnData.concatenate(*adatals, batch_key='sample', join='outer', batch_categories= sampleIDs[0].astype("category"))

for i in list(range(len(adatals))):
    print(i)
    # add back the spatial coordinates as separate embeddings
    adata_merge.obsm['X_spatial_'+list(adatals[i].uns["spatial"])[0]] = np.zeros(adata_merge.obsm['spatial'].shape)
    adata_merge.obsm['X_spatial_'+list(adatals[i].uns["spatial"])[0]][np.where(adata_merge.obs['sample']==list(adatals[i].uns["spatial"])[0])] = adatals[i].obsm['spatial']

adata_merge.uns['spatial'] = dict()
for i in list(range(len(adatals))):
    adata_merge.uns['spatial']["spatial_"+list(adatals[i].uns["spatial"])[0]] = adatals[i].uns['spatial'][list(adatals[i].uns["spatial"])[0]]
    
# QC metric
adata_merge.var["mt"] = adata_merge.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(adata_merge, qc_vars=["mt"], inplace=True)

# QC plots
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata_merge.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata_merge.obs["total_counts"][adata_merge.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata_merge.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata_merge.obs["n_genes_by_counts"][adata_merge.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])

# filtering, turn this off to keep all spots for visualization also cutoffs are case-by-case based on the QC plots
#sc.pp.filter_cells(adata, min_counts=5000)
#sc.pp.filter_cells(adata, max_counts=35000)
#adata = adata[adata.obs["pct_counts_mt"] < 20]
#print(f"#cells after MT filter: {adata.n_obs}")
#sc.pp.filter_genes(adata, min_cells=10)

# normalization, log1p transformation and select HVGs
sc.pp.normalize_total(adata_merge, inplace=True)
sc.pp.log1p(adata_merge)
sc.pp.highly_variable_genes(adata_merge, flavor="seurat", n_top_genes=2000)

# PCA, UMAP and clustering by leiden
sc.pp.pca(adata_merge)
sc.pp.neighbors(adata_merge)
sc.tl.umap(adata_merge)
sc.tl.leiden(adata_merge, key_added="clusters")

# collect sample names
sampleNames = list()
for f in list(adata_merge.obsm):
    if "spatial_" in f: # search for the pattern
        library_id=f.replace("X_spatial_","") # parse the string and get the sample id
        #library_id=library_id.replace("V1_","")
        sampleNames.append(library_id)
        
from PIL import Image
spatial=adata_merge.uns["spatial"]
dim=''
import math
if dim=='':
    height = math.ceil(math.sqrt(len(samples)))
    width = math.ceil(len(samples)/height)
else:
    width,height = dim.split('x')

idx = 0
size=700
#creates a new empty image, RGB mode, and size 1400 by 1400.
new_im = Image.new('RGB', (size*width,size*height))
for i in range(0,size*width,size):
    for j in range(0,size*height,size):
        # load the image from the object
        #im = Image.fromarray((spatial["spatial_V1_"+samples[idx]]["images"]["lowres"]* 255).round().astype(np.uint8)) # found a solution to covert float32 to unit8
        im = Image.fromarray((spatial["spatial_"+sampleNames[idx]]["images"]["lowres"]* 255).round().astype(np.uint8)) # found a solution to covert float32 to unit8
        # paste images together
        new_im.paste(im, (j,i))
        print(idx)
        idx = idx+1
        if idx>=len(sampleNames):
            break

# fake a adata.uns by providing merged lowres image and scale factors 1
adata_merge.uns['spatial']['spatial_Merged'] = copy.deepcopy(adata_merge.uns['spatial'][list(adata_merge.uns['spatial'])[0]])
adata_merge.uns['spatial']['spatial_Merged']['images']["hires"] = np.asarray(new_im)
adata_merge.uns['spatial']['spatial_Merged']['images']["lowres"] = np.asarray(new_im)
adata_merge.uns['spatial']['spatial_Merged']['scalefactors']['tissue_lowres_scalef'] = 1
adata_merge.uns['spatial']['spatial_Merged']['scalefactors']['tissue_hires_scalef'] = 1

# add back the spatial coordinates as separate embeddings
idx = 0
adata_merge.obsm['X_spatial_Merged'] = adata_merge.obsm['spatial']
for i in range(0,size*width,size):
    for j in range(0,size*height,size):
        #library_id='spatial_V1_'+samples[idx] # parse the string and get the sample id
        library_id='spatial_'+sampleNames[idx] # parse the string and get the sample id
        print(library_id)
        tissue_lowres_scalef = spatial[library_id]['scalefactors']['tissue_lowres_scalef']
        adata_merge.obsm['X_spatial_Merged'][np.where(adata_merge.obs['sample']==sampleNames[idx])] = copy.deepcopy(adatals[idx].obsm['spatial'])
        adata_merge.obsm['X_spatial_Merged'][np.where(adata_merge.obs['sample']==sampleNames[idx]),1] = adatals[idx].obsm['spatial'][:,1]*tissue_lowres_scalef - i
        adata_merge.obsm['X_spatial_Merged'][np.where(adata_merge.obs['sample']==sampleNames[idx]),0] = adatals[idx].obsm['spatial'][:,0]*tissue_lowres_scalef + j
        idx = idx+1
        if idx>=len(sampleNames):
            break
outputfile = '10X_data.h5ad'
adata_merge.write_h5ad(outputfile)
```


### NanoString CosMx data

This is the sample script for the spatial demo data set adapted from [cellxgene_VIP public cosMx liver notebook](https://github.com/interactivereport/cellxgene_VIP/blob/master/notebook/public_cosMx_liver.ipynb). Please go to [cellxgene_VIP public cosMx liver notebook](https://github.com/interactivereport/cellxgene_VIP/blob/master/notebook/public_cosMx_liver.ipynb) to check the intermediate outputs.


#### Download input demo data
Download the public cosMx data from [NanoString website](https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/). The seurat object contains the expression, cell meta and transcript coordinates(some may contains cell boundaries). If the cell boundaries are not included in the seurat, they can be estimated from the transcript coordinates by cellPoly. The following will use [Human Liver (RNA)](https://nanostring.com/resources/seurat-object-cosmx-smi-human-liver-ffpe-dataset/) as an example.

#### Export Expression, cell meta, transcript coordinate and Reduction UMAP from seurat object using R

```r
require(Seurat)
require(CellPoly)
require(dplyr)
D <- readRDS('LiverDataReleaseSeurat_newUMAP.RDS')
saveSlideFOV <- function(slide,fov){
    selC <- rownames(D@meta.data)[D@meta.data$fov==fov&D@meta.data$slide_ID_numeric==slide]
    prefix <- paste0("Slide",slide,"_FOV",gsub(" ","0",format(fov,width=3)))
    data.table::fwrite(data.frame(gene=rownames(D@assays$RNA),D@assays$RNA[,selC]),paste0(prefix,'_exp.csv'))
    data.table::fwrite(data.frame(cID=selC,D@meta.data[selC,]),paste0(prefix,'_meta.csv'))
    data.table::fwrite(data.frame(cID=selC,D@reductions$approximateumap@cell.embeddings[selC,]),paste0(prefix,'_UMAP.csv'))
    data.table::fwrite(data.frame(cID=selC,D@reductions$approximateUMAP_bySlide@cell.embeddings[selC,]),paste0(prefix,'_slideUMAP.csv'))
    data.table::fwrite(D@misc$transcriptCoords[D@misc$transcriptCoords$fov==fov&D@misc$transcriptCoords$slideID==slide,] %>%
        dplyr::rename(cell_ID=cell_id),paste0(prefix,'_tx.csv'))
}

saveSlideFOV(1,2)
saveSlideFOV(1,3)
```

#### Extract cell boundaries from cell ID labeled tif file using python

```python
import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
def extractCellBoundaryPolygon(fov):
    print(fov)
    img = plt.imread(('NormalLiverFiles/CellStatsDir/FOV%*d/CellLabels_F%*d.tif'%(3,fov,3,fov)).replace(" ","0"))
    contours = []
    nStep = int(np.max(img)/10)
    for i in range(1,np.max(img)+1):
        if i%nStep==0:
            print("%d0%%"%(i/nStep),end=" ")
        A = img.copy()
        A[A!=i]=0
        _, B = cv2.threshold(A, i-0.5, 255, cv2.THRESH_BINARY)
        A,_ = cv2.findContours(B.astype(np.uint8),cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        A = pd.DataFrame(A[0].reshape(len(A[0]),2),columns=['x_FOV_px','y_FOV_px'])
        A['cell_ID'] = 'c_1_%d_%d'%(fov,i)
        contours.append(A.copy())
    print()
    return pd.concat(contours)

extractCellBoundaryPolygon(2).to_csv('Slide1_FOV002_cellPolygons.csv',index=False)
extractCellBoundaryPolygon(3).to_csv('Slide1_FOV003_cellPolygons.csv',index=False)
```

##### The function to read one FOV data

```python
import os, re
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
cellID_col = 'cell_ID'
suf_px={'local':'_local_px','global':'_global_px'}
sample_col = 'library_id'
cosMx_uns = 'cosMx'
cell_uns = 'cell_polygons'
transc_uns = 'tx_loc'
img_uns = 'images'
def readOneCapture(strExp, # the path to the expression matrix of a FOV
                   strMeta, # the path to the cell meta table matrix of a FOV
                   strReduc, # a list of the paths to the expression matrix of a FOV
                   strTargetCoord, # the path to the target coordinate of a FOV
                   strCellBoundaries, # the path to the cell boundary polygons of a FOV
                   strImg, # the image associated with the FOV, like histology please make sure the scale/rotation matching the coordinates
                   sID, #the unique name of the FOV
                   cell_column='cell_id',
                   lpx='_FOV_px', # the column suffix which indicate the location on each FOV image in cell meta table
                   gpx=4256 # the column suffix which indicate the location on the global of all FOV images in cell meta table. If a numeric, a fake global will be created vertically stack all FOVs
                  ):
    # exp
    print(" exp")
    gExp = pd.read_csv(strExp,index_col=0).T
    gExp.columns.name=None
    # meta
    print(" meta")
    cInfo = pd.read_csv(strMeta,index_col=0)
    cID = cInfo.index.intersection(gExp.index)
    # create anndata
    print(" create AnnData")
    D = ad.AnnData(X=gExp.loc[cID,:],obs=cInfo.loc[cID,:])
    del gExp
    # add histology image
    D.uns[cosMx_uns]={sID:{}}
    D.uns[cosMx_uns][sID][img_uns]=plt.imread(strImg)
    # add px coordinates
    print(" add reduction")
    D.obsm['X_FOV_local'] = cInfo.loc[cID,[_ for _ in cInfo.columns if _.endswith(lpx)]].to_numpy('float32')
    if isinstance(gpx, (int, float)):
        FOV_global = cInfo.loc[cID,[_ for _ in cInfo.columns if _.endswith(lpx)]]#
        FOV_global.iloc[:,0] += gpx
        FOV_global.iloc[:,1] = D.uns[cosMx_uns][sID][img_uns].shape[1]-FOV_global.iloc[:,1]-1 # seems like the y(height) is flipped from the image and coordinates of transcript/cell boundaries
        D.obsm['X_FOV_global'] = FOV_global.to_numpy('float32')   
    else:
        D.obsm['X_FOV_global'] = cInfo.loc[cID,[_ for _ in cInfo.columns if _.endswith(lpx)]].to_numpy('float32')    
    for oneF in strReduc:
        reducName = re.sub(".csv","",os.path.basename(oneF).split("_")[-1])
        D.obsm['X_'+reducName] = pd.read_csv(oneF,index_col=0).loc[cID,].to_numpy('float32')
    # add cell boundaries
    print(" add cell boundaries")
    X = pd.read_csv(strCellBoundaries)
    X.columns = [re.sub(cell_column,cellID_col,re.sub(lpx,suf_px['local'],_)) for _ in X.columns]
    if isinstance(gpx, (int, float)):
        X['x'+suf_px['global']] = X['x'+suf_px['local']] + gpx
        X['y'+suf_px['global']] = X['y'+suf_px['local']]
    else:
        X.columns = [re.sub(gpx,suf_px['global'],_) for _ in X.columns]
    D.uns[cosMx_uns][sID][cell_uns] = X    
    # add transcript coordinates
    print(" add transcript coordinates")
    X = pd.read_csv(strTargetCoord)
    X.columns = [re.sub(cell_column,cellID_col,re.sub(lpx,suf_px['local'],_)) for _ in X.columns]
    if isinstance(gpx, (int, float)):
        X['x'+suf_px['global']] = X['x'+suf_px['local']] + gpx
        X['y'+suf_px['global']] = X['y'+suf_px['local']]
    else:
        X.columns = [re.sub(gpx,suf_px['global'],_) for _ in X.columns]
    D.uns[cosMx_uns][sID][transc_uns] = X
    return D
```

##### Use two FOVs from human liver data as an example

```python
samples = {}
for i in [2,3]:
    samples['Normal%d'%i] = {}
    for one in ['exp','meta','cellPolygons','tx','reduct','img']:
        if one == 'reduct':
            samples['Normal%d'%i][one] = ['Slide1_FOV00%d_UMAP.csv'%i,'Slide1_FOV00%d_slideUMAP.csv'%i]
        elif one == 'img':
            samples['Normal%d'%i][one] = 'NormalLiverFiles/CellStatsDir/CellComposite/CellComposite_F00%d.jpg'%(i) # this can be replaced with the right histology image
        else:#
            samples['Normal%d'%i][one] = 'Slide1_FOV00%d_%s.csv'%(i,one)

adatas = []
w = 0
for sID in samples.keys():
    D = readOneCapture(samples[sID]['exp'],
                       samples[sID]['meta'],
                       samples[sID]['reduct'],
                       samples[sID]['tx'],
                       samples[sID]['cellPolygons'],
                       samples[sID]['img'],
                       sID,
                       gpx=w
                      )
    # If no global coordinates, a fake one align all of the input vertically
    w += D.uns[cosMx_uns][sID][img_uns].shape[1]
    adatas.append(D)
adata = ad.AnnData.concatenate(*adatas,
            join="outer",
            batch_categories=list(samples.keys()),
            batch_key=sample_col,
            index_unique=None,
            uns_merge='unique')
adata.uns[cosMx_uns]['keys'] = {'%s_px'%i:{j:'%s%s'%(j,suf_px[i]) for j in ['x','y']} for i in suf_px}
adata.uns[cosMx_uns]['keys']['cell'] = cell_uns
adata.uns[cosMx_uns]['keys']['tx_loc'] = transc_uns
adata.uns[cosMx_uns]['keys']['img'] = img_uns
adata.uns[cosMx_uns]['keys']['tx_col'] = 'target'
adata.write("human_liver.h5ad")
```

#### Check the coordinates among cell boundaries, transcript locations and the image are aligned

```python
import random
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
#adata = ad.read_h5ad('human_liver.h5ad')
sID= 'Normal3'
def showImg(img):
    plt.figure(figsize=(12, 12))
    plt.imshow(img,cmap='gray')
    plt.axis('off')
    plt.show()
    plt.close()
def bresenham_line(pt1,pt2):
    x1,y1 = pt1
    x2,y2 = pt2
    points = []
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    sx = 1 if x1 < x2 else -1
    sy = 1 if y1 < y2 else -1
    error = dx - dy
    while x1 != x2 or y1 != y2:
        points.append((x1, y1))
        e2 = 2 * error
        if e2 > -dy:
            error -= dy
            x1 += sx
        if e2 < dx:
            error += dx
            y1 += sy
    points.append((x2, y2))
    return points
def loc_img(pt,imgH):
    x,y = pt
    return((round(y),round(x)))#imgH-y-1
def square_integer(ct,halfS):
    x_center, y_center = ct
    coordinates_in_square = []
    for x in range(x_center - halfS, x_center + halfS + 1):
        for y in range(y_center - halfS, y_center + halfS + 1):
            coordinates_in_square.append((x, y))
    return coordinates_in_square
cosMxKey='cosMx'
keys = adata.uns[cosMxKey]['keys']
cosMxArray={}
cood_pair={}
imgC = adata.uns[cosMxKey][sID][keys['img']].copy()
cell_color= 'ffffff' #'000000'
rgb = np.array(tuple(int(cell_color[i:i+2], 16) for i in (0, 2, 4)),dtype='uint8')
cellB = adata.uns[cosMxKey][sID][keys['cell']]
for i in cellB.cell_ID.unique():
    oneC = cellB[cellB.cell_ID==i].reset_index(drop=True)
    N = oneC.shape[0]
    for j in range(N):
        p = bresenham_line(loc_img((oneC[keys['local_px']['x']][j%N],oneC[keys['local_px']['y']][j%N]),imgC.shape[1]),
                           loc_img((oneC[keys['local_px']['x']][(j+1)%N],oneC[keys['local_px']['y']][(j+1)%N]),imgC.shape[1]))
        for a in p:
            imgC[a]=rgb
showImg(imgC)
# cell boundries with transcript within the same randomly selected 10 cell
imgC = np.ones((adata.uns[cosMxKey][sID][keys['img']].shape[0],adata.uns[cosMxKey][sID][keys['img']].shape[1],3),dtype='uint8')*255
cell_color="000000"
rgb = np.array(tuple(int(cell_color[i:i+2], 16) for i in (0, 2, 4)),dtype='uint8')
cID = random.sample(list(cellB.cell_ID.unique()),10)
txLoc = adata.uns['cosMx'][sID][keys['tx_loc']]
for i in cID:
    oneC = cellB[cellB.cell_ID==i].reset_index(drop=True)
    N = oneC.shape[0]
    for j in range(N):
        p = bresenham_line(loc_img((oneC[keys['local_px']['x']][j%N],oneC[keys['local_px']['y']][j%N]),imgC.shape[1]),
                           loc_img((oneC[keys['local_px']['x']][(j+1)%N],oneC[keys['local_px']['y']][(j+1)%N]),imgC.shape[1]))
        for a in p:
            imgC[a]=rgb
col = 'aa0000'
rgb = np.array(tuple(int(col[i:i+2], 16) for i in (0, 2, 4)),dtype='uint8')
for i in txLoc.index[txLoc[cellID_col].isin(cID)].tolist():
        for a in square_integer(loc_img((txLoc[keys['local_px']['x']][i],txLoc[keys['local_px']['y']][i]),imgC.shape[1]),4):
            imgC[a] = rgb
showImg(imgC) 
```
