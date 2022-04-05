
# Methods

## Client-side	Integration	by	a	jsPanel	Window	(VIP) {-}
Following section in config.sh file.



```bash
<script src="static/jquery.min.js"></script>
<link href='static/jspanel/dist/jspanel.css' rel='stylesheet'>
<script src='static/jspanel/dist/jspanel.js'></script>
<script src='static/jspanel/dist/extensions/modal/jspanel.modal.js'></script>
<script src='static/jspanel/dist/extensions/tooltip/jspanel.tooltip.js'></script>
<script src='static/jspanel/dist/extensions/hint/jspanel.hint.js'></script>
<script src='static/jspanel/dist/extensions/layout/jspanel.layout.js'></script>
<script src='static/jspanel/dist/extensions/contextmenu/jspanel.contextmenu.js'></script>
<script src='static/jspanel/dist/extensions/dock/jspanel.dock.js'></script>
```

```js
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
      url: 'static/interface.html',
      done: function (panel) {
            setInnerHTML(panel.content, this.responseText);
      }
    },
    onwindowresize: function(event, panel) {
      var jptop = parseInt(this.currentData.top);
      var jpleft = parseInt(this.currentData.left);
      
      if (jptop<-10 || window.innerHeight-jptop<10 || window.innerWidth-jpleft<10 ||
      jpleft+parseInt(this.currentData.width)<10) {
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
```

```bash
</script>
EOF
insertL=$(sed -e 's/[&\\/]/\\&/g; s/|/\\|/g; s/$/\\/;' -e '$s/\\$//' <<<"$insertL")
sed -i "s|<div id=\"root\"></div>|$insertL\n&|" "cellxgene/client/index_template.html"
```



All functional VIP HTML and JavaScript code will be in “interface.html” that is independent of cellxgene code bases.

## Server-side	Integration {-}
Following section in config.sh file.

```bash
echo '
from server.app.biogenInterface import route  
@webbp.route("/biogen", methods=["POST"]) 
def biogen():
  return route(request.data,current_app.app_config)' >> cellxgene/server/app/app.py
.
.
.

strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
strweb="${strPath}/server/common/web/static/."
echo $strweb
cp interface.html $strweb
cp jquery.min.js $strweb
cp color_map.png $strweb

cp -R DataTables $strweb
cp -R jspanel $strweb

cp cellxgene/server/test/decode_fbs.py $strPath/server/app/.
cp VIPInterface.py $strPath/server/app/.
```



## Communication	between	VIP	and	cellxgene	web	GUI {-}
Cellxgene client utilizes React Redux that is the official React binding for Redux. It lets your React components read data from a Redux store, and dispatch actions to the store to update data.

So, this line of code is appended to the end of client/src/reducers/index.js of cellxgene source code to expose the store to the browser.



```js
window.store = store;
```


<script type="text/javascript">
window.store = store;
</script>


By doing this, Redux store holding client data and user selections are visible to VIP to access variables and dispatch actions to control cellxgene user interface. For example,

- Unselect / select a feature. GUI is refreshed automatically after dispatching.



```js
window.store.dispatch({type: "categorical metadata filter deselect", metadataField: "louvain", categoryIndex: 5})
window.store.dispatch({type: "categorical metadata filter select", metadataField: "louvain", categoryIndex: 5})
```

- Get state of just finished action and synchronize gene input and cell selections from main window to VIP if corresponding action was performed.



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

## Diffxpy	Integration {-}
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


## Create	h5ad	file	from	Seurat	object
First, export the following from Seurat object in R: **expression matrix (assume normalized), metadata and coordinates (pca, tsne, umap) as separate txt files.**

Next in Python, create an AnnData object from 10x (scanpy.read_h5ad function) as a starting point. Then replace the expression matrix, meta data and coordinates as following, a h5ad file would be generated.

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

# read in expression matrix as numpy.ndarray as following:
exp_mat = np.loadtxt(fname =”expression matrix .txt")
adata.X = exp_mat

# convert dense matrix into sparse matrix to save storage space and memory usage
adata.X = csc_matrix(adata.X)_matrix

# add short description and initial graph settings. “|” and “by” are delimiters for VIP to parse the initial settings. Please follow the same rule for your own h5ad files.
adata.obs['>Description'] = ['Human brain snRNAseq 46k cells (MS Nature 2019 Schirmer et al.); data normalized, log transformed and scaled UMI; platform - 10X v2 chemistry | embedding by umap; color by cell_type']*adata.n_obs

# Then last step to save h5ad:
adata.write_h5ad("final output.h5ad")
```




