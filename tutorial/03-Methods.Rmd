# Methods

## Client-side	Integration	by	a	jsPanel	Window	(VIP) {-}
Following section in config.sh file.

NEED sh code

All functional VIP HTML and JavaScript code will be in “interface.html” that is independent of cellxgene code bases.

## Server-side	Integration {-}
Following section in config.sh file.
(NEED sh code)

## Communication	between	VIP	and	cellxgene	web	GUI {-}
Cellxgene client utilizes React Redux that is the official React binding for Redux. It lets your React components read data from a Redux store, and dispatch actions to the store to update data.

So, this line of code is appended to the end of client/src/reducers/index.js of cellxgene source code to expose the store to the browser.

(NEED sh code)

By doing this, Redux store holding client data and user selections are visible to VIP to access variables and dispatch actions to control cellxgene user interface. For example,

- Unselect / select a feature. GUI is refreshed automatically after dispatching.

(NEED sh code)

- Get state of just finished action and synchronize gene input and cell selections from main window to VIP if corresponding action was performed.

(NEED sh code)

## Diffxpy	Integration {-}
This is the sample pseudocode, please see VIPInterface.py for actual implementation.

(PY CODE)

## Create	h5ad	file	from	Seurat	object
First, export the following from Seurat object in R: expression matrix (assume normalized), metadata and coordinates (pca, tsne, umap) as separate txt files.

Next in Python, create an AnnData object from 10x (scanpy.read_h5ad function) as a starting point. Then replace the expression matrix, meta data and coordinates as following, a h5ad file would be generated.

(PY CODE)

