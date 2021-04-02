import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import copy
import sys, getopt
import numpy as np
from PIL import Image, ImageOps

def parseArgs(argv):
    inputfolder = ''
    outputfile = ''
    path = ''
    try:
        opts, args = getopt.getopt(argv,"hi:p:o:",["input=","path=","output="])
    except getopt.GetoptError:
        print('merge_multiple_ST_image_and_spatial_coordinates.py -i <inputfolder pattern> -p <path> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('merge_multiple_ST_image_and_spatial_coordinates.py -i <inputfolder pattern> -p <path> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfolder = arg
        elif opt in ("-p", "--path"):
            path = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
    print('Input pattern is "',inputfolder,'"')
    print('path file is "',path,'"')
    print('Output file is "',outputfile,'"')
    return(inputfolder, path, outputfile)

(inputfolder, path, outputfile) = parseArgs(sys.argv[1:])

#sc.logging.print_versions()

import glob
samples = glob.glob(path+"/"+inputfolder)
print(samples)

def read_each(i):
    adata = sc.read_visium(i)
    adata.var_names_make_unique()
    # flip Y axis to show correctly in cellxgene VIP
    adata.obsm['spatial'][:,1] = -adata.obsm['spatial'][:,1]
    return(adata)

adatals = [read_each(i) for i in samples]

import anndata
adata_merge = sc.AnnData.concatenate(*adatals, batch_key='sample', join='outer')

for i in list(range(len(adatals))):
    print(i)
    # add back the spatial coordinates as separate embeddings
    adata_merge.obsm['X_spatial_'+list(adatals[i].uns["spatial"])[0]] = np.zeros(adata_merge.obsm['spatial'].shape)
    adata_merge.obsm['X_spatial_'+list(adatals[i].uns["spatial"])[0]][np.where(adata_merge.obs['sample']==str(i))] = adatals[i].obsm['spatial']

adata_merge.uns['spatial'] = dict()
for i in list(range(len(adatals))):
    adata_merge.uns['spatial']["spatial_"+list(adatals[i].uns["spatial"])[0]] = adatals[i].uns['spatial'][list(adatals[i].uns["spatial"])[0]]

# plot the images and coordinates
#spatial=adata_merge.uns["spatial"]
#figC = 1
#for f in list(adata_merge.obsm):
#    if "spatial_" in f: # search for the pattern
#        print(f)
#        library_id=f.replace("X_","") # parse the string and get the sample id
#        print(library_id)
#        tissue_hires_scalef = spatial[library_id]['scalefactors']['tissue_hires_scalef']
#        matplotlib.pyplot.figure()
#        matplotlib.pyplot.imshow(spatial[library_id]["images"]["hires"])
#        spatial1=adata_merge.obsm['X_'+library_id]
#        matplotlib.pyplot.scatter(spatial1[:,0]*tissue_hires_scalef,spatial1[:,1]*tissue_hires_scalef)
#        figC=figC+1


# QC metric
adata_merge.var["mt"] = adata_merge.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata_merge, qc_vars=["mt"], inplace=True)

# QC plots
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(adata_merge.obs["total_counts"], kde=False, ax=axs[0])
sns.distplot(adata_merge.obs["total_counts"][adata_merge.obs["total_counts"] < 10000], kde=False, bins=40, ax=axs[1])
sns.distplot(adata_merge.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.distplot(adata_merge.obs["n_genes_by_counts"][adata_merge.obs["n_genes_by_counts"] < 4000], kde=False, bins=60, ax=axs[3])

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
samples = list()
for f in list(adata_merge.obsm):
    if "spatial_" in f: # search for the pattern
        library_id=f.replace("X_spatial_","") # parse the string and get the sample id
        library_id=library_id.replace("V1_","")
        samples.append(library_id)

from PIL import Image
spatial=adata_merge.uns["spatial"]
idx = 0
size=700
#creates a new empty image, RGB mode, and size 1400 by 1400.
new_im = Image.new('RGB', (size*2,size*2))
for i in range(0,size*2,size):
    for j in range(0,size*2,size):
        # load the image from the object
        im = Image.fromarray((spatial["spatial_V1_"+samples[idx]]["images"]["lowres"]* 255).round().astype(np.uint8)) # found a solution to covert float32 to unit8
        # paste images together
        new_im.paste(im, (j,i))
        idx = idx+1
#new_im


# fake a adata.uns by providing merged lowres image and scale factors 1
adata_merge.uns['spatial']['spatial_Merged'] = copy.deepcopy(adata_merge.uns['spatial'][list(adata_merge.uns['spatial'])[0]])
adata_merge.uns['spatial']['spatial_Merged']['images']["hires"] = np.asarray(new_im)
adata_merge.uns['spatial']['spatial_Merged']['images']["lowres"] = np.asarray(new_im)
adata_merge.uns['spatial']['spatial_Merged']['scalefactors']['tissue_lowres_scalef'] = 1
adata_merge.uns['spatial']['spatial_Merged']['scalefactors']['tissue_hires_scalef'] = 1


# add back the spatial coordinates as separate embeddings
idx = 0
adata_merge.obsm['X_spatial_Merged'] = adata_merge.obsm['spatial']
for i in range(0,size*2,size):
    for j in range(0,size*2,size):
        library_id='spatial_V1_'+samples[idx] # parse the string and get the sample id
        print(library_id)
        tissue_lowres_scalef = spatial[library_id]['scalefactors']['tissue_lowres_scalef']
        adata_merge.obsm['X_spatial_Merged'][np.where(adata_merge.obs['sample']==str(idx))] = copy.deepcopy(adatals[idx].obsm['spatial'])
        adata_merge.obsm['X_spatial_Merged'][np.where(adata_merge.obs['sample']==str(idx)),1] = adatals[idx].obsm['spatial'][:,1]*tissue_lowres_scalef - i
        adata_merge.obsm['X_spatial_Merged'][np.where(adata_merge.obs['sample']==str(idx)),0] = adatals[idx].obsm['spatial'][:,0]*tissue_lowres_scalef + j
        idx = idx+1    


## plotting out how the image and coordinates look
## note, in python the image and coordinates are NOT overlapping, but when load in VIP, they would. Reason Y coordinate values are negative
#library_id="spatial_Merged" # parse the string and get the sample id
#matplotlib.pyplot.figure()
#matplotlib.pyplot.imshow(spatial[library_id]['images']['lowres'])
#spatial1=adata_merge.obsm['X_'+library_id]
##matplotlib.pyplot.scatter(spatial1[:,0]*tissue_hires_scalef,2000-spatial1[:,1]*tissue_hires_scalef)
#matplotlib.pyplot.scatter(spatial1[:,0],spatial1[:,1])
#
adata_merge.write_h5ad(outputfile)
