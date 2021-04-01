import sys, getopt
import json
import matplotlib
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import copy
from PIL import Image, ImageOps

sc.logging.print_versions()

def parseArgs(argv):
    inputfile = ''
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:j:o:",["input=","json=","output="])
    except getopt.GetoptError:
        print('st_h5ad_image_operation.py -i <inputfile> -j <jsonfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('st_h5ad_image_operation.py -i <inputfile> -j <jsonfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-i", "--input"):
            inputfile = arg
        elif opt in ("-j", "--json"):
            jsonfile = arg
        elif opt in ("-o", "--output"):
            outputfile = arg
    print('Input file is "', inputfile)
    print('json file is "', jsonfile)
    print('Output file is "', outputfile)
    return(inputfile, jsonfile, outputfile)

def image_size(i):
    x_max = 0 # width of the merged image
    y_max = 0 # height of the merged image
    for t in i['layout']:
#         print(t)
        x_max = max(x_max, t['x'])
        y_max = max(y_max, t['y'])
#     print(x_max, y_max)
    return(x_max+1, y_max+1)

# flip x axis, new x is -x
def flipx_img(im):
    im_new=ImageOps.mirror(im)
    return(im_new)

# flip y axis, new y is -y
def flipy_img(im):
    im_new=ImageOps.flip(im)
    return(im_new)

# rotate image N degrees
def rotate_img(im, deg):
    im_new=im.rotate(-deg) # this function rotate counter clockwise, so using negative is making it clockwise
    return(im_new)

# flip x axis, new x is -x
def flipx_sp(i):
    i_new = copy.deepcopy(i)
    i_new[:,0] = -i_new[:,0]+600 # adding 600 is to make coordinates positive
    return(i_new)

# flip y axis, nw y is -y
def flipy_sp(i):
    i_new = copy.deepcopy(i)
    i_new[:,1] = -i_new[:,1]-600
    return(i_new)

# rotate sp N degrees
# only accept 90, 180 and 270 as options
def rotate_sp(i, deg):
    if (deg==90):
        # y <= -x and x <= y = 600
        i_new_tmp = copy.deepcopy(i)
        i_new = copy.deepcopy(i)
        i_new[:,0] = i_new_tmp[:,1] + 600
        i_new[:,1] = -i_new_tmp[:,0]
    if (deg==180):
        # x <= -x+600 and y <= -y-600 
        i_new_tmp = copy.deepcopy(i)
        i_new = copy.deepcopy(i)
        i_new[:,0] = -i_new_tmp[:,0] + 600
        i_new[:,1] = -i_new_tmp[:,1] - 600
    if (deg==270):
        # y <= -600 + x and x <= -y 
        i_new_tmp = copy.deepcopy(i)
        i_new = copy.deepcopy(i)
        i_new[:,0] = -i_new_tmp[:,1]
        i_new[:,1] = i_new_tmp[:,0] - 600
    return(i_new)

def standardizeImg(im):
    w,h = im.size
    delta_w = 600 - w
    delta_h = 600 - h
    padding = (0, 0, delta_w, delta_h)
    im_new = ImageOps.expand(im, padding)
    return(im_new)

def sp_ops(entry, adata):
    name = entry['name']
    # identify which sample this is
    sampleID = adata.obs['sample'][np.where((adata.obsm['X_spatial_'+name] != (0, 0)).all(axis=1))[0][0]]
    # extract coordinates
    coord = adata.obsm['X_spatial_'+name][np.where((adata.obsm['X_spatial_'+name] != (0, 0)).all(axis=1))]
    # find scale factor
    tissue_lowres_scalef = adata.uns['spatial']['spatial_'+name]['scalefactors']['tissue_lowres_scalef']
    coord1 = copy.deepcopy(coord)
    coord1[:,1] = coord[:,1]*tissue_lowres_scalef
    coord1[:,0] = coord[:,0]*tissue_lowres_scalef
    # coord1 = rotate_sp(coord1,270)
    img = Image.fromarray((adata.uns['spatial']["spatial_"+name]["images"]["lowres"]* 255).round().astype(np.uint8)) # found a solution to covert float32 to unit8
    w, h = img.size
    print(w, h)
    # all images pad to 600x600 if not 
    img1 = standardizeImg(img)
    if entry['flipx']==0 and entry['flipy']==0 and entry['rotate']==0: # no operation is needed
        img1=img1
        coord1=coord1
    if entry['flipx']==1:
        img1=flipx_img(img1)
        coord1=flipx_sp(coord1)
    if entry['flipy']==1:
        img1=flipy_img(img1)
        coord1=flipy_sp(coord1)
    if entry['rotate']!=0:
        img1=rotate_img(img1,entry['rotate'])
        coord1=rotate_sp(coord1,entry['rotate'])
    return(img1, coord1)


(inputfile, jsonfile, outputfile) = parseArgs(sys.argv[1:])
with open(jsonfile) as f:
    operation = json.load(f)
print(operation)

adata=sc.read_h5ad(inputfile)
imgs=list()
coords=list()
for i in operation['layout']:
    (img2, coord2) = sp_ops(i, adata)
    #matplotlib.pyplot.figure()
    #matplotlib.pyplot.imshow(img2)
    #matplotlib.pyplot.scatter(coord2[:,0],coord2[:,1])
    imgs.append(img2)
    coords.append(coord2)

# stich images
idx=0
size=640
dim = image_size(operation)
new_im = Image.new('RGB', (size*dim[0],size*dim[1]))
for i in operation['layout']:
    print(i)
    im = imgs[idx]
    width, height = im.size
    new_im.paste(im, (i['x']*size,i['y']*size))
    idx = idx+1

# fake a adata.uns by providing merged lowres image and scale factors 1
adata.uns['spatial']['spatial_Merged_customized'] = copy.deepcopy(adata.uns['spatial'][list(adata.uns['spatial'])[0]])
adata.uns['spatial']['spatial_Merged_customized']['images']["hires"] = np.asarray(new_im)
adata.uns['spatial']['spatial_Merged_customized']['images']["lowres"] = np.asarray(new_im)
adata.uns['spatial']['spatial_Merged_customized']['scalefactors']['tissue_lowres_scalef'] = 1
adata.uns['spatial']['spatial_Merged_customized']['scalefactors']['tissue_hires_scalef'] = 1

# stich coords
idx = 0
adata.obsm['X_spatial_Merged_customized'] = copy.deepcopy(adata.obsm['spatial'])
for i in operation['layout']:
    coord = coords[idx]
    name=operation['layout'][idx]['name']
    sampleID = adata.obs['sample'][np.where((adata.obsm['X_spatial_'+name] != (0, 0)).all(axis=1))[0][0]]
    adata.obsm['X_spatial_Merged_customized'][np.where(adata.obs['sample']==str(sampleID)),1] = coord[:,1] - i['y']*size # y axis values are all negative so need to substract
    adata.obsm['X_spatial_Merged_customized'][np.where(adata.obs['sample']==str(sampleID)),0] = coord[:,0] + i['x']*size
    idx = idx+1


adata.write_h5ad(outputfile)
