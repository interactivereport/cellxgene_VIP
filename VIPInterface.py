import requests
import json
import traceback
import sqlite3
#import server.app.decode_fbs as decode_fbs
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import diffxpy.api as de
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rcParams
import plotly.graph_objects as go
import plotly.io as plotIO
import base64
import math
import scipy
from io import BytesIO
import sys
import time
import os
import re
import glob
import subprocess
import gc #,psutil
import errno
from pathlib import Path

strExePath = os.path.dirname(os.path.abspath(__file__))

import pprint
ppr = pprint.PrettyPrinter(depth=6,width=500)

import server.common.compute.diffexp_generic as diffDefault
import pickle
from pyarrow import feather

sys.setrecursionlimit(10000)
sc.settings.verbosity = 2
rcParams.update({'figure.autolayout': True})

api_version = "/api/v0.2"

import threading
jobLock = threading.Lock()
def getLock(lock):
    while not lock.acquire():
        time.sleep(1.0)
def freeLock(lock):
  lock.release()

def route(data,appConfig):
  data = initialization(data,appConfig)
  #ppr.pprint(appConfig.server_config.single_dataset__datapath)
  #ppr.pprint(data)
  try:
    getLock(jobLock)
    setTimeStamp(data)
    taskRes = distributeTask(data["method"])(data)
    freeLock(jobLock)
    gc.collect()
    #ppr.pprint("memory usage: rss (%dM) and vms (%dM)"%(int(psutil.Process().memory_info().rss / 1024 **2),
    #                                                    int(psutil.Process().memory_info().vms / 1024 **2)))
    return taskRes
  except Exception as e:
    freeLock(jobLock)
    return 'ERROR @server: '+traceback.format_exc()

import server.app.app as app
class getAdapter(object):
    def __init__(self, adapter):
        self.adapter = adapter
    def __enter__(self):
        return self.adapter
    def __exit__(self):
        return

def initialization(data,appConfig):
  # obtain the server host information
  data = json.loads(str(data,encoding='utf-8'))

  # update the environment information
  data.update(VIPenv)

  # updatting the hosting data information
  # v1.1.1 no multi datasets option
  # if appConfig.is_multi_dataset():
  #  data["url_dataroot"]=appConfig.server_config.multi_dataset__dataroot['d']['base_url']
  #  data['h5ad']=os.path.join(appConfig.server_config.multi_dataset__dataroot['d']['dataroot'], data["dataset"])
  # else:
  #  data["url_dataroot"]=None
  #  data["dataset"]=None
  #  data['h5ad']=appConfig.server_config.single_dataset__datapath
  data["url_dataroot"]=None
  data["dataset"]=None
  data['h5ad']=appConfig.server_config.single_dataset__datapath
  # setting the plotting options
  if 'figOpt' in data.keys():
    setFigureOpt(data['figOpt'])

  # get the var (gene) and obv index
  #ppr.pprint(dir(appConfig.server_config.data_adaptor))
  #with app.get_data_adaptor() as scD:
  data['data_adapter'] = appConfig.server_config.data_adaptor
  scD=appConfig.server_config.data_adaptor
  data['obs_index'] = scD.get_schema()["annotations"]["obs"]["index"]
  data['var_index'] = scD.get_schema()["annotations"]["var"]["index"]
  return data

def setFigureOpt(opt):
  sc.set_figure_params(dpi_save=int(opt['dpi']),fontsize= float(opt['fontsize']),vector_friendly=(opt['vectorFriendly'] == 'Yes'),transparent=(opt['transparent'] == 'Yes'),color_map=opt['colorMap'])
  rcParams.update({'savefig.format':opt['img']})

def getObs(data):
  selC = list(data['cells'].values())
  cNames = ["cell%d" %i for i in selC]
  ## obtain the category annotation
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    selAnno = [data['obs_index']]+data['grp']
    dAnno = list(scD.get_obs_keys())
    anno = []
    sel = list(set(selAnno)&set(dAnno))
    if len(sel)>0:
      tmp = scD.data.obs.iloc[selC,:].loc[:,sel].astype('str')
      tmp.index = cNames
      anno += [tmp]
    sel = list(set(selAnno)-set(dAnno))
    if len(sel)>0:
      annotations = scD.dataset_config.user_annotations
      if annotations:
        labels = annotations.read_labels(scD)
        tmp = labels.loc[list(scD.data.obs.loc[selC,data['obs_index']]),sel]
        tmp.index = cNames
        anno += [tmp]
    obs = pd.concat(anno,axis=1)
  #ppr.pprint(obs)
  ## update the annotation Abbreviation
  combUpdate = cleanAbbr(data)
  if 'abb' in data.keys():
    for i in data['grp']:
      obs[i] = obs[i].map(data['abb'][i])
  return combUpdate, obs

def getObsNum(data):
  selC = list(data['cells'].values())
  cNames = ["cell%d" %i for i in selC]
  ## obtain the category annotation
  obs = pd.DataFrame()
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    selAnno = data['grpNum']
    dAnno = list(scD.get_obs_keys())
    sel = list(set(selAnno)&set(dAnno))
    if len(sel)>0:
      obs = scD.data.obs.loc[selC,sel]
      obs.index = cNames
  return obs

def getVar(data):
  ## obtain the gene annotation
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    gInfo = scD.data.var
  gInfo.index = list(gInfo[data['var_index']])
  gInfo = gInfo.drop([data['var_index']],axis=1)
  return gInfo

def collapseGeneSet(data,expr,gNames,cNames,fSparse):
  Y = expr
  if 'geneGrpColl' in data.keys() and not data['geneGrpColl']=='No' and 'geneGrp' in data.keys() and len(data['geneGrp'])>0:
    data['grpLoc'] = []
    data['grpID'] = []
    if fSparse:
      Y = pd.DataFrame.sparse.from_spmatrix(Y,columns=gNames,index=cNames)
    for aN in data['geneGrp'].keys():
      if data['geneGrpColl']=='mean':
        Y = pd.concat([Y,Y[data['geneGrp'][aN]].mean(axis=1).rename(aN)],axis=1,sort=False)
      if data['geneGrpColl']=='median':
        Y = pd.concat([Y,Y[data['geneGrp'][aN]].median(axis=1).rename(aN)],axis=1,sort=False)
      for gene in data['geneGrp'][aN]:
        if gene in data['genes']:
          data['genes'].remove(gene)
      data['genes'] += [aN]
    gNames = list(Y.columns)
  return Y,gNames

def createData(data):
  selC = list(data['cells'].values())
  cNames = ["cell%d" %i for i in selC]

  ## onbtain the expression matrix
  gNames = []
  expr = []
  fSparse = False
  X = []
  if 'genes' in data.keys():
    #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
    scD=data['data_adapter']
    if scD is not None:
      if not type(scD.data.X) is np.ndarray:
        fSparse = True
      if len(data['genes'])>0:
        fullG = list(scD.data.var[data['var_index']])
        selG = sorted([fullG.index(i) for i in data['genes']]) #when data loaded backed, incremental is required
        X = scD.data.X[:,selG]
        gNames = [fullG[i] for i in selG] #data['genes']
      else:
        X = scD.data.X
        gNames = list(scD.data.var[data['var_index']])
    if 'figOpt' in data.keys() and data['figOpt']['scale'] == 'Yes':
      X = sc.pp.scale(X,zero_center=(data['figOpt']['scaleZero'] == 'Yes'),max_value=(float(data['figOpt']['scaleMax']) if data['figOpt']['clipValue']=='Yes' else None))
    X = X[selC]
  if fSparse:
    expr = X
  else:
    expr = pd.DataFrame(X,columns=gNames,index=cNames)

  expr,gNames = collapseGeneSet(data,expr,gNames,cNames,fSparse)
  #ppr.pprint("finished expression ...")
  ## obtain the embedding
  embed = {}
  if 'layout' in data.keys():
    layout = data['layout']
    if isinstance(layout,str):
      layout = [layout]
    if len(layout)>0:
      for one in layout:
        #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
        scD=data['data_adapter']
        if scD is not None:
          embed['X_%s'%one] = pd.DataFrame(scD.data.obsm['X_%s'%one][selC][:,[0,1]],columns=['%s1'%one,'%s2'%one],index=cNames)
  #ppr.pprint("finished layout ...")
  ## obtain the category annotation
  combUpdate, obs = getObs(data)

  ## create a custom annotation category and remove cells which are not in the selected annotation
  if combUpdate and len(data['grp'])>1:
    newGrp = 'Custom_combine'
    combineGrp = list(data['combine'].keys());
    obs[newGrp] = obs[combineGrp[0]]
    for i in combineGrp:
      if not i==combineGrp[0]:
        obs[newGrp] += ":"+obs[i]
    selC = ~obs[newGrp].str.contains("Other").to_numpy()
    expr = expr[selC]
    for i in embed.keys():
      embed[i] = embed[i][selC]
    obs = obs[selC].astype('category')
    obs[newGrp].cat.set_categories(data['combineOrder'],inplace=True)
    data['grp'] = [newGrp]

  obs = obs.astype('category')
  ## empty selection
  if expr.shape[0]==0 or expr.shape[1]==0:
    return []
  #ppr.pprint("finished obv ...")

  return sc.AnnData(expr,obs,var=pd.DataFrame([],index=gNames),obsm={layout:embed[layout].to_numpy() for layout in embed.keys()})

def cleanAbbr(data):
  updated = False
  if 'abb' in data.keys() and 'combine' in data.keys():
    if len(data['combine'])>0:
      updated = True
      for cate in data['abb'].keys():
        if cate in data['combine'].keys():
          for anName in data['abb'][cate].keys():
            if not anName in data['combine'][cate]:
              data['abb'][cate][anName] = "Other";
            else:
              if not data['abb'][cate][anName]==anName:
                data['combineOrder'] = [one.replace(anName,data['abb'][cate][anName]) for one in data['combineOrder']]
        else:
          data['abb'][cate] = {key:"Other" for key in data['abb'][cate].keys()}
  return updated

def errorTask(data):
  raise ValueError('Error task!')

def distributeTask(aTask):
  return {
    'SGV':SGV,
    'SGVcompare':SGVcompare,
    'PGV':PGV,
    'PGVcompare':PGVcompare,
    'VIOdata':VIOdata,
    'HEATplot':pHeatmap,
    'HEATdata':HeatData,
    'GD':GD,
    'DEG':DEG,
    'DOT':DOT,
    'DOTdata':DOTdata,
    'EMBED':EMBED,
    'TRAK':TRACK,
    'DUAL':DUAL,
    'MARK': MARK,
    'MINX':MINX,
    'DENS':DENS,
    'DENS2D':DENS2D,
    'SANK':SANK,
    'STACBAR':STACBAR,
    'HELLO':HELLO,
    'CLI':CLI,
    'preDEGname':getPreDEGname,
    'preDEGvolcano':getPreDEGvolcano,
    'preDEGmulti':getPreDEGbubble,
    'preDEGfgsea':getPreDEGfgsea,
    'mergeMeta': mergeMeta,
    'isMeta': isMeta,
    'testVIPready':testVIPready,
    'Description':getDesp,
    'GSEAgs':getGSEA,
    'SPATIAL':SPATIAL,
    'saveTest':saveTest,
    'getBWinfo':getBWinfo,
    'GSP':GSP,
    'plotBW':plotBW,
    'checkCosMx':getCosMx,
    'plotCOSMX':plotCosMx,
    'checkVisium':getVisium,
    'plotVisium':plotVisium
  }.get(aTask,errorTask)

def HELLO(data):
  return 'Hi, connected.'

def iostreamFig(fig):
  #getLock(iosLock)
  figD = BytesIO()
  #ppr.pprint('io located at %d'%int(str(figD).split(" ")[3].replace(">",""),0))
  fig.savefig(figD,bbox_inches="tight")
  #ppr.pprint(sys.getsizeof(figD))
  #ppr.pprint('io located at %d'%int(str(figD).split(" ")[3].replace(">",""),0))
  imgD = base64.encodebytes(figD.getvalue()).decode("utf-8")
  figD.close()
  #ppr.pprint("saved Fig")
  #freeLock(iosLock)
  if 'matplotlib' in str(type(fig)):
    plt.close(fig)#'all'
  return imgD

def Msg(msg):
  fig = plt.figure(figsize=(5,2))
  plt.text(0,0.5,msg)
  ax = plt.gca()
  ax.axis('off')
  return iostreamFig(fig)

def SPATIAL(data):
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    #ppr.pprint(vars(scD.data.uns["spatial"]))
    spatial=scD.data.uns["spatial"]
    if (data['embedding'] == "get_spatial_list"):
      return json.dumps({'list':list(spatial)})
    library_id=list(spatial)[0]
    if (data['embedding'] in list(spatial)):
      library_id=data['embedding']

    height, width, depth = spatial[library_id]["images"][data['resolution']].shape

    embedding = 'X_'+data['embedding']
    spatialxy = scD.data.obsm[embedding]
    tissue_scalef = spatial[library_id]['scalefactors']['tissue_' + data['resolution'] + '_scalef']
    i = data['spots']['spoti_i']
    x = 0
    y = 1
    # from original embedding to (0,1) coordinate system (cellxgene embedding)
    scalex = (data['spots']['spot0_x'] - data['spots']['spoti_x']) / (spatialxy[0][x] - spatialxy[i][x])
    scaley = (data['spots']['spot0_y'] - data['spots']['spoti_y']) / (spatialxy[0][y] - spatialxy[i][y])

    # image is in (-1,0,1) coordinate system, so multiplied by 2
    translatex = (spatialxy[i][x]*scalex - data['spots']['spoti_x']) * 2
    translatey = (spatialxy[i][y]*scaley - data['spots']['spoti_y']) * 2
    scale = 1/tissue_scalef * scalex * 2
    # Addtional translate in Y due to flipping of the image if needed
    ppr.pprint(scalex)
    ppr.pprint(scaley)
    ppr.pprint(translatex)
    ppr.pprint(translatey)

    # from (-1,0,1) (image layer) to (0,1) coordinate system (cellxgene embedding). Overlapping (0,0) origins of both.
    translatex = -(1+translatex)
    if (translatey > -0.1):
      flip = True
      translatey = -(1+translatey) + height*scale
    else:
      flip = False
      translatey = -(1+translatey)

    returnD = [{'translatex':translatex,'translatey':translatey,'scale':scale}]

    dpi=100
    figsize = width / float(dpi), height / float(dpi)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    if (flip):
      ax.imshow(np.flipud(spatial[library_id]["images"][data['resolution']]), interpolation='nearest')
    else:
      ax.imshow(spatial[library_id]["images"][data['resolution']], interpolation='nearest')

    figD = BytesIO()
    plt.savefig(figD, dpi=dpi)
    ppr.pprint(sys.getsizeof(figD))
    imgD = base64.encodebytes(figD.getvalue()).decode("utf-8")
    figD.close()
    plt.close(fig)
  return json.dumps([returnD, imgD])

def MINX(data):
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    minV = min(scD.data.X[0])
  return '%.1f'%minV

def geneFiltering(adata,cutoff,opt):
  ## 1. remove cells if the max expression of all genes is lower than the cutoff
  if opt==1:
    #sT = time.time()
    #ix = adata.to_df().apply(lambda x: max(x)>float(cutoff),axis=1)
    #ppr.pprint(time.time()-sT)
    #sT=time.time()
    df = adata.to_df()
    ix = df[df>float(cutoff)].count(axis=1)>0
    #ppr.pprint(time.time()-sT)
    #sT = time.time()
    #ix = pd.DataFrame((adata.X>float(cutoff)).sum(1)>0,index=list(adata.obs.index)).iloc[:,0]
    #ppr.pprint(time.time()-sT)

    adata = adata[ix,]
  ## 2. Set all expression level smaller than the cutoff to be NaN not for plotting without removing any cells
  elif opt==2:
    def cutoff(x):
        return x if x>float(cutoff) else None
    X = adata.to_df()
    X=X.applymap(cutoff)
    adata = sc.AnnData(X,adata.obs)
  return adata

def SGV(data):
  # figure width and heights depends on number of unique categories
  # characters of category names, gene number
  #ppr.pprint("SGV: creating data ...")
  adata = createData(data)
  #ppr.pprint("SGV: data created ...")
  adata = geneFiltering(adata,data['cutoff'],1)
  if len(adata)==0:
    raise ValueError('No cells in the condition!')
  a = list(set(list(adata.obs[data['grp'][0]])))
  ncharA = max([len(x) for x in a])
  w = len(a)/4+1
  h = ncharA/6+2.5
  ro = math.acos(10/max([15,ncharA]))/math.pi*180
  ##
  fig = plt.figure(figsize=[w,h])
  if data.get('dotsize') is None or float(data['dotsize'])==0:
    sc.pl.violin(adata,data['genes'],groupby=data['grp'][0],ax=fig.gca(),show=False,stripplot=False)
  else:
    sc.pl.violin(adata,data['genes'],groupby=data['grp'][0],ax=fig.gca(),show=False,size=float(data['dotsize']))
  #sc.pl.violin(adata,data['genes'],groupby=data['grp'][0],ax=fig.gca(),show=False,size=0.5)
  fig.autofmt_xdate(bottom=0.2,rotation=ro,ha='right')
  return iostreamFig(fig)

def SGVcompare(data):
  adata = createData(data)
  #adata = geneFiltering(adata,data['cutoff'],1)
  if len(adata)==0:
    raise ValueError('No cells in the condition!')

  # plot in R
  strF = ('%s/SGV%f.csv' % (data["CLItmp"],time.time()))
  X=pd.concat([adata.to_df(),adata.obs[data['grp']]],axis=1,sort=False)
  X[X.iloc[:,0]>=float(data['cellCutoff'])].to_csv(strF,index=False)
  #strCMD = " ".join(["%s/Rscript"%data['Rpath'],strExePath+'/violin.R',strF,str(data['cutoff']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']])
  #ppr.pprint(" ".join([strExePath+'/violin.R',strF,str(data['cutoff']),str(data['dotsize']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']]))
  res = subprocess.run([strExePath+'/violin.R',strF,str(data['cutoff']),str(data['dotsize']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']],capture_output=True)#
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in R: "+res.stderr.decode('utf-8'))

  return img

def VIOdata(data):
  adata = createData(data)
  adata = geneFiltering(adata,data['cutoff'],1)
  if len(adata)==0:
    raise ValueError('No cells in the condition!')
  return pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
def updateGene(data):
  grpID = []
  grpLoc=[]
  allG = []
  if 'geneGrp' in data.keys():
    for aN in data['geneGrp'].keys():
      grpLoc += [(len(allG),len(allG)+len(data['geneGrp'][aN])-1)]
      allG += data['geneGrp'][aN]
      grpID += [aN]

  data['genes'] = unique(allG+data['genes'])
  data['grpLoc'] = grpLoc
  data['grpID'] = grpID

def PGV(data):
  # figure width and heights depends on number of unique categories
  # characters of category names, gene number #pecam1 pdpn
  updateGene(data)
  #ppr.pprint("PGV: creating data ...")
  adata = createData(data)
  #ppr.pprint("PGV: data created ...")
  adata = geneFiltering(adata,data['cutoff'],1)
  if adata.shape[0]==0 or adata.shape[1]==0:
    return Msg('No cells in the condition!')
  a = list(set(list(adata.obs[data['grp'][0]])))
  ncharA = max([len(x) for x in a])
  w = max([3,ncharA/8])+len(data['genes'])/2+1.5
  h = len(a)+0.5
  swapAx = False
  ##
  if data['by']=='Columns':
    a = w
    w = h
    h = a
    swapAx = True
  if 'split_show' in data['figOpt']['scanpybranch']: #.dev140+ge9cbc5f
    vp = sc.pl.stacked_violin(adata,data['genes'],groupby=data['grp'][0],return_fig=True,figsize=(w,h),swap_axes=swapAx,var_group_positions=data['grpLoc'],var_group_labels=data['grpID'])
    vp.add_totals().style(yticklabels=True, cmap=data['color']).show()
    #vp.add_totals().show()
    fig = vp#plt.gcf()
  else:
    fig = plt.figure(figsize=[w,h])
    axes = sc.pl.stacked_violin(adata,data['genes'],groupby=data['grp'][0],show=False,ax=fig.gca(),swap_axes=swapAx,
                                var_group_positions=data['grpLoc'],var_group_labels=data['grpID'])
  return iostreamFig(fig)

def PGVcompare(data):
  adata = createData(data)
  #adata = geneFiltering(adata,data['cutoff'],1)
  sc.pp.filter_cells(adata,min_counts=float(data['cutoff']))
  if len(adata)==0:
    raise ValueError('No cells in the condition!')
  strF = ('%s/PGVcompare%f.csv' % (data["CLItmp"],time.time()))
  X=pd.concat([adata.to_df(),adata.obs[data['grp']]],axis=1,sort=False).to_csv(strF,index_label="cellID")

  # plot in R
  strCMD = " ".join(["Rscript",strExePath+'/complex_vlnplot_multiple.R',strF,','.join(data['genes']),data['grp'][0],data['grp'][1],str(data['width']),str(data['height']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']])
  res = subprocess.run(strCMD,shell=True,capture_output=True)
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in R: "+res.stderr.decode('utf-8'))
  return img

def silentRM(strF):
  try:
    os.remove(strF)
  except OSError as e:
    if e.errno != errno.ENOENT:
      raise

def pHeatmap(data):
  # figure width is depends on the number of categories was choose to show
  # and the character length of each category term
  # if the number of element in a category is smaller than 10, "Set1" or "Set3" is choosen
  # if the number of element in a category is between 10 and 20, default is choosen
  # if the number of element in a category is larger than 20, husl is choosen
  #Xsep = createData(data,True)
  #adata = sc.AnnData(Xsep['expr'],Xsep['obs'])
  #sT = time.time()

  adata = createData(data)
  if len(data['grpNum'])>0:
    adata.obs = pd.concat([adata.obs,getObsNum(data)],axis=1)
  data['grp'] += data['addGrp']
  #Xdata = pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()
  #ppr.pprint('HEAT data reading cost %f seconds' % (time.time()-sT) )
  #sT = time.time()


  h = 8
  w = len(data['genes'])/3+0.3 + 2
  heatCol=data['color']
  Zscore=None
  heatCenter=None
  colTitle="Expression"

  if data['plotMethod']=='sns':
    exprOrder = True
    if data['order'][0]!="Expression":
      exprOrder = False;
      adata = adata[adata.obs.sort_values(data['order']).index,]
      #s = adata.obs[data['order']]
      #ix = sorted(range(len(s)), key=lambda k: s[k])
      #adata = adata[ix,]
    colCounter = 0
    colName =['Set1','Set3']
    grpCol = list()
    grpLegend = list()
    grpWd = list()
    grpLen = list()
    for gID in data['grp']:
      grp = adata.obs[gID]
      Ugrp = grp.unique()
      if len(Ugrp)<10:
        lut = dict(zip(Ugrp,sns.color_palette(colName[colCounter%2],len(Ugrp)).as_hex()))
        colCounter += 1
      elif len(Ugrp)<20:
        lut = dict(zip(Ugrp,sns.color_palette(n_colors=len(Ugrp)).as_hex()))
      else:
        lut = dict(zip(Ugrp,sns.color_palette("husl",len(Ugrp)).as_hex()))
      grpCol.append(grp.map(lut))
      grpLegend.append([mpatches.Patch(color=v,label=k) for k,v in lut.items()])
      grpWd.append(max([len(x) for x in Ugrp]))#0.02*fW*max([len(x) for x in Ugrp])
      grpLen.append(len(Ugrp)+2)
    if data['norm']=='zscore':
      Zscore=1
      #heatCol="vlag"
      heatCenter=0
      colTitle="Z-score"
    #ppr.pprint('HEAT data preparing cost %f seconds' % (time.time()-sT) )
    #sT = time.time()

    try:
      g = sns.clustermap(adata.to_df(),
                       method="ward",row_cluster=exprOrder,z_score=Zscore,cmap=heatCol,center=heatCenter,
                       row_colors=pd.concat(grpCol,axis=1).astype('str'),yticklabels=False,xticklabels=True,
                       figsize=(w,h),colors_ratio=0.05,
                       cbar_pos=(.3, .95, .55, .02),
                       cbar_kws={"orientation": "horizontal","label": colTitle,"shrink": 0.5})
    except Exception as e:
      return 'ERROR: Z score calculation failed for 0 standard diviation. '+traceback.format_exc() # 'ERROR @server: {}, {}'.format(type(e),str(e))


    #ppr.pprint('HEAT plotting cost %f seconds' % (time.time()-sT) )
    #sT = time.time()
    g.ax_col_dendrogram.set_visible(False)
    #g.ax_row_dendrogram.set_visible(False)
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    grpW = [1.02]
    grpH = [1.2]
    cumulaN = 0
    cumulaMax = 0
    characterW=1/40 # a character is 1/40 of heatmap width
    characterH=1/40 # a character is 1/40 of heatmap height
    for i in sorted(range(len(grpLen)),key=lambda k:grpLen[k]):#range(5):#
      cumulaN += grpLen[i]
      if cumulaN>(10+1/characterH):
          grpW.append(grpW[-1]+cumulaMax)
          grpH = [1.2]
          cumulaN =0
          cumulaMax=0
      leg = g.ax_heatmap.legend(handles=grpLegend[i],frameon=True,title=data['grp'][i],loc="upper left",
                                bbox_to_anchor=(grpW[-1],grpH[-1]),fontsize=5)#grpW[i],0.5,0.3
      #leg = g.ax_heatmap.legend(handles=grpLegend[0],frameon=True,title=data['grp'][0],loc="upper left",
      #                          bbox_to_anchor=(1.02,1-i*0.25),fontsize=5)#grpW[i],0.5,0.
      cumulaMax = max([cumulaMax,grpWd[i]*characterW])
      grpH.append(grpH[-1]-grpLen[i]*characterH)
      leg.get_title().set_fontsize(6)#min(grpSize)+2
      g.ax_heatmap.add_artist(leg)
    #ppr.pprint('HEAT post plotting cost %f seconds' % (time.time()-sT) )
    return iostreamFig(g)#json.dumps([iostreamFig(g),Xdata])#)#
  elif data['plotMethod']=='cHeatmap':
    characterW=1/40
    grpW = list()
    grpN = list()
    for gID in data['grp']:
      grp = adata.obs[gID]
      Ugrp = grp.unique()
      grpN.append(1+len(Ugrp))
      grpW.append(max([len(x) for x in Ugrp]))
    legendCol = math.ceil(sum(grpN)/40)
    w=w+characterW*sum(sorted(grpW,reverse=True)[0:legendCol])
    strF = ('%s/HEAT%f.csv' % (data["CLItmp"],time.time()))
    D = adata.to_df()
    if data['norm']=='zscore':
      for one in D.columns:
        D[one] = (D[one]-D[one].mean())/D[one].std()
      colTitle="Z-score"
    D = pd.concat([D,adata.obs],axis=1,sort=False)
    D.to_csv(strF,index=False)

    if len(data['gAnno'])>0:
      gAnno = pd.DataFrame(data['gAnno'])
      gAnno.to_csv(strF.replace("HEAT","HEATgene"),index=False)
    elif data['gAnnoDef']:
      gAnno = pd.read_csv("%s/proteinatlas_protein_class.csv"%strExePath,index_col=0,header=0)
      gNames = {i:i.upper() for i in data['genes']}
      gAnno = gAnno.loc[list(set(gNames.values()) & set(gAnno.index)),:].iloc[:,1:]
      for one in [i for i in gNames.values() if i not in gAnno.index]:
        gAnno.loc[one,:] = np.nan
      gAnno = gAnno.loc[gNames.values(),:]
      gAnno.index=gNames.keys()
      gAnno.fillna("N",inplace=True)
      gAnno.to_csv(strF.replace("HEAT","HEATgene"))
      #ppr.pprint(gAnno)
    ## plot in R
    cmd = "%s/complexHeatmap.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s '%s'"%(strExePath,strF,','.join(data['genes']),colTitle,','.join(data['order']),str(data['width']),str(data['height']),heatCol,
      str(data['legendRow']),str(data['fontadj']),
      data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),str(data['swapAxes']),data['figOpt']['vectorFriendly'],data['Rlib'])

    #ppp = pprint.PrettyPrinter(depth=6,width=300)
    #ppp.pprint(cmd)
    #return SyntaxError("in R: ")
    res = subprocess.run(cmd,check=True,shell=True,capture_output=True)#
    img = res.stdout.decode('utf-8')
    os.remove(strF)
    #silentRM(strF.replace("HEAT","HEATgene"))
    if 'Error' in res.stderr.decode('utf-8'):
      raise SyntaxError("in R: "+res.stderr.decode('utf-8'))
    return img
  else:
    raise ValueError('Unknown heatmap plotting method (%s)!'%data['plotMethod'])

def HeatData(data):
  adata = createData(data)
  Xdata = pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()
  return Xdata

def GD(data):
  adata = None;
  for one in data['cells'].keys():
    #sT = time.time()
    oneD = data.copy()
    oneD.update({'cells':data['cells'][one],
            'genes':[],
            'grp':[]})
    D = createData(oneD)
    #ppr.pprint("one grp aquire data cost %f seconds" % (time.time()-sT))
    D.obs['cellGrp'] = one
    if adata is None:
      adata = D
    else:
      #sT =time.time()
      adata = adata.concatenate(D)
      #ppr.pprint("Concatenate data cost %f seconds" % (time.time()-sT))
  if adata is None:
    return Msg("No cells were satisfied the condition!")

  ##
  adata.obs.astype('category')
  cutOff = 'geneN_cutoff'+data['cutoff']
  #sT = time.time()
  #adata.obs[cutOff] = adata.to_df().apply(lambda x: sum(x>float(data['cutoff'])),axis=1)
  #ppr.pprint(time.time()-sT)
  #sT = time.time()
  #df = adata.to_df()
  #adata.obs[cutOff] = df[df>float(data['cutoff'])].count(axis=1)
  #ppr.pprint(time.time()-sT)
  sT = time.time()
  adata.obs[cutOff] = (adata.X >float(data['cutoff'])).sum(1)
  ppr.pprint(time.time()-sT)
  ##
  w = 3
  if len(data['cells'])>1:
    w += 3
  fig = plt.figure(figsize=[w,4])
  sc.pl.violin(adata,cutOff,groupby='cellGrp',ax=fig.gca(),show=False,rotation=0,size=2)
  return iostreamFig(fig)

def getGSEA(data):
  strGSEA = '%s/gsea/'%strExePath
  return json.dumps(sorted([os.path.basename(i).replace(".symbols.gmt","") for i in glob.glob(strGSEA+"*.symbols.gmt")]))

def DEG(data):
  adata = None
  genes = data['genes']
  data['genes'] = []
  comGrp = 'cellGrp'
  if 'combine' in data.keys():
    if data['DEmethod']=='default':
      combUpdate, obs = getObs(data)
      if combUpdate and len(data['grp'])>1:
        obs[comGrp] = obs[data['grp'][0]]
        for i in data['grp']:
          if i!=data['grp'][0]:
            obs[comGrp] += ":"+obs[i]
      mask = [obs[comGrp].isin([data['comGrp'][i]]) for i in [0,1]]
    else:
      data['figOpt']['scale'] = 'No'
      adata = createData(data)
      comGrp = data['grp'][0]
      adata = adata[adata.obs[comGrp].isin(data['comGrp'])]
  else:
    mask = [pd.Series(range(data['cellN'])).isin(data['cells'][one].values()) for one in data['comGrp']]
    for one in data['comGrp']:
      oneD = data.copy()
      oneD['cells'] = data['cells'][one]
      oneD['genes'] = []
      oneD['grp'] = []
      oneD['figOpt']['scale']='No'
      #oneD = {'cells':data['cells'][one],
      #        'genes':[],
      #        'grp':[],
      #        'figOpt':{'scale':'No'},
      #        'url':data['url']}

      D = createData(oneD)
      D.obs[comGrp] = one
      if adata is None:
        adata = D
      else:
        adata = adata.concatenate(D)

  if data['DEmethod']=='default':
    if sum(mask[0]==True)<10 or sum(mask[1]==True)<10:
      raise ValueError('Less than 10 cells in a group!')
    #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
    scD=data['data_adapter']
    if scD is not None:
      #res = diffDefault.diffexp_ttest(scD,mask[0].to_numpy(),mask[1].to_numpy(),scD.data.shape[1])# shape[cells as rows, genes as columns]
      res = scD.compute_diffexp_ttest(mask[0].to_numpy(),mask[1].to_numpy(),scD.data.shape[1]-1,0.01)
      gNames = list(scD.data.var[data['var_index']])
    deg = pd.DataFrame(res['positive']+res['negative'],columns=['gID','log2fc','pval','qval'])
    deg = deg.drop_duplicates(ignore_index=True)
    gName = pd.Series([gNames[i] for i in deg['gID']],name='gene')
    deg = pd.concat([deg,gName],axis=1).loc[:,['gene','log2fc','pval','qval']]
  else:
    if not 'AnnData' in str(type(adata)):
      raise ValueError('No data extracted by user selection')
    adata.obs.astype('category')
    nm = None
    if data['DEmethod']=='wald':
      nm = 'nb'
    if data['DEmethod']=='wald':
        res = de.test.wald(adata,formula_loc="~1+"+comGrp,factor_loc_totest=comGrp)
    elif data['DEmethod']=='t-test':
        res = de.test.t_test(adata,grouping=comGrp)
    elif data['DEmethod']=='rank':
        res = de.test.rank_test(adata,grouping=comGrp)
    else:
        raise ValueError('Unknown DE methods:'+data['DEmethod'])
    #res = de.test.two_sample(adata,comGrp,test=data['DEmethod'],noise_model=nm)
    deg = res.summary()
    deg['log2fc'] = -1 * deg['log2fc']
  deg = deg.sort_values(by=['qval']).loc[:,['gene','log2fc','pval','qval']]
  #del adata
  ## plot in R
  #strF = ('/tmp/DEG%f.csv' % time.time())
  strF = ('%s/DEG%f.csv' % (data["CLItmp"],time.time()))
  deg.to_csv(strF,index=False)
  #ppr.pprint([strExePath+'/volcano.R',strF,'"%s"'%';'.join(genes),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),str(data['logFC']),data['comGrp'][1],data['comGrp'][0]])
  res = subprocess.run([strExePath+'/volcano.R',strF,';'.join(genes),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),str(data['logFC']),data['comGrp'][1],data['comGrp'][0],str(data['sigFDR']),str(data['sigFC']),str(data['labelSize']),str(data['dotSize']),str(data['ymin']),str(data['ymax']),data['figOpt']['vectorFriendly'],data['Rlib']],capture_output=True)#
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in volcano.R: "+res.stderr.decode('utf-8'))
  img = res.stdout.decode('utf-8')

  # GSEA
  GSEAimg=""
  GSEAtable=pd.DataFrame()
  if data['gsea']['enable']:
    res = subprocess.run([strExePath+'/fgsea.R',
                          strF,
                          '%s/gsea/%s.symbols.gmt'%(strExePath,data['gsea']['gs']),
                          str(data['gsea']['gsMin']),
                          str(data['gsea']['gsMax']),
                          str(data['gsea']['padj']),
                          data['gsea']['up'],
                          data['gsea']['dn'],
                          str(data['gsea']['collapse']),
                          data['figOpt']['img'],
                          str(data['figOpt']['fontsize']),
                          str(data['figOpt']['dpi']),
                          data['Rlib']],capture_output=True)#
    if 'Error' in res.stderr.decode('utf-8'):
        raise SyntaxError("in fgsea.R: "+res.stderr.decode('utf-8'))
    GSEAimg = res.stdout.decode('utf-8')
    GSEAtable = pd.read_csv(strF)
    GSEAtable['leadingEdge'] = GSEAtable['leadingEdge'].apply(lambda x:'|'.join(x.split('|')[:10]))

  os.remove(strF)
  #####
  gInfo = getVar(data)
  deg.index = deg['gene']
  deg = pd.concat([deg,gInfo],axis=1,sort=False)
  #return deg.to_csv()

  if not data['topN']=='All':
    deg = deg.iloc[range(int(data['topN'])),]
  #deg.loc[:,'log2fc'] = deg.loc[:,'log2fc'].apply(lambda x: '%.2f'%x)
  #deg.loc[:,'pval'] = deg.loc[:,'pval'].apply(lambda x: '%.4E'%x)
  #deg.loc[:,'qval'] = deg.loc[:,'qval'].apply(lambda x: '%.4E'%x)
  #ppr.pprint(GSEAtable)
  #ppr.pprint(GSEAtable.sort_values('pval'))
  return json.dumps([deg.to_csv(index=False),img,GSEAtable.to_csv(index=False),GSEAimg])#json.dumps([deg.values.tolist(),img])

def specificity_score(adata=None, ctype_col:str=None, ctypes:list=None, ctype_sets:list=None, glist:list=None):
  if not (adata and ctype_col):
    raise ValueError('No adata or No ctype_col_name provided!')
  
  if not (ctypes and ctype_sets):
    ctypes = adata.obs[ctype_col].unique().tolist()
    ctype_dict = {_: set([_]) for _ in ctypes}
  else:
    if len(ctypes)!=len(ctype_sets):
      raise ValueError("Please check your cell type partition!")
    if len(ctypes)>1:
      sumset = ctype_sets[0]
      for nextset in ctype_sets[1:]:
        if sumset.intersection(nextset):
          raise ValueError("Self-defined cell type partitions have overlap!!")
        else:
          sumset = sumset.union(nextset)    
    ctype_dict = {ctypes[i]: ctype_sets[i] for i in range(len(ctypes))}
    
  if glist and len(glist):
    adata = adata[:, adata.var_names.isin(set(glist))]
    if adata.shape[-1]==0:
      raise ValueError("No gene found! Pls check your gene list!")
    
  adata_dict = {ctype:adata[adata.obs[ctype_col].isin(ctype_dict[ctype])] for ctype in ctypes}
  mean_exp_dict = {ctype: np.squeeze(np.asarray(adata_dict[ctype].X.sum(axis=0))) / adata_dict[ctype].n_obs for ctype in ctypes}
  df = pd.DataFrame(data=mean_exp_dict, index=adata_dict[ctypes[0]].var_names)
  df['all'] = df.sum(axis=1)
  for col in df.columns:
    df[col] = df[col]/df['all']
    
  # assign nan for genes not found
  # gene in glist order
  if glist and len(glist): df = df.reindex(glist)   
  return df
def restoreX(d):
  if d.max()<50:
    d1 = np.exp(d)
    d2 = np.exp2(d)
    if abs(d1-d1.round()).sum()<abs(d2-d2.round()).sum():
      d = d1
      ppr.pprint("...exp is taken")
    else:
      d = d2
      ppr.pprint("...exp2 is taken")
    if d.min()>0.1 and (d==d.min()).sum()>100:
      d = d-d.min()
      ppr.pprint("...min value %f is taken"%d.min())
  return(d)
def GSP(data):
  data['figOpt']['scale']='No'
  D = createData(data)
  if scipy.sparse.issparse(D.X):
    D.X.data = restoreX(D.X.data)
  else:
    D.X = restoreX(D.X)
  grp = data['grp'][0]
  X = specificity_score(adata=D,ctype_col=grp)
  if 'all' in X.columns:
    Xshow=X.drop('all',axis=1)
  else:
    Xshow = X
  # plot depends on the gene number
  if X.shape[0]==1:
    XT = Xshow.T
    gene = XT.columns[0]
    XT[grp] = XT.index
    fig = plt.figure()
    #ax = XT.plot.bar(x=grp,y=gene)
    ax = sns.barplot(data=XT,x=grp,y=gene)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)#rotation=45, horizontalalignment='right'
    img = iostreamFig(fig)
  elif X.shape[0]<20:
    fig = plt.figure(figsize=[0.8*Xshow.shape[1],0.5*Xshow.shape[0]+0.1*max([len(x) for x in Xshow.columns])])
    ax = sns.heatmap(Xshow,annot=True,fmt=".1f",cmap=sns.cubehelix_palette(as_cmap=True))
    img = iostreamFig(fig)
  else:
    g = sns.clustermap(
      Xshow.sample(min(1000,X.shape[0])),
      method="ward",col_cluster=False,
      yticklabels=False,xticklabels=True)
    img = iostreamFig(g)
  X.insert(0,"gene",X.index)
  return json.dumps([X.to_csv(index=False),img])

def DOT(data):
  #ppr.pprint("DOT, starting ...")
  updateGene(data)
  # Dot plot, The dotplot visualization provides a compact way of showing per group, the fraction of cells expressing a gene (dot size) and the mean expression of the gene in those cell (color scale). The use of the dotplot is only meaningful when the counts matrix contains zeros representing no gene counts. dotplot visualization does not work for scaled or corrected matrices in which zero counts had been replaced by other values, see http://scanpy-tutorials.readthedocs.io/en/multiomics/visualizing-marker-genes.html
  data['figOpt']['scale'] = 'No';
  #ppr.pprint("DOT: creating data ...")
  adata = createData(data)
  #ppr.pprint("DOT: data created!")
  if len(adata)==0:
    return Msg('No cells in the condition!')
  #return adata
  grp = adata.obs[data['grp'][0]].unique()
  if len(grp)<10:
      col = np.array(sns.color_palette('Set1',len(grp)).as_hex())
  elif len(grp)<20:
      col = np.array(sns.color_palette(n_colors=len(grp)).as_hex())
  else:
      col = np.array(sns.color_palette("husl",len(grp)).as_hex())
  adata.uns[data['grp'][0]+'_colors'] = col
  #ppr.pprint(sc.__version__)
  if 'split_show' in data['figOpt']['scanpybranch']:#.dev140+ge9cbc5f
    dp = sc.pl.dotplot(adata,data['genes'],groupby=data['grp'][0],expression_cutoff=float(data['cutoff']),mean_only_expressed=(data['mean_only_expressed'] == 'Yes'),
                       var_group_positions=data['grpLoc'],var_group_labels=data['grpID'],
                       return_fig=True)#
    dp = dp.add_totals(size=1.2).legend(show_size_legend=True,width=float(data['legendW'])).style(cmap=data['color'], dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5)
    dp.show()
    fig = dp.get_axes()['mainplot_ax'].figure
  else:
    sc.pl.dotplot(adata,data['genes'],groupby=data['grp'][0],show=False,expression_cutoff=float(data['cutoff']),mean_only_expressed=(data['mean_only_expressed'] == 'Yes'),var_group_positions=data['grpLoc'],var_group_labels=data['grpID'], color_map=data['color'])
    fig = plt.gcf()
  #ppr.pprint(adata)

  return iostreamFig(fig)
def DOTdata(data):
  updateGene(data)
  data['figOpt']['scale'] = 'No';
  adata = createData(data)
  if len(adata)==0:
    raise ValueError('No cells in the condition!')
  return pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()

def EMBED(data):
  adata = createData(data)
  if len(data['grpNum'])>0:
    adata.obs = pd.concat([adata.obs,getObsNum(data)],axis=1)
  subSize = 4
  ncol = int(data['ncol'])
  ngrp = len(data['grp'])
  ngrpNum = len(data['grpNum'])
  ngene = len(data['genes'])
  nrow = ngrp+math.ceil(ngrpNum/ncol)+math.ceil(ngene/ncol)
  if 'splitGrp' in data.keys():
    splitName = list(adata.obs[data['splitGrp']].unique())
    nsplitRow = math.ceil(len(splitName)/ncol)
    nrow = ngrp+math.ceil(ngrpNum/ncol)+ngene*nsplitRow
  step =11
  grpCol = {gID:math.ceil(len(list(adata.obs[gID].unique()))/step) for gID in data['grp']}

  rcParams['figure.constrained_layout.use'] = False
  fig = plt.figure(figsize=(ncol*subSize,subSize*nrow))
  gs = fig.add_gridspec(nrow,ncol,wspace=0.2)
  for i in range(ngrp):
    grpName = adata.obs[data['grp'][i]].value_counts().to_dict()
    grpPalette = None
    plotOrder = None
    dotSize = None
    if len(grpName)==2 and max(grpName.values())/min(grpName.values())>10:
      grpPalette = {max(grpName,key=grpName.get):'#c0c0c030',min(grpName,key=grpName.get):'#de2d26ff'}
      plotOrder = min(grpName,key=grpName.get) #list(grpPalette.keys()) #
      grpPalette = [grpPalette[k] for k in list(adata.obs[data['grp'][i]].cat.categories)]
      dotSize = adata.obs.apply(lambda x: 360000/adata.shape[1] if x['HIVcell']==plotOrder else 120000/adata.shape[1],axis=1).tolist()
    ax = sc.pl.embedding(adata,data['layout'],color=data['grp'][i],ax=fig.add_subplot(gs[i,0]),show=False,palette=grpPalette,groups=plotOrder,size=dotSize)
    if grpCol[data['grp'][i]]>1:
      ax.legend(ncol=grpCol[data['grp'][i]],loc=6,bbox_to_anchor=(1,0.5),frameon=False)
    ax.set_xlabel('%s1'%data['layout'])
    ax.set_ylabel('%s2'%data['layout'])

  for i in range(ngrpNum):
    x = int(i/ncol)+ngrp
    y = i % ncol
    ax = sc.pl.embedding(adata,data['layout'],color=data['grpNum'][i],ax=fig.add_subplot(gs[x,y]),show=False)#,wspace=0.25
    ax.set_xlabel('%s1'%data['layout'])
    ax.set_ylabel('%s2'%data['layout'])

  if 'splitGrp' in data.keys():
    vMax = adata.to_df().apply(lambda x: max(x))
    vMin = adata.to_df().apply(lambda x: min(x))
    dotSize = 120000 / adata.n_obs
    for i in range(ngene):
      for j in range(len(splitName)):
        x = ngrp + math.ceil(ngrpNum/ncol) + i*nsplitRow+int(j/ncol)
        y = j % ncol
        ax = sc.pl.embedding(adata,data['layout'],ax=fig.add_subplot(gs[x,y]),show=False)#color=data['genes'][i],wspace=0.25,
        ax = sc.pl.embedding(adata[adata.obs[data['splitGrp']]==splitName[j]],data['layout'],color=data['genes'][i],
                vmin=vMin[data['genes'][i]],vmax=vMax[data['genes'][i]],ax=ax,show=False,
                size=dotSize,title='{} in {}'.format(data['genes'][i],splitName[j]))
        ax.set_xlabel('%s1'%data['layout'])
        ax.set_ylabel('%s2'%data['layout'])
  else:
    for i in range(ngene):
      x = int(i/ncol)+ngrp+math.ceil(ngrpNum/ncol)
      y = i % ncol
      ax = sc.pl.embedding(adata,data['layout'],color=data['genes'][i],ax=fig.add_subplot(gs[x,y]),show=False)
      ax.set_xlabel('%s1'%data['layout'])
      ax.set_ylabel('%s2'%data['layout'])

  return iostreamFig(fig)

def TRACK(data):
  updateGene(data)
  adata = createData(data)
  if len(adata)==0:
    return Msg('No cells in the condition!')
  w = math.log2(adata.n_obs)
  h = adata.n_vars/2

  ## a bug in scanpy reported: https://github.com/theislab/scanpy/issues/1265, if resolved the following code is not needed
  if len(data['grpLoc'])>0 and data['grpLoc'][len(data['grpLoc'])-1][1] < (len(data['genes'])-1):
    data['grpLoc'] += [(data['grpLoc'][len(data['grpLoc'])-1][1]+1,len(data['genes'])-1)]
    data['grpID'] += ['others']
  ##############
  #ppr.pprint(data['grpLoc'])
  #ppr.pprint(data['grpID'])

  ax = sc.pl.tracksplot(adata,data['genes'],groupby=data['grp'][0],figsize=(w,h),
                        var_group_positions=data['grpLoc'],var_group_labels=data['grpID'],
                        show=False)
  fig=ax['track_axes'][0].figure
  return iostreamFig(fig)

def cut(x,cutoff,anno):
  iC = x[x>cutoff].count()
  if iC ==0:
    return "None"
  elif iC==2:
    return "Both"
  elif x[0]>cutoff:
    return anno[0]
  elif x[1]>cutoff:
    return anno[1]
  return "ERROR"
def dualExp(df,cutoff,anno):
  label = ['None']+list(anno)+['Both']
  a = df.iloc[:,0]>cutoff
  b = df.iloc[:,1]>cutoff
  return pd.Series([label[i] for i in list(a+2*b)],index=df.index,dtype='category')

def DUAL(data):
  adata = createData(data)
  adata.obs['Expressed'] = dualExp(adata.to_df(),float(data['cutoff']),adata.var_names)
  sT = time.time()
  pCol = {"None":"#AAAAAA44","Both":"#EDDF01AA",data['genes'][0]:"#1CAF82AA",data['genes'][1]:"#FA2202AA"}
  adata.uns["Expressed_colors"]=[pCol[i] for i in adata.obs['Expressed'].cat.categories]

  rcParams['figure.figsize'] = 4.5, 4
  fig = sc.pl.embedding(adata,data['layout'],color='Expressed',return_fig=True,show=False,legend_fontsize="small")
  plt.xlabel('%s1'%data['layout'])
  plt.ylabel('%s2'%data['layout'])
  rcParams['figure.figsize'] = 4, 4
  return iostreamFig(fig)

def MARK(data):
  adata = createData(data)
  if len(adata)==0:
    return Msg('No cells in the condition!')
  ## remove the annotation whose cell counts are smaller than 2 to avoid division by zero
  vCount = adata.obs[data["grp"][0]].value_counts()
  keepG = [key for key,val in vCount.items() if val>2]
  adata = adata[adata.obs[data["grp"][0]].isin(keepG),:]

  if len(adata.obs[data['grp'][0]].unique())<3:
    return 'ERROR @server: {}'.format('Less than 3 groups in selected cells! Please use DEG for 2 groups')
    #return json.dumps([[['name','scores'],['None','0']],Msg('Less than 3 groups in selected cells!Please use DEG for 2 groups')])

  sc.tl.rank_genes_groups(adata,groupby=data["grp"][0],n_genes=int(data['geneN']),method=data['markMethod'])#
  #ppr.pprint(int(data['geneN']))
  sc.pl.rank_genes_groups(adata,n_genes=int(data['geneN']),ncols=min([3,len(adata.obs[data['grp'][0]].unique())]),show=False)
  fig =plt.gcf()

  gScore = adata.uns['rank_genes_groups']
  #ppr.pprint(gScore)
  pKeys = [i for i in ['names','scores','logfoldchanges','pvals','pvals_adj'] if i in gScore.keys()]
  scoreM = [pKeys+['Group']]
  for i in gScore['scores'].dtype.names:
    for j in range(len(gScore['scores'][i])):
      one = []
      for k in pKeys:
        if k=='logfoldchanges':
          one += ['%.2f' % gScore[k][i][j]]
        elif k in ['pvals','pvals_adj']:
          one += ['%.4E' % gScore[k][i][j]]
        elif k=='scores':
          one += ['%.4f' % gScore[k][i][j]]
        else:
          one += [gScore[k][i][j]]
      scoreM += [one+[i]]
  return json.dumps([scoreM,iostreamFig(fig)])

def DENS(data):
  #sT = time.time()
  adata = createData(data)
  #ppr.pprint("read data cost: %f seconds" % (time.time()-sT))
  #sT = time.time()
  adata.obs['None'] = pd.Categorical(['all']*adata.shape[0])
  bw=float(data['bw'])
  sGrp = data['category'][0]
  cGrp = data['category'][1]

  defaultFontsize = 16
  if 'figOpt' in data.keys():
    defaultFontsize = float(data['figOpt']['fontsize'])
  subSize = 4
  #split = list(adata.obs[sGrp].unique())
  split = sorted(list(adata.obs[sGrp].cat.categories))
  genes = sorted(list(adata.var.index))
  #colGrp = list(adata.obs[cGrp].unique())
  colGrp = sorted(list(adata.obs[cGrp].cat.categories))
  legendCol = math.ceil(len(colGrp)/(len(split)*11))
  fig = plt.figure(figsize=(len(genes)*subSize,len(split)*(subSize-1)))
  plt.xlabel("Expression",labelpad=20,fontsize=defaultFontsize+1)
  #plt.ylabel(sGrp,labelpad=50,fontsize=defaultFontsize+1)
  plt.xticks([])
  plt.yticks([])
  plt.box(on=None)

  #plt.xlabel("Expression")
  #plt.ylabel(sGrp)
  gs = fig.add_gridspec(len(split),len(genes),wspace=0.2)#
  #dataT = 0
  #plotT = 0
  for i in range(len(split)):
    #resT = time.time()
    Dobs = adata[adata.obs[sGrp]==split[i]].obs[cGrp]
    D = adata[adata.obs[sGrp]==split[i]].to_df()
    #dataT += (time.time()-resT)
    for j in range(len(genes)):
      ax = fig.add_subplot(gs[i,j])
      #resT = time.time()
      for one in colGrp:
        if sum(Dobs==one)<1:
          sns.kdeplot([0],label=one)
        else:
          sns.kdeplot(D[Dobs==one][genes[j]].to_numpy(),bw_method=bw,label=one)

      ax.set_ylabel("",fontsize=defaultFontsize)
      if i==0:
        ax.set_title(genes[j],fontsize=defaultFontsize+2)
      if j==0:
        ax.set_ylabel(split[i],fontsize=defaultFontsize)
      if i==0 and j==(len(genes)-1):
        ax.legend(prop={'size': 10},title = cGrp,loc=2,bbox_to_anchor=(1,1),ncol=legendCol,frameon=False)#
      else:
        leg = ax.get_legend()
        if not leg==None:
          leg.remove()
  #fig.text(0.6,0.09,"Expression",ha='center')
  #ppr.pprint("plotting data cost: %f seconds" % dataT)
  #ppr.pprint("plotting plot cost: %f seconds" % plotT)
  #ppr.pprint("plotting total cost: %f seconds" % (time.time()-sT))
  return iostreamFig(fig)

def SANK(data):
  updateGene(data)
  if len(data['genes'])==0:
    tmp, D = getObs(data)
    D = D.apply(lambda x:x.apply(lambda y:x.name+":"+y))
  else:
    adata = createData(data)
    D = pd.concat([adata.obs.apply(lambda x:x.apply(lambda y:x.name+":"+y)),
                   adata.to_df().apply(lambda x:pd.cut(x,int(data['sankBin'])).apply(lambda y:x.name+":"+'%.1f_%.1f'%(y.left,y.right)))],
                  axis=1,sort=False)
  D = D.astype('str').astype('category')
  if data['obs_index'] in D.columns:
    del D[data['obs_index']]

  colName =['Set1','Set3','viridis']
  labels = []
  cols = []
  colindex = 0
  for gID in D.columns:
    gNames = list(D[gID].unique())
    labels += gNames
    if len(gNames) <10:
      cols += sns.color_palette(colName[colindex%2],len(gNames)).as_hex()
      colindex += 1
    else:
      cols += sns.color_palette(colName[2],len(gNames)).as_hex()

  sIDs =[]
  dIDs =[]
  v=[]
  Dnames = data['sankOrder']#list(D.columns)
  #maxGrp = 0
  #ppr.pprint(Dnames)
  for i in range(len(Dnames)-1):
    oneName = Dnames[i:i+2]
    #maxGrp = max(maxGrp,len(D[oneName[0]].unique()))
    summaryOne = D.groupby(oneName,observed=False).size().reset_index(name='Count')
    summaryOne=summaryOne[summaryOne['Count']>0]
    sIDs += list(summaryOne[oneName[0]].apply(lambda x: labels.index(x)))
    dIDs += list(summaryOne[oneName[1]].apply(lambda x: labels.index(x)))
    v += list(summaryOne['Count'])

  data_trace = dict(
    type='sankey',
    domain=dict(x=[0,1],y=[0,1]),
    orientation='h',
    valueformat = ".0f",
    node = dict(
      pad = 10,
      thickness = 15,
      line = dict(
        color = "black",
        width = 0.5
      ),
      label =  labels,
      color =  cols
    ),
    link = dict(
      source = sIDs,
      target = dIDs,
      value = v
    )
  )
  ## if the image is requested
  if 'imgSave' in data.keys():
    layout = dict(
      font = dict(size=int(data['figOpt']['fontsize'])),
      height= int(data['imgH']),
      width = int(data['imgW'])*D.shape[1]
    )
    fig = go.Figure(data=[go.Sankey(data_trace)],layout=layout)
    img = plotIO.to_image(fig,data['imgSave'])
    return base64.encodebytes(img).decode('utf-8')

  layout = dict(
    font = dict(size=int(data['figOpt']['fontsize'])),
    height= int(data['imgH']),
    width = int(data['imgW'])*D.shape[1],
    updatemenus= [
            dict(
                y=0.9,
                buttons=[
                    dict(
                        label='Thick',
                        method='restyle',
                        args=['node.thickness', 15]
                    ),
                    dict(
                        label='Thin',
                        method='restyle',
                        args=['node.thickness', 8]
                    )
                ]
            ),
            dict(
                y=0.8,
                buttons=[
                    dict(
                        label='Small gap',
                        method='restyle',
                        args=['node.pad', 15]
                    ),
                    dict(
                        label='Large gap',
                        method='restyle',
                        args=['node.pad', 20]
                    )
                ]
            ),
            dict(
                y=0.7,
                buttons=[
                    dict(
                        label='Snap',
                        method='restyle',
                        args=['arrangement', 'snap']
                    ),
                    dict(
                        label='Perpendicular',
                        method='restyle',
                        args=['arrangement', 'perpendicular']
                    ),
                    dict(
                        label='Freeform',
                        method='restyle',
                        args=['arrangement', 'freeform']
                    ),
                    dict(
                        label='Fixed',
                        method='restyle',
                        args=['arrangement', 'fixed']
                    )
                ]
            ),
            dict(
                y=0.6,
                buttons=[
                    dict(
                        label='Horizontal',
                        method='restyle',
                        args=['orientation','h']#{,'height':700,'width':250*D.shape[1]}
                    ),
                    dict(
                        label='Vertical',
                        method='restyle',
                        args=['orientation','v']#{'orientation': 'v','height':250*D.shape[1],'width':700}
                    )
                ]

            )
        ]
  )
  fig = go.Figure(data=[go.Sankey(data_trace)],layout=layout)
  div = plotIO.to_html(fig)
  return div#[div.find('<div>'):(div.find('</div>')+6)]

def DENS2D(data):
  adata = createData(data)

  ## plot in R
  strF = ('%s/DENS2D%f.csv' % (data["CLItmp"],time.time()))
  adata.to_df().to_csv(strF)#
  res = subprocess.run([strExePath+'/Density2D.R',strF,data['figOpt']['img'],str(data['cutoff']),str(data['bandwidth']),data['figOpt']['colorMap'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']],capture_output=True)#
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in R: "+res.stderr.decode('utf-8'))

  return img

def toInt(x):
  if x.shape[0]==0:
    return 0
  return int(x.iloc[0])

def STACBAR(data):
  if len(data['genes'])==0:
    tmp, D = getObs(data)
    D = D.apply(lambda x:x.apply(lambda y:y))
  else:
    adata = createData(data)

    D = pd.concat([adata.obs.apply(lambda x:x.apply(lambda y:y)),
                   adata.to_df().apply(lambda x:pd.cut(x,int(data['Nbin'])).apply(lambda y:'%s:%.1f_%.1f'%(x.name,y.left,y.right)))],
                  axis=1,sort=False)
  D = D.astype('str').astype('category')
  if data['obs_index'] in D.columns:
    del D[data['obs_index']]
  cellN = D.groupby(list(D.columns),observed=False).size().reset_index(name="Count")

  strCol = data['colorBy']
  tmp = list(D.columns)
  tmp.remove(strCol)
  strX = tmp[0]
  returnD = [{'name':i,
              'sales':[{'year':j,#.replace(strX+':',''),
                        'profit':toInt(cellN[(cellN[strCol]==i) & (cellN[strX]==j)]['Count'])}
                        for j in cellN[strX].unique()]}
              for i in cellN[strCol].unique()]
  return json.dumps(returnD)

def CLI(data):
  strPath = data["CLItmp"]+('/CLI%f' % time.time())
  script = data['script']
  del data['script']

  adata = createData(data)

  strData = strPath + '.h5ad'
  adata.write(strData)
  #with open(strData,'wb') as f:
  #pickle.dump(adata,f)
  ppr.pprint(len(re.findall(r'```',script)))
  if (len(re.findall(r'```',script)) >0):
    strScript = strPath + '.Rmd'
    with open(strScript,'w') as f:
     f.writelines(['---\noutput:\n  html_document:\n    code_folding: hide\n---\n\n```{r}\nstrPath <- "%s"\n```\n\n'%strPath])
     f.write(script)
    ppr.pprint(subprocess.run('which Rscript',capture_output=True,shell=True).stdout.decode('utf-8'))
    res = subprocess.run('Rscript -e \'rmarkdown::render("%s", output_file="%s.html")\''%(strScript,strPath),capture_output=True,shell=True)
    if (os.path.exists('%s.html'%strPath)):
      with open('%s.html'%strPath,'r') as file:
        html = file.read()
    else:
      html = ''
    ppr.pprint(res.stdout.decode('utf-8'))
    ppr.pprint(res.stderr.decode('utf-8'))
  else:
    strScript = strPath + '.py'
    with open(strScript,'w') as f:
      f.writelines(['%load_ext rpy2.ipython\n','from anndata import read_h5ad\n','adata=read_h5ad("%s")\n'%strData, 'strPath="%s"\n\n'%strPath])
      #f.writelines(['%load_ext rpy2.ipython\n','import pickle\n','with open("%s","rb") as f:\n'%strData,'  adata=pickle.load(f)\n','strPath="%s"\n\n'%strPath])
      f.writelines(['%%R\n','strPath="%s"\n\n'%strPath])
      f.write(script)
    ppr.pprint(subprocess.run('which Rscript',capture_output=True,shell=True).stdout.decode('utf-8'))
    ppr.pprint(subprocess.run('which pandoc',capture_output=True,shell=True).stdout.decode('utf-8'))
    # (time consuming) ppr.pprint(subprocess.run("Rscript -e 'reticulate::py_config()'",capture_output=True,shell=True).stdout.decode('utf-8'))
    try:
      res = subprocess.run('jupytext --to notebook --output - %s | jupyter nbconvert --ExecutePreprocessor.timeout=1800 --to html --template classic --execute --stdin --stdout'%strScript,capture_output=True,shell=True)
    except:
      res = subprocess.run('jupytext --to notebook --output - %s | jupyter nbconvert --ExecutePreprocessor.timeout=1800 --to html --execute --stdin --stdout'%strScript,capture_output=True,shell=True)
    html = res.stdout.decode('utf-8')
  #    h,s,e = html.partition('<div class="cell border-box-sizing code_cell rendered">')
  #    h1,s,e = e.partition('<div class="cell border-box-sizing code_cell rendered">') ## remove the first cell
  #    h1,s,e = e.partition('<div class="cell border-box-sizing code_cell rendered">') ## remove the second cell
  #    html = h+s+e
  if 'Error' in res.stderr.decode('utf-8'):
     html = 'ERROR @server:\nstderr:\n' + re.sub(r"\x1b[^m]*m", "", res.stderr.decode('utf-8')) + '\nstdout:\n' + res.stdout.decode('utf-8')
  for f in glob.glob(strPath+"*"):
    try:
      os.remove(f)
    except:
      continue

  return html

def getDesp(data):
  strF = re.sub("h5ad$","txt",data["h5ad"])
  if not os.path.isfile(strF):
    return ""
  txt = ""
  with open(strF,'r') as fp:
    for line in fp:
      txt = "%s<br>%s"%(txt,line)
  return txt

def setTimeStamp(data):
  if data["h5ad"].startswith("http"):
      return
  strF = re.sub("h5ad$","timestamp",data["h5ad"])
  Path(strF).touch()

def getPreDEGname(data):
  strF = re.sub("h5ad$","db",data["h5ad"])
  if not os.path.isfile(strF):
    #ppr.pprint(strF+" is NOT found!")
    return ""
  conn = sqlite3.connect(strF)
  df = pd.read_sql_query("select DISTINCT contrast,tags from DEG;", conn)
  conn.close()

  return json.dumps(list(df['contrast']+"::"+df['tags']))

def getPreDEGvolcano(data):
  strF = re.sub("h5ad$","db",data["h5ad"])
  comGrp = data["compSel"].split("::")

  conn = sqlite3.connect(strF)
  df = pd.read_sql_query("select gene,log2fc,pval,qval from DEG where contrast=? and tags=?;", conn,params=comGrp)
  conn.close()
  deg = df.sort_values(by=['qval'])
  data["comGrp"] = comGrp[0].split(".vs.")

  ## plot in R
  strF = ('%s/DEG%f.csv' % (data["CLItmp"],time.time()))
  deg.to_csv(strF,index=False)
  #ppr.pprint([strExePath+'/volcano.R',strF,';'.join(data['genes']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),str(data['logFC']),data['comGrp'][1],data['comGrp'][0],str(data['sigFDR']),str(data['sigFC']),str(data['labelSize']),str(data['dotSize']),str(data['ymin']),str(data['ymax']),data['Rlib']])
  res = subprocess.run([strExePath+'/volcano.R',strF,';'.join(data['genes']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),str(data['logFC']),data['comGrp'][1],data['comGrp'][0],str(data['sigFDR']),str(data['sigFC']),str(data['labelSize']),str(data['dotSize']),str(data['ymin']),str(data['ymax']),data['figOpt']['vectorFriendly'],data['Rlib']],capture_output=True)#
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in R: "+res.stderr.decode('utf-8'))
  #####
  gInfo = getVar(data)
  deg.index = deg['gene']
  deg = pd.concat([deg,gInfo],axis=1,join='inner',sort=False)
  #return deg.to_csv()

  if not data['topN']=='All':
    deg = deg.iloc[range(min(deg.shape[0],int(data['topN']))),]
  #deg.loc[:,'log2fc'] = deg.loc[:,'log2fc'].apply(lambda x: '%.2f'%x)
  #deg.loc[:,'pval'] = deg.loc[:,'pval'].apply(lambda x: '%.4E'%x)
  #deg.loc[:,'qval'] = deg.loc[:,'qval'].apply(lambda x: '%.4E'%x)

  return json.dumps([deg.to_csv(index=False),img])#json.dumps([deg.values.tolist(),img])

def getPreDEGbubble(data):
  #data={'compSel':['MS.vs.Control::EN.L4','MS.vs.Control::Endo.cells','MS.vs.Control::EN.PYR'],'genes':['RASGEF1B','SLC26A3','UNC5C','AHI1','CD9']}
  sql = "select gene,log2fc,pval,qval,contrast || '::' || tags as tag from DEG where tag in ({comp}) and gene in ({gList}) order by case tag {oList} end;".format(
    comp=','.join(['?']*len(data['compSel'])),
    gList=','.join(['?']*len(data['genes'])),
    oList=' '.join(['WHEN ? THEN %d'%i for i in range(len(data['compSel']))]))

  strF = re.sub("h5ad$","db",data["h5ad"])
  conn = sqlite3.connect(strF)
  deg = pd.read_sql_query(sql,conn,params=data['compSel']+data['genes']+data['compSel'])
  conn.close()
  if deg.shape[0]==0:
    raise ValueError("No data for selected genes ("+", ".join(data['genes'])+") in selected comparison ("+", ".join(data['compSel'])+")!")

  ## add selected genes which is not in the database back to the dataframe as NA
  addG = [[i,np.nan,np.nan,np.nan,data['compSel'][0]] for i in data['genes'] if i not in list(deg.gene.unique())]
  if len(addG)>0:
    deg = pd.concat([deg,pd.DataFrame(addG,columns=deg.columns)])
  ## add selected comparison which is not in the database back to the dataframe as NA
  addComp = [[data['genes'][0],np.nan,np.nan,np.nan,i] for i in data['compSel'] if i not in list(deg.tag.unique())]
  if len(addComp)>0:
    deg = pd.concat([deg,pd.DataFrame(addComp,columns=deg.columns)])
  #ppr.pprint(deg)
  ## plot in R
  strF = ('%s/DEG%f.csv' % (data["CLItmp"],time.time()))
  deg.to_csv(strF,index=False)
  #ppr.pprint(' '.join([strExePath+'/bubbleMap.R',strF,data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['scale'],data['Rlib']]))
  res = subprocess.run([strExePath+'/bubbleMap.R',strF,data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['scale'],data['Rlib']],capture_output=True)#
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in R: "+res.stderr.decode('utf-8'))

  #RASGEF1B SLC26A3 UNC5C AHI1 CD9
  return img

def getPreDEGfgsea(data):
  strF = re.sub("h5ad$","db",data["h5ad"])
  comGrp = data["compSel"].split("::")
  conn = sqlite3.connect(strF)
  df = pd.read_sql_query("select gene,log2fc,pval,qval from DEG where contrast=? and tags=?;", conn,params=comGrp)
  conn.close()
  deg = df.sort_values(by=['qval'])
  data["comGrp"] = comGrp[0].split(".vs.")
  
  strF = ('%s/DEG%f.csv' % (data["CLItmp"],time.time()))
  deg.to_csv(strF,index=False)
  GSEAimg=""
  GSEAtable=pd.DataFrame()
  res = subprocess.run([strExePath+'/fgsea.R',
                        strF,
                        '%s/gsea/%s.symbols.gmt'%(strExePath,data['gs']),
                        str(data['gsMin']),
                        str(data['gsMax']),
                        str(data['padj']),
                        data['up'],
                        data['dn'],
                        str(data['collapse']),
                        data['figOpt']['img'],
                        str(data['figOpt']['fontsize']),
                        str(data['figOpt']['dpi']),
                        data['Rlib']],capture_output=True)#
  if 'Error' in res.stderr.decode('utf-8'):
    raise SyntaxError("in fgsea.R: "+res.stderr.decode('utf-8'))
  GSEAimg = res.stdout.decode('utf-8')
  GSEAtable = pd.read_csv(strF)
  GSEAtable['leadingEdge'] = GSEAtable['leadingEdge'].apply(lambda x:'|'.join(x.split('|')[:10]))
  os.remove(strF)
  return json.dumps([GSEAtable.to_csv(index=False),GSEAimg])

def getEnv():
  config = {'CLItmp':'/tmp','Rpath':'','Rlib':'','METAtmp':'/tmp','METAurl':'','METAmax':1e4}
  strEnv = '%s/vip.env'%strExePath
  if os.path.isfile(strEnv):
    with open(strEnv,'r') as fp:
      for line in fp:
        one = line.strip().replace("\t", "").replace(" ", "").split("=")
        if not len(one)==2:
          continue
        config[one[0]]=one[1]
  #ppr.pprint(config)
  if len(config['Rpath'])>3:
    os.stat("%s/Rscript"%config['Rpath'])
    os.environ['PATH'] = config['Rpath']+os.pathsep+os.environ['PATH']
  return config
try:
  VIPenv = getEnv()
except Exception as e:
  ppr.pprint("The specified R path is incorrect, please check or remove from vip.env!")
  raise e

def mergeMeta(data):
  selC = list(data['cells'].values())
  ## obtain the category annotation
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    if not 'cellN' in scD.data.obs:
      raise ValueError('This is not a metacell data!')
    obs = scD.data.obs.loc[selC,[data['obs_index'],'cellN']]
  ppr.pprint(obs)
  ppr.pprint(obs['cellN'].sum())
  if obs['cellN'].sum() > int(data['METAmax']):
    raise ValueError('The selected meta cells include more than maximum %d cells!'% int(data['METAmax']))
  strPath = re.sub(".h5ad$","",data["h5ad"])
  selCells = []
  for i in obs[data['obs_index']]:
    strOne = strPath+"/"+i+".h5ad"
    if os.path.exists(strOne):
      selCells += [ad.read(strOne)]
  strOut = data['METAtmp']+"/"+os.path.basename(strPath)+"_"+data['metaPostfix']+".h5ad"
  ad.concat(selCells).write(strOut)
  return data['METAurl']+"/d/"+os.path.basename(strOut)+"/"

def isMeta(data):
  #with app.get_data_adaptor(url_dataroot=data['url_dataroot'],dataset=data['dataset']) as scD:
  scD=data['data_adapter']
  if scD is not None:
    if not 'cellN' in scD.data.obs:
      return "FALSE"
  strPath = re.sub(".h5ad$","",data["h5ad"])
  if not os.path.exists(strPath):
    return "FALSE"
  return "TRUE"

def getBWinfo(data):
    BWinfo = {"BWfile":[],"BWannotation":[],"BWlink":[],"BWpeak":[],"BWcluster":[]}
    strD = re.sub(".h5ad$","/",data["h5ad"])
    if os.path.isdir(strD):
        for one in os.listdir(strD):
            if one.endswith("bw"):#not re.search("bw$",one)==None:
                BWinfo["BWfile"].append(one)
            elif one=="annotation.rds":
                BWinfo["BWannotation"]="annotation.rds"
            elif one=="peaks.rds":
                BWinfo["BWpeak"]="peaks.rds"
            elif one=="links.rds":
                BWinfo["BWlink"]="links.rds"
            elif one=="bw.cluster":
                BWinfo["BWcluster"]=pd.read_csv(strD+one,sep="\t",header=0) #"bw.cluster" .to_csv(index=False)
        if len(BWinfo["BWcluster"])>0:
            BWinfo["BWcluster"]=BWinfo["BWcluster"][BWinfo["BWcluster"]['Wig'].isin(BWinfo["BWfile"])].to_csv(index=False)
    return json.dumps(BWinfo)

def plotBW(data):
    strD = re.sub(".h5ad$","/",data["h5ad"])
    strCSV = ('%s/BW%f.csv' % (data["CLItmp"],time.time()))
    ## select all cells
    strType = strD + 'bw.cluster'
    if os.path.isfile(strType) and len(data['genes'])>0:
        clusterInfo = pd.read_csv(strType,sep="\t",header=0,index_col=0)
        clusterInfo = clusterInfo.loc[data['bw'],:]
        grp = []
        scD=data['data_adapter']
        if scD is not None:
          dAnno = list(scD.get_obs_keys())
          grp = [i for i in clusterInfo.columns if i in dAnno]
        if len(grp)>0:
            data['grp'] = grp
            adata = createData(data)
            if len(adata)>0:
                selC=adata.obs[grp[0]].isin(list(clusterInfo[grp[0]]))
                for oneG in grp:
                    selC = selC | adata.obs[oneG].isin(list(clusterInfo[oneG]))
                adata = adata[selC,:]
                pd.concat([adata.obs,adata.to_df()], axis=1, sort=False).to_csv(strCSV)
    ## plot in R
    #strCMD = ' '.join([strExePath+'/browserPlot.R',strD,data['region'],','.join(data['bw']),str(data['exUP']),str(data['exDN']),strCSV,str(data['cutoff']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']])
    #ppr.pprint(strCMD)
    res = subprocess.run([strExePath+'/browserPlot.R',strD,data['region'],','.join(data['bw']),str(data['exUP']),str(data['exDN']),strCSV,str(data['cutoff']),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),data['Rlib']],capture_output=True)#
    img = res.stdout.decode('utf-8')
    if len(adata)>0:
        os.remove(strCSV)
    if 'Error' in res.stderr.decode('utf-8'):
        raise SyntaxError("in R: "+res.stderr.decode('utf-8'))

    return img

def getCosMx(data):
    scD=data['data_adapter']
    cosMxKey='cosMx'
    cosMxID=[]
    if scD is not None and cosMxKey in scD.data.uns_keys():
        cosMxID = list(scD.data.uns[cosMxKey].keys())
        cosMxID.remove('keys')
    return json.dumps(cosMxID)
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
    return((round(y),round(x)))#imgH-y
def square_integer(ct,halfS):
    x_center, y_center = ct
    coordinates_in_square = []
    for x in range(x_center - halfS, x_center + halfS + 1):
        for y in range(y_center - halfS, y_center + halfS + 1):
            coordinates_in_square.append((x, y))
    return coordinates_in_square
def npArray2jpg(img):
    buff = BytesIO()
    plt.imsave(buff,img,format='jpg')
    img_base64=base64.b64encode(buff.getvalue()).decode("utf-8")
    return img_base64
def mergeGlobal(imgs,cood):
    #https://stackoverflow.com/questions/20546182/how-to-perform-coordinates-affine-transformation-using-python-part-2?answertab=votes#tab-top
    pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
    unpad = lambda x: x[:,:-1]
    imgMin={}
    imgMax=[]
    for sID in imgs:
        X = pad(cood[sID]['local'])
        Y = pad(cood[sID]['global'])
        A, res, rank, s = np.linalg.lstsq(X, Y,rcond=None)
        transform = lambda x: unpad(np.dot(pad(x), A))
        imgMin[sID]=transform(np.array([[0,0]]))
        imgMax.append(transform(np.array([[imgs[sID].shape[0],imgs[sID].shape[1]]])))
    Xmin,Ymin = np.amin(np.concatenate(list(imgMin.values()),axis=0),axis=0)
    Xmax,Ymax = np.amax(np.concatenate(imgMax,axis=0),axis=0)
    if Xmax-Xmin>65000 or Ymax-Ymin>65000:
      raise ValueError('Selected area is too large, try to select adjacent FOVs for global!')
    # flip the coordinates
    img = np.zeros((math.ceil(Ymax-Ymin),math.ceil(Xmax-Xmin),3),dtype='uint8')
    for sID in imgs:
        ix,iy = loc_img(imgMin[sID][0]-(Xmin,Ymin),img.shape[0]-imgs[sID].shape[0])
        img[ix:(ix+imgs[sID].shape[0]),iy:(iy+imgs[sID].shape[1])] = imgs[sID]
    return img
def plotCosMx(data):
    scD=data['data_adapter']
    cosMxKey='cosMx'
    keys = scD.data.uns[cosMxKey]['keys']
    cosMxArray={}
    cood_pair={}
    cood_N=10
    for sID in data['sIDs']:
        if data['histology']:
            imgC = scD.data.uns[cosMxKey][sID][keys['img']].copy()
        else:
            imgC = np.ones((scD.data.uns[cosMxKey][sID][keys['img']].shape[0],scD.data.uns[cosMxKey][sID][keys['img']].shape[1],3),dtype='uint8')*255
        if data['cell']:
            cell_color=data['cell'].lstrip('#')
            rgb = np.array(tuple(int(cell_color[i:i+2], 16) for i in (0, 2, 4)),dtype='uint8')
            cellB = scD.data.uns[cosMxKey][sID][keys['cell']]
            for i in cellB.cell_ID.unique():
                oneC = cellB[cellB.cell_ID==i].reset_index(drop=True)
                N = oneC.shape[0]
                for j in range(N):
                    p = bresenham_line(loc_img((oneC[keys['local_px']['x']][j%N],oneC[keys['local_px']['y']][j%N]),imgC.shape[1]),
                            loc_img((oneC[keys['local_px']['x']][(j+1)%N],oneC[keys['local_px']['y']][(j+1)%N]),imgC.shape[1]))
                    for a in p:
                        imgC[a]=rgb
        if len(data['genes'])>0:
            txLoc = scD.data.uns['cosMx'][sID][keys['tx_loc']]
            for one in data['genes']:
                col = data['genes'][one].lstrip('#')
                rgb = np.array(tuple(int(col[i:i+2], 16) for i in (0, 2, 4)),dtype='uint8')
                for i in txLoc.index[txLoc[keys['tx_col']]==one].tolist():
                    for a in square_integer(loc_img((txLoc[keys['local_px']['x']][i],txLoc[keys['local_px']['y']][i]),imgC.shape[1]),int(data['geneSize'])-1):
                        imgC[a] = rgb
        cosMxArray[sID]=imgC
        txLoc_sub=scD.data.uns['cosMx'][sID][keys['tx_loc']].sample(n=cood_N)
        cood_pair[sID]={'global':txLoc_sub[keys['global_px'].values()].to_numpy(),'local':txLoc_sub[keys['local_px'].values()].to_numpy()}
    cosMxImg = {}
    if data['cood']=='global' and len(cosMxArray)>1:
        cosMxImg['Global'] = npArray2jpg(mergeGlobal(cosMxArray,cood_pair))
    else:
        for sID in cosMxArray:
            cosMxImg[sID]=npArray2jpg(cosMxArray[sID])
    return json.dumps(cosMxImg)

def getVisium(data):
    scD=data['data_adapter']
    k='visium'
    visiumID=[]
    if scD is not None and k in scD.data.uns_keys():
        adata = scD.data
        keys = adata.uns[k]['keys']
        visiumID = list(set(adata.obs[keys['slide_column']].unique()) & set(adata.uns['spatial'].keys()))
    return json.dumps(visiumID)
def adjustAlpha(x,adjA):
    if max(x)<=0:
        return adjA
    x = x/max(x) #(alpha-min(alpha))/(max(alpha)-min(alpha))
    a = 1-adjA
    b = 20**a
    alpha =(1+a**b)/(1+(a/x)**b)
    return alpha
def plotVisiumOne(adata,sid,col,ax,fig,alpha=1,cmap='viridis',dotsize=4):
    keys = adata.uns['visium']['keys']
    img = adata.uns['spatial'][sid]
    for v in keys['img']:
        img = img[v]
    scaler = adata.uns['spatial'][sid]
    for v in keys['scale']:
        scaler = scaler[v]
    subD = adata[adata.obs[keys['slide_column']]==sid,]
    a=ax.imshow(img)
    x = subD.obsm[keys['coordinates']][:,0]*scaler
    y = subD.obsm[keys['coordinates']][:,1]*scaler
    g_column = 'name_0'
    if 'feature_name' in subD.var.columns:
        g_column = 'feature_name'
    df = sc.get.obs_df(subD,[col],gene_symbols=g_column)
    if pd.api.types.is_numeric_dtype(df[col]):
        a=ax.scatter(x,y,
            c=df[col].to_numpy(),cmap=cmap,s=dotsize,
            alpha=adjustAlpha(df[col].to_numpy(),alpha))
        a=fig.colorbar(a,ax=ax)
    else:
        try:
            grps = sorted(adata.obs[col].unique().tolist(),key=int)
        except:
            grps = sorted(adata.obs[col].unique().tolist())
        if len(grps)<10:
            colors = dict(zip(grps,sns.color_palette('Set1',n_colors=len(grps)).as_hex()))
        else:
            colors = dict(zip(grps,sns.color_palette('husl',n_colors=len(grps)).as_hex()))
        a=ax.scatter(x,y,
            color=[colors[_] for _ in df[col]],
            s=dotsize,alpha=alpha)
        legend_handles = [Line2D([],[],marker='.',color=colors[i],linestyle='None') for i in colors.keys()]
        a=ax.legend(handles=legend_handles,labels=colors.keys(),ncols=math.ceil(len(grps)/11),
            bbox_to_anchor=(1.05, 1),loc='upper left')
    a=ax.set_title("%s: %s"%(sid,col))
def plotVisium(data):
    subSize = float(data['subsize'])
    adata = data['data_adapter'].data
    sel = data['grpNum'] + data['genes']+data['grp']
    keys = adata.uns['visium']['keys']
    ncol = int(data['ncol'])
    nSel = len(sel)
    slides = data['sIDs'] #adata.obs[keys['slide_column']].unique()
    nSlides = len(slides) #adata.obs[keys['slide_column']].nunique()
    if data['sortID']:
        nrowSection = math.ceil(nSel/ncol)
        nrow = nSlides*nrowSection
        fig, axs = plt.subplots(nrow, ncol, figsize=(ncol*subSize,nrow*subSize))
        axs_flat = axs.flatten()
        for i in range(nSlides):
            for j in range(nSel):
                plotVisiumOne(adata,slides[i],sel[j],
                    axs_flat[i*nrowSection*ncol+j],fig,
                    alpha=float(data['alpha']),
                    cmap=data['figOpt']['colorMap'],
                    dotsize=float(data['dotsize']))
    else:
        nrowSection = math.ceil(nSlides/ncol)
        nrow = nSel*nrowSection
        fig, axs = plt.subplots(nrow, ncol, figsize=(ncol*subSize,nrow*subSize))
        axs_flat = axs.flatten()
        for j in range(nSel):
            for i in range(nSlides):
                plotVisiumOne(adata,slides[i],sel[j],
                    axs_flat[j*nrowSection*ncol+i],fig,
                    alpha=float(data['alpha']),
                    cmap=data['figOpt']['colorMap'],
                    dotsize=float(data['dotsize']))
    for ax in axs_flat:
        a=ax.axis('off')
    return iostreamFig(fig)

#make sure the h5ad file full name is listed in vip.env as a variable 'testVIP';
def testVIPready(data):
  strH5ad = os.path.basename(data["h5ad"])
  if 'testVIP' in data and strH5ad==data["testVIP"]:
    both = True
    for one in [re.sub("h5ad$","info.txt",strH5ad),re.sub("h5ad$","img.txt",strH5ad)]:
      both = both and os.path.exists(strExePath+"/../common/web/static/testVIP/"+one)
    if both:
      return "SHOW"
    else:
      return "TRUE"
  return "FALSE"

def saveTest(data):
    strPath = strExePath+"/../common/web/static/testVIP/"
    if not os.path.exists(strPath):
        os.makedirs(strPath)
    strH5ad = os.path.basename(data["h5ad"])

    if len(data['info'])>100:
        #ppr.pprint(strPath+re.sub("h5ad$","info.txt",strH5ad))
        with open(strPath+re.sub("h5ad$","info.txt",strH5ad),'w') as f:
            f.write(data['info'])
    if len(data['img'])>100:
        with open(strPath+re.sub("h5ad$","img.txt",strH5ad),'w') as f:
            f.write(data['img'])
    return 'success'
