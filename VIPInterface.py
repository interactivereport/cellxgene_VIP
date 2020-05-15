import requests
import json
import server.app.decode_fbs as decode_fbs
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rcParams
import base64
import math
from io import BytesIO
import sys

import pprint
ppr = pprint.PrettyPrinter(depth=6)

import server.app.app as app
import pickle

sys.setrecursionlimit(10000)
sc.set_figure_params(dpi=150, color_map='viridis')
sc.settings.verbosity = 2
rcParams.update({'figure.autolayout': True})

api_version = "/api/v0.2"

def route(data,appConfig=None):
  if appConfig is None:
    data["url"] = f'http://127.0.0.1:8888/{api_version}'
  else:
    data = json.loads(str(data,encoding='utf-8'))
    data["url"] = f'http://{appConfig.server__host}:{appConfig.server__port}/{api_version}'

  return distributeTask(data["method"])(data)
  
def subData(data):
  selC = list(data['cells'].values())
  cName = ["cell%d" %i for i in selC]
  
  ## onbtain the expression matrix
  gNames = []
  X = []
  if 'genes' in data.keys():
    if len(data['genes'])>0:
      with app.get_data_adaptor() as scD:
        fullG = list(scD.data.var['name_0'])
        selG = [fullG.index(i) for i in data['genes']]
        X = scD.data.X[selC][:,selG]
        gNames = data['genes']
    else:
      with app.get_data_adaptor() as scD:
        X = scD.data.X[selC]
        gNames = list(scD.data.var["name_0"])
  expr = pd.DataFrame(X,columns=gNames,index=cName)

  ## obtain the embedding
  strEmbed = 'umap'
  if 'layout' in data.keys():## tsne or umap
    strEmbed = data['layout']
  with app.get_data_adaptor() as scD:
    embed = pd.DataFrame(scD.data.obsm['X_%s'%strEmbed][selC],columns=['%s1'%strEmbed,'%s2'%strEmbed],index=cName)

  ## obtain the category annotation
  with app.get_data_adaptor() as scD:
    obs = scD.data.obs.loc[selC,['name_0']+data['grp']].astype('str')
  obs.index = cName

  ## update the annotation Abbreviation
  combUpdate = cleanAbbr(data)
  if 'abb' in data.keys():
    for i in data['grp']:
      obs[i] = obs[i].map(data['abb'][i])

  ## create a custom annotation category and remove cells which are not in the selected annotation
  if combUpdate and len(data['grp'])>1:
    newGrp = 'Custom_combine'
    obs[newGrp] = obs[data['grp'][0]]
    for i in data['grp']:
      if i!=data['grp'][0]:
        obs[newGrp] += "_"+obs[i]
    expr = expr[~obs[newGrp].str.contains("Other")]
    embed = embed[~obs[newGrp].str.contains("Other")]
    obs = obs[~obs[newGrp].str.contains("Other")]
    data['grp'] = [newGrp]
    
  obs = obs.astype('category')
  ## empty selection
  if expr.shape[0]==0 or expr.shape[1]==0:
    return []
  return sc.AnnData(expr,obs,obsm={'X_%s'%strEmbed:embed.to_numpy()})

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
          data['abb'][cate] = {key:"Other" for key in data['abb'][cate].keys()}
      
  return updated

def createData(data,seperate=False):
  #print("CreateData")
  return subData(data)
  headers = {'content-type':'application/json'}
  if not 'genes' in data:
    data['genes'] = []
  fil = json.dumps({'filter':{'var':{'annotation_value':[{'name':'name_0','values':data['genes']}]}}})
  res = requests.put('%s/data/var' % data["url"],fil,headers=headers)
  expr = decode_fbs.decode_matrix_FBS(res.content)
  res = requests.get('%s/annotations/var' % data["url"],params={'annotation-name':'name_0'})
  gNames = decode_fbs.decode_matrix_FBS(res.content)['columns'][0]

  cNames = ["cell%d" % x for x in data['cells'].values()]
  expr = pd.DataFrame([[expr['columns'][i][x] for x in data['cells'].values()] for i in range(len(expr['columns']))],
                    index=[gNames[x] for x in expr['col_idx']],columns=cNames).T
  expr.dropna()
  cNames = expr.index
  obsL = [cNames]
  combUpdate = cleanAbbr(data)
  for one in data['grp']:
    res = requests.get('%s/annotations/obs' % data["url"],params={'annotation-name':one})
    grp = decode_fbs.decode_matrix_FBS(res.content)
    if 'abb' in data.keys():
      subGrp = [data['abb'][one][str(grp['columns'][0][i])] for i in data['cells'].values()]
    else:
      subGrp = [str(grp['columns'][0][i]) for i in data['cells'].values()]
    obsL += [subGrp]
  obs = pd.DataFrame(obsL,index=['name_0']+data['grp'],columns=cNames).T.astype('category')
  strGN = [i for i in data['grp'] if 'genes' in i]
  
  if len(strGN)>0:
    obs[strGN] = obs[strGN].apply(pd.to_numeric)
  if combUpdate and len(data['grp'])>1:
    newGrp = 'Custom_combine'
    obs[newGrp] = obs[data['grp'][0]]
    for i in data['grp']:
      if i!=data['grp'][0]:
        obs[newGrp] += "_"+obs[i]
    expr = expr[~obs[newGrp].str.contains("Other")]
    obs = obs[~obs[newGrp].str.contains("Other")]
    data['grp'] = [newGrp]
  if seperate:
    return {'expr':expr,'obs':obs}
  if 'layout' in data.keys():## tsne or umap
    res = requests.get('%s/layout/obs' % data["url"],params={'layout-name':data['layout']})
    embed= decode_fbs.decode_matrix_FBS(res.content)
    embed = pd.DataFrame([[embed['columns'][i][x] for x in data['cells'].values()] for i in range(len(embed['columns']))],
                    index=embed['col_idx'],columns=cNames).T.to_numpy()
    adata = sc.AnnData(expr,obs,obsm={'X_%s'%data['layout']:embed})
    return adata

  adata = sc.AnnData(expr,obs)
  return adata
  
def errorTask(data):
  return "Error task!"
  
def distributeTask(aTask):
  return {
    'SGV':SGV,
    'PGV':PGV,
    'HEAT':pHeatmap,
    'GD':GD,
    'DEG':DEG,
    'DOT':DOT,
    'EMBED':EMBED,
    'TRAK':TRACK,
    'DUAL':DUAL,
    'MARK': MARK,
    'MINX':MINX
  }.get(aTask,errorTask)

def iostreamFig(fig):
  figD = BytesIO()
  fig.savefig(figD,format='png',bbox_inches="tight")
  imgD = base64.encodestring(figD.getvalue()).decode("utf-8")
  figD.close()
  plt.close('all')
  return imgD
  
def Msg(msg):
  fig = plt.figure(figsize=(5,2))
  plt.text(0,0.5,msg)
  ax = plt.gca()
  ax.axis('off')
  return iostreamFig(fig)

def MINX(data):
  with app.get_data_adaptor() as scD:
    minV = min(scD.data.X[0])
  return '%.1f'%minV

def SGV(data):
  # figure width and heights depends on number of unique categories
  # characters of category names, gene number
  adata = createData(data)
  if len(adata)==0:
    return Msg('No cells in the condition!')

  a = list(set(list(adata.obs[data['grp'][0]])))
  ncharA = max([len(x) for x in a])
  w = len(a)/4+1
  h = ncharA/6+2.5
  ro = math.acos(10/max([15,ncharA]))/math.pi*180
  ##
  
  fig = plt.figure(figsize=[w,h])
  sc.pl.violin(adata,data['genes'],groupby=data['grp'][0],ax=fig.gca(),show=False)
  fig.autofmt_xdate(bottom=0.2,rotation=ro,ha='right')
  return iostreamFig(fig)

def PGV(data):
  # figure width and heights depends on number of unique categories
  # characters of category names, gene number
  adata = createData(data)
  if len(adata)==0:
    return Msg('No cells in the condition!')
  a = list(set(list(adata.obs[data['grp'][0]])))
  ncharA = max([len(x) for x in a])
  w = ncharA/8+len(data['genes'])/2+0.5
  h = len(a)+0.5
  swapAx = False
  ##
  if data['by']=='Columns':
    a = w
    w = h
    h = a
    swapAx = True

  fig = plt.figure(figsize=[w,h])
  axes = sc.pl.stacked_violin(adata,data['genes'],groupby=data['grp'][0],show=False,ax=fig.gca(),swap_axes=swapAx)
  for ax in axes:
    ax.tick_params(axis='y',which='major',labelsize=7)
  return iostreamFig(fig)
  
def pHeatmap(data):
  # figure width is depends on the number of categories was choose to show
  # and the character length of each category term
  # if the number of element in a category is smaller than 10, "Set1" or "Set3" is choosen
  # if the number of element in a category is between 10 and 20, default is choosen
  # if the number of element in a category is larger than 20, husl is choosen
  #Xsep = createData(data,True)
  #adata = sc.AnnData(Xsep['expr'],Xsep['obs'])
  adata = createData(data)
  Xdata = pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()
  
  exprOrder = True
  if data['order']!="Expression":
    exprOrder = False;
    s = adata.obs[data['order']]
    ix = sorted(range(len(s)), key=lambda k: s[k])
    adata = adata[ix,]
  colCounter = 0
  colName =['Set1','Set3']
  grpCol = list()
  grpLegend = list()
  grpWd = list()
  grpLen = list()
  h = 8
  w = len(data['genes'])/3+0.3
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

  w += 2 
  Zscore=None
  heatCol=None
  heatCenter=None
  colTitle="Expression"
  if data['norm']=='zscore':
    Zscore=1
    heatCol="vlag"
    heatCenter=0
    colTitle="Z score"

  g = sns.clustermap(pd.DataFrame(adata.X,index=list(adata.obs.index),columns=list(adata.var.index)),
                     method="ward",row_cluster=exprOrder,z_score=Zscore,cmap=heatCol,center=heatCenter,
                     row_colors=pd.concat(grpCol,axis=1).astype('str'),yticklabels=False,xticklabels=True,
                     figsize=(w,h),colors_ratio=0.05,
                     cbar_pos=(.3, .95, .55, .02),
                     cbar_kws={"orientation": "horizontal","label": colTitle,"shrink": 0.5})
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

  return json.dumps([iostreamFig(g),Xdata])#)#

def GD(data):
  adata = None;
  for one in data['cells'].keys():
    oneD = {'cells':data['cells'][one],
            'genes':[],
            'grp':[]}
    D = createData(oneD)
    D.obs['cellGrp'] = one
    if adata is None:
      adata = D
    else:
      adata = adata.concatenate(D)
  if adata is None:
    return ""
  ##
  adata.obs.astype('category')
  cutOff = 'geneN'+data['cutoff']
  adata.obs[cutOff] = adata.to_df().apply(lambda x: sum(x>float(data['cutoff'])),axis=1)
  ##
  w = 3
  if len(data['cells'])>1:
    w += 3
  fig = plt.figure(figsize=[w,4])
  sc.pl.violin(adata,cutOff,groupby='cellGrp',ax=fig.gca(),show=False,rotation=0,size=2)
  return iostreamFig(fig)

def DEG(data):
  res = requests.get('%s/annotations/var' % data["url"],params={'annotation-name':'name_0'})
  gNames = decode_fbs.decode_matrix_FBS(res.content)['columns'][0]
  fil = json.dumps({"mode":"topN","count":int(data['topN']),
                  "set1":{"filter":{"obs":{"index":list(data['cells']['group1'].values())}}},
                  "set2":{"filter":{"obs":{"index":list(data['cells']['group2'].values())}}}})
  res = requests.post('%s/diffexp/obs' % data["url"],fil,headers={'content-type':'application/json'})
  data = json.loads(res.content)
  diff = [[gNames[data[i][0]]]+['%.2f' % data[i][1]]+['%.4E' % data[i][j] for j in range(2,len(data[i]))] for i in range(len(data))]
  return json.dumps(diff)
  
def DOT(data):
  adata = createData(data)
  if len(adata)==0:
    return Msg('No cells in the condition!')
  #return adata
  a = list(set(list(adata.obs[data['grp'][0]])))
  ncharA = max([len(x) for x in a])
  w = ncharA/8+len(data['genes'])/2+0.5
  h = len(a)/4+0.5
  
  #fig = 
  sc.pl.dotplot(adata,data['geneGrp'],groupby=data['grp'][0],figsize=(w,h),show=False,expression_cutoff=float(data['cutoff']))#,show=False,ax=fig.gca()
  fig = plt.gcf()
  return iostreamFig(fig)
  
def EMBED(data):
  adata = createData(data)
  subSize = 4
  ncol = int(data['ncol'])
  ngrp = len(data['grp'])
  ngene = len(data['genes'])
  nrow = ngrp+math.ceil(ngene/ncol)
  
  step =11
  grpCol = {gID:math.ceil(len(list(adata.obs[gID].unique()))/step) for gID in data['grp']}
  
  fig = plt.figure(figsize=(ncol*subSize,subSize*nrow))
  gs = fig.add_gridspec(nrow,ncol)
  for i in range(ngrp):
      ax = sc.pl.umap(adata,color=data['grp'][i],ax=fig.add_subplot(gs[i,0]),show=False)
      if grpCol[data['grp'][i]]>3:
          ax.legend(ncol=grpCol[data['grp'][i]],loc=6,bbox_to_anchor=(1,0.5),frameon=False)
  for i in range(ngene):
      x = int(i/ncol)+ngrp
      y = i % ncol
      sc.pl.umap(adata,color=data['genes'][i],ax=fig.add_subplot(gs[x,y]),show=False)

  return iostreamFig(fig)
  
def TRACK(data):
  adata = createData(data)
  if len(adata)==0:
    return Msg('No cells in the condition!')
  w = math.log2(adata.n_obs)
  h = adata.n_vars/2
  
  ax = sc.pl.tracksplot(adata,data['geneGrp'],groupby=data['grp'][0],figsize=(w,h),show=False)
  fig=ax[0].figure
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

def DUAL(data):
  adata = createData(data)
  adata.obs['Expressed'] = adata.to_df().apply(cut,axis=1,args=(float(data['cutoff']),adata.var_names)).astype('category')
  pCol = {"None":"#AAAAAA44","Both":"#EDDF01AA",data['genes'][0]:"#1CAF82AA",data['genes'][1]:"#FA2202AA"}
  adata.uns["Expressed_colors"]=[pCol[i] for i in adata.obs['Expressed'].cat.categories]
  
  rcParams['figure.figsize'] = 4.5, 4
  if data['layout']=='umap':
    fig = sc.pl.umap(adata,color='Expressed',return_fig=True,show=False,legend_fontsize="small")
  if data['layout']=='tsne':
    fig = sc.pl.tsne(adata,color='Expressed',return_fig=True,show=False,legend_fontsize="small")
  
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
    
  sc.tl.rank_genes_groups(adata,groupby=data["grp"][0],n_genes=int(data['geneN']),method=data['markMethod'])
  sc.pl.rank_genes_groups(adata,n_genes=int(data['geneN']),ncols=3,show=False)
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

def version():
  print("1.0.2")
  ## 1.0.2: April 27, 2020
  ## 1. add Spinning button
  ## 2. Selection on both groups of selected cells
  ## ------------
  ## 1.0.1: April 26, 2020
  ## 1. Removed “close” & “Max” according to baohong's code;
  ## 2. Added footerToolbar according to baohong' code;
  ## 3. Relocated the “refresh” button;
  ## 4. Fixed a bug to handle any numerical group information;
  ## 5. Added the group element number into group information at Heatmap page.
  ## ---------------
  ## 1.0.2: May 6, 2020
  ## 1. Panel violin add an option to swap the axes
  ## 2. Provide the user to add the annotation abbreviation, as well as a customized categoryby combining annotation across existing categories.
  ## 3. Add dot plot as expression level and cell proportion for specified gene sets 
  ## ---------------
  ## 1.0.3: May 10, 2020
  ## 1. Keep the selection when refresh is clicked;
  ## 2. Multi-tsne/umap, embedding plots for genes; for annotations as well;
  ## 3. Batch adding genes with verifications as well as adding gene sets;
  ## 4. Download heatmap data including meta info as csv;
  ## 5. Uncheck All features and check all features with the dispatch method;
  ## 6. Updated using "_" to separate the combined groups;
  ## the following is required reinstall cellxgene
  ## 7. Change "PLOTTING PANEL" to "Visulization in Plugin";
  ## 8. Change biogenInterface to VIPInterface.py, and change the ajax call to VIP instead of biogen
  ## -------------------
  ## 1.0.4: May 13, 2020
  ## 1. Adding genes are case insensitive
  ## 2. update the data obtaining from ajax API to direct call by server.app.app.get_data_adaptor method
  ## 3. Trackplot;
  ## 4. Two gene embedding plots;
  ## ------------------
  ## 1.0.5: May 15, 2020
	## 1. Used “annotation” instead of group consistently across VIP;
	## 2. Removed "',' separated" from adding gene groups;
	## 3. Set the minimal value for gene expression in dual gene and dot plots
	## 4. Added a button to create Combine Annotation which can only be changed in "Combine & Abbr" tab
	## 5. Initialized the panel with full load main page by detecting the category number in window.store is NOT increased in 0.5 second interval
	## 6. Updated the gene detection plot with cut-off instead of pre-calculation
  ## 7. Change the location of the menu button from the top to the left side due to too many buttons;	
	## 8. Visualize marker genes with download;
	## 9. Fixed the tSNE/UMAP annotation legend is too wide by each annotation is in separated row while the number of gene plots per row is specified by user.
	##------------------------------
	## 1.0.6: not done yet
	## 1. DEG (it seems like I have difficulty to install/import diffxpy and I submitted an issue: https://github.com/theislab/diffxpy/issues/158);
	## 2. Save all user current information into file and load from a file;

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
