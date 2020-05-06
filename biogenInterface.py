import requests
import json
import server.app.decode_fbs as decode_fbs
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rcParams
import base64
import math
from io import BytesIO
import sys

import pprint
ppr = pprint.PrettyPrinter(depth=6)

#sys.setrecursionlimit(10000)
sc.set_figure_params(dpi=150, color_map='viridis')
sc.settings.verbosity = 2
rcParams.update({'figure.autolayout': True})

api_version = "/api/v0.2"

#urlExpr = 'http://127.0.0.1:8888/api/v0.2/data/var'
#urlGrp = 'http://127.0.0.1:8888/api/v0.2/annotations/obs'
#urlGen = 'http://127.0.0.1:8888/api/v0.2/annotations/var'
#urlDiff = "http://127.0.0.1:8888/api/v0.2/diffexp/obs"

def route(data,appConfig=None):#data,
  #adata = createData(data)
  #return(adata)
  #url = f"http://{appConfig.server__host}:{appConfig.server__port}|{appConfig.single_dataset__datapath}"
  if appConfig is None:
    data["url"] = f'http://127.0.0.1:8888/{api_version}'
  else:
    data = json.loads(str(data,encoding='utf-8'))
    data["url"] = f'http://{appConfig.server__host}:{appConfig.server__port}/{api_version}'

  return distributeTask(data["method"])(data)
  
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
  
def createData(data):
  #print("CreateData")
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
  obs = pd.DataFrame(obsL,index=['name_0']+data['grp'],columns=cNames).T
  strGN = [i for i in data['grp'] if 'genes' in i]
  if len(strGN)>0:
    obs[strGN] = obs[strGN].apply(pd.to_numeric)
  if combUpdate and len(data['grp'])>1:
    newGrp = 'Custom_combine'
    obs[newGrp] = obs[data['grp'][0]]
    for i in data['grp']:
      if i!=data['grp'][0]:
        obs[newGrp] += "|"+obs[i]
    expr = expr[~obs[newGrp].str.contains("Other")]
    obs = obs[~obs[newGrp].str.contains("Other")]
    data['grp'] = [newGrp]
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
    'DOT':DOT
  }.get(aTask,errorTask)

def iostreamFig(fig):
  figD = BytesIO()
  fig.savefig(figD,format='png',bbox_inches="tight")
  imgD = base64.encodestring(figD.getvalue()).decode("utf-8")
  figD.close()
  plt.close('all')
  return imgD

def SGV(data):
  # figure width and heights depends on number of unique categories
  # characters of category names, gene number
  adata = createData(data)
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
  adata = createData(data)
  #return adata
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
  #grpSize = list()
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
      #fSize = int(min([15,math.sqrt(len(grp))/len(Ugrp)+4]))
      #grpSize.append(fSize)
      #fW = int(max([len(x) for x in Ugrp])*fSize/15)/10
      #w += 1
      grpWd.append(max([len(x) for x in Ugrp]))#0.02*fW*max([len(x) for x in Ugrp])
      grpLen.append(len(Ugrp)+2)
  #grpW = [1.05]
  #for x in grpLen:
  #  oneW = 1.2*max([0.2,x/max(grpLen)])
  #  grpW.append(grpW[-1]+0.5)#oneW
    #w += 1#oneW*1.5
  #ix = sorted(range(len(grpLen)),key=lambda k:grpLen[k])

  w += 2 
  Zscore=None
  heatCol=None
  heatCenter=None
  colTitle="Expression"
  if data['norm']=='zscore':
    Zscore=1
    heatCol="vlag"
    heatCenter=0
    colTitle="column Z score"
  g = sns.clustermap(pd.DataFrame(adata.X,index=list(adata.obs['name_0']),columns=list(adata.var.index)),
                     method="ward",row_cluster=exprOrder,z_score=Zscore,cmap=heatCol,center=heatCenter,
                     row_colors=pd.concat(grpCol,axis=1),yticklabels=False,xticklabels=True,
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

  return iostreamFig(g)

def GD(data):
  adata = None;
  for one in data['cells'].keys():
    oneD = {'cells':data['cells'][one],
            'grp':data['grp'],
            'url':data['url']}
    D = createData(oneD)
    D.obs['grp'] = one
    if adata is None:
      adata = D
    else:
      adata = adata.concatenate(D)
  if adata is None:
    return ""
  ##
  w = 3
  if len(data['cells'])>1:
    w += 3
  fig = plt.figure(figsize=[w,4])
  sc.pl.violin(adata,data['grp'],groupby='grp',ax=fig.gca(),show=False,rotation=0,size=2)
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
  #return adata
  a = list(set(list(adata.obs[data['grp'][0]])))
  ncharA = max([len(x) for x in a])
  w = ncharA/8+len(data['genes'])/2+0.5
  h = len(a)/4+0.5
  
  fig = sc.pl.dotplot(adata,data['geneGrp'],groupby=data['grp'][0],figsize=(w,h),show=False,expression_cutoff=float(data['cutoff']))#,show=False,ax=fig.gca()
  
  return iostreamFig(fig)
#Fxyd3,Lgals7
#Acta2,Krt5
#Cd79a,Rac2,Coro1a
  
  
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
  ## 3. Add dot plot as expression level and cell propotion for sepecified gene sets 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
