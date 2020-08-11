import requests
import json
import traceback
import server.app.decode_fbs as decode_fbs
import scanpy as sc
import pandas as pd
import numpy as np
import diffxpy.api as de
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib import rcParams
import plotly.graph_objects as go
import plotly.io as plotIO
import base64
import math
from io import BytesIO
import sys
import time
import os
import subprocess
strExePath = os.path.dirname(os.path.abspath(__file__))

import pprint
ppr = pprint.PrettyPrinter(depth=6)

import server.app.app as app
import server.compute.diffexp_generic as diffDefault
import pickle

sys.setrecursionlimit(10000)
sc.settings.verbosity = 2
rcParams.update({'figure.autolayout': True})

api_version = "/api/v0.2"

def route(data,appConfig=None):
  if appConfig is None:
    data["url"] = f'http://127.0.0.1:8888/{api_version}'
  else:
    data = json.loads(str(data,encoding='utf-8'))
    data["url"] = f'http://localhost:{appConfig.server__port}/{api_version}'#{appConfig.server__host}
  #ppr.pprint(data['figOpt']) 
  if 'figOpt' in data.keys():
    setFigureOpt(data['figOpt'])
  try:
    return distributeTask(data["method"])(data)
  except Exception as e:
    return 'ERROR @server: '+traceback.format_exc() # 'ERROR @server: {}, {}'.format(type(e),str(e))
  #return distributeTask(data["method"])(data)

def setFigureOpt(opt):
  sc.set_figure_params(dpi_save=int(opt['dpi']),fontsize= float(opt['fontsize']),vector_friendly=(opt['vectorFriendly'] == 'Yes'),transparent=(opt['transparent'] == 'Yes'),color_map=opt['colorMap'])
  rcParams.update({'savefig.format':opt['img']})

def getObs(data):
  selC = list(data['cells'].values())
  cNames = ["cell%d" %i for i in selC]
  ## obtain the category annotation
  with app.get_data_adaptor() as scD:
    obs = scD.data.obs.loc[selC,['name_0']+data['grp']].astype('str')
  obs.index = cNames
  ## update the annotation Abbreviation
  combUpdate = cleanAbbr(data)
  if 'abb' in data.keys():
    for i in data['grp']:
      obs[i] = obs[i].map(data['abb'][i])
  return combUpdate, obs
  
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

def subData(data):
  selC = list(data['cells'].values())
  cNames = ["cell%d" %i for i in selC]
  
  ## onbtain the expression matrix
  gNames = []
  expr = []
  if True:
    fSparse = False
    X = []
    if 'genes' in data.keys():
      with app.get_data_adaptor() as scD:
        if not type(scD.data.X) is np.ndarray:
          fSparse = True
        if len(data['genes'])>0:
          fullG = list(scD.data.var['name_0'])
          selG = [fullG.index(i) for i in data['genes']]
          X = scD.data.X[:,selG]
          gNames = data['genes']
        else:
          X = scD.data.X
          gNames = list(scD.data.var["name_0"])
      if 'figOpt' in data.keys() and data['figOpt']['scale'] == 'Yes':
        X = sc.pp.scale(X,zero_center=(data['figOpt']['scaleZero'] == 'Yes'),max_value=(float(data['figOpt']['scaleMax']) if data['figOpt']['clipValue']=='Yes' else None))
      X = X[selC]
    if fSparse:
      expr = X
    else:
      expr = pd.DataFrame(X,columns=gNames,index=cNames)
  else:
    fSparse = False
    if 'genes' in data.keys():
      with app.get_data_adaptor() as scD:
        if not type(scD.data.X) is np.ndarray:
          fSparse = True
        if len(data['genes'])>0:
          fullG = list(scD.data.var['name_0'])
          selG = [fullG.index(i) for i in data['genes']]
          X = scD.data.X[selC][:,selG]
          gNames = data['genes']
          if fSparse:
            expr = pd.DataFrame.sparse.from_spmatrix(X,index=cNames,columns=gNames)
          else:
            expr = pd.DataFrame(X,columns=gNames,index=cNames) 
        else:
          X = scD.data.X[selC]
          gNames = list(scD.data.var["name_0"])
          if fSparse:
            expr = X
          else:
            expr = pd.DataFrame(X,columns=gNames,index=cNames) 
  expr,gNames = collapseGeneSet(data,expr,gNames,cNames,fSparse)
  
  ## obtain the embedding
  if False:
    strEmbed = 'umap'
    #embed = pd.DataFrame([[0 for x in range(len(cNames))] for i in range(2)],
    #                      index=['%s1'%strEmbed,'%s2'%strEmbed],columns=cNames).T
    embed = pd.DataFrame([],index=cNames)
    if 'layout' in data.keys():## tsne or umap
      strEmbed = data['layout']
      with app.get_data_adaptor() as scD:
        embed = pd.DataFrame(scD.data.obsm['X_%s'%strEmbed][selC][:,[0,1]],columns=['%s1'%strEmbed,'%s2'%strEmbed],index=cNames)
      strEmbed = 'umap'

  embed = {}
  if 'layout' in data.keys():
    layout = data['layout']
    if isinstance(layout,str):
      layout = [layout]
    if len(layout)>0:
      for one in layout:
        with app.get_data_adaptor() as scD:
          embed['X_%s'%one] = pd.DataFrame(scD.data.obsm['X_%s'%one][selC][:,[0,1]],columns=['%s1'%one,'%s2'%one],index=cNames)

  ## obtain the category annotation
  combUpdate, obs = getObs(data)

  ## create a custom annotation category and remove cells which are not in the selected annotation
  if combUpdate and len(data['grp'])>1:
    newGrp = 'Custom_combine'
    obs[newGrp] = obs[data['grp'][0]]
    for i in data['grp']:
      if i!=data['grp'][0]:
        obs[newGrp] += ":"+obs[i]
    selC = ~obs[newGrp].str.contains("Other").to_numpy()
    expr = expr[selC]
    for i in embed.keys():
      embed[i] = embed[i][selC]
    obs = obs[selC]
    data['grp'] = [newGrp]
    
  obs = obs.astype('category')
  ## empty selection
  if expr.shape[0]==0 or expr.shape[1]==0:
    return []
    
  return sc.AnnData(expr,obs,var=pd.DataFrame([],index=gNames),obsm={layout:embed[layout].to_numpy() for layout in embed.keys()})
  #return sc.AnnData(expr,obs,var=pd.DataFrame([],index=gNames),obsm={'X_%s'%strEmbed:embed.to_numpy()})

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
  return subData(data)
  #print("CreateData")
  with app.get_data_adaptor() as scD:
    if (type(scD.data.X) is np.ndarray):
      return subData(data)

  headers = {'content-type':'application/json'}
  # obtain the expression
  res = requests.get('%s/annotations/var' % data["url"],params={'annotation-name':'name_0'})
  gNames = decode_fbs.decode_matrix_FBS(res.content)['columns'][0]
  if not 'genes' in data:
    data['genes'] = gNames[0]
  if len(data['genes'])==0:# obtain all genes
    data['genes'] = gNames
    
  fil = json.dumps({'filter':{'var':{'annotation_value':[{'name':'name_0','values':data['genes']}]}}})
  res = requests.put('%s/data/var' % data["url"],fil,headers=headers)
  expr = decode_fbs.decode_matrix_FBS(res.content)
  cNames = ["cell%d" % x for x in data['cells'].values()]
  expr = pd.DataFrame([[expr['columns'][i][x] for x in data['cells'].values()] for i in range(len(expr['columns']))],
                        index=[gNames[x] for x in expr['col_idx']],columns=cNames).T

  ## obtain the embedding
  strEmbed = 'umap'
  embed = pd.DataFrame([[0 for x in range(len(cNames))] for i in range(2)],
                        index=['%s1'%strEmbed,'%s2'%strEmbed],columns=cNames).T
  if 'layout' in data.keys():## tsne or umap
    strEmbed = data['layout']
    res = requests.get('%s/layout/obs' % data["url"],params={'layout-name':strEmbed})
    embed= decode_fbs.decode_matrix_FBS(res.content)
    embed = pd.DataFrame([[embed['columns'][i][x] for x in data['cells'].values()] for i in range(len(embed['columns']))],
                          index=embed['col_idx'],columns=cNames).T

  # obtain the meta grouping
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

def errorTask(data):
  raise ValueError('Error task!')
  
def distributeTask(aTask):
  return {
    'SGV':SGV,
    'PGV':PGV,
    'VIOdata':VIOdata,
    'HEATplot':pHeatmap,
    'HEATdata':HeatData,
    'GD':GD,
    'DEG':DEG,
    'DOT':DOT,
    'EMBED':EMBED,
    'TRAK':TRACK,
    'DUAL':DUAL,
    'MARK': MARK,
    'MINX':MINX,
    'DENS':DENS,
    'DENS2D':DENS2D,
    'SANK':SANK,
    'STACBAR':STACBAR,
    'CLI':CLI
  }.get(aTask,errorTask)

def iostreamFig(fig):
  figD = BytesIO()
  fig.savefig(figD,bbox_inches="tight")
  imgD = base64.encodebytes(figD.getvalue()).decode("utf-8")
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

  adata = createData(data)
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
  sc.pl.violin(adata,data['genes'],groupby=data['grp'][0],ax=fig.gca(),show=False)
  fig.autofmt_xdate(bottom=0.2,rotation=ro,ha='right')
  return iostreamFig(fig)

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
  # characters of category names, gene number
  updateGene(data)
  adata = createData(data)
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
    fig = plt.gcf()
  else:
    fig = plt.figure(figsize=[w,h])
    axes = sc.pl.stacked_violin(adata,data['genes'],groupby=data['grp'][0],show=False,ax=fig.gca(),swap_axes=swapAx,
                                var_group_positions=data['grpLoc'],var_group_labels=data['grpID'])
  return iostreamFig(fig)

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
  #Xdata = pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()
  #ppr.pprint('HEAT data reading cost %f seconds' % (time.time()-sT) )
  #sT = time.time()
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
    colTitle="Z-score"
  #ppr.pprint('HEAT data preparing cost %f seconds' % (time.time()-sT) )
  #sT = time.time()
  g = sns.clustermap(adata.to_df(),
                     method="ward",row_cluster=exprOrder,z_score=Zscore,cmap=heatCol,center=heatCenter,
                     row_colors=pd.concat(grpCol,axis=1).astype('str'),yticklabels=False,xticklabels=True,
                     figsize=(w,h),colors_ratio=0.05,
                     cbar_pos=(.3, .95, .55, .02),
                     cbar_kws={"orientation": "horizontal","label": colTitle,"shrink": 0.5})
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

def HeatData(data):
  adata = createData(data)
  Xdata = pd.concat([adata.to_df(),adata.obs], axis=1, sort=False).to_csv()
  return Xdata
  
def GD(data):
  adata = None;
  for one in data['cells'].keys():
    sT = time.time()
    oneD = {'cells':data['cells'][one],
            'genes':[],
            'grp':[],
            'figOpt':data['figOpt'],
            'url':data['url']}
    D = createData(oneD)
    ppr.pprint("one grp aquire data cost %f seconds" % (time.time()-sT))
    D.obs['cellGrp'] = one
    if adata is None:
      adata = D
    else:
      sT =time.time()
      adata = adata.concatenate(D)
      ppr.pprint("Concatenate data cost %f seconds" % (time.time()-sT))
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

def DEG(data):
  
  adata = None;
  genes = data['genes']
  data['genes'] = []
  comGrp = 'cellGrp'
  if 'combine' in data.keys():
    data['figOpt']['scale'] = 'false'
    adata = createData(data)
    comGrp = data['grp'][0]
    mask = [adata.obs[comGrp].isin(data['comGrp'][i]) for i in [0,1]]
    adata = adata[adata.obs[comGrp].isin(data['comGrp'])]
  else:
    mask = [pd.Series(range(data['cellN'])).isin(data['cells'][one].values()) for one in data['cells'].keys()]
    for one in data['cells'].keys():
      oneD = {'cells':data['cells'][one],
              'genes':[],
              'grp':[],
              'figOpt':{'scale':'false'},
              'url':data['url']}
      D = createData(oneD)
      D.obs[comGrp] = one
      if adata is None:
        adata = D
      else:
        adata = adata.concatenate(D)
  #ppr.pprint(adata) 
  #with open("adata.pkl",'wb') as f:
  #  pickle.dump(adata,f)
  if not 'AnnData' in str(type(adata)):
    raise ValueError('No data extracted by user selection')
  if data['DEmethod']=='t-test':
    with app.get_data_adaptor() as scD:
      res = diffDefault.diffexp_ttest(scD,mask[0].to_numpy(),mask[1].to_numpy(),scD.data.shape[0])
      gNames = list(scD.data.var["name_0"])
    deg = pd.DataFrame(res,columns=['gID','log2fc','pval','qval'])
    gName = pd.Series([gNames[i] for i in deg['gID']],name='gene')
    deg = pd.concat([deg,gName],axis=1).loc[:,['gene','log2fc','pval','qval']]
  else:
    adata.obs.astype('category')
    nm = None
    if data['DEmethod']=='wald': 
      nm = 'nb'
    res = de.test.two_sample(adata,comGrp,test=data['DEmethod'],noise_model=nm)
    deg = res.summary()
    deg = deg.sort_values(by=['qval']).loc[:,['gene','log2fc','pval','qval']]
  ## plot in R
  strF = ('/tmp/DEG%f.csv' % time.time())
  deg.to_csv(strF,index=False)
  res = subprocess.run([strExePath+'/volcano.R',strF,';'.join(genes),data['figOpt']['img'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi']),str(data['logFC'])],capture_output=True)#
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  #####
  deg = deg.iloc[range(int(data['topN'])),]
  deg.loc[:,'log2fc'] = deg.loc[:,'log2fc'].apply(lambda x: '%.2f'%x)
  deg.loc[:,'pval'] = deg.loc[:,'pval'].apply(lambda x: '%.4E'%x)
  deg.loc[:,'qval'] = deg.loc[:,'qval'].apply(lambda x: '%.4E'%x)
  
  return json.dumps([deg.values.tolist(),img])

def DOT(data):
  updateGene(data)
  adata = createData(data)
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

  return iostreamFig(fig)

def EMBED(data):
  adata = createData(data)
  subSize = 4
  ncol = int(data['ncol'])
  ngrp = len(data['grp'])
  ngene = len(data['genes'])
  nrow = ngrp+math.ceil(ngene/ncol)
  if 'splitGrp' in data.keys():
    splitName = list(adata.obs[data['splitGrp']].unique())
    nsplitRow = math.ceil(len(splitName)/ncol)
    nrow = ngrp+ngene*nsplitRow
  
  step =11
  grpCol = {gID:math.ceil(len(list(adata.obs[gID].unique()))/step) for gID in data['grp']}
  
  rcParams['figure.constrained_layout.use'] = False
  fig = plt.figure(figsize=(ncol*subSize,subSize*nrow))
  gs = fig.add_gridspec(nrow,ncol,wspace=0.2)
  for i in range(ngrp):
      ax = sc.pl.embedding(adata,data['layout'],color=data['grp'][i],ax=fig.add_subplot(gs[i,0]),show=False)#,wspace=0.25
      if grpCol[data['grp'][i]]>1:
          ax.legend(ncol=grpCol[data['grp'][i]],loc=6,bbox_to_anchor=(1,0.5),frameon=False)
      ax.set_xlabel('%s1'%data['layout'])
      ax.set_ylabel('%s2'%data['layout'])

  if 'splitGrp' in data.keys():
    vMax = adata.to_df().apply(lambda x: max(x))
    vMin = adata.to_df().apply(lambda x: min(x))
    dotSize = 120000 / adata.n_obs
    for i in range(ngene):
      for j in range(len(splitName)):
        x = ngrp + i*nsplitRow+int(j/ncol)
        y = j % ncol
        ax = sc.pl.embedding(adata,data['layout'],ax=fig.add_subplot(gs[x,y]),show=False)#color=data['genes'][i],wspace=0.25,
        ax = sc.pl.embedding(adata[adata.obs[data['splitGrp']]==splitName[j]],data['layout'],color=data['genes'][i],
                vmin=vMin[data['genes'][i]],vmax=vMax[data['genes'][i]],ax=ax,show=False,
                size=dotSize,title='{} in {}'.format(data['genes'][i],splitName[j]))
        ax.set_xlabel('%s1'%data['layout'])
        ax.set_ylabel('%s2'%data['layout'])
  else:
    for i in range(ngene):
        x = int(i/ncol)+ngrp
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
  #if len(data['grpLoc'])>0 and data['grpLoc'][len(data['grpLoc'])-1][1] < (len(data['genes'])-1):
  #  data['grpLoc'] += [(data['grpLoc'][len(data['grpLoc'])-1][1]+1,len(data['genes'])-1)]
  #  data['grpID'] += ['others']
  ##############
  ppr.pprint(data['genes'])
  
  ax = sc.pl.tracksplot(adata,data['genes'],groupby=data['grp'][0],figsize=(w,h),
                        var_group_positions=data['grpLoc'],var_group_labels=data['grpID'],
                        show=False)
  ppr.pprint("test")
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
  #sT = time.time()
  adata = createData(data)
  #ppr.pprint('DUAL data reading cost %f seconds' % (time.time()-sT) )
  #sT = time.time()
  #adata.obs['Expressed'] = adata.to_df().apply(cut,axis=1,args=(float(data['cutoff']),adata.var_names)).astype('category')
  adata.obs['Expressed'] = dualExp(adata.to_df(),float(data['cutoff']),adata.var_names)
  #ppr.pprint(data['cutoff'])
  #ppr.pprint(adata.obs['Expressed'].cat.categories)
  #ppr.pprint(adata.obs['Expressed'].value_counts())
  
  #ppr.pprint('DUAL filtering cost %f seconds' % (time.time()-sT) )
  #sT = time.time()
  pCol = {"None":"#AAAAAA44","Both":"#EDDF01AA",data['genes'][0]:"#1CAF82AA",data['genes'][1]:"#FA2202AA"}
  adata.uns["Expressed_colors"]=[pCol[i] for i in adata.obs['Expressed'].cat.categories]
  
  rcParams['figure.figsize'] = 4.5, 4
  fig = sc.pl.embedding(adata,data['layout'],color='Expressed',return_fig=True,show=False,legend_fontsize="small")
  plt.xlabel('%s1'%data['layout'])
  plt.ylabel('%s2'%data['layout'])
  rcParams['figure.figsize'] = 4, 4
  #ppr.pprint('DUAL plotting cost %f seconds' % (time.time()-sT) )
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
    
  sc.tl.rank_genes_groups(adata,groupby=data["grp"][0],n_genes=2,method=data['markMethod'])#int(data['geneN'])
  ppr.pprint(int(data['geneN']))
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
  adata.obs['None'] = 'all'
  bw=float(data['bw'])
  sGrp = data['category'][0]
  cGrp = data['category'][1]
  
  defaultFontsize = 16
  if 'figOpt' in data.keys():
    defaultFontsize = float(data['figOpt']['fontsize'])
  subSize = 4
  split = list(adata.obs[sGrp].unique())
  genes = list(adata.var.index)
  colGrp = list(adata.obs[cGrp].unique())
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
          sns.kdeplot(D[Dobs==one][genes[j]].to_numpy(),bw=bw,label=one)

      if i==0:
        ax.set_title(genes[j],fontsize=defaultFontsize+2)
      if j==0:
        ax.set_ylabel(split[i],fontsize=defaultFontsize)
      if i==0 and j==(len(genes)-1):
        ax.legend(prop={'size': 10},title = cGrp,loc=2,bbox_to_anchor=(1,1),ncol=legendCol,frameon=False)#
      else:
        ax.get_legend().remove()
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
  if 'name_0' in D.columns:
    del D['name_0']
  
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
    summaryOne = D.groupby(oneName).size().reset_index(name='Count')
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
  strF = ('/tmp/DEG%f.csv' % time.time())
  adata.to_df().to_csv(strF)#
  res = subprocess.run([strExePath+'/Density2D.R',strF,data['figOpt']['img'],str(data['cutoff']),str(data['bandwidth']),data['figOpt']['colorMap'],str(data['figOpt']['fontsize']),str(data['figOpt']['dpi'])],capture_output=True)#
  if 'Error' in res.stderr.decode('utf-8'):
    raise ValueError(res.stderr.decode('utf-8'))
  img = res.stdout.decode('utf-8')
  os.remove(strF)
  
  return img

def toInt(x):
  if len(x)==0:
    return 0
  return int(x)
  
def STACBAR(data):
  if len(data['genes'])==0:
    tmp, D = getObs(data)
    D = D.apply(lambda x:x.apply(lambda y:y))
  else:
    adata = createData(data)
    D = pd.concat([adata.obs.apply(lambda x:x.apply(lambda y:y)),
                   adata.to_df().apply(lambda x:pd.cut(x,int(data['Nbin'])).apply(lambda y:'%.1f_%.1f'%(y.left,y.right)))],
                  axis=1,sort=False)
  D = D.astype('str').astype('category')
  if 'name_0' in D.columns:
    del D['name_0']
  cellN = D.groupby(list(D.columns)).size().reset_index(name="Count")
  
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

def ajaxData(strData):
  with open(strData,'rb') as f:
    data=pickle.load(f)
  url = 'http://localhost:8888/api/v0.2'#data['url']#
  genes = data['genes']#['BTK','SALL1']
  #import random
  #cix = random.sample(range(40000),10000)
  cells = data['cells']#{str(x):cix[x] for x in range(len(cix))}
  layout= data['layout']#['umap_harmony','umap_liger']
  grps = data['grp']#['cell_type','diagnosis']
  
  ## 
  headers = {'content-type':'application/json'}
  #### obtain the expression -------------
  res = requests.get('%s/annotations/var' % url,params={'annotation-name':'name_0'})
  gNames = decode_fbs.decode_matrix_FBS(res.content)['columns'][0]
  if len(genes)==0:
    genes=gNames
  fil = json.dumps({'filter':{'var':{'annotation_value':[{'name':'name_0','values':genes}]}}})
  res = requests.put('%s/data/var' % url,fil,headers=headers)    
  expr = decode_fbs.decode_matrix_FBS(res.content)  
  cNames = ["cell%d" % x for x in cells.values()]
  expr = pd.DataFrame([[expr['columns'][i][x] for x in cells.values()] for i in range(len(expr['columns']))],
                          index=[gNames[x] for x in expr['col_idx']],columns=cNames).T
  
  #### obtain the layout -------------
  layX = {}
  if len(layout)>0:
    for one in layout:
      res = requests.get('%s/layout/obs' % url,params={'layout-name':one})
      embed= decode_fbs.decode_matrix_FBS(res.content)
      embed = pd.DataFrame([[embed['columns'][i][x] for x in cells.values()] for i in range(len(embed['columns']))],
                            index=embed['col_idx'],columns=cNames).T
      layX['X_%s'%one]=embed.to_numpy()
  
  ## obtain the meta ---------
  obsL = [cNames]
  for one in grps:
    res = requests.get('%s/annotations/obs' % url,params={'annotation-name':one})
    grp = decode_fbs.decode_matrix_FBS(res.content)
    subGrp = [str(grp['columns'][0][i]) for i in cells.values()]
    obsL += [subGrp]
  obs = pd.DataFrame(obsL,index=['name_0']+grps,columns=cNames).T
  obs = obs.astype('category')
  
  adata = sc.AnnData(expr,obs,obsm=layX)
  return adata

def CLI(data):
  strPath = ('/tmp/CLI%f' % time.time())
  script = data['script']
  del data['script']
  
  adata = createData(data)

  strData = strPath + '.pkl'
  with open(strData,'wb') as f:
    pickle.dump(adata,f)
  
  strScript = strPath + '.py'
  #addedScript=['import os','os.chdir("%s")'%strExePath,'import VIPInterface as vip','adata=vip.ajaxData("%s")'%strData]
  with open(strScript,'w') as f:
    f.writelines(['import pickle\n','with open("%s","rb") as f:\n'%strData,'  adata=pickle.load(f)\n\n'])
    f.write(script)
  
  res = subprocess.run('jupytext --to notebook --output - %s | jupyter nbconvert --ExecutePreprocessor.timeout=600 --to html --execute --stdin --stdout'%strScript,capture_output=True,shell=True)
  if 'Error' in res.stderr.decode('utf-8'):
    raise ValueError(res.stderr.decode('utf-8'))
  html = res.stdout.decode('utf-8')
  h,s,e = html.partition('<div class="cell border-box-sizing code_cell rendered">')
  h1,s,e = e.partition('<div class="cell border-box-sizing code_cell rendered">')
  html = h+s+e
  os.remove(strData)
  os.remove(strScript)
  return html
  
def version():
  print("1.0.8")
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
	## 6. Added a gene expression cut-off in the gene detection plot instead of pre-calculated in the meta data
  ## 7. Change the location of the menu button from the top to the left side due to too many buttons;	
	## 8. Visualize marker genes with download;
	## 9. Fixed the tSNE/UMAP when annotation legend is too wide with each annotation is in a separated row while the number of gene plots per row is specified by user.
	##------------------------------
	## 1.0.6: May 17, 2020
	## 1. Update the scanpy to incorporate a fancier dotplot and updated the stacked violin plots;
	## 2. Disabled the panel button until all categorical information ajaxed over, buttons will be automatically enabled after;
	## 3. Adjust the image panel not to be limited by the VIP window size;
	## 4. Fixed a bug of ignoring the tsne for embeding plot;
	## 5. Adjust the embeding plot to improve the legend for annotation
	## 6. Added handing ajax error to return control to the users
	## 7. Added a new tab for figure setting, update the python to respond to those setting
	## 8. Added a spliter between side menu buttons and content using more efficient JS code
	## 9. Other UI bug fix pointed by Baohong
	## 10. Add different DEG methods from diffxpy;
	## 11. Save all user current information into local file and load from a local file;
  ## --------------------------
  ## 1.0.7: May 26, 2020
  ## 1. Performance (time) was significantly improved for several plots (gene detection, violin, stack violin, tSNE/UMAP,...) on large data half million cells
  ## 2. Separated the heatmap data downloading from heatmap plotting, which improved the time for heatmap plotting
  ## 3. Added gene expression density plots splitted by one annoation and colored by one annotation
  ## -------------------------
  ## 1.0.8: May 28, 2020
  ## 1. Optimize the legend for the density plots;
  ## 2. Add the display on cell numbers for custom combined annotations;
  ## 3. Add DEG option on custom combined annotations;
  ## 4. Add the python error return to the user interface.
  ## -------------------------
  ## 1.0.9: June 3, 2020
  ## 1. Add the annotation split for gene express in tsne/umap plot 
  ## 2. Add the gene sets selection for stack violin
  ## 3. Add the gene selection for Dot plot and track plot
  ## -------------------------
  ## 1.0.10: Jun 8, 2020
  ## 1. Add Sankey diagram
  ## -----------------------
  ## 1.0.11: June 14, 2020
  ## 1. volcano plot added for DEG
  
  
  
  
  
  
  
  
  
  
  
  
  
