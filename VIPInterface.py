import requests
import json
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
import base64
import math
from io import BytesIO
import sys

import pprint
ppr = pprint.PrettyPrinter(depth=6)

import server.app.app as app
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
    data["url"] = f'http://{appConfig.server__host}:{appConfig.server__port}/{api_version}'
  
  if 'figOpt' in data.keys():
    setFigureOpt(data['figOpt'])
  
  return distributeTask(data["method"])(data)

def setFigureOpt(opt):
  sc.set_figure_params(dpi_save=int(opt['dpi']),fontsize= int(opt['fontsize']),vector_friendly=(opt['vectorFriendly']=='true'),transparent=(opt['transparent']=='true'),color_map=opt['colorMap'])
  rcParams.update({'savefig.format':opt['img']})

def subData(data):
  selC = list(data['cells'].values())
  cNames = ["cell%d" %i for i in selC]
  
  fSparse = False
  ## onbtain the expression matrix
  gNames = []
  X = []
  if 'genes' in data.keys():
    with app.get_data_adaptor() as scD:
      if not type(scD.data.X) is np.ndarray:
        fSparse = True
      if len(data['genes'])>0:
        fullG = list(scD.data.var['name_0'])
        selG = [fullG.index(i) for i in data['genes']]
        X = scD.data.X[selC][:,selG]
        gNames = data['genes']
      else:
        X = scD.data.X[selC]
        gNames = list(scD.data.var["name_0"])
  if fSparse:
    expr = X
  else:
    expr = pd.DataFrame(X,columns=gNames,index=cNames)
  #ppr.pprint("Finished expression...")
  
  ## obtain the embedding
  strEmbed = 'umap'
  embed = pd.DataFrame([[0 for x in range(len(cNames))] for i in range(2)],
                        index=['%s1'%strEmbed,'%s2'%strEmbed],columns=cNames).T
  if 'layout' in data.keys():## tsne or umap
    strEmbed = data['layout']
    with app.get_data_adaptor() as scD:
      embed = pd.DataFrame(scD.data.obsm['X_%s'%strEmbed][selC],columns=['%s1'%strEmbed,'%s2'%strEmbed],index=cNames)

  ## obtain the category annotation
  with app.get_data_adaptor() as scD:
    obs = scD.data.obs.loc[selC,['name_0']+data['grp']].astype('str')
  obs.index = cNames

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
    selC = ~obs[newGrp].str.contains("Other").to_numpy()
    expr = expr[selC]
    embed = embed[selC]
    obs = obs[selC]
    data['grp'] = [newGrp]
    
  obs = obs.astype('category')
  ## empty selection
  if expr.shape[0]==0 or expr.shape[1]==0:
    return []
  return sc.AnnData(expr,obs,var=pd.DataFrame([],index=gNames),obsm={'X_%s'%strEmbed:embed.to_numpy()})

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
    ix = adata.to_df().apply(lambda x: max(x)>float(cutoff),axis=1)
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
  adata = geneFiltering(adata,data['cutoff'],1)
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

  if '1.4.7' in sc.__version__:
    vp = sc.pl.stacked_violin(adata,data['genes'],groupby=data['grp'][0],return_fig=True,figsize=(w,h),swap_axes=swapAx)
    vp.add_totals().style(yticklabels=True).show()
    fig = plt.gcf()
  else:
    fig = plt.figure(figsize=[w,h])
    axes = sc.pl.stacked_violin(adata,data['genes'],groupby=data['grp'][0],show=False,ax=fig.gca(),swap_axes=swapAx)

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
    colTitle="column Z score"

  g = sns.clustermap(adata.to_df(),
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
  ppr.pprint("finished legend")
  return json.dumps([iostreamFig(g),Xdata])#)#

def GD(data):
  adata = None;
  for one in data['cells'].keys():
    oneD = {'cells':data['cells'][one],
            'genes':[],
            'grp':[],
            'url':data['url']}
    D = createData(oneD)
    D.obs['cellGrp'] = one
    if adata is None:
      adata = D
    else:
      adata = adata.concatenate(D)
  if adata is None:
    return Msg("No cells were satisfied the condition!")
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
  adata = None;
  for one in data['cells'].keys():
    oneD = {'cells':data['cells'][one],
            'genes':[],
            'grp':[],
            'url':data['url']}
    D = createData(oneD)
    D.obs['cellGrp'] = one
    if adata is None:
      adata = D
    else:
      adata = adata.concatenate(D)
  if adata is None:
    return Msg("No cells were satisfied the condition!")
  
  adata.obs.astype('category')
  nm = None
  if data['DEmethod']=='wald': 
    nm = 'nb'
  res = de.test.two_sample(adata,'cellGrp',test=data['DEmethod'],noise_model=nm)
  deg = res.summary()
  deg = deg.sort_values(by=['qval']).loc[:,['gene','log2fc','pval','qval']]
  deg = deg.iloc[range(int(data['topN'])),]
  deg.loc[:,'log2fc'] = deg.loc[:,'log2fc'].apply(lambda x: '%.2f'%x)
  deg.loc[:,'pval'] = deg.loc[:,'pval'].apply(lambda x: '%.4E'%x)
  deg.loc[:,'qval'] = deg.loc[:,'qval'].apply(lambda x: '%.4E'%x)
  
  return json.dumps(deg.values.tolist())

def DOT(data):
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
  
  if '1.4.7' in sc.__version__:
    dp = sc.pl.dotplot(adata,data['geneGrp'],groupby=data['grp'][0],expression_cutoff=float(data['cutoff']),return_fig=True)#
    dp = dp.add_totals(size=1.2).legend(show_size_legend=True).style(cmap='Blues', dot_edge_color='black', dot_edge_lw=1, size_exponent=1.5)
    dp.show()
    fig = dp.get_axes()['mainplot_ax'].figure
  else:
    sc.pl.dotplot(adata,data['geneGrp'],groupby=data['grp'][0],figsize=(w,h),show=False,expression_cutoff=float(data['cutoff']))
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
  
  rcParams['figure.constrained_layout.use'] = False
  fig = plt.figure(figsize=(ncol*subSize,subSize*nrow))
  gs = fig.add_gridspec(nrow,ncol,wspace=0.2)
  for i in range(ngrp):
      ax = getattr(sc.pl,data['layout'])(adata=adata,color=data['grp'][i],ax=fig.add_subplot(gs[i,0]),wspace=0.25,show=False)
      if grpCol[data['grp'][i]]>1:
          ax.legend(ncol=grpCol[data['grp'][i]],loc=6,bbox_to_anchor=(1,0.5),frameon=False)
  for i in range(ngene):
      x = int(i/ncol)+ngrp
      y = i % ncol
      getattr(sc.pl,data['layout'])(adata,color=data['genes'][i],ax=fig.add_subplot(gs[x,y]),wspace=0.25,show=False)

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
  fig = getattr(sc.pl,data['layout'])(adata,color='Expressed',return_fig=True,show=False,legend_fontsize="small")
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
  
  if len(adata.obs[data['grp'][0]].unique())<2:
    return json.dumps([[['name','scores'],['None','0']],Msg('Less than 2 groups in selected cells!')])
  
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
  print("1.0.6")
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
	## 1.0.6: not done yet
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

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
