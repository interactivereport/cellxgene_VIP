import sys,json,re,time,warnings,math,os,contextlib,textwrap,traceback,distinctipy,resource
from datetime import timedelta
import pandas as pd
import seaborn as sns
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import fastcluster as fc
from scipy.cluster import hierarchy
from difflib import SequenceMatcher
import PyComplexHeatmap as pch
warnings.simplefilter("ignore", UserWarning)
verbose=False
def main():
  if len(sys.argv)==1:
    data = json.load(sys.stdin)
  else:
    with open(sys.argv[1],'r') as f:
      data = json.load(f)
  try:
    taskRes = distributeTask(data['plot'])(data)
  except Exception as e:
    msgPlot(traceback.format_exc(),data)
  if verbose:
    print("Final main memory %.2fG<br>"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))

def errorTask(data):
  msgPlot('Error plot task (unknown %s)!'%data['plot'],data)
def errorCheck(data):
  if data['plot'] in ["violin","dotplot","heatmap"]:
    if len(data["genes"])<1:
      msgPlot('Error: No matched gene!',data)
    if len(data["groups"])<1:
      msgPlot('Error: No matched annotation!',data)
  if data['plot']=="embedding":
    if len(data["reductions"])<1:
      msgPlot('Error: No matched embedding keys or genes or annotations)!',data)
    if len(data["genes"])<1 and len(data["groups"])<1:
      msgPlot('Error: No matched genes and annotations)!',data)
  if data['plot']=="stackbar":
    if len(data["groups"])<2:
      msgPlot('Error: At least two matched annotations are required!',data)
def msgPlot(msg,data):
  fig = plt.figure(figsize=(4,3))
  if msg.startswith("Traceback"):
    a= plt.text(0.01,0.5,msg,fontsize=4,horizontalalignment="left",verticalalignment="center")
  else:
    a= plt.text(0.5,0.5,textwrap.fill(msg,35),fontsize=14,horizontalalignment="center",verticalalignment="center")
  a= plt.axis("off")
  toHTML(fig,data)
  exit()
def distributeTask(aTask):
  return {
    'violin': complexViolin,
    'dotplot': twofactorDotplot,
    'embedding': reductionPlot,
    'stackbar':stackBar,
    'heatmap':complexHeatmap
  }.get(aTask,errorTask)
def is_numeric(var):
  try:
    float(var)  # Convert to float (handles integers and floats)
    return True
  except ValueError:
    return False
def isOptionDefined(data,k):
  return (data['options'].get(k) is not None and (is_numeric(data['options'][k]) or len(data['options'][k])>0))
def get_n_distinct_colors(n,lightness=0.5,saturation=0.9,cName=None):
  cmap = plt.get_cmap("Set3" if cName is None else cName)
  if n<len(cmap.colors):
    return([cmap(i) for i in range(n)])
  else:
    return distinctipy.get_colors(n,[cmap(i)[:3] for i in range(len(cmap.colors))])
  #return [colorsys.hls_to_rgb(i/n, lightness, saturation) for i in range(n)]
def toHTML(fig,data):
  st = time.time()
  imgD = iostreamFig(fig,data['options']['img_format'])
  imgID=""
  if len(data['options']['img_id'])>0:
    imgID='id="%s" '%data['options']['img_id']
  imgFormat = re.sub("svg","svg+xml",data['options']['img_format'])
  if data['options']['img_html']:
    if verbose:
      print("toHtml: %s<br>"%str(timedelta(seconds=time.time()-st)))
    print('<html><body><img %s src="data:image/%s;base64,%s" width="100%%" height="auto"/></body></html>'%(imgID,imgFormat,imgD))
  else:
    print('data:image/%s;base64,%s'%(imgFormat,imgD))
def iostreamFig(fig,img_format):
  figD = BytesIO()
  fig.savefig(figD,bbox_inches="tight",
    format=img_format)
  imgD = base64.encodebytes(figD.getvalue()).decode("utf-8")
  figD.close()
  if 'matplotlib' in str(type(fig)):
    plt.close(fig)#'all'
  return imgD
def getData(data,dataframe=True):
  st = time.time()
  errorCheck(data)
  D = ad.read_h5ad(data['h5ad'],backed='r')
  if verbose:
    print("Read: %s<br>"%str(timedelta(seconds=time.time()-st)))
    st = time.time()
  data['options']["img_format"] = data['options']["img_format"] if data['options'].get("img_format") in ['png','svg'] else "png"
  data['options']["img_width"]=6 if not isOptionDefined(data,"img_width") else data['options']['img_width']
  data['options']["img_height"]=4 if not isOptionDefined(data,"img_height") else data['options']['img_height']
  data['options']['cutoff']=0 if not isOptionDefined(data,"cutoff") else data['options']['cutoff']
  data['options']['titlefontsize']=6 if not isOptionDefined(data,"titlefontsize") else data['options']['titlefontsize']
  # checking existing genes/annotations/reduction keys
  if len(data['var_col'])>0 and data['var_col'] in D.var.columns:
    genes = {(D.var_names[D.var[data['var_col']]==k][0]):k for k in data['genes']}
  else:
    g = list(D.var_names[D.var_names.str.lower().isin([s.lower() for s in data['genes']])])
    genes = dict(zip(g,g))
  data["genes"] = list(genes.values())
  data["groups"] = {k:data["groups"][k] for k in data["groups"] if k in D.obs.columns}
  # only needs when plotting embedding
  reduc=[]
  if data['plot']=='embedding':
    reducName = []
    for one in data['reductions']:
      s=0.5
      selK=None
      for k in D.obsm.keys():
        if SequenceMatcher(None,one.lower(),k.lower()).ratio()>s:
          s=SequenceMatcher(None,one.lower(),k.lower()).ratio()
          selK=k
      if selK is not None and not selK in reducName:
        reducName += [selK]
        reduc += [(selK,0),(selK,1)]
    data['reductions'] = reducName
  errorCheck(data)
  if verbose:
    print("Init: %s<br>"%str(timedelta(seconds=time.time()-st)))
    st = time.time()
  #filter cells by annotation selections
  selC = [True] * D.shape[0]
  for one in data["groups"]:
    if len(data["groups"][one])>0:
      delGrp = [re.sub("^-","",_) for _ in data['groups'][one] if _.startswith('-')]
      if len(delGrp)>0:
        selC = selC & ~D.obs[one].isin(delGrp)
      else:
        selC = selC & D.obs[one].isin(data["groups"][one])
  if verbose:
    print("Filter: %s<br>"%str(timedelta(seconds=time.time()-st)))
    st = time.time()
  if dataframe:
    df = sc.get.obs_df(D,list(genes.keys())+list(data["groups"].keys()))[selC].rename(columns=genes)
    for k in data["groups"]:
      df[k] = df[k].astype(str).astype('category')
    if len(reduc)>0:
      df = df.merge(sc.get.obs_df(D,obsm_keys=reduc),how="left",left_index=True,right_index=True)
    if verbose:
      print("Get data: %s<br>"%str(timedelta(seconds=time.time()-st)))
    return df
  if verbose:
    print("Get data: %s<br>"%str(timedelta(seconds=time.time()-st)))
  return D,selC

def complexViolin(data):
  df = getData(data)
  st=time.time()
  w=data['options']["img_width"]
  h=data['options']["img_height"]
  genes=data['genes']
  grps=list(data['groups'].keys())
  gN = len(genes)
  fig, axes = plt.subplots(gN, 1, figsize=(w, h*gN), sharey='row')
  if gN==1:
    axes = [axes]
  for i in range(gN):
    subDF = df
    strTitle = "Total of %d cells" %df.shape[0]
    if data['options']['cutoff']>0:
      subDF = df[(df[genes[i]]>data['options']['cutoff']).values]
      strTitle="%d out of selected %d cells passed the expression filter %.2f"%(subDF.shape[0],df.shape[0],data['options']['cutoff'])
    if subDF.shape[0]<5:
        msgPlot("Less than 5 cells are satisfied with cutoff %.3f"%data['options']['cutoff'],data)
    sns.violinplot(x=grps[0],y=genes[i],ax=axes[i],
      data=subDF,cut=0,
      palette="bright" if not isOptionDefined(data,"palette") else data['options']['palette'],
      #fill=False,inner_kws={"alpha":0.5}, seaborn v0.13.0
      hue=None if len(grps)<2 else grps[1])
    if isOptionDefined(data,"dotsize"):
      dotColor='#000' if not isOptionDefined(data,"dotcolor") else data['options']['dotcolor']
      sns.stripplot(x=grps[0],y=genes[i],ax=axes[i],legend=False,
        data=subDF,size=data['options']["dotsize"],
        palette=[dotColor] if len(grps)<2 else [dotColor]*df[grps[1]].nunique(),
        dodge=False if len(grps)<2 else True,
        hue=None if len(grps)<2 else grps[1])
    axes[i].set_title(strTitle,loc="left",fontdict={'fontsize':data['options']['titlefontsize']})
    if i<(len(genes)-1):
      axes[i].get_xaxis().set_visible(False)
    else:
      plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    if len(grps)>1:
      if i==0:
        axes[i].legend(loc='lower right', bbox_to_anchor=(1, 1),
          ncol=1 if len(grps)<2 else df[grps[1]].nunique())
      else:
        axes[i].get_legend().remove()
  if verbose:
    print("complexViolin: %s<br>"%str(timedelta(seconds=time.time()-st)))
  return(toHTML(fig,data))
def twofactorDotplot(data):
  #D,selC = getData(data,False)
  #cellN = pd.DataFrame(selC).sum()[0]
  df = getData(data)
  st = time.time()
  cellN = df.shape[0]
  if verbose:
    print("Dotplot Backed peak memory %.2fG %s<br>"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
      str(timedelta(seconds=time.time()-st))))
    st = time.time()
  w=data['options']["img_width"]
  h=data['options']["img_height"]
  grps=list(data['groups'].keys())
  genes=data['genes']
  D=ad.AnnData(X=df[genes],obs=df[grps])

  strGrp=grps[0]
  if len(grps)>1:
    strGrp = "_".join(grps[:2])
    D.obs[strGrp] = (D.obs[grps[0]].astype(str)+"_"+ D.obs[grps[1]].astype(str)).astype('category') #D.obs.apply(lambda x: "_".join(x[grps[:2]]),axis=1)
  if verbose:
    print("Dotplot merge groups peak memory %.2fG %s<br>"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
      str(timedelta(seconds=time.time()-st))))
    st = time.time()
  strTitle = "%d selected cells"%cellN
  if data['options']['cutoff']>0:
    strTitle = "%d selected cells with expression cutoff %.2f"%(cellN,data['options']['cutoff'])
  dp=sc.pl.dotplot(D,genes,groupby=strGrp,figsize=(w,h),
    expression_cutoff=data['options']['cutoff'],mean_only_expressed=True,
    return_fig=True)
  if verbose:
    print("Dotplot Plot1 peak memory %.2fG %s<br>"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
      str(timedelta(seconds=time.time()-st))))
    st = time.time()
  dp = (dp.add_totals(size=1.2).
    legend(show_size_legend=True). #,width=float(data['legendW'])
    style(cmap="Reds" if not isOptionDefined(data,'color_map') else data['options']['color_map'],
      dot_edge_color='black', dot_edge_lw=0.5, size_exponent=1.5))
  if verbose:
    print("Dotplot Plot2 peak memory %.2fG %s<br>"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
      str(timedelta(seconds=time.time()-st))))
    st = time.time()
  fig = dp.show(True)['mainplot_ax'].figure
  if len(grps)>1:
    #n = df[selC][grps[1]].nunique()
    n = df[grps[1]].nunique()
    for i in range(df[grps[0]].nunique()):#[selC]
      if i==0:
        fig.axes[0].set_title(strTitle,loc="left",fontdict={'fontsize':data['options']['titlefontsize']})
      else:
        fig.axes[0].axhline(y=i*n,color="#0002",linestyle="--")
  if verbose:
    print("Dotplot Add line peak memory %.2fG %s <br>"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
      str(timedelta(seconds=time.time()-st))))
  return(toHTML(fig,data))
def reductionPlot(data):
  df = getData(data)
  st=time.time()
  w=data['options']["img_width"]
  h=data['options']["img_height"]
  grps=list(data['groups'].keys())
  genes=data['genes']
  obsm={}
  for one in data['reductions']:
    obsm[one]=df[["%s-0"%one,"%s-1"%one]].to_numpy()
  dotsize=120000/df.shape[0]
  subSize = 4
  groupN = len(grps)
  geneN = len(genes)
  ncol = 4 if groupN<2 else df[grps[1]].nunique()
  nrow = (groupN + geneN) if groupN>1 else (groupN+math.ceil(geneN/ncol))
  fig = plt.figure(figsize=(ncol*subSize,subSize*nrow))
  gs = fig.add_gridspec(nrow,ncol,wspace=0.2)
  oneReduc = re.sub("^X_","",data['reductions'][0])
  D=ad.AnnData(X=None if geneN==0 else df[genes],obs=None if groupN==0 else df[grps],obsm=obsm)
  for i in range(groupN):
    ix = groupN-i-1
    ax = sc.pl.embedding(D,oneReduc,color=grps[ix],ax=fig.add_subplot(gs[i,0]),
      show=False,
      palette=None if not isOptionDefined(data,'palette') else data['options']["palette"])
    ax.legend(ncol=math.ceil(df[grps[ix]].nunique()/10),loc=6,bbox_to_anchor=(1,0.5),
      frameon=False,fontsize=8-df[grps[ix]].nunique()/20)
    ax.set_xlabel('%s 1'%oneReduc)
    ax.set_ylabel('%s 2'%oneReduc)
  if groupN<2:
    for i in range(geneN):
      x = int(i/ncol)+groupN
      y = i % ncol
      ax = sc.pl.embedding(D,oneReduc,color=genes[i],ax=fig.add_subplot(gs[x,y]),show=False,size=dotsize)
      ax.set_xlabel('%s 1'%oneReduc)
      ax.set_ylabel('%s 2'%oneReduc)
  else:
    splitNames = list(df[grps[1]].unique())
    for i in range(geneN):
      for j in range(len(splitNames)):
        x = groupN + i
        y = j
        ax = sc.pl.embedding(D,oneReduc,ax=fig.add_subplot(gs[x,y]),show=False,size=dotsize)
        ax = sc.pl.embedding(D[D.obs[grps[1]]==splitNames[j]],oneReduc,color=genes[i],
          color_map="viridis" if not isOptionDefined(data,'color_map') else data['options']["color_map"],
          vmin=df[genes[i]].min(),vmax=df[genes[i]].max(),ax=ax,show=False,
          size=dotsize,
          title='{} in {}'.format(genes[i],splitNames[j]))
        ax.set_xlabel('%s 1'%oneReduc)
        ax.set_ylabel('%s 2'%oneReduc)
  fig.suptitle("%d selected cells"%df.shape[0],x=0.9,y=0.9,ha="right",va="top",
    fontsize=data['options']['titlefontsize'])
  if verbose:
    print("reductionPlot: %s<br>"%str(timedelta(seconds=time.time()-st)))
  return(toHTML(fig,data))
def stackBar(data):
  df = getData(data)
  st=time.time()
  strTitle = "%d selected cells"%df.shape[0]
  w=data['options']["img_width"]
  h=data['options']["img_height"]
  grps=list(data['groups'].keys())
  x = list(df[grps[1]].unique()) if len(data["groups"][grps[1]])==0 else data["groups"][grps[1]]
  df = (df[grps[:2]].value_counts().to_frame("count").reset_index().
    pivot_table(index=grps[0],columns=grps[1],values="count"))
  fig = plt.figure(figsize=(w,h))
  if data["options"].get("yscale") is not None and data["options"]["yscale"]=="proportion":
    df = df.apply(lambda x: x/x.sum())
    plt.ylabel("Proportion")
  else:
    plt.ylabel("Count")
  plt.xlabel(grps[1])
  color=get_n_distinct_colors(df.shape[0],cName=data['options']["palette"] if isOptionDefined(data,"palette") else None)
  for i in range(df.shape[0]):
    plt.bar(x,df.iloc[i,:][x],color=color[i],
      bottom=df.iloc[:i,:][x].sum())
  plt.legend(df.index,loc=4,bbox_to_anchor=(1,1),
    ncol=math.ceil(df.shape[0]/10),
    fontsize=max(2,6-df.shape[0]/20))
  fig.axes[0].set_title(strTitle,
    loc="left",fontdict={'fontsize':data['options']['titlefontsize']})
  if verbose:
    print("stackBar: %s<br>"%str(timedelta(seconds=time.time()-st)))
  return(toHTML(plt,data))
def complexHeatmap(data):
  df = getData(data)
  st=time.time()
  w=data['options']["img_width"]
  h=data['options']["img_height"]
  grps=list(data['groups'].keys())
  genes = data["genes"]
  selN = df.shape[0]
  df = df[df[genes].apply(lambda x: max(x)>data["options"]["cutoff"],axis=1)]
  if df.shape[0]<5:
    msgPlot("Less than 5 cells are satisfied with cutoff %.3f"%data['options']['cutoff'],data)
  heat_scale=None
  heat_title="Expression"
  if data["options"].get("heat_scale") is not None and data["options"]["heat_scale"]=="z-score":
    heat_scale=0
    heat_title="Row Z-score"
  if data["options"].get("cell_order") is None or data["options"]["cell_order"]=="groups":
    df = df.sort_values(list(data["groups"].keys()))
  elif data["options"]["cell_order"]=="expression":
    ix = hierarchy.leaves_list(fc.linkage_vector(df[genes],method="ward"))
    df = df.iloc[ix,]
  colors={_:dict(zip(df[_].unique(),get_n_distinct_colors(df[_].nunique(),
    cName=data["options"]["palette"] if isOptionDefined(data,"palette") else None))) for _ in grps}
  fig = plt.figure(figsize=(w, h))
  df.to_csv("test.csv")
  with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
    left_anno = pch.HeatmapAnnotation(df[grps],colors=colors,axis=0)
    left_anno.plot_annotations()
    plt.close()
    plt.rc('legend',fontsize=8 if not isOptionDefined(data,"heat_legend_fontsize") else data['options']['heat_legend_fontsize'])
    cm = pch.ClusterMapPlotter(
        data=df[genes],z_score=heat_scale,
        label=heat_title,cmap="jet" if not isOptionDefined(data,"color_map") else data['options']["color_map"],
        left_annotation=left_anno,
        show_rownames=False,show_colnames=True,
        row_dendrogram=False,col_dendrogram=False,
        col_cluster=False,row_cluster=False,
        #row_cluster_method="complete",col_cluster_method="complete",
        rasterized=True,legend=True,legend_anchor='ax_heatmap',
        verbose=0)
  #print(len(fig.axes))
  fig.axes[1].set_title("%d of %d selected cells passed expression threshold %.2f"%(df.shape[0],selN,data["options"]["cutoff"]),
    loc="left",fontdict={'fontsize':data['options']['titlefontsize']})
  if verbose:
    print("complexHeatmap: %s<br>"%str(timedelta(seconds=time.time()-st)))
  return(toHTML(plt,data))
main()
# cat ../testVIP/violin.json | python -u plotH5ad.py
# python -u ./plotH5ad.py ../testVIP/violin.json


