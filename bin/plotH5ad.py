import sys,json,re,time
import pandas as pd
import seaborn as sns
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import base64
from io import BytesIO

def main():
  if len(sys.argv)==1:
    data = json.load(sys.stdin)
  else:
    with open(sys.argv[1],'r') as f:
      data = json.load(f)
  taskRes = distributeTask(data['plot'])(data)

def errorTask(data):
  raise ValueError('Error task!')
def distributeTask(aTask):
  return {
    'violin': complexViolin,
  }.get(aTask,errorTask)
def toHTML(fig,data):
  st = time.time()
  imgD = iostreamFig(fig,data['options']['img_format'])
  imgID=""
  if len(data['options']['img_id'])>0:
    imgID='id="%s" '%data['options']['img_id']
  imgFormat = re.sub("svg","svg+xml",data['options']['img_format'])
  if data['options']['img_html']:
    print("toHtml: %.2f"%(time.time()-st))
    print('<html><body><img %s src="data:image/%s;base64,%s" width="100%%" height="auto"/></body></html>'%(imgID,imgFormat,imgD))
  else:
    print('<img %s src="data:image/%s;base64,%s" width="100%%" height="auto"/>'%(imgID,imgFormat,imgD))
def iostreamFig(fig,img_format):
  figD = BytesIO()
  fig.savefig(figD,bbox_inches="tight",
    format=img_format)
  imgD = base64.encodebytes(figD.getvalue()).decode("utf-8")
  figD.close()
  if 'matplotlib' in str(type(fig)):
    plt.close(fig)#'all'
  return imgD
def getData(data):
  D = ad.read_h5ad(data['h5ad'],backed='r')
  if len(data['var_col'])>0 and data['var_col'] in D.var.columns:
    D.var_names = list(D.var[data['var_col']])
  data["genes"] = list(D.var_names[D.var_names.str.lower().isin([s.lower() for s in data['genes']])])
  data['options']["img_format"] = data['options']["img_format"] if data['options'].get("img_format") in ['png','svg'] else "png"
  data['options']["img_width"]=6 if data['options'].get('img_width') is None else data['options'].get('img_width')
  data['options']["img_height"]=4 if data['options'].get('img_height') is None else data['options'].get('img_height')
  return D
def complexViolin(data):
  st=time.time()
  recordT = {}
  w=data['options']["img_width"]
  h=data['options']["img_height"]
  genes=data['genes']
  grps=list(data['groups'].keys())
  gN = len(genes)
  if gN<1:
    raise ValueError('Missing genes!')
  recordT["init"]=time.time()-st
  D = getData(data)
  recordT["Get obs"]=time.time()-st
  df = sc.get.obs_df(D,genes+grps)
  recordT["Create data"]=time.time()-st
  for one in grps:
    if len(data["groups"][one])>0:
      delGrp = [re.sub("^-","",_) for _ in data['groups'][one] if _.startswith('-')]
      if len(delGrp)>0:
        df = df[~df[one].isin(data["groups"][one])]
      else:
        df = df[df[one].isin(data["groups"][one])]
      df[one] = df[one].cat.remove_unused_categories()
  recordT["Filter"]=time.time()-st
  
  fig, axes = plt.subplots(gN, 1, figsize=(w, h*gN), sharey='row')
  if gN==1:
    axes = [axes]
  for i in range(gN):
    subDF = df
    if data['options'].get("cutoff") is not None and data['options']['cutoff']>0:
      subDF = df[(df[genes[i]]>data['options']['cutoff']).values]
    sns.violinplot(x=grps[0],y=genes[i],ax=axes[i],
      data=subDF,cut=0,
      hue=None if len(grps)<2 else grps[1])
    if data['options'].get("dotsize") is not None and data['options']["dotsize"]>0:
      sns.stripplot(x=grps[0],y=genes[i],ax=axes[i],
        palette=['#000'],legend=False,
        data=subDF,size=data['options']["dotsize"],
        dodge=False if len(grps)<2 else True,
        hue=None if len(grps)<2 else grps[1])
    axes[i].set_title("%d cells"%df.shape[0],loc="left")
    if i<(len(genes)-1):
      axes[i].get_xaxis().set_visible(False)
    else:
      plt.setp(axes[i].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    if len(grps)>1:
      if i==0:
        axes[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=1 if len(grps)<2 else df[grps[1]].nunique())
      else:
        axes[i].get_legend().remove()
  recordT["Plot"]=time.time()-st
  if data['options']['img_html']:
    print(pd.DataFrame(recordT,index=["Time"]).transpose())
  #plt.savefig('f.pdf',bbox_inches="tight")
  return(toHTML(fig,data))

main()
# cat ../testVIP/violin.json | python -u plotH5ad.py
# python -u ./plotH5ad.py ../testVIP/violin.json
