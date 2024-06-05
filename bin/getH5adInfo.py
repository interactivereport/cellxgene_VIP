import sys, json, math, resource, itertools, time, h5py
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import vstack
from scipy.sparse import csc_matrix
from scipy.sparse import hstack

def baseInfo(strH5ad):
  D = ad.read_h5ad(strH5ad,backed='r')
  ## check CSC
  CSC=False
  if hasattr(D.X,'format_str'):
    if D.X.format_str=="csc": #format_str when backed="r"
      CSC=True
  obs = D.obs.select_dtypes('category')
  obs = obs[obs.columns[obs.nunique()<500]]
  # check the expressed genes
  info = {'version':"2024.05",
    'cellN':D.shape[0],'geneN':D.shape[1],
    'layout':[x.replace("X_","") for x in D.obsm.keys()],
    'csc':CSC,
    'annotation':obs.apply(lambda x: x.value_counts().to_dict()).to_dict()}
  
  print(json.dumps(info))
def geneMax(strH5ad):
  #st = time.time()
  D = ad.read_h5ad(strH5ad,backed='r')
  nC,nG = D.shape
  gName = np.array(D.var.index)
  stepMax=[]
  with h5py.File(strH5ad,'r') as f:
    X = f['X']
    maxM = 1e9
    stepN = math.ceil(nC*nG/maxM)
    if hasattr(X,'keys') and "indptr" in X.keys() and X['indptr'].shape[0]==(nC+1): #csr
      # read all genes across each batch of cells
      k,m = divmod(nC,stepN)
      steps = [0]+[(i+1)*k+min(i+1, m) for i in range(stepN)]
      indptr=X['indptr'][()]
      for i in range(stepN):
        iStart = indptr[steps[i]]
        iEnd = indptr[steps[i+1]]
        subX=csr_matrix((X['data'][iStart:iEnd],X['indices'][iStart:iEnd],indptr[steps[i]:(1+steps[i+1])]-indptr[steps[i]]),shape=(steps[i+1]-steps[i],nG))
        stepMax += [subX.max(0)]
      stepMax = vstack(stepMax).max(0).toarray()[0]
    elif hasattr(X,'keys') and "indptr" in X.keys() and X['indptr'].shape[0]==(nG+1): #csc
      k,m = divmod(nG,stepN)
      steps = [0]+[(i+1)*k+min(i+1, m) for i in range(stepN)]
      indptr=X['indptr'][()]
      stepMax=[]
      for i in range(stepN):
        iStart = indptr[steps[i]]
        iEnd = indptr[steps[i+1]]
        subX=csc_matrix((X['data'][iStart:iEnd],X['indices'][iStart:iEnd],indptr[steps[i]:(1+steps[i+1])]-indptr[steps[i]]),shape=(nC,steps[i+1]-steps[i]))
        stepMax += [subX.max(0)]
      stepMax = hstack(stepMax).toarray()[0]
    elif hasattr(X,'shape') and X.shape[1]==nG:
      stepMax=X[()].max(axis=0)
    else:
      print("Error: only support X in dense and CSR or CSC sparse matrix")
      return
  genes = dict(zip(gName.astype(str),np.round(stepMax,2).astype("float64")))
  print(json.dumps(genes))
  #print("\tTotal %d cells\tPeak memory %.2fG with %.1f seconds"%(nC,resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,time.time()-st))
def verify(strH5ad):
  info={}
  D = None
  try:
    D = ad.read_h5ad(strH5ad,backed='r')
    info["format"]={"Name":"h5ad Format","Pass":"True","Message":""}
  except:
    info["format"]={"Name":"h5ad Format","Pass":"False","Message":"Cannot be opened by AnnData v%s"%ad.__version__}
  if D is not None:
    if len(D.obsm)==0:
      info["obsm"]={"Name":"embedding","Pass":"False","Message":"There is no embedding in obsm"}
    else:
      info["obsm"]={"Name":"embedding","Pass":"True","Message":""}
    if len([_ for _ in D.obsm.keys() if _.startswith("X_")])==0:
      info["obsm"]={"Name":"embedding names","Pass":"False","Message":"There is no embedding name starts with 'X_'"}
    else:
      info["obsm"]={"Name":"embedding names","Pass":"True","Message":""}
  print(json.dumps(info))
def main():
  if len(sys.argv)<3:
    print("ERROR: path to a h5ad file is required!")
    exit()
  if sys.argv[1].lower()=="basic":
    baseInfo(sys.argv[2])
  elif sys.argv[1].lower()=="gene":
    geneMax(sys.argv[2])
  elif sys.argv[1].lower()=="verify":
    verify(sys.argv[2])
  else:
    print("Error: unknown task %s"%sys.argv[1])
if __name__ == "__main__":
  main()
