import os,resource,sys
import anndata as ad
from scipy.sparse import csc_matrix

def main():
  if len(sys.argv)>2:
    slimH5ad(sys.argv[1],sys.argv[2])
    print("Final peak memory %.2fG"%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
  else:
    raise Exception("Missing input file and/or path to the ouput file")
  
def slimH5ad(strIn,strOut):
  if not os.path.isfile(strIn):
    raise Exception("Error: The input h5ad does NOT exist: %s!"%strIn)
  if os.path.isfile(strOut):
    print("Warning: The output h5ad exists and will be overwritten: %s!"%strOut)
  print("Loading %s"%strIn)
  D = ad.read_h5ad(strIn)
  print("Slimming ...")
  D.uns = dict()
  D.raw = None
  for one in [_ for _ in D.obsm.keys() if 'spatial' in _]:
    del D.obsm[one]
  D.X = csc_matrix(D.X)
  print("Saving %s"%strOut)
  D.write(strOut)

main()
