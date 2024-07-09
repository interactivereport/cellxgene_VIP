import scanpy as sc
import pandas as pd
import glob 

file_name = 'name_of_annData_File' # e.g. Malaria_Data.h5ad

adata = sc.read_h5ad(file_name) 

LineageCounter = 1

for x in glob.glob('Lineage*.csv'):
    lineage = pd.read_csv(x)
    lin_title = "Lineage_" + str(LineageCounter)
    LineageCounter += 1

    lineage_dictionary = {
        'dim1': lineage[lineage.columns[0]].tolist(), 
        'dim2': lineage[lineage.columns[1]].tolist() 
        }

    adata.uns[lin_title] = lineage_dictionary

adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

adata.write_h5ad(file_name)
