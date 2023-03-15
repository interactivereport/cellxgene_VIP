import scanpy as sc
import pandas as pd

cowthei = sc.read_h5ad("cowthei_example.h5ad")

#Read in Host-Parasite Gene Lists

parasite_genes = pd.read_csv("parasite_genes.csv")
parasite_genes = parasite_genes["x"]

host_genes = pd.read_csv("host_genes.csv")
host_genes = host_genes["x"]

#Add Gene Lists to Object

cowthei.uns["parasite_genes"] = parasite_genes.values
cowthei.uns["host_genes"] = host_genes.values

#Add Features metadata

cowthei.var["features"] = cowthei.var_names.values

# Save AnnData object

cowthei.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

cowthei.write_h5ad("cowthei_test.h5ad")