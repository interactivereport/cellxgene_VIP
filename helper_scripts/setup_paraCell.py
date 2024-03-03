import scanpy as sc

#Load AnnData file

file_name = '' # e.g. Malaria_Data.h5ad
adata = sc.read_h5ad(file_name) 

#Format Dictionary
#Guidance on this can be found here: https://github.com/sii-cell-atlas/paraCell/wiki/Setting-Up-paraCell

paraCell_setup = {

"Database": "",

"Initial Embedding":"",

"Initial Coloring":"",

"includeGeneSearch":"",

"include_goSearch":"",

"isHP": "",

"includePseudo": "",

"pseudoEmbed": "",

"includeTradeSeq": "",

"tsAnnot": "",

"includeDescription":"",

"Description":"",

"Original_Paper":"",

"Author(s)":""

}

#Save Dictionary to the uns slot of the AnnData Object

adata.uns["paraCell_setup"] = paraCell_setup

#Save updated AnnData Object

adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

adata.write_h5ad(file_name)

