'''
Prepare scripts from notebooks/example_notebook.ipynb
'''
import pickle
import json 
import scanpy as sc
import scipy 
import numpy as np
import os
import argparse
from pprint import pprint

# get boolean argument use_cell_types, use_highly_variable, use_weights, lam, kappa, rho, delta
parser = argparse.ArgumentParser()
parser.add_argument('--use_cell_types', '-u', type=str, required=True)
parser.add_argument('--use_highly_variable', '-v', type=str, default=True)
parser.add_argument('--use_weights', '-w', type=str, default=True)
parser.add_argument('--lam', '-l', type=float, default=0.1)
parser.add_argument('--kappa', '-k', type=float, default=0.00001)
parser.add_argument('--rho', '-r', type=float, default=0.00001)
parser.add_argument('--delta', '-d', type=float, default=0.001)
args = parser.parse_args()

args.use_cell_types = args.use_cell_types == 'True'
args.use_highly_variable = args.use_highly_variable == 'True'
args.use_weights = args.use_weights == 'True'

print("Arguments:")
for arg in vars(args):
    print('  - ', arg, getattr(args, arg))

if args.use_cell_types:
    cell_type_key = None
else:
    cell_type_key = 'cell_type_annotations'

def get_factor_celltypes(adata, obs_key, cellscore_obsm_key = 'SPECTRA_cell_scores'):
    import pandas as pd
    cell_scores_df = pd.DataFrame(adata.obsm[cellscore_obsm_key])
    cell_scores_df['celltype'] = list(adata.obs[obs_key])
    
    global_factors_series = (cell_scores_df.groupby('celltype').mean() !=0).all()
    global_factors = [factor for factor in global_factors_series.index if global_factors_series[factor]]
    specific_cell_scores = (cell_scores_df.groupby('celltype').mean()).T[~global_factors_series].T
    specific_factors = {}
    
    for i in set(cell_scores_df['celltype']):
        specific_factors[i]=[factor for factor in specific_cell_scores.loc[i].index if specific_cell_scores.loc[i,factor]]
    
    factors_inv = {}
    for i,v in specific_factors.items():
        for factor in v:
            factors_inv[factor]=i
    
    for factor in global_factors:
        factors_inv[factor]= 'global'
            
    return factors_inv

def check_gene_set_dictionary(adata, annotations, obs_key='cell_type_annotations',global_key='global', return_dict = True):
    adata_labels  = list(set(adata.obs[obs_key]))+['global']#cell type labels in adata object
    annotation_labels = list(annotations.keys())
    if set(annotation_labels)==set(adata_labels):
        print('Cell type labels in gene set annotation dictionary and AnnData object are identical')
        dict_keys_OK = True
    if len(annotation_labels)<len(adata_labels):
        print('The following labels are missing in the gene set annotation dictionary:',set(adata_labels)-set(annotation_labels))
        dict_keys_OK = False
    if len(adata_labels)<len(annotation_labels):
        print('The following labels are missing in the AnnData object:',set(annotation_labels)-set(adata_labels))
        dict_keys_OK = False
        
    Counter = 0
    annotations_new = {}
    for k,v in annotations.items():
        annotations_new[k] = {}
        for k2,v2 in v.items():
            annotations_new[k][k2]= [x for x in v2 if x in adata.var_names]
            length = len(v2)
            if length<3:
                print('gene set',k2,'for cell type',k,'is of length',length)
                Counter = Counter+1
            
    if Counter > 0:
        print(Counter,'gene sets are too small. Gene sets must contain at least 3 genes')
    elif Counter == 0 and dict_keys_OK:
        print('Your gene set annotation dictionary is correctly formatted.')
    if return_dict:
        return annotations_new


print('Read annotations from json file')
dict_path = os.getcwd() + '/data/Spectra_dict.json'
with open(dict_path, "rb") as file:
    annotations = json.load(file)

print('Read data from h5ad file')
adata_path = os.getcwd() + '/data/sample_data.h5ad'
adata = sc.read(adata_path)

print('Check gene set annotation dictionary')
annotations = check_gene_set_dictionary(adata, annotations, obs_key='cell_type_annotations',global_key='global')

# store adata and annotations in data/temp
if not os.path.exists(os.getcwd() + '/data/temp'):
    os.makedirs(os.getcwd() + '/data/temp')

adata.write(os.getcwd() + '/data/temp/adata.h5ad')
with open(os.getcwd() + '/data/temp/Spectra_dict.json', 'w') as fp:
    json.dump(annotations, fp)

print('adata.h5ad is stored in:', os.getcwd() + '/data/temp/adata.h5ad')
print('annotations.json is stored in:', os.getcwd() + '/data/temp/Spectra_dict.json')

######### from est_spectra() function - prepare for faster iterations #########
print('Prepare model')
# prepare model


if args.use_cell_types:
    gene_set_dictionary = annotations
    L = {}
    for key in gene_set_dictionary.keys(): 
        length = len(list(gene_set_dictionary[key].values()))
        L[key] = length + 1 
else:
    gene_set_dictionary = annotations['global']
    length = len(list(gene_set_dictionary.values()))
    L = length 

#create vocab list from gene_set_dictionary
lst = []
if args.use_cell_types:
    for key in gene_set_dictionary:
        for key2 in gene_set_dictionary[key]:
            gene_list = gene_set_dictionary[key][key2] 
            lst += gene_list
else:
    for key in gene_set_dictionary:
        gene_list = gene_set_dictionary[key]
        lst += gene_list

#lst contains all of the genes that are in the gene sets --> convert to boolean array 
bools = [] 
for gene in adata.var_names:
    if gene in lst:
        bools.append(True)
    else: 
        bools.append(False)
bools = np.array(bools)

if args.use_highly_variable:
    idx_to_use = bools | adata.var.highly_variable #take intersection of highly variable and gene set genes (todo: add option to change this at some point)
    X = adata.X[:,idx_to_use] 
    vocab = adata.var_names[idx_to_use]
    adata.var["spectra_vocab"] = idx_to_use
else: 
    X = adata.X
    vocab = adata.var_names 

if cell_type_key is not None:
    labels = adata.obs[cell_type_key].values
else:
    labels = None 
if type(X) == scipy.sparse.csr.csr_matrix:
    X = np.array(X.todense())
word2id = dict((v, idx) for idx, v in enumerate(vocab))
id2word = dict((idx, v) for idx, v in enumerate(vocab))

params = dict(
    X = X, 
    labels = labels,
    L = L, 
    vocab = vocab, 
    gs_dict = gene_set_dictionary, 
    use_weights = args.use_weights,
    lam = args.lam, 
    delta=args.delta,
    kappa = args.kappa, 
    rho = args.rho, 
    use_cell_types = args.use_cell_types
)
# save params as pickle
if args.use_cell_types:
    pickle_path = os.getcwd() + '/data/temp/params_cell_types.pkl'
else:
    pickle_path = os.getcwd() + '/data/temp/params.pkl'

with open(pickle_path, 'wb') as f:
    pickle.dump(params, f)

print('params.pkl is stored in:', pickle_path)