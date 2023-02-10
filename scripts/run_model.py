import os
import json
import pickle
import scanpy as sc
from spectra import spectra as spc

# Get argument for number of epochs, use_cell_types
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--epochs', '-e', type=int, default=1)
parser.add_argument('--use_cell_types', '-u', type=str, required=True)
args = parser.parse_args()

args.use_cell_types = args.use_cell_types == 'True'

print("Arguments:")
for arg in vars(args):
    print('  - ', arg, getattr(args, arg))

print("Loading data...")
if args.use_cell_types:
    with open(os.getcwd() + '/data/temp/params_cell_types.pkl', 'rb') as f:
        params = pickle.load(f)
else:
    with open(os.getcwd() + '/data/temp/params.pkl', 'rb') as f:
        params = pickle.load(f)


print("Running Spectra...")
print("labels:", params['labels'])
params.pop("cell_type_key")
spectra = spc.SPECTRA_Model(**params)

print("Training model...")
params_train = {k: v for k, v in params.items() if k in ['X', 'labels', 'lr_schedule', 'num_epochs']}
params_train['num_epochs'] = args.epochs
spectra.train(**params_train)
