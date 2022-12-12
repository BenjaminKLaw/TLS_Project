'''
Author: BL
Last_date_modified: 2022.08.02
Purpose: To run the AM-DNA-097 (TLS_120h_Rep1) dataset through Cassiopeia preprocessing and tree reconstruction.

Returns:
    Creates itol plots 
'''

# Plot Tree
#######################################################################

import cassiopeia as cas
import pandas as pd
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt
import seaborn as sns
from ete3 import Tree
from typing import Tuple
import scipy

name = 'AM-DNA-097'
output_dir = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/2_lineage_reconstruction/'
tree_loc = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/2_lineage_reconstruction/AM-DNA-097_hybrid_newick_noMutationlessEdges.txt'
meta_loc = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/2_lineage_reconstruction/AM-DNA-097_metadata.txt'
character_matrix_loc = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/2_lineage_reconstruction/AM-DNA-097_character_matrix.txt'
priors_loc = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/2_lineage_reconstruction/AM-DNA-097_priors.pickle'
clusterColorsFile = "/Genomics/chanlab/mchan/Adriano/TLS/TLS_TLSCL/20211102_clusterColorsTLSCL.p"

# Plot in itol
##################################################################

# Load in tree and tree metadata
tree = Tree(tree_loc, format=1)
tree_meta = pd.read_csv(meta_loc, sep='\t')
character_matrix = pd.read_csv(character_matrix_loc, sep='\t', index_col = 0)

with open(priors_loc, 'rb') as f:
    priors = pickle.load(f)
        
with open(clusterColorsFile,'rb') as fp:
    colorDict = pickle.load(fp)

# Create a cassiopeia Tree        
test_tree = cas.data.CassiopeiaTree(character_matrix = character_matrix, priors = priors, tree = tree)
    
missing_proportion = (character_matrix == -1).sum(axis=0) / character_matrix.shape[0]
uncut_proportion = (character_matrix == 0).sum(axis=0) / character_matrix.shape[0]
n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis=0)

character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 'uncut_prop', 'n_unique_states']).T
test_tree.cell_meta = tree_meta
test_tree.character_meta = character_meta

test_tree.cell_meta['cluster'] = test_tree.cell_meta['cell_state'].str[2:-2]
test_tree.cell_meta.set_index("cellBC",inplace = True)

color_order = ['#FFFF00', '#0000cd',
             '#30D5C8', '#228B22',
             '#023020','#ff0000',
             '#ba55d3', '#800080',
             '#9400d3', '#FF00FF',
             '#FFC0CB', '#808080',
             '#483d8b', '#6a5acd']

cas.pl.upload_and_export_itol(test_tree, 'TLS_120h_1_noMutationLess_final', api_key = 'FH5LLzi6mypVzMx1oYxyoA',
                              project_name = 'AM-DNA-097', export_filepath = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/3_tree_graphs/TLS_120h_1_noMutationLess_final.pdf', 
                              meta_data = ['cluster'], palette = color_order)

cas.pl.upload_and_export_itol(test_tree, 'TLS_120h_1_noMutationLess_indels_final', api_key = 'FH5LLzi6mypVzMx1oYxyoA',
                              project_name = 'AM-DNA-097', export_filepath = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/3_lineage_stats/TLS_120h_1_noMutationLess_indels_final.pdf', 
                              meta_data = ['cluster'], palette = color_order, allele_table = allele_table)
