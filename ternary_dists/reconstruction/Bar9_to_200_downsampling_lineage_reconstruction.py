'''
Author: BL
Last_date_modified: 2023.05.23
Purpose: To run the TLS1 subsampled trees that match the TLS M barcode sizes. Run on gen-comp2

Returns:
    Lineage Reconstruction - A filtered lineage table, a newick file containing a tree created by Cassiopeia's Hybrid Solver, tree metadata, and logfiles from the ILP portion of the solver
'''
import sys
import os       
import pandas as pd
import numpy as np
import cassiopeia as cas
import gurobipy
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

# Lineage Reconstruction
#######################################################################
cell_state_loc = '/Genomics/chanlab/blaw/TLS/LineageTracer/scRNAseq/TLS_TLSCL_cellBC_cellState.tsv'
best_seed = 222
n_threads = 10

# Remove any rows with NaN alleles, since they disrupt Cassiopeia
def remove_NaN(table):
    '''
    input:
        table - an allele_table that contains an r1, r2, and r3 cute site column
    output:
        returns a new table with any rows that contain NaN removed
    '''
    table = table[table['r1'].isnull() == False]
    table = table[table['r2'].isnull() == False]
    table = table[table['r3'].isnull() == False]
    
    return table

def hybrid_tree_pipeline(table, name, total_priors, output_loc, best_seed = 222):
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(table, allele_rep_thresh = 1, mutation_priors = total_priors)  

    character_matrix.to_csv(output_loc + '{}_downsampling_character_matrix.txt'.format(name), sep='\t')

    with open(output_loc + '{}_downsampling_priors.pickle'.format(name), 'wb') as handle:
        pickle.dump(priors, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open(output_loc + '{}_downsampling_state_2_indel.pickle'.format(name), 'wb') as handle:
        pickle.dump(state_2_indel, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    
    cas_tree = cas.data.CassiopeiaTree(character_matrix = character_matrix, priors = priors)
    
    missing_proportion = (character_matrix == -1).sum(axis = 0) / character_matrix.shape[0]
    uncut_proportion = (character_matrix == 0).sum(axis = 0) / character_matrix.shape[0]
    n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis = 0)

    cell_meta = table.groupby('cellBC').agg({"intBC": 'nunique', 'UMI': 'sum', 'cell_state': 'unique'})
    character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 'uncut_prop', 'n_unique_states']).T

    cas_tree.cell_meta = cell_meta
    cas_tree.character_meta = character_meta
    cas_tree.cell_meta.to_csv(output_loc + '{}_downsampling_metadata.txt'.format(name), sep = '\t')
    
    cas_tree.cell_meta['cell_state'] = cas_tree.cell_meta['cell_state'].astype('str')
    
    top = cas.solver.VanillaGreedySolver()
    bottom = cas.solver.ILPSolver(convergence_time_limit=3600, maximum_potential_graph_layer_size=1000, weighted=True, seed=best_seed)
    
    hybrid_solver = cas.solver.HybridSolver(top_solver=top, bottom_solver=bottom, cell_cutoff=200, threads=10)
    hybrid_solver.solve(cas_tree, logfile = output_loc + 'logs/{}_hybrid.log'.format(name))
    
    cas_tree.collapse_mutationless_edges(infer_ancestral_characters = True)
    
    noMutationlessEdgesFile = open(output_loc + '{}_downsampling_hybrid_newick_noMutationlessEdges.txt'.format(name), 'wt')
    noMutationlessEdgesFile.write(cas_tree.get_newick())
    noMutationlessEdgesFile.close()

### Calculate the total indels using allele_tables from TLS1, TLS2, and TLSCL combined
TLSCL_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-258/2_add_metadata/allele_table_multiSeq+lenti.txt'
TLS1_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-097/1_preprocessing/allele_table.txt'
TLS2_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-098/1_preprocessing/allele_table.txt'

TLSCL_allele = pd.read_csv(TLSCL_filepath, sep='\t')
TLS1_allele = pd.read_csv(TLS1_filepath, sep='\t')
TLS2_allele = pd.read_csv(TLS2_filepath, sep='\t')

# Combine all three datasets
TLSCL_allele = TLSCL_allele[(TLSCL_allele['finalCalls'] != 'Doublet') & (TLSCL_allele['finalCalls'] != 'Negative')]
TLS1_allele['finalCalls'] = 'Bar25'
TLS2_allele['finalCalls'] = 'Bar26'
total_table = pd.concat([TLSCL_allele, TLS1_allele, TLS2_allele])
total_table = remove_NaN(total_table)

# Calculate the priors using each intBC / Multiseq Barcode as a different indel instance
total_priors = cas.pp.compute_empirical_indel_priors(total_table, grouping_variables=['intBC', 'finalCalls'])

barcode = 'Bar9'

for i in range(30):
    name = '{}_{}_to_200'.format(i, barcode)
    temp_allele = pd.read_csv('/Genomics/chanlab/blaw/TLS/sandbox/ternary_studies/downsampling_to_200/{}_to_200/allele_tables/{}_downsample_allele_table.txt'.format(barcode, name), sep = '\t', index_col = 0)

    temp_allele = remove_NaN(temp_allele)

    hybrid_tree_pipeline(temp_allele, name, total_priors, '/Genomics/chanlab/blaw/TLS/sandbox/ternary_studies/downsampling_to_200/{}_to_200/trees/'.format(barcode))