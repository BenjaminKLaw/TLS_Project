'''
Author: BL
Last_date_modified: 2023.06.16
Purpose: To run the lineage reconstruction of the explant multiseq and timepoint allele tables

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
output = '/Genomics/chanlab/blaw/TLS/sandbox/explant/3_lineage_reconstruction/'
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

def hybrid_tree_pipeline(table, name, total_priors, best_seed, output_loc):
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(table, allele_rep_thresh = 1, mutation_priors = total_priors)  

    character_matrix.to_csv(output_loc + name + '_character_matrix.txt', sep = '\t')

    with open(output_loc + name + '_priors.pickle', 'wb') as handle:
        pickle.dump(priors, handle, protocol = pickle.HIGHEST_PROTOCOL)
    with open(output_loc + name + '_state_2_indel.pickle', 'wb') as handle:
        pickle.dump(state_2_indel, handle, protocol = pickle.HIGHEST_PROTOCOL)

    cas_tree = cas.data.CassiopeiaTree(character_matrix = character_matrix, priors = priors)
    print(cas_tree.n_cell, cas_tree.n_character, '\n')
    
    missing_proportion = (character_matrix == -1).sum(axis = 0) / character_matrix.shape[0]
    uncut_proportion = (character_matrix == 0).sum(axis = 0) / character_matrix.shape[0]
    n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis = 0)

    cell_meta = table.groupby('cellBC').agg({"intBC": 'nunique', 'UMI': 'sum', 'cell_state': 'unique', 'Timepoint': 'unique', 'orig.ident':'unique'})
    character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 'uncut_prop', 'n_unique_states']).T

    cas_tree.cell_meta = cell_meta
    cas_tree.character_meta = character_meta
    cas_tree.cell_meta.to_csv(output_loc + name + '_metadata.txt', sep = '\t')
    
    cas_tree.cell_meta['cell_state'] = cas_tree.cell_meta['cell_state'].astype('str')
    
    top = cas.solver.VanillaGreedySolver()
    bottom = cas.solver.ILPSolver(convergence_time_limit=3600, maximum_potential_graph_layer_size=10000, weighted=True, seed=best_seed)
    
    hybrid_solver = cas.solver.HybridSolver(top_solver=top, bottom_solver=bottom, cell_cutoff=200, threads=10)
    hybrid_solver.solve(cas_tree, logfile = output_loc + 'logs/'+ name + '_hybrid.log')
    
    newickFile = open(output_loc + name + '_hybrid_newick.txt', 'wt')
    newickFile.write(cas_tree.get_newick())
    newickFile.close()
    
    cas_tree.collapse_mutationless_edges(infer_ancestral_characters = True)
    
    noMutationlessEdgesFile = open(output_loc + name + '_hybrid_newick_noMutationlessEdges.txt', 'wt')
    noMutationlessEdgesFile.write(cas_tree.get_newick())
    noMutationlessEdgesFile.close()

tableFiltered = pd.read_csv('/Genomics/chanlab/blaw/TLS/sandbox/explant/3_lineage_reconstruction/allele_table_filtered.txt', sep = '\t', index_col = 0)

### Calculate the total indels using allele_tables from TLS1, TLS2, and TLSCL combined
TLSCL_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-258/2_add_metadata/allele_table_multiSeq+lenti.txt'
TLS1_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-097/1_preprocessing/allele_table.txt'
TLS2_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-098/1_preprocessing/allele_table.txt'
explant_filepath = '/Genomics/chanlab/blaw/TLS/sandbox/explant/2_add_metadata/merged_allele_table_with_metadata.txt'

TLSCL_allele = pd.read_csv(TLSCL_filepath, sep='\t')
TLS1_allele = pd.read_csv(TLS1_filepath, sep='\t')
TLS2_allele = pd.read_csv(TLS2_filepath, sep='\t')
explant_allele = pd.read_csv(explant_filepath, sep = '\t')

# Combine all three datasets
TLSCL_allele = TLSCL_allele[(TLSCL_allele['finalCalls'] != 'Doublet') & (TLSCL_allele['finalCalls'] != 'Negative')]
TLS1_allele['finalCalls'] = 'Bar25'
TLS2_allele['finalCalls'] = 'Bar26'
explant_allele['finalCalls'] = ['explant_' + i for i in explant_allele['finalCalls']]
total_table = pd.concat([TLSCL_allele, TLS1_allele, TLS2_allele, explant_allele])
total_table = remove_NaN(total_table)

# Calculate the priors using each intBC / Multiseq Barcode as a different indel instance
total_priors = cas.pp.compute_empirical_indel_priors(total_table, grouping_variables=['intBC', 'finalCalls'])

barcodes = ['Bar4']
timepoints = [['144h']]

for barcode in barcodes:
    for time in timepoints:
        temp_tableFiltered = tableFiltered[(tableFiltered['finalCalls'] == barcode) & (tableFiltered['Timepoint'].isin(time))].copy()
        
        temp_tableFiltered = remove_NaN(temp_tableFiltered)
        
        #temp_tableFiltered = remove_NaN(temp_tableFiltered)
        if len(time) > 1:
            timepoint = '{}_{}'.format(time[0][:-1], time[1][:-1])
        else:
            timepoint = '{}'.format(time[0][:-1])
        # I will use the best see for TLS1 for now
        hybrid_tree_pipeline(table = temp_tableFiltered,
                             name = '{}_{}'.format(barcode, timepoint), 
                             total_priors = total_priors,
                             best_seed = 222,
                             output_loc = '/Genomics/chanlab/blaw/TLS/sandbox/explant/3_lineage_reconstruction/{}/{}/'.format(barcode, timepoint))