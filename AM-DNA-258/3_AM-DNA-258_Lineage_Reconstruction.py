'''
Author: BL
Last_date_modified: 2022.07.21
Purpose: To run the AM-DNA-258 (TLS_TLSCL) dataset through tree reconstruction.

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
name = 'AM-DNA-258'
output_dir = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-258/lineage/3_lineage_reconstruction/'
alleleTable_filepath = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-258/lineage/2_add_metadata/allele_table_multiSeq+lenti.txt'
seed_table_loc = '/Genomics/chanlab/blaw/TLS/sandbox/AM-DNA-258/top_seeds_tree_statistics.txt'
cell_state_loc = '/Genomics/chanlab/blaw/TLS/LineageTracer/scRNAseq/TLS_TLSCL_cellBC_cellState.tsv'
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

def make_folders(output_dir, name):
    '''
    input: a name for the folders
    makes a folder structure to store the TLS_TLSCL dataset in the data folders
    '''
    if not os.path.exists(output_dir + name):
        os.makedirs(output_dir + name)
        os.makedirs(output_dir + name + '/logs')
    elif not os.path.exists(output_dir + name + '/logs'):
        os.makedirs(output_dir + name + '/logs')
    else:
        print('the folder name is already made')

def ilp_tree_pipeline(table, name, indel_priors, best_seed, output):
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(table, allele_rep_thresh = 1,
                                                                                             mutation_priors = indel_priors)  

    character_matrix.to_csv(output + name + '_character_matrix.txt', sep = '\t')

    with open(output + name + '_priors.pickle', 'wb') as handle:
        pickle.dump(priors, handle, protocol = pickle.HIGHEST_PROTOCOL)
    with open(output + name + '_state_2_indel.pickle', 'wb') as handle:
        pickle.dump(state_2_indel, handle, protocol = pickle.HIGHEST_PROTOCOL)
    
    cas_tree = cas.data.CassiopeiaTree(character_matrix = character_matrix, priors = priors)
    print(cas_tree.n_cell, cas_tree.n_character, '\n')
    
    missing_proportion = (character_matrix == -1).sum(axis = 0) / character_matrix.shape[0]
    uncut_proportion = (character_matrix == 0).sum(axis = 0) / character_matrix.shape[0]
    n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis = 0)

    cell_meta = table.groupby('cellBC').agg({"intBC": 'nunique', 'UMI': 'sum', 'cell_state': 'unique', 'StF_Group':'unique'})
    character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 
                                                                                                    'uncut_prop', 'n_unique_states']).T

    cas_tree.cell_meta = cell_meta
    cas_tree.character_meta = character_meta
    cas_tree.cell_meta.to_csv(output + name + '_metadata.txt', sep = '\t')
    cas_tree.cell_meta['cell_state'] = cas_tree.cell_meta['cell_state'].astype('str')
    
    ilp_solver = cas.solver.ILPSolver(convergence_time_limit=3600, maximum_potential_graph_layer_size=1000, weighted=True, seed=best_seed)
    ilp_solver.solve(cas_tree, logfile=output + 'logs/' + name + '_logfile.txt')
        
    newfile = open(output + name + '_ilp_newick.txt', 'wt')
    newfile.write(cas_tree.get_newick())
    newfile.close()

    cas_tree.collapse_mutationless_edges(infer_ancestral_characters = True)
    
    noMutationlessEdgesFile = open(output + name + '_ilp_newick_noMutationlessEdges.txt', 'wt')
    noMutationlessEdgesFile.write(cas_tree.get_newick())
    noMutationlessEdgesFile.close()
    
def hybrid_tree_pipeline(table, name, indel_priors, best_seed, output):
    character_matrix, priors, state_2_indel = cas.pp.convert_alleletable_to_character_matrix(table, allele_rep_thresh = 1,
                                                                                             mutation_priors = indel_priors)  

    character_matrix.to_csv(output + name + '_character_matrix.txt', sep = '\t')

    with open(output + name + '_priors.pickle', 'wb') as handle:
        pickle.dump(priors, handle, protocol = pickle.HIGHEST_PROTOCOL)
    with open(output + name + '_state_2_indel.pickle', 'wb') as handle:
        pickle.dump(state_2_indel, handle, protocol = pickle.HIGHEST_PROTOCOL)

    cas_tree = cas.data.CassiopeiaTree(character_matrix = character_matrix, priors = priors)
    print(cas_tree.n_cell, cas_tree.n_character, '\n')
    
    missing_proportion = (character_matrix == -1).sum(axis = 0) / character_matrix.shape[0]
    uncut_proportion = (character_matrix == 0).sum(axis = 0) / character_matrix.shape[0]
    n_unique_states = character_matrix.apply(lambda x: len(np.unique(x[(x != 0) & (x != -1)])), axis = 0)

    cell_meta = table.groupby('cellBC').agg({"intBC": 'nunique', 'UMI': 'sum', 'cell_state': 'unique', 'StF_Group':'unique'})
    character_meta = pd.DataFrame([missing_proportion, uncut_proportion, n_unique_states], index = ['missing_prop', 
                                                                                                    'uncut_prop', 'n_unique_states']).T

    cas_tree.cell_meta = cell_meta
    cas_tree.character_meta = character_meta
    cas_tree.cell_meta.to_csv(output + name + '_metadata.txt', sep = '\t')
    cas_tree.cell_meta['cell_state'] = cas_tree.cell_meta['cell_state'].astype('str')
    
    top = cas.solver.VanillaGreedySolver()
    bottom = cas.solver.ILPSolver(convergence_time_limit=3600, maximum_potential_graph_layer_size=1000, weighted=True, seed=best_seed)
    
    hybrid_solver = cas.solver.HybridSolver(top_solver=top, bottom_solver=bottom, cell_cutoff=200, threads=10)
    hybrid_solver.solve(cas_tree, logfile = output + 'logs/'+ name + '_hybrid.log')
    
    newickFile = open(output + name + '_hybrid_newick.txt', 'wt')
    newickFile.write(cas_tree.get_newick())
    newickFile.close()
    
    cas_tree.collapse_mutationless_edges(infer_ancestral_characters = True)
    
    noMutationlessEdgesFile = open(output + name + '_hybrid_newick_noMutationlessEdges.txt', 'wt')
    noMutationlessEdgesFile.write(cas_tree.get_newick())
    noMutationlessEdgesFile.close()

# Read in the allele table and cell state metadata
allele_table = pd.read_csv(alleleTable_filepath, sep='\t')
cell_state = pd.read_csv(cell_state_loc, sep='\t')

tableMerge = allele_table.merge(cell_state, how='left', left_on='cellBC', right_on='cellBC')

# Remove all cells that are not edited (contain no indels)
testSet = set(tableMerge['cellBC'])
for index, row in tableMerge.iterrows():
    if row['allele'] != '[None][None][None]':
        if row['cellBC'] in testSet:
            testSet.remove(row['cellBC'])

tableFiltered = tableMerge[~tableMerge['cellBC'].isin(testSet)]
tableFiltered = remove_NaN(tableFiltered)
tableMerge.to_csv(output_dir + 'allele_table_total.txt', sep='\t', index = False)
tableFiltered.to_csv(output_dir + 'allele_table_filtered.txt', sep='\t', index = False)

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

# Read in the table of best random seeds
seed_table = pd.read_csv(seed_table_loc, index_col = False)
                    
best_seeds = dict(zip(seed_table.Barcode, seed_table.Seed))
                    
for b, s in best_seeds.items():
    name = b
    best_seed = s
    table = tableFiltered[tableFiltered['finalCalls'] == name]

    make_folders(output_dir, name)
    test = cas.pp.convert_alleletable_to_lineage_profile(table)
    test.to_csv(output_dir + name + '/' + name + '_lineage_table.txt', sep='\t')
    
    if name == 'Bar3':
        hybrid_tree_pipeline(table, name, total_priors, best_seed, output_dir + name + '/')
    else:
        ilp_tree_pipeline(table, name, total_priors, best_seed, output_dir + name + '/')