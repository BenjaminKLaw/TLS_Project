'''
Author: BL
Last_date_modified: 2022.07.13
Purpose: To run the AM-DNA-097 (TLS_120h_Rep1) dataset through Cassiopeia preprocessing and tree reconstruction.

Returns:
    Preprocessing - QC, Filtered BAM files, lineage tables, and allele tables
    Lineage Reconstruction - A filtered lineage table, a newick file containing a tree created by Cassiopeia's Hybrid Solver, tree metadata, and logfiles from the ILP portion of the solver
'''

# Preprocessing
#######################################################################

import sys
import os       
import pandas as pd
import numpy as np
import cassiopeia as cas
import gurobipy
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

name = 'AM-DNA-097'
genome_bam = '/Genomics/chanlab/blaw/TLS/raw_data/AM-DNA-097.bam'
output_dir = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/1_preprocessing/'
reference_filepath = '/Genomics/chanlab/mchan/PCT/v2/ppct60_piggybac_ucoe-ef1a-mch-pct61-ts1-1_amplicon.fa'
BC_whitelist = '/Genomics/chanlab/blaw/TLS/LineageTracer/annotation/List_of_10_Tracer_BCs_IDs.txt'
allow_allele_conflicts = False
n_threads = 8

# Performs initial filtering on the BAM file. Each base of the UMI must pass the PHRED score to pass the filtering.
bam_fp = cas.pp.filter_bam(genome_bam, output_directory = output_dir, quality_threshold = 10, n_threads = n_threads)

# Collapse the UMI
umi_table = cas.pp.collapse_umis(bam_fp, output_directory = output_dir, max_hq_mismatches = 3, max_indels = 2, method = 'likelihood', n_threads = n_threads)

# This function performs UMI and cellBC filtering based off the reads per UMI and UMIs per cell and assigns the most abundant sequence to each UMI
umi_table_resolved = cas.pp.resolve_umi_sequence(umi_table, output_directory = output_dir, min_umi_per_cell = 1, min_avg_reads_per_umi = 2.0, plot = False)
umi_table_resolved.to_csv(output_dir + 'umi_table_resolved.txt', sep='\t', index = False)

# This function takes the cellBC-UMI mapped file from umi_table_resolved and aligns each to a reference sequence
umi_table_aligned = cas.pp.align_sequences(umi_table_resolved, ref_filepath = reference_filepath, gap_open_penalty = 20, gap_extend_penalty = 1, n_threads = n_threads)
umi_table_aligned.to_csv(output_dir + 'umi_table_aligned.txt', sep = '\t', index = False)

umi_table = umi_table_aligned.copy()
umi_table = umi_table[["readName","AlignmentScore","CIGAR","QueryBegin","ReferenceBegin","Seq","UMI","cellBC","readCount"]]
umi_table.reset_index(inplace=True)

# This function compares the indel cutsites to the CIGAR strings of each alignment and produces a dataFrame mapping for each
umi_table_alleles = cas.pp.call_alleles(umi_table, ref_filepath = reference_filepath, barcode_interval = (21, 35), cutsite_locations = [113, 167, 221], cutsite_width = 12, context = True, context_size = 0)
umi_table_alleles.to_csv(output_dir + 'umi_table_alleles.txt', sep='\t', index = False)

umi_table_correct_intbc_whitelist = cas.pp.error_correct_intbcs_to_whitelist(umi_table_alleles, whitelist = BC_whitelist, intbc_dist_thresh = 1)
umi_table_correct_intbc_whitelist.to_csv(output_dir + 'umi_table_correct_intbc_whitelist.txt', sep='\t', index = False)

umi_table_error_correct = cas.pp.error_correct_umis(umi_table_correct_intbc_whitelist, max_umi_distance = 2, allow_allele_conflicts = allow_allele_conflicts, n_threads = n_threads)

umi_table_error_correct.to_csv(output_dir + 'umi_table_error_correct.txt', sep='\t', index = False)

umi_table_filtered = cas.pp.filter_molecule_table(umi_table_error_correct,
                                                 output_directory = output_dir,
                                                  min_umi_per_cell = 1,
                                                  min_avg_reads_per_umi = 2.0,
                                                  min_reads_per_umi = 0,
                                                  intbc_prop_thresh = 0.5,
                                                  intbc_umi_thresh = 10,
                                                  intbc_dist_thresh = 1,
                                                  doublet_threshold = 0.35,
                                                  allow_allele_conflicts = allow_allele_conflicts,
                                                  plot=False)

umi_table_filtered.to_csv(output_dir + 'umi_table_filtered.txt', sep='\t', index = False)

allele_table = cas.pp.call_lineage_groups(umi_table_filtered,
    output_directory = output_dir,
    min_umi_per_cell = 1,
    min_avg_reads_per_umi = 1,
    min_cluster_prop = 0.005,
    min_intbc_thresh = 0,
    inter_doublet_threshold = 0.35,
    kinship_thresh = 0.25,
    plot = False)
allele_table.to_csv(output_dir + 'allele_table.txt', sep='\t', index = False)

