'''
Author: BL
Last_date_modified: 2022.07.21
Purpose: To run the AM-DNA-258 (TLS_TLSCL) dataset through adding lenti and cell state metadata.

Returns:
    Plots over the multiSeq BC and lenti BC coverage
    A new allele table with the added multiSeq and Lenti metadata
'''

# QC and Add Metadata
#######################################################################

import sys
import os
import pandas as pd
import numpy as np
import cassiopeia as cas
import matplotlib as plt
import scipy

name = 'AM-DNA-258'
output_dir = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-258/lineage/2_add_metadata/'
allele_table_loc = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-258/lineage/1_preprocessing/allele_table.txt'
lenti_loc = '/Genomics/chanlab/blaw/TLS/LineageTracer/processing/AM-DNA-194/'
multiSeq_BC_loc = '/Genomics/chanlab/blaw/TLS/LineageTracer/processing/AM-DNA-197/MultiSeqBC.tsv'

# Read in the data
allele_table = pd.read_csv(allele_table_loc, sep='\t',
                          usecols = ['cellBC', 'intBC', 'r1', 'r2', 'r3', 'allele', 'lineageGrp', 'readCount', 'UMI'])

# Read in the multiSeq barcode metadata
multiSeq_BC = pd.read_csv(multiSeq_BC_loc, sep='\t')

multiSeq_BC['cellBC'] = multiSeq_BC['cellBC'].astype(str) + '-1'

# Add the multiSeq metadata and rename columns
multiSeqMerge = allele_table.merge(multiSeq_BC, left_on='cellBC', right_on='cellBC')

multiSeqMerge['finalCalls'] = multiSeqMerge['final.calls']
multiSeq_BC['finalCalls'] = multiSeq_BC['final.calls']
multiSeqMerge.drop(axis='columns', labels='final.calls')
multiSeq_BC.drop(axis='columns', labels='final.calls')

# Plot the multiSeq label coverage of all the UMIs as a pieplot
ax = multiSeqMerge.finalCalls.value_counts().plot(kind='pie', figsize=(10,10))
ax.figure.savefig(output_dir + 'multiSeqBC_UMI_Coverage_pie.png', dpi=900)

# Plot the multiSeq label coverage of all the UMIs as a barplot
ax = multiSeqMerge.finalCalls.value_counts().plot(kind='barh', figsize=(10,10))
ax.figure.savefig(output_dir + 'multiSeqBC_UMI_Coverage_bar.png', dpi=900)

# Plot the multiSeq label coverage of all the cellBCs as a pieplot
ax = multiSeq_BC.finalCalls.value_counts().plot(kind='pie', figsize=(10,10))
ax.figure.savefig(output_dir + 'multiSeqBC_cellBC_Coverage_pie.png', dpi=900)

# Plot the multiSeq label coverage of all the cellBCs as a barplot
ax = multiSeq_BC.finalCalls.value_counts().plot(kind='barh', figsize=(10,10))
ax.figure.savefig(output_dir + 'multiSeqBC_cellBC_Coverage_bar.png', dpi=900)

# Save the allele table with the added multiSeq barcode metadata
multiSeqMerge.to_csv(output_dir + 'allele_table_multiSeq.txt', sep='\t')
print('Saved allele_table_multiSeq.txt')


# Read in the lenti barcode metadata files
lentiBarcodes = ['Bar2', 'Bar3', 'Bar4', 'Bar5', 'Bar6', 'Bar7', 'Bar8', 'Bar9', 'Bar10', 'Bar11', 'Bar12', 'Bar13', 'Bar14',
            'Bar15', 'Bar16', 'Bar18', 'Bar19', 'Bar20', 'Bar21', 'Bar22', 'Bar23', 'Bar24']

lenti_total = pd.read_csv(lenti_loc + 'AM-DNA-194_Bar1_cellBC_to_StFgroup.txt', sep='\t')
for barcode in lentiBarcodes:
    lentiFile = pd.read_csv(lenti_loc + 'AM-DNA-194_' + barcode + '_cellBC_to_StFgroup.txt', sep='\t')
    lenti_total = lenti_total.append(lentiFile)

# Merge the lenti information with the multiSeq allele table
lentiMerge = multiSeqMerge.merge(lenti_total, how='left', left_on='cellBC', right_on='cellBC')

lentiMerge['StF_Group'] = lentiMerge['group']
lentiMerge['StF_Group'] = lentiMerge[['StF_Group']].fillna('empty')

# Plot the lenti barcode coverage of all the cellBC as a pieplot
ax = lentiMerge.StF_Group.value_counts().plot(kind='pie', figsize=(10,10))
ax.figure.savefig(output_dir + 'LentiBC_cellBC_Coverage_pie.png', dpi=900)

# Plot the lenti barcode coverage of all the cellBC as a barplot
ax = lentiMerge.StF_Group.value_counts().plot(kind='barh', figsize=(10,10))
ax.figure.savefig(output_dir + '/Lenti_Barcode_Coverage_bar.png', dpi=900)

# Plot the lenti barcode coverage excluding the cells with no lenti as a pieplot
ax = lentiMerge.group.value_counts().plot(kind='pie', figsize=(10,10))
ax.figure.savefig(output_dir + 'Lenti_Barcode_Coverage_Without_Null_pie.png', dpi=900)

# Plot the lenti barcode coverage excluding the cells with no lenti as a barplot
ax = lentiMerge.group.value_counts().plot(kind='barh', figsize=(10,10))
ax.figure.savefig(output_dir + '/Lenti_Barcode_Coverage_Without_Null_bar.png', dpi=900)

lentiMerge.to_csv(output_dir + 'allele_table_multiSeq+lenti.txt', sep='\t', index = False)
print('Saved allele_table_multiSeq+lenti.txt')
