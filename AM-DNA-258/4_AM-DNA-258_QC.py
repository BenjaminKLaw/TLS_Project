'''
Author: BL
Last_date_modified: 2022.07.14
Purpose: To run the AM-DNA-098 (TLS_120h_Rep2) dataset through Cassiopeia preprocessing and tree reconstruction.

Returns:
    QC plots showing the
'''

# QC
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
output_dir = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/1_preprocessing/QC/'
umi_table_loc = '/Genomics/chanlab/blaw/TLS/data/AM-DNA-097/lineage/1_preprocessing/umi_table_error_correct.txt'

# Read in data
umi_table = pd.read_csv(umi_table_loc, sep='\t')

# Create tables of the UMI counts grouped by cellBC, cellBC / intBC, and cellBC / UMI
umis_per_cellBC = umi_table.groupby("cellBC", sort=False).size().values
umis_per_intBC = umi_table.groupby(["cellBC", "intBC"], sort=False).size().values
reads_per_umi = umi_table.groupby(['cellBC', 'UMI'])['readCount'].sum()

# Calculate histrograms for distribution plots
reads_per_umi_hist = np.histogram(reads_per_umi, bins = 100)
umis_per_cellBC_hist = np.histogram(umis_per_cellBC, bins = 50)
umis_per_intBC_hist = np.histogram(umis_per_intBC, bins = 50)

# Create a table of the number of UMI counts per cellBC / intBC pair
umis_per_intBC_df = pd.DataFrame({'count' : umi_table.groupby([ "cellBC", "intBC"] ).size()}).reset_index()

# Calculate the average number of UMIs and the percent cell coverage per intBC
mean_UMI = []
percent_cells = []

intBCs = np.unique(umis_per_intBC_df['intBC'])
total_cells = len(np.unique(umis_per_intBC_df['cellBC']))

for intBC in intBCs:
    temp = umis_per_intBC_df[umis_per_intBC_df['intBC'] == intBC]
    
    test_cells = len(np.unique(temp['cellBC']))
    
    mean_UMI.append(temp['count'].mean())
    
    percent_cells.append(100 * test_cells / total_cells)

# Plot the distribution of the number of reads per UMI
fig = plt.subplots()
plt.bar(x = reads_per_umi_hist[1][:-1], height = reads_per_umi_hist[0])
plt.title('TLS1 Reads per UMI')
plt.ylabel('UMI Count')
plt.xlabel('Number of Reads')
#plt.yscale('log')
plt.savefig(output_dir + 'reads_per_UMI.png', dpi = 300)
plt.close()

# Plot the distributino of the number of UMIs per cellBC
fig = plt.subplots()
x_values = range(0, len(umis_per_cellBC))
plt.plot(x_values, np.flip(np.sort(umis_per_cellBC)), '-')
plt.title('TLS1 UMI per CellBC')
plt.ylabel('Number of UMI')
plt.xlabel('Rank Order')
plt.xscale('log')
plt.yscale('log')
plt.savefig(output_dir + 'UMI_per_cellBC.png', dpi = 300)
plt.close()

# Plot the distribution of the number of UMIs per intBC
fig = plt.subplots()
plt.bar(x = umis_per_intBC_hist[1][:-1], height = umis_per_intBC_hist[0], width = 10)
plt.title('TLS1 UMI per intBC')
plt.ylabel('intBC Count')
plt.xlabel('Number of UMIs')
#plt.yscale('log')
plt.savefig(output_dir + 'UMI_per_intBC.png', dpi = 300)
plt.close()

# Plot the distribution of UMI counts for each intBC
fig, ax = plt.subplots(figsize = (10, 4))
sns.violinplot(ax=ax, data = umis_per_intBC_df, x = 'intBC', y = 'count', scale = 'count', cut = 0, color = 'skyblue')
plt.title('TLS1 UMIs per intBC')
ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
ax.set_ylabel('intBC Count')
plt.savefig(output_dir + 'UMI_per_Target_Site.png', bbox_inches='tight', dpi = 300)
plt.close()

# Plot the average UMI per intBC against the percent of cells per intBC
fig = plt.subplots()
to_label = ['TAACTTTTGAGACA', 'TCTTCAATAGTTTT', 'TGAGTGTAACACTG']

plt.scatter(x = mean_UMI, y = percent_cells, s = 10)
for i, label in enumerate(intBCs):
    if label in to_label:
        plt.annotate(label, (mean_UMI[i], percent_cells[i]))

plt.title('TLS1 Mean UMI vs % Cells')
plt.ylabel('% of Cells')
plt.xlabel('Mean UMIs per intBC')
#plt.yscale('log')
plt.savefig(output_dir + 'Mean_UMI_vs_Percent_Cells.png', dpi = 300)
plt.close()