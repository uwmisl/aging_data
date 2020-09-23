import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.colors
import numpy as np
import pickle
import os
from scipy import stats

# This script constructs a table that shows what percent of p-values for all storage conditions is significant, by
# trimer

file_type = 'f8'  # can either be for f8 or f15
top_percentile = 5  # change this if you want to calculate the top x percentile of coverages

columns = ['ACA', 'ACG', 'ACT', 'AGA', 'AGC', 'AGT', 'ATA', 'ATC', 'ATG', 'CAC', 'CAG', 'CAT', 'CGA', 'CGC', 'CGT',
           'CTA', 'CTC', 'CTG', 'GAC', 'GAG', 'GAT', 'GCA', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG', 'TAC', 'TAG', 'TAT',
           'TCA', 'TCG', 'TCT', 'TGA', 'TGC', 'TGT', 'GC:AT']

infile = open('3mer_analysis/3mer_frequency_percent_significant_f8.pickle', 'rb')
trimer_significant_tally_f8 = pickle.load(infile)
infile.close()

infile = open('3mer_analysis/3mer_frequency_percent_significant_f15.pickle', 'rb')
trimer_significant_tally_f15 = pickle.load(infile)
infile.close()

pickle_file = ''
if file_type == 'f8':
    pickle_file = trimer_significant_tally_f8
if file_type == 'f15':
    pickle_file = trimer_significant_tally_f15


entries = np.array([])
i = 0
for item in columns:
    percent_significant = int(pickle_file[item]['significant'] / (0.0 + pickle_file[item]['significant'] + pickle_file[item]['not significant']) * 100)
    entries = np.append(entries, percent_significant)

print entries
# make columns and rows the correct size
# print 'Entries: ' + str(entries)
row_labels = ['Percent significant']
entries = np.reshape(entries, (1, len(columns)))
# np.reshape(len(columns), len(row_labels))
# fig, ax = plt.subplots()

fig, ax = plt.subplots(1, 1)
ax.axis('tight')
ax.axis('off')
the_table = ax.table(cellText=entries, colLabels=columns, loc='center', fontsize=20)  # , cellLoc='center')
the_table.set_fontsize(20)
plt.subplots_adjust(right=1, left=0)

ax.set_title('Percent significance of p-values by trimer ' + file_type)  # , fontsize=16, weight='bold')
# table = plt.table(cellText=entries, cellLoc='center', colLabels=columns, rowLabels=rows)

plt.savefig('Plots/3merAnalysis/percent_significant_table_' + file_type + '.png', dpi=1000)

