import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.colors as mcolors
import numpy as np
import pickle
import os
from scipy import stats

# this script analyzes the frequency of all trimer possibilities, (e.g. ATC, GTC, ...) to determine if there is a
# connection between coverage and certain patterns of trimers. It uses a t-test to determine whether the difference in
# primer patterns of the strands that have zero coverage and stay missing compared to the top 5% coverage strands is
# significant. Also calculates the significance of the GC:AT ratio to see if GC content significantly differs among
# different file coverages
# ALSO pickles a file that counts the total number of significant vs not significant p-values for every trimer over
# ALL files
# Uses subsampled pickled data

# Calculates maps for all conditions of specified file (id8 or id15)

# Displayed as a heat map with columns as trimer patterns and rows as time points 0-5 for temps of 65, 75, and 85
# significance is any p-value below or equal to 0.05
# ALSO creates a table of counts of what percentage of time points were significant for each trimer over all
# temperatures and conditions

# Also pickles the dict for all of the p-values

conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

columns = ['ACA', 'ACG', 'ACT', 'AGA', 'AGC', 'AGT', 'ATA', 'ATC', 'ATG', 'CAC', 'CAG', 'CAT', 'CGA', 'CGC', 'CGT',
           'CTA', 'CTC', 'CTG', 'GAC', 'GAG', 'GAT', 'GCA', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG', 'TAC', 'TAG', 'TAT',
           'TCA', 'TCG', 'TCT', 'TGA', 'TGC', 'TGT', 'GC:AT']

folders = ['TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5']  # need to iterate through these folders to get all conditions

temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp

big_trimer_tally = {}
for trimer in columns:
    big_trimer_tally[trimer] = {}
    big_trimer_tally[trimer]['significant'] = 0.0
    big_trimer_tally[trimer]['not significant'] = 0.0

for current_condition in conditions:
    condition = current_condition  # change this to calculate tables for different conditions

    file_type = 'f15'  # can either be for f8 or f15
    top_percentile = 5  # change this if you want to calculate the top x percentile of coverages

    # you'll have to change the infile location to your own path
    infile = open('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/3mer_frequency_by_strand_f8.pickle', 'rb')
    trimer_frequency_f8 = pickle.load(infile)
    infile.close()

    infile = open('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/3mer_frequency_by_strand_f15.pickle', 'rb')
    trimer_frequency_f15 = pickle.load(infile)
    infile.close()

    infile = open('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/stay_missing_strands_f8.pickle', 'rb')
    stay_missing_f8 = pickle.load(infile)
    infile.close()

    infile = open('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/stay_missing_strands_f15.pickle', 'rb')
    stay_missing_f15 = pickle.load(infile)
    infile.close()

    frequency_file = ''  # this will change depending on whether the coverage pickle file to analyze is for id8 or id15
    # it will be equal to either trimer_frequency_f8 or trimer_frequency_f15

    # will be of certain condition, then: '65: {TP0: {ACA: p-value 1, p-value 2, ...}, ACG: ...}, TP1: {}}
    p_value_dict = {'65': {'TP0': {}, 'TP1': {}, 'TP2': {}, 'TP3': {}, 'TP4': {}, 'TP5': {}}, '75': {'TP0': {},
                           'TP1': {}, 'TP2': {}, 'TP3': {}, 'TP4': {},
                           'TP5': {}}, '85': {'TP0': {}, 'TP1': {},
                           'TP2': {}, 'TP3': {}, 'TP4': {}, 'TP5': {}}}

    # instantiate p-value dict so that each 'TPX' has all trimers, with a numpy array corresponding to a given trimer
    for temp in p_value_dict.keys():
        for time_point in p_value_dict.get(temp).keys():
            for trimer in columns:
                p_value_dict[temp][time_point][trimer] = np.array([])

    # if pickle file contains certain id number in name, then it should inspect that text
    id15_file = 'ID15_sequences.txt'
    strand_file = ''  # this is the file to check trimer occurrences of
    front_primer = ''
    end_primer = ''

    # chooses the correct 3mer frequency pickle to use based on the specified file type (f8 or f15)
    if file_type == 'f8':
        frequency_file = trimer_frequency_f8
        missing_strands_file = stay_missing_f8
        file_id = 'f8'
    else:
        frequency_file = trimer_frequency_f15
        missing_strands_file = stay_missing_f15
        file_id = 'f15'

    # loop that will calculate p-values for every 3mer for every file
    directory = ''
    for temp in temp_list:
        for folder in folders:
            # you will have to change directory information to your own path
            directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/subsampled_pickled_coverage_data'
            for pickle_file in os.listdir(directory):
                # have to do '_' + temp + '_' because the run 65 and temp 65 within the name could be confused
                if condition in pickle_file and '_' + temp + '_' in pickle_file and file_type in pickle_file and \
                        folder in pickle_file:
                    # the value of the high coverage at the given percentile will change for each file
                    infile = open(directory + '/' + pickle_file, 'rb')
                    file_conditions = pickle.load(infile)
                    high_coverage = np.percentile(file_conditions.values(), 100 - top_percentile)
                    infile.close()
                    # after frequency list for each 3mer created, compare that to array of stayed missing strands
                    for item in columns:
                        strand_frequencies_high = []  # for each 3mer for each pickle file, compared with stay missing
                        # strands frequencies
                        strand_frequencies_low = []  # just using this for now before we find the strands that stay
                        # missing
                        for strand in file_conditions.keys():
                            if file_conditions[strand] >= high_coverage:
                                # then add frequency from 3mer frequency pickle file to list
                                strand_frequencies_high.append(frequency_file[strand][item])
                            if strand in missing_strands_file[condition][temp][folder]:
                                strand_frequencies_low.append(frequency_file[strand][item])
                        # print strand_frequencies_high
                        # print 'P-value: ' + str(p_value)
                        # print p_value_dict[temp][folder][item]
                        # if not enough to compare, just fill in table as n/a with -1
                        if len(strand_frequencies_low) < 5:
                            p_value_dict[temp][folder][item] = np.append(p_value_dict[temp][folder][item], -1)
                        else:  # only calculate p-value when there are sufficient frequencies to compare
                            # do t-test and append p-value to file HERE, result returns two things,
                            # only care about p-value
                            other, p_value = stats.ttest_ind(strand_frequencies_high, strand_frequencies_low)
                            p_value_dict[temp][folder][item] = np.append(p_value_dict[temp][folder][item], p_value)

    # pickles the p-value dict
    '''
    file_name = '3mer_analysis/p_values_' + file_type + '.pickle'
    with open(file_name, 'wb') as handle:
        pickle.dump(p_value_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    '''

    # Make a heat map based on the calculated p-values in the dict, columns are trimers in columns array
    rows = ['65: TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5', '75: TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5', '85: TP0',
            'TP1', 'TP2', 'TP3', 'TP4', 'TP5']

    entries = np.array([])
    i = 0
    for temp in temp_list:
        for folder in folders:
            column_entries = np.array([])
            for item in columns:
                # if entry is empty or contains -1, add -1 at entry. Otherwise, take average p_value
                if -1 not in p_value_dict[temp][folder][item] and len(p_value_dict[temp][folder][item]) > 0:
                    column_entries = np.append(column_entries, np.mean(p_value_dict[temp][folder][item]))
                else:
                    column_entries = np.append(column_entries, np.float('nan'))  # this will be displayed as different
                    # color
            entries = np.append(entries, column_entries)
            i += 1

    entries = np.reshape(entries, (len(rows), len(columns)))
    # print 'Entries: ' + str(entries)

    fig, ax = plt.subplots()

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(columns)))
    ax.set_yticks(np.arange(len(rows)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(columns)
    ax.set_yticklabels(rows)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor", fontsize='small')

    # Loop over data dimensions and create text annotations.
    '''
    for i in range(len(rows)):
        for j in range(len(columns)):
            text = ax.text(j, i, entries[i, j],
                           ha="center", va="center", color="w")
    '''

    # method from stackOverflow:


    # trying to create custom colormap
    cim = plt.imread("C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/3mer_analysis_colorgradient.png")
    cim = cim[cim.shape[0]//2, 100:620, :]
    # print cim

    cmap = mcolors.ListedColormap(cim)


    # norm = mcolors.Normalize(vmin=0, vmax=0.1)

    plt.imshow(entries, vmin=0, vmax=0.1, cmap=cmap)
    plt.colorbar(ticks=[0.0, 0.05, 0.1])


    if condition is '_DNAStable_':
        condition = 'DNAStable'
    ax.set_title("P-values of " + condition + ' (' + file_type + ')')
    # ax.set_facecolor('grey')  # should set background to grey so it isn't confused with white
    fig.tight_layout()
    plt.suptitle(condition, y=0.95, x=0.45)

    plt.savefig('/yourpath/outputfolder' + condition + '_' + file_type + '.svg')  # , dpi=1000)

    

# MAKE TABLE FOR BIG TALLIES!

    # for condition in conditions:
    # if condition == current_condition:
    for temp in temp_list:
        for folder in folders:
            for trimer in columns:
                if -1 not in p_value_dict[temp][folder][trimer] and len(p_value_dict[temp][folder][trimer]) > 0:
                    if np.mean(p_value_dict[temp][folder][trimer]) <= 0.05:
                        big_trimer_tally[trimer]['significant'] += 1
                    else:
                        big_trimer_tally[trimer]['not significant'] += 1

print 'All conditions in one' + file_type + ':'
print big_trimer_tally

# pickle the big_trimer_tally
file_name = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/3mer_frequency_percent_significant_' + file_type + '.pickle'
with open(file_name, 'wb') as handle:
    pickle.dump(big_trimer_tally, handle, protocol=pickle.HIGHEST_PROTOCOL)
