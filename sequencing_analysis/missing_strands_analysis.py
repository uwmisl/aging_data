import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import pickle
import os
import csv


# Total missing does NOT necessarily have to grow since it is not based on triplicates. So there could be
# more with zero reads in the beginning that do not necessarily stay missing

# this script constructs a table that shows for every time point and temperature of a specific storage condition,
# how many strands stay missing vs how many strands reappear in later time points and how many strands out of all
# strands are missing


file_type = 'f15' # can also be set to 'f8'
# condition = '_DNAStable_'

conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

folders = ['TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5']  # need to iterate through these folders to get all conditions

temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp

columns = ['Total percent missing', 'Percent missing that reappear', 'Percent that stay missing']

for condition in conditions:
    missing_dict = {}
    missing_strands_file = {}

    for temp in temp_list:
        missing_dict[temp] = {}
        for folder in folders:
            missing_dict[temp][folder] = {}
            for column in columns:
                missing_dict[temp][folder][column] = 0.0

    #infile = open('3mer_analysis/stay_missing_strands_f8.pickle', 'rb')
    # you will need to set infile paths to your specific configuration
    infile = open('/Users/leeorg/Documents/GitHub/aging_sequencing_analysis/3mer_analysis/stay_missing_strands_f8.pickle', 'rb')

    stay_missing_f8 = pickle.load(infile)
    infile.close()

    #infile = open('3mer_analysis/stay_missing_strands_f15.pickle', 'rb')
    infile = open('/Users/leeorg/Documents/GitHub/aging_sequencing_analysis/3mer_analysis/stay_missing_strands_f15.pickle', 'rb')

    stay_missing_f15 = pickle.load(infile)
    infile.close()

    file_size = 0.0  # should correspond with total number of strands for given file type

    if file_type == 'f8':
        file_size = 21601.0
        missing_strands_file = stay_missing_f8
    if file_type == 'f15':
        file_size = 7373.0
        missing_strands_file = stay_missing_f15

    # should be percentages that are missing since triplicates affects the raw count
    for temp in temp_list:
        for folder in folders:
            directory = 'subsampled_pickled_coverage_data'
            total_missing = set()
            reappear_missing = 0.0  # keeps track of the count of strands that are missing but reappear later
            triplicate_count = 0.0  # will determine how much the stay_missing and reappear_missing needs to
            # be multiplied by. If it stays 0 after, need to put in NaN for not existing
            for pickle_file in os.listdir(directory):
                # have to do '_' + temp + '_' because the run 65 and temp 65 within the name could be confused
                if condition in pickle_file and '_' + temp + '_' in pickle_file and file_type in pickle_file and \
                        folder in pickle_file:
                    triplicate_count += 1  # make sure to show that there was a file found at a specific place
                    # the value of the high coverage at the given percentile will change for each file
                    infile = open(directory + '/' + pickle_file, 'rb')
                    file_conditions = pickle.load(infile)
                    for strand in file_conditions.keys():
                        if file_conditions[strand] == 0:
                            if strand not in total_missing:
                                total_missing.add(strand)
            if triplicate_count == 0 or folder is 'TP5':  # if not read, set all columns to NaN
                for column in columns:
                    missing_dict[temp][folder][column] = float('nan')
            else:
                # use triplicate count to average percent total of missing
                stay_missing_count = len(missing_strands_file[condition][temp][folder]) * 1.0
                missing_dict[temp][folder][columns[0]] = len(total_missing) * 1.0 / file_size * 100.0
                # count of missing by finding average of total missing subtracted by what stayed missing
                # then take percentage
                print '\n' + folder
                print 'total missing: ' + str(len(total_missing))
                print 'stay missing: ' + str(stay_missing_count)
                print 'triplicate count: ' + str(triplicate_count)
                missing_dict[temp][folder][columns[1]] = (len(total_missing) * 1.0 - stay_missing_count) / len(total_missing) * 100.0
                # calculate stay missing strands percentage based on pickle file of missing strands
                missing_dict[temp][folder][columns[2]] = stay_missing_count / file_size * 100.0

    # print condition + ':\n' + str(missing_dict)

    # print str(missing_strands_file[condition]['75'])

    # put into a table
    rows = ['65: TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5', '75: TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5', '85: TP0',
            'TP1', 'TP2', 'TP3', 'TP4', 'TP5']

    entries = np.array([])
    i = 0
    for temp in temp_list:
        for folder in folders:
            column_entries = np.array([])
            for item in columns:
                entries = np.append(entries, round(missing_dict[temp][folder][item], 3))

    # make columns and rows the correct size
    entries = np.reshape(entries, (len(rows), len(columns)))
    # print 'Entries: ' + str(entries)
    # print
    # print

    # write data to csv file
    csvname = condition + " " + file_type + " " + "missing sequences analysis.csv"
    with open(csvname, 'wb') as csvfile:
        datawriter = csv.writer(csvfile, delimiter=',')
        datawriter.writerow(['Temperature and Time Point', 'Total Percent Missing', 'Percent Missing that Reappear',
                             'Percent that Stay Missing'])
        label = -1
        tp_temp_labels = ['65C: TP0','TP1','TP2','TP3','TP4','TP5','75C: TP0','TP1','TP2','TP3','TP4','TP5','85C: TP0','TP1','TP2','TP3','TP4','TP5']
        for e in entries:
            label += 1
            # skip time point 5 due to lack of general data, since we're leaving TP5 out of all analysis
            if tp_temp_labels[label] != 'TP5':
                datarow = [tp_temp_labels[label]]
                for datapoint in e:
                    datarow.append(datapoint)
                datawriter.writerow(datarow)

    # optional plotting
    # # fig, ax = plt.subplots()
    #
    # fig, ax = plt.subplots(1, 1)
    # ax.axis('tight')
    # ax.axis('off')
    # the_table = ax.table(cellText=entries, colLabels=columns, rowLabels=rows, loc='center', cellLoc='center')
    # plt.subplots_adjust(right=0.8)
    #
    # ax.set_title(condition + ' ' + file_type + ' Missing Strands Analysis')  # , fontsize=16, weight='bold')
    # # table = plt.table(cellText=entries, cellLoc='center', colLabels=columns, rowLabels=rows)
    #
    # #plt.savefig('Plots/MissingStrands/' + 'Subsampled_' + condition + '_' + file_type + '.svg', dpi=1000)
    # plt.savefig('Plots/MissingStrands/' + 'updated_Subsampled_' + condition + '_' + file_type + 'normaltp0.png', dpi=1000)
