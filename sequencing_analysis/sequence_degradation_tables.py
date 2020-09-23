import os
import pickle
import csv

# This script attempts so show a small display of the randomness of sequence degradation
# It will do so by showing the percentage overlap of missing sequences for the triplicates for all conditions
# and time points through a table in the form of a csv file

conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

folders = ['TP0', 'TP1', 'TP2', 'TP3', 'TP4']  # need to iterate through these folders to get all conditions

temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp

files = ['f8', 'f15']


# create a dictionary for all conditions for a certain time point. Format will be something like the following,
# where the innermost data is a list of sets corresponding to which sequences are missing.
# Gentegra:{65: f8:{TP1: 234, 4533, 31, 134}, {8372, 23, 8013}..., f15: ....}, 75: TP1..., TP2,}, Imagene: ...
missing_sequence_set_dict = {}
for condition in conditions:
    missing_sequence_set_dict[condition] = {}
    for temp in temp_list:
        missing_sequence_set_dict[condition][temp] = {}
        for file_type in files:
            # start of list of sets
            missing_sequence_set_dict[condition][temp][file_type] = {}
            for time_point in folders:
                missing_sequence_set_dict[condition][temp][file_type][time_point] = []



# iterate through and tally all the missing sequences in sets separate for each triplicate
directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/subsampled_pickled_coverage_data'
# directory = '/Users/leeorg/Documents/GitHub/aging_sequencing_analysis/subsampled-pickled_coverage_data'
for pickle_file in os.listdir(directory):
    for condition in conditions:
        for temp in temp_list:
            for file_type in files:
                for time_point in folders:
                    # have to do '_' + temp + '_' because the run 65 and temp 65 within the name could be confused
                    if condition in pickle_file and '_' + temp + '_' in pickle_file and file_type in pickle_file and \
                            time_point in pickle_file:
                        # the value of the high coverage at the given percentile will change for each file
                        infile = open(directory + '/' + pickle_file, 'rb')
                        file_conditions = pickle.load(infile)
                        infile.close()
                        # add a set to the inner list in the dictionary
                        missing_sequence_set = set()
                        # search for missing strands. If missing, add the set
                        for strand in file_conditions.keys():
                            if file_conditions[strand] == 0:
                                missing_sequence_set.add(strand)
                        missing_sequence_set_dict[condition][temp][file_type][time_point].append(missing_sequence_set)


# iterate again. If triplicate size 3: do 3-way venn, if triplicate size 2: do 2-way venn, else: no venn at all
percent_overlap_dict = {}
for condition in conditions:
    percent_overlap_dict[condition] = {}
    for temp in temp_list:
        percent_overlap_dict[condition][temp] = {}
        for file_type in files:
            percent_overlap_dict[condition][temp][file_type] = {}
            for time_point in folders:
                all_sets = missing_sequence_set_dict[condition][temp][file_type][time_point]
                if len(all_sets) == 3:
                    # create and save plot right here
                    union = all_sets[0].union(all_sets[1].union(all_sets[2]))
                    intersection = all_sets[0].intersection(all_sets[1].intersection(all_sets[2]))
                    percent_overlap = float(len(intersection)) / float(len(union)) * 100.0
                    percent_overlap_dict[condition][temp][file_type][time_point] = int(round(percent_overlap))
                    print "Percent overlapping: " + str(percent_overlap)
                elif len(all_sets) == 2:
                    # create and save plot right here
                    union = all_sets[0].union(all_sets[1])
                    intersection = all_sets[0].intersection(all_sets[1])
                    percent_overlap = float(len(intersection)) / float(len(union)) * 100.0
                    percent_overlap_dict[condition][temp][file_type][time_point] = int(round(percent_overlap))
                    print "Percent overlapping: " + str(percent_overlap)


# iterate through again and construct csv file
# each table will have f8 and f15 in two columns
# the rows are timepoint and temp. N/A in entry if not more than 1 triplicate remaining
for condition in conditions:
    fieldnames = ['Condition', 'TP0', 'TP1', 'TP2', 'TP3', 'TP4']  # top of csv file
    name_of_csv = "Plots/TableTriplicates/%s_analyzed.csv" % condition
    with open(name_of_csv, 'w') as csvfile:  # create top of list
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for file_id in files:
            for temp in temp_list:
                # row will be the percent overlapping among triplicates for all time points
                datadict = percent_overlap_dict[condition][temp][file_id]
                datadict['Condition'] = temp + ' ' + file_id
                writer.writerow(datadict)  # add rows here
