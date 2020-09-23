import pickle
import os
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

# NOTE: this uses the SUBSAMPLED pickle files

# will produce a dictionary of sets for missing strands for different time points for 65, 75, and 85 for specified
# id number. It will be pickled for later use in the 3mer analysis to table script

# DEFINITION: to be deemed "missing", a strand must not reappear from the specified time point on. If strand 1,2,3
# are part of time point three 65 missing set, then they will never have a coverage greater than zero in any 65
# files for any time point past either. In addition, the strands that stay missing for time point zero are NOT
# included in the sets of missing strands for TP1-5. This is because it is assumed that these strands are lost in
# pre-processing rather than through the degradation process. For tp1-5, the size of the sets of missing strands
# will either stay the same or grow, because the missing strands in the previous time points (excluding tp0) are
# included in the following sets of missing strands.
# There are different sets of missing strands for not only file id but also by temperature and storage condition.

# FORMAT: A set containing -1 means that there is no sequencing for the specified condition at the specified time point


file_type = 'f15'  # can be either f8 for id8 or f15 for id15
folders = ['TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5']  # need to iterate through these folders to get all conditions
temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp

# NOTE TO SELF: there is currently inconsistent naming among the pickle files so these condition names will NOT
# collect all these files correctly (as of 7/18)
conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

# will have temps as keys and then within that, time points as keys with set a of strands for each
# e.g. {65: {TP0: {212, 22, 1343, 50}, ...}, ...}
stay_missing_dictionary = {}

# instantiate dictionary
for condition in conditions:
    stay_missing_dictionary[condition] = {}
    for temp in temp_list:
        stay_missing_dictionary[condition][temp] = {}
        for folder in folders:
            # = set()
            stay_missing_dictionary[condition][temp][folder] = set()  # creating a set within


# first, do one for time point zero where you find strands that are probably lost in pre-processing
# check that it is zero coverage for all the files
for condition in conditions:
    for temp in temp_list:
        directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/subsampled_pickled_coverage_data'  # /TP0'
        # need to figure out what strands stay missing (how to deal with triplicates???) <---- IMPORTANT
        for pickle_file in os.listdir(directory):  # a certain time point directory
            if 'TP0' in pickle_file and temp in pickle_file and condition in pickle_file and file_type in pickle_file:
                infile = open(directory + '/' + pickle_file, 'rb')
                file_conditions = pickle.load(infile)
                infile.close()
                # find all zero coverage strands in given pickle file
                for strand in file_conditions.keys():
                    if file_conditions[strand] == 0:
                        stay_missing_dictionary[condition][temp]['TP0'].add(strand)
                        # stay_missing_dictionary[condition][temp]['TP0'] = stay_missing_dictionary[condition][temp]
                        # ['TP0'].append(strand)

# now go through and remove from TP0 any strand that ended up having coverage later (since it was not lost in
# pre-processing)
# TP0 will now be set to strands that never show up
for condition in conditions:
    for temp in temp_list:
        for folder in folders:
            directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/subsampled_pickled_coverage_data'  # /'  + folder
            if folder != 'TP0':  # obviously ignore time point since we are comparing following to this set
                # need to figure out what strands stay missing (how to deal with triplicates???) <---- IMPORTANT
                for pickle_file in os.listdir(directory):  # a certain time point directory
                    if folder in pickle_file and temp in pickle_file and condition in pickle_file and file_type in \
                            pickle_file:
                        infile = open(directory + '/' + pickle_file, 'rb')
                        file_conditions = pickle.load(infile)
                        infile.close()
                        # find all zero coverage strands in given pickle file
                        for strand in file_conditions.keys():
                            if file_conditions[strand] != 0 and strand in stay_missing_dictionary[condition][temp]['TP0']:
                                # if strand not in stay_missing_dictionary[condition][temp][folder]:
                                # remove strand if it pops up again
                                stay_missing_dictionary[condition][temp]['TP0'].remove(strand)

# print stay_missing_dictionary['GenTegra']['85']['TP0']

# now check all of zeros, add all except the ones already in time point zero
for condition in conditions:
    for temp in temp_list:
        i = 0
        for folder in folders:
            condition_found = False
            directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/subsampled_pickled_coverage_data'  # /' + folder
            # need to figure out what strands stay missing (how to deal with triplicates???) <---- IMPORTANT
            for pickle_file in os.listdir(directory):  # a certain time point directory
                if folder in pickle_file and temp in pickle_file and condition in pickle_file and file_type in \
                        pickle_file:
                    condition_found = True
                    infile = open(directory + '/' + pickle_file, 'rb')
                    file_conditions = pickle.load(infile)
                    infile.close()
                    # find all zero coverage strands in given pickle file except for those already part of TP0 set
                    for strand in file_conditions.keys():
                        if file_conditions[strand] == 0 and strand not in \
                                stay_missing_dictionary[condition][temp]['TP0']:
                            stay_missing_dictionary[condition][temp][folder].add(strand)
            # if folders does not contain a certain condition, just set it as copy of missing strands as last condition
            i += 1
            if condition_found is False:
                # add something that could not indicate a strand to show it's nonexistent
                stay_missing_dictionary[condition][temp][folders[i - 1]].add(-1)

# for TP1 and on, if zero coverage strand is not included in the time points beyond, remove them (find intersection of
# sets)
# iterate through the stay missing dictionary
for condition in conditions:
    for temp in temp_list:
        i = 1
        while i < len(folders):
            j = i + 1
            # skip the removal if a condition at certain time point does not exist
            if -1 not in stay_missing_dictionary[condition][temp][folders[i]]:
                while j < len(folders):
                    # set current set as intersection of those beyond (so that it will not include strands that reappear
                    # in sequencing
                    if -1 not in stay_missing_dictionary[condition][temp][folders[j]]:
                        updated_set = stay_missing_dictionary[condition][temp][folders[i]].copy()
                        for strand in stay_missing_dictionary[condition][temp][folders[i]]:
                            if strand not in stay_missing_dictionary[condition][temp][folders[j]]:
                                # print "Removed"
                                updated_set.remove(strand)
                        stay_missing_dictionary[condition][temp][folders[i]] = updated_set.copy()
                    j += 1
            i += 1

# print stay_missing_dictionary['Imagene']

# pickle the outer_dictionary object
# Store data (serialize)
file_name = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/3mer_analysis/stay_missing_strands_' + file_type + '.pickle'
with open(file_name, 'wb') as handle:
    pickle.dump(stay_missing_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

