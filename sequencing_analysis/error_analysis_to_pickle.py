import pickle
import csv


# This file converts the master file sheet to a Pickle file. The CSV file is converted to a dictionary
# of dictionaries with file id's as keys, whose values are dict's from a metric to a float {e.g. 'Average': 2.342}

# create DictReader that can loop through data in CSV
input_file = csv.DictReader(open("master_error_analysis_for_pickle.csv"))  # for some reason seems to print to stdout

outer_dictionary = {}  # for dict with keys of file id's, values will be another dictionary

adjusted_long_id = ""

for row in input_file:
    # print row
    adjusted_long_id = ''
    inner_dictionary = {}
    for value in row:  # if 'Condition' make key of outer_dict. Otherwise, key/value of inner_dict
        if value != 'Condition' and value != 'Reference' and value != '\xef\xbb\xbfRun':
            inner_dictionary[value] = float(row.get(value))  # set value to float
        if value == 'Reference':
            inner_dictionary[value] = row.get(value)
        if value == '\xef\xbb\xbfRun':  # Run column added extra characters in CSV, so remove that from key name
            inner_dictionary['Run'] = int(row.get(value))
        # string for long id must have file id appended to support unique keys:
        adjusted_long_id = row['Condition'] + "_" + row['Reference']
    # make condition a name based on the id name of the outer dictionary
    condition_title = adjusted_long_id[adjusted_long_id.index('_') + 1:]
    condition_title = condition_title[0: condition_title.index('_')]
    inner_dictionary['Condition'] = condition_title
    # print "Condition title:"
    # print condition_title
    outer_dictionary[adjusted_long_id] = inner_dictionary  # associate filename with dict of metrics

# test to see if constructed correctly
for i in outer_dictionary.items():
    print i

# testing to see if files constructed correctly
# print "Number of files: " + str(len(outer_dictionary))

# pickle the outer_dictionary object
# Store data (serialize)
with open('master_file_metrics.pickle', 'wb') as handle:
    pickle.dump(outer_dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
