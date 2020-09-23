import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.lines as mlines
import numpy as np
import pickle
import os
import seaborn as sns

# This script attempts to visually get an idea of an approximate threshold for average file coverage (e.g. the
# average coverage over all strands for a specific run) by plotting average file coverage versus the percent of
# strands missing. To do so, Time Point zero will be plotted since all conditions can be included in a single
# plot since this is before the strands are affected by their various storage conditions.

# Much of this logic may be similar to missing_strands_analysis.py

sns.set()

file_type = 'f15' # can also be set to 'f8'
file_to_color = {'f8': 'red', 'f15': 'blue'}

conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

# excluding TP0 since it has its own plot since percent missing is calculated differently than that of following
# time points
folders = ['TP1', 'TP2', 'TP3', 'TP4', 'TP5']  # need to iterate through these folders to get all conditions

temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp

columns = ['Total percent missing', 'Percent missing that reappear', 'Percent that stay missing']

fig, ax = plt.subplots()

average_coverage_missing_dict_f8 = {}
average_coverage_missing_dict_f15 = {}

file_size = 0.0

# you'll need to set your own path here
directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/pickled_coverage_data'
for pickle_file in os.listdir(directory):
    if 'TP0' in pickle_file:  # 'TP0' in pickle_file and
        if 'f8' in pickle_file:
            file_size = 21601.0
        if 'f15' in pickle_file:
            file_size = 7373.0
        infile = open(directory + '/' + pickle_file, 'rb')
        file_conditions = pickle.load(infile)
        missing_in_file = 0
        for strand in file_conditions.keys():
            if file_conditions[strand] == 0:
                missing_in_file += 1.0
        # print missing_in_file
        total_percent_missing = missing_in_file / file_size * 100
        average_file_coverage = int(np.mean(file_conditions.values()))
        if 'f8' in pickle_file:
            average_coverage_missing_dict_f8[average_file_coverage] = total_percent_missing
        if 'f15' in pickle_file:
            average_coverage_missing_dict_f15[average_file_coverage] = total_percent_missing
        if average_file_coverage > 3 and average_file_coverage < 15:
            print "Average coverage: " + str(average_file_coverage) + ", Percent missing: " + str(total_percent_missing)
        # update the large condition with key of average coverage and value of percent missing

# creating a plot:
x_objects_f8 = (average_coverage_missing_dict_f8.keys())
x_objects_f15 = (average_coverage_missing_dict_f15.keys())

# want the y-axis to be logarithmic
# plt.axis([0, 150, 0, 100])
plt.plot(x_objects_f8, average_coverage_missing_dict_f8.values(), 'ro', color='red', alpha=0.7)
plt.plot(x_objects_f15, average_coverage_missing_dict_f15.values(), 'ro', color='blue', alpha=0.7)
plt.ylabel('Total Percent Missing', fontsize=12)
plt.xlabel('Average File Coverage', fontsize=12)

plt.yscale('log')

# annotating the least amount of coverage that is plotted for determining coverage threshold (ignoring average
# coverage outliers near zero)
ax.annotate('(14, 2.075)', (7, 2.7), fontsize='x-small')

# the following code shared between id8 and id15
legend_elements = [mlines.Line2D([], [], color='red', marker='.', markersize=15, label='id8'),
                   mlines.Line2D([], [], color='blue', marker='.', markersize=15, label='id15')]

plt.figlegend(handles=legend_elements, loc='upper right', fontsize='x-small')

# if deletion there, then insertion and substitution are also guaranteed to be there
plt.title('Average Coverage vs Total Percent of Strands Missing')

# not sure why, but it's making me write the entire path or else there's an error
plt.savefig('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/Plots/CoverageThresholdPlots/CoverageThresholdTP0Log.svg', dpi=1000)
plt.close()
