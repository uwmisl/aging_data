import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.lines as mlines
import numpy as np
import pickle
import os
import seaborn as sns

# This script attempts to visually justify why it is acceptable to exclude files in the 9 subplots of ins/del/sub
# can be excluded, by plotting average coverage vs. insertion/deletion/substitution rate

# Uses mater_file_metrics.pickle since this is the only pickle associating average coverage with run (since we only
# want TP0, i.e. Run54,

sns.set()

file_type = 'f15' # can also be 'f8'
file_to_color = {'f8': 'red', 'f15': 'blue'}

conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

# excluding TP0 since it has its own plot since percent missing is calculated differently than that of following
# time points
folders = ['TP1', 'TP2', 'TP3', 'TP4', 'TP5']  # need to iterate through these folders to get all conditions

temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp

columns = ['Total percent missing', 'Percent missing that reappear', 'Percent that stay missing']

fig, ax = plt.subplots()

substitution_rate_dict_f8 = {}
substitution_rate_dict_f15 = {}

# this pickle contains ins/del/sub rate info, as well as time point info based on file name and average coverage info
# you'll need to set infile to your own path
infile = open('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/master_file_metrics.pickle', 'rb')
file_conditions = pickle.load(infile)


directory = 'C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/pickled_coverage_data'
for file_key in file_conditions:  # run through all files in pickle to check if TP0
    inner_dict = file_conditions[file_key]
    if 'Run54' in file_key:  # Run54 is equivalent to TP0
        for data in file_key:
            # print inner_dict['Average']
            average_file_coverage = float(inner_dict['Average'])
            if 'id8' in file_key:
                substitution_rate_dict_f8[average_file_coverage] = float(inner_dict['Deletion Average'])
            if 'id15' in file_key:
                substitution_rate_dict_f15[average_file_coverage] = float(inner_dict['Deletion Average'])
            # update the large condition with key of average coverage and value of percent missing

# creating a plot:
x_objects_f8 = (substitution_rate_dict_f8.keys())
x_objects_f15 = (substitution_rate_dict_f15.keys())

# want the y-axis to be logarithmic
# plt.axis([0, 150, 0, 100])
plt.plot(x_objects_f8, substitution_rate_dict_f8.values(), 'ro', color='red', alpha=0.7)
plt.plot(x_objects_f15, substitution_rate_dict_f15.values(), 'ro', color='blue', alpha=0.7)
plt.ylabel('Deletion Rate (Percent)', fontsize=12)
plt.xlabel('Average File Coverage', fontsize=12)

# plt.yscale('log')


# the following code shared between id8 and id15
legend_elements = [mlines.Line2D([], [], color='red', marker='.', markersize=15, label='id8'),
                   mlines.Line2D([], [], color='blue', marker='.', markersize=15, label='id15')]

plt.figlegend(handles=legend_elements, loc='upper right', fontsize='x-small')

# if deletion there, then insertion and substitution are also guaranteed to be there
plt.title('Average Coverage vs Deletion Rate')

# not sure why, but it's making me write the entire path or else there's an error
plt.savefig('C:/Users/Rachel/PycharmProjects/aging_sequencing_analysis/Plots/CoverageThresholdPlots/CoverageThresholdJustificationDeletion.svg', dpi=1000)
plt.close()
