import matplotlib.pyplot as plt
plt.rcdefaults()
import numpy as np
import pickle
import matplotlib.lines as mlines
import seaborn as sns

# order of time points goes 54-->57-->60-->58 or 59-->62

# this script produces bar plots for averages of insertion/deletion/substitution triplicates among all different runs

# will generate 9 subplots based: first row will be del/sub/ins rates for 65, then next row 75, then 85
# this will be done by plotting the points of the average rate per triplicate as well as error bars

# NOTE: Similar to generate_rates_subplots.py except that all conditions are plotted together!
sns.set()
coverage_threshold = 14  # will exclude plots whose average coverage is below this

conditions_altered = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', 'DNAStable',
              'GenTegra', 'Dried']

# colors in hex for plotting different conditions, according to patterns from other plots
colors = ['#7e3737', '#86d64e', '#cd9fba', '#8f44c8', '#598140', '#d5553a', '#8bd9a8', '#7877c8', '#d4c24d']
# colors = ['#a361c7', '#6ab64c', '#c65c8a', '#4fb79b', '#cb4f42', '#6587cd', '#b9aa46', '#607a38', '#c78143']

# instantiated before loop so that axes super-titles can be added using fig but not drawing over previous plots
fig, ax = plt.subplots()

z = int(0)
for current_condition in conditions_altered:
    # opens the pickle file
    infile = open('master_file_metrics.pickle', 'rb')
    file_conditions = pickle.load(infile)
    infile.close()

    # print "File conditions:"
    # print file_conditions

    # generate title based on condition, condition 2-4 are for variations in naming among id files in Pickle
    condition_title = current_condition
    condition = current_condition

    id_num = 'id8'

    # variables will be a dictionary with a key of condition, time points found should be the same for ins/del/sub
    # dictionary for y's for Insertion and specific id (del and sub would have separate ones) as such :
    # {Gentegra:{65:[tp0:[a,b,c],tp1:[a,b],tp3[c],etc...},{75:},{85:}, Other condition: etc...}

    # will keep track of whether things are triplicates by having a list for
    # afterwards, plot points (with error bars for standard deviation)

    # new condition/temp list that will have 9 elements of inner lists
    all_conditions_dict_id8 = {}
    all_conditions_dict_id15 = {}
    std_dev_conditions_id8 = {}
    std_dev_conditions_id15 = {}
    mean_conditions_id8 = {}
    mean_conditions_id15 = {}

    # we will iterate over these lists to fill in info for each subplot in the all_conditions_dict_id8
    temp_list = ['65', '75', '85']
    rate_list = ['Insertion Average', 'Deletion Average', 'Substitution Average']

    # every subplot has its own set of keys with temp then ins/del/sub as inner dict
    for i in temp_list:
        all_conditions_dict_id8[i] = {}
        std_dev_conditions_id8[i] = {}
        mean_conditions_id8[i] = {}
        all_conditions_dict_id15[i] = {}
        std_dev_conditions_id15[i] = {}
        mean_conditions_id15[i] = {}
        for j in rate_list:
            all_conditions_dict_id8[i][j] = [None] * 5
            std_dev_conditions_id8[i][j] = [None] * 5
            mean_conditions_id8[i][j] = [None] * 5
            all_conditions_dict_id15[i][j] = [None] * 5
            std_dev_conditions_id15[i][j] = [None] * 5
            mean_conditions_id15[i][j] = [None] * 5
            # create numpy array for every time point, which may or may not have a triplicate
            for k in range(0, 5):
                all_conditions_dict_id8[i][j][k] = np.array([])
                all_conditions_dict_id15[i][j][k] = np.array([])

    # Iterating over keys to add error rate info for each file (if it exists)
    for a in file_conditions.keys():
        # check if key contains condition substring and go by temp
        # keep following line on its own line to easily change condition names for variance in names of pdfs
        if file_conditions[a]['Condition'] == condition and file_conditions[a]['Average'] >= coverage_threshold:
            # and 'Beads' not in a and 'PCR' not in a:
            print a
            for b in temp_list:
                for c in rate_list:
                    if file_conditions[a]['Reference'] == 'id8':
                        # if temp part of name
                        if b in a:
                            if file_conditions[a]['Run'] == 54:
                                all_conditions_dict_id8[b][c][0] = np.append(all_conditions_dict_id8[b][c][0],
                                                                             file_conditions[a][c])
                            if file_conditions[a]['Run'] == 57:
                                all_conditions_dict_id8[b][c][1] = np.append(all_conditions_dict_id8[b][c][1],
                                                                             file_conditions[a][c])
                            if file_conditions[a]['Run'] == 60:
                                all_conditions_dict_id8[b][c][2] = np.append(all_conditions_dict_id8[b][c][2],
                                                                             file_conditions[a][c])
                            if file_conditions[a]['Run'] == 58 or file_conditions[a]['Run'] == 59:
                                all_conditions_dict_id8[b][c][3] = np.append(all_conditions_dict_id8[b][c][3],
                                                                             file_conditions[a][c])
                            if file_conditions[a]['Run'] == 62:
                                all_conditions_dict_id8[b][c][4] = np.append(all_conditions_dict_id8[b][c][4],
                                                                             file_conditions[a][c])
                    if file_conditions[a]['Reference'] == 'id15':
                        # if temp part of name
                        if b in a:
                            if file_conditions[a]['Run'] == 54:
                                all_conditions_dict_id15[b][c][0] = np.append(all_conditions_dict_id15[b][c][0],
                                                                              file_conditions[a][c])
                            if file_conditions[a]['Run'] == 57:
                                all_conditions_dict_id15[b][c][1] = np.append(all_conditions_dict_id15[b][c][1],
                                                                              file_conditions[a][c])
                            if file_conditions[a]['Run'] == 60:
                                all_conditions_dict_id15[b][c][2] = np.append(all_conditions_dict_id15[b][c][2],
                                                                              file_conditions[a][c])
                            if file_conditions[a]['Run'] == 58 or file_conditions[a]['Run'] == 59:
                                all_conditions_dict_id15[b][c][3] = np.append(all_conditions_dict_id15[b][c][3],
                                                                              file_conditions[a][c])
                            if file_conditions[a]['Run'] == 62:
                                all_conditions_dict_id15[b][c][4] = np.append(all_conditions_dict_id15[b][c][4],
                                                                              file_conditions[a][c])

    # the keys for all 9 conditions will be the same for the std dev as the all_conditions_dict_id8
    # set the mean dynamically. Calculate mean and std if size is > 0
    # if NaN, numpy just doesn't plot it
    for b in temp_list:
        for c in rate_list:
            # for each time point, zero through 5
            for d in range(5):
                mean_id8 = np.NaN
                std_id8 = np.NaN
                mean_id15 = np.NaN
                std_id15 = np.NaN
                if all_conditions_dict_id8[b][c][d].size > 0:
                    # print(b)
                    mean_id8 = np.mean(all_conditions_dict_id8[b][c][d])
                    std_id8 = np.std(all_conditions_dict_id8[b][c][d])
                if all_conditions_dict_id15[b][c][d].size > 0:
                    mean_id15 = np.mean(all_conditions_dict_id15[b][c][d])
                    std_id15 = np.std(all_conditions_dict_id15[b][c][d])
                mean_conditions_id8[b][c][d] = mean_id8
                std_dev_conditions_id8[b][c][d] = std_id8  # error needs to be same size as columns
                mean_conditions_id15[b][c][d] = mean_id15
                std_dev_conditions_id15[b][c][d] = std_id15  # error needs to be same size as columns

    objects = ('0', '1', '2', '3', '4')
    y_pos = np.arange(len(objects))

    # make 9 subplots with the given information.
    i = 1
    for b in temp_list:
        for c in rate_list:
            cur_plot = plt.subplot(3, 3, i)
            # set dif x and y based on temp or condition
            # keep the ticks but remove the labels
            cur_plot.tick_params(axis='x', labelbottom=False, bottom=True, labelsize=9)
            cur_plot.tick_params(axis='y', labelleft=False, left=True, labelsize=9)
            if i < 4:
                cur_plot.set_title(c, fontsize='small')
            if i > 6:
                cur_plot.tick_params(bottom=True, labelbottom=True)
                cur_plot.set_xticks([0, 1, 2, 3, 4, 5])
            if i % 3 == 1:
                cur_plot.set_ylabel(str(b + 'C'), fontsize='small')
                cur_plot.set_yticks([0, 1, 2])
                cur_plot.tick_params(axis='y', labelleft=True)
            cur_plot.plot(objects, mean_conditions_id8[b][c], 'ro', alpha=0.9, color=colors[z], marker='^')
            cur_plot.plot(objects, mean_conditions_id15[b][c], 'ro', alpha=0.9, color=colors[z])
            cur_plot.set_xlim(-0.5, 4.5)
            cur_plot.set_ylim(-0.5, 2.5)
            cur_plot.errorbar(objects, mean_conditions_id8[b][c],
                              std_dev_conditions_id8[b][c],
                              color=colors[z], ms=5, mew=4, ls='', alpha=0.9, elinewidth=3)
            cur_plot.errorbar(objects, mean_conditions_id15[b][c],
                              std_dev_conditions_id15[b][c], color=colors[z],
                              ms=5, mew=4, ls='', alpha=0.9, elinewidth=3)
            i += 1
    z += 1  # needed to index the colors


# the following code shared between id8 and id15
# common for both subplots
legend_elements = [mlines.Line2D([], [], color='black', marker='^', markersize=10, label='id8'),
                   mlines.Line2D([], [], color='black', marker='.', markersize=10, label='id15'),
                   mlines.Line2D([], [], color=colors[0], markersize=15, label=conditions_altered[0]),
                   mlines.Line2D([], [], color=colors[1], markersize=15, label=conditions_altered[1]),
                   mlines.Line2D([], [], color=colors[2], markersize=15, label=conditions_altered[2]),
                   mlines.Line2D([], [], color=colors[3], markersize=15, label=conditions_altered[3]),
                   mlines.Line2D([], [], color=colors[4], markersize=15, label=conditions_altered[4]),
                   mlines.Line2D([], [], color=colors[5], markersize=15, label=conditions_altered[5]),
                   mlines.Line2D([], [], color=colors[6], markersize=15, label=conditions_altered[6]),
                   mlines.Line2D([], [], color=colors[7], markersize=15, label=conditions_altered[7]),
                   mlines.Line2D([], [], color=colors[8], markersize=15, label=conditions_altered[8])]

# fig, ax = plt.subplots()
plt.figlegend(handles=legend_elements, loc='upper right', fontsize='x-small')
plt.subplots_adjust(right=0.77, bottom=0.11, hspace=0.07, left=0.1, wspace=0.07)

plt.suptitle('All Storage Conditions', x=0.44, y=0.97, fontsize='x-large')
fig.text(0.43, 0.02, 'Time points', ha='center', fontsize='medium')
fig.text(0.01, 0.5, 'Percent', va='center', rotation='vertical', fontsize='medium')
plt.savefig('Plots/ComparingRatesSubplots/NineSubplots/GiantSubplots.svg', dpi=1000)
