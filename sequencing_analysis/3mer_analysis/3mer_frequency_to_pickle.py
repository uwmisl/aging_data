import pickle
import linecache

# this script will pickle the 3mer frequency of every strand for every file for either id8 or id15
# so that the calculations for the 3mer to table t-tests take much less time and are simpler to code
# It will also calculate the ratio of GC:AT for the GC content script.
file_type = 'f15'  # can be either f8 for id8 or f15 for id15
strand_file = ''
front_primer = ''
end_primer = ''


columns = ['ACA', 'ACG', 'ACT', 'AGA', 'AGC', 'AGT', 'ATA', 'ATC', 'ATG', 'CAC', 'CAG', 'CAT', 'CGA', 'CGC', 'CGT',
           'CTA', 'CTC', 'CTG', 'GAC', 'GAG', 'GAT', 'GCA', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG', 'TAC', 'TAG', 'TAT',
           'TCA', 'TCG', 'TCT', 'TGA', 'TGC', 'TGT']

conditions = ['Imagene', 'Sugars', 'Trehalose', 'Filterpaper', 'BeadsDNAStable', 'DNAStablePCR', '_DNAStable_',
              'GenTegra', 'Dried']

folders = ['TP0', 'TP1', 'TP2', 'TP3', 'TP4', 'TP5']  # need to iterate through these folders to get all conditions

temp_list = ['65', '75', '85']  # need to calculate trimers for correct temp


# by strand #, then a dict of 3mer frequencies within the strand # {strand 1: {ACA: 0.02, ACG: 0.04, ...}}
frequency_dict = {}

# Pickle file will be of a dictionary type (strand then frequency of each trimer):
# {'1': 'ACA': 0.33, 'ACG': 0.002, ....} '2': {...}}
# if file_type == 'f8' or file_type == 'id8':
#     strand_file = 'this file is not available to the public'
#     front_primer = 'TTCGTTCGTCGTTGATTGGT'
#     end_primer = 'AAACGGAGCCATGAGTTTGT'
if file_type == 'f15' or file_type == 'id15':
    strand_file = 'ID15_sequences.txt'
    front_primer = 'ACATTCCGTGCCATTGGATT'
    end_primer = 'TCGGCAAATCGTTCCACAAA'

i = 1
for line in open(strand_file):
    strand = linecache.getline(strand_file, i)
    # add 18 because it should start two before end of primer, which is of length 20
    start_index = strand.index(front_primer) + 18
    end_index = strand.index(end_primer) + 2
    strand = strand[start_index:end_index]
    # find frequency of every trimer in each strand
    coverage_dict = {}  # temporary dict to count occurrences of all trimers to then find percentage
    coverage_total = 0.0
    inner_dict = {}  # to add within frequency dict
    for item in columns:
        # add to occurrences for the high coverage files
        occurrence = strand.count(item)
        coverage_dict[item] = occurrence
        coverage_total += occurrence
    for item in columns:  # need second loop to find percentages to create inner_dict
        inner_dict[item] = coverage_dict[item] / coverage_total
    # Calculate GC:AT once per strand
    total_gc = strand.count('G') + strand.count('C') + 0.0
    total_at = strand.count('A') + strand.count('T') + 0.0
    if total_at == 0:  # avoid divide by zero error, though probably will never be necessary
        inner_dict['GC:AT'] = float('nan')
    else:
        inner_dict['GC:AT'] = total_gc / total_at
    frequency_dict[i] = inner_dict
    i = i + 1

print frequency_dict

# pickle the outer_dictionary object
# Store data (serialize)
file_name = '3mer_analysis/3mer_frequency_by_strand_' + file_type + '.pickle'
with open(file_name, 'wb') as handle:
    pickle.dump(frequency_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
