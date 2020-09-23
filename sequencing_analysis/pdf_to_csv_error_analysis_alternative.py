#import pdf2txt
import os
import csv
import sys

# This is exactly like pdf_to_csv_error_analysis,
# except it changes the pattern of finding deletion and
# insertion columns for the minority of pdf files that were
# interpreted slightly differently

##### CAUTION
##### THIS SCRIPT WILL DELETE ALL .txt FILES IN WORKING DIRECTORY

# takes the first argument after the program name in the terminal, this is the directory name with the relevant info
directory_of_interest = sys.argv[1]
# can also specify directory here
# directory_of_interest = "run54"

# makes a list of all .pdf files in the below specified directory within current directory
all_files = os.listdir(directory_of_interest)

# only looks at files with .pdf in the filename
for filename in all_files:
    if ".pdf" not in filename:
        all_files.remove(filename)

# dumps the .pdf data into .txt format...ugly...but workable!
for file_to_analyze in all_files:
    # changes .pdf to .txt
    output_name = file_to_analyze.split('.pdf')[0] + '.txt'
    print "analyzing file %s" % (file_to_analyze)
    os.system(
        "python pdfminer-20140328/tools/pdf2txt.py -o %s %s/%s" % (output_name, directory_of_interest, file_to_analyze))

# gathers appropriate data from all .txt files and puts it into a dictionary (datadict), then deletes the .txt when done
datadict = {}

all_analyzed_files = os.listdir(os.curdir)
for filename in all_analyzed_files:
    ref_line = 0
    reference = 0
    avg = 0
    min = 0
    max = 0
    total_oligos = 0
    prcnt_no_reads = 0
    insertion_avg = 0
    ins_a = 0
    ins_c = 0
    ins_t = 0
    ins_g = 0
    sub_avg = 0
    sub_a_c = 0
    sub_a_g = 0
    sub_a_t = 0
    sub_c_a = 0
    sub_c_g = 0
    sub_c_t = 0
    sub_g_a = 0
    sub_g_c = 0
    sub_g_t = 0
    sub_t_a = 0
    sub_t_c = 0
    sub_t_g = 0
    deletion_avg = 0
    del_a = 0
    del_c = 0
    del_t = 0
    del_g = 0
    if ".txt" in filename:
        file = open("%s" % (filename), "rt")
        line_count = 0
        for line in file:
            line_count += 1
            if "Reference" in line:
                ref_line = line_count
                reference = line_count + 2
                avg = line_count + 4
                min = line_count + 6
                max = line_count + 8
                total_oligos = line_count + 10
                prcnt_no_reads = line_count + 12
            if "Substitution Average" in line:
                sub_avg = line_count + 2
                sub_a_c = line_count + 28
                sub_a_g = line_count + 30
                sub_a_t = line_count + 32
                sub_c_a = line_count + 34
                sub_c_g = line_count + 36
                sub_c_t = line_count + 38
                sub_g_a = line_count + 40
                sub_g_c = line_count + 42
                sub_g_t = line_count + 44
                sub_t_a = line_count + 46
                sub_t_c = line_count + 48
                sub_t_g = line_count + 50
            if "Insertion Average" in line:
                insertion_avg = line_count + 12
                ins_a = line_count + 14
                ins_c = line_count + 16
                ins_t = line_count + 18
                ins_g = line_count + 20
            if "Deletion Average" in line:
                deletion_avg = line_count + 2
                del_a = line_count + 12
                del_c = line_count + 14
                del_t = line_count + 16
                del_g = line_count + 18
            if line_count == reference:
                highest_key = filename + line[:-1]
                datadict[highest_key] = {}  ##had trouble making csv play nice with defaultdict, sorry it's verbose
                datadict[highest_key]['Condition'] = filename.split('.txt')[0]
                datadict[highest_key]['Reference'] = line[:-1]  # to get rid of \n
            if line_count == avg:
                datadict[highest_key]['Average'] = line[:-1]  # to get rid of \n
            if line_count == min:
                datadict[highest_key]['Min'] = line[:-1]  # to get rid of \n
            if line_count == max:
                datadict[highest_key]['Max'] = line[:-1]  # to get rid of \n
            if line_count == total_oligos:
                datadict[highest_key]['Total Oligos'] = line[:-1]  # to get rid of \n
            if line_count == prcnt_no_reads:
                datadict[highest_key]['Prcnt 0 Reads'] = line[:-1]  # to get rid of \n
            if line_count == insertion_avg:
                datadict[highest_key]['Insertion Average'] = line[:-1]  # to get rid of \n
            if line_count == ins_a:
                datadict[highest_key]['%InsA'] = line[:-1]  # to get rid of \n
            if line_count == ins_c:
                datadict[highest_key]['%InsC'] = line[:-1]  # to get rid of \n
            if line_count == ins_t:
                datadict[highest_key]['%InsT'] = line[:-1]  # to get rid of \n
            if line_count == ins_g:
                datadict[highest_key]['%InsG'] = line[:-1]  # to get rid of \n
            if line_count == deletion_avg:
                datadict[highest_key]['Deletion Average'] = line[:-1]  # to get rid of \n
            if line_count == del_a:
                datadict[highest_key]['%DelA'] = line[:-1]  # to get rid of \n
            if line_count == del_c:
                datadict[highest_key]['%DelC'] = line[:-1]  # to get rid of \n
            if line_count == del_t:
                datadict[highest_key]['%DelT'] = line[:-1]  # to get rid of \n
            if line_count == del_g:
                datadict[highest_key]['%DelG'] = line[:-1]  # to get rid of \n
            if line_count == sub_avg:
                datadict[highest_key]['Substitution Average'] = line[:-1]  # to get rid of \n
            if line_count == sub_a_c:
                datadict[highest_key]['%Sub AtoC'] = line[:-1]  # to get rid of \n
            if line_count == sub_a_g:
                datadict[highest_key]['%Sub AtoG'] = line[:-1]  # to get rid of \n
            if line_count == sub_a_t:
                datadict[highest_key]['%Sub AtoT'] = line[:-1]  # to get rid of \n
            if line_count == sub_c_a:
                datadict[highest_key]['%Sub CtoA'] = line[:-1]  # to get rid of \n
            if line_count == sub_c_g:
                datadict[highest_key]['%Sub CtoG'] = line[:-1]  # to get rid of \n
            if line_count == sub_c_t:
                datadict[highest_key]['%Sub CtoT'] = line[:-1]  # to get rid of \n
            if line_count == sub_g_a:
                datadict[highest_key]['%Sub GtoA'] = line[:-1]  # to get rid of \n
            if line_count == sub_g_c:
                datadict[highest_key]['%Sub GtoC'] = line[:-1]  # to get rid of \n
            if line_count == sub_g_t:
                datadict[highest_key]['%Sub GtoT'] = line[:-1]  # to get rid of \n
            if line_count == sub_t_a:
                datadict[highest_key]['%Sub TtoA'] = line[:-1]  # to get rid of \n
            if line_count == sub_t_c:
                datadict[highest_key]['%Sub TtoC'] = line[:-1]  # to get rid of \n
            if line_count == sub_t_g:
                datadict[highest_key]['%Sub TtoG'] = line[:-1]  # to get rid of \n
        # comment out the next line of code for debugging extraction
        #os.remove(filename)

# makes the .csv with relevant data

### CHANGE .CSV DATA NAME HERE##
name_of_csv = "%s_analyzed_alternative.csv" % (directory_of_interest)
with open(name_of_csv, 'w') as csvfile:
    fieldnames = ['Condition', 'Reference', 'Average', 'Min', 'Max', 'Total Oligos', 'Prcnt 0 Reads',
                  'Insertion Average', '%InsA', '%InsC', '%InsT', '%InsG', 'Deletion Average', '%DelA',
                  '%DelC', '%DelT', '%DelG', 'Substitution Average', '%Sub AtoC', '%Sub AtoG', '%Sub AtoT',
                  '%Sub CtoA', '%Sub CtoG', '%Sub CtoT', '%Sub GtoA', '%Sub GtoC', '%Sub GtoT', '%Sub TtoA',
                  '%Sub TtoC', '%Sub TtoG']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for key in datadict:
        writer.writerow(datadict[key])
