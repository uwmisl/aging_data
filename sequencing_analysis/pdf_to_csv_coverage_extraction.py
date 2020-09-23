import pdf2txt
import os
import csv
import sys

##### CAUTION
##### THIS SCRIPT WILL DELETE ALL .txt FILES IN WORKING DIRECTORY

# takes the first argument after the program name in the terminal, this is the directory name with the relevant info
directory_of_interest = sys.argv[1]

# or you can run the script on a directory of your chosing here
# directory_of_interest = "run57"

# makes a list of all .pdf files in the below specified directory within current directory
all_files = os.listdir(directory_of_interest)

#only looks at files with .pdf in the filename
for filename in all_files:
    if ".pdf" not in filename:
        all_files.remove(filename)

#dumps the .pdf data into .txt format...ugly...but workable!
for file_to_analyze in all_files:
    #changes .pdf to .txt
    output_name = file_to_analyze.split('.pdf')[0] + '.txt'
    print "analyzing file %s" %(file_to_analyze)
    os.system("python pdfminer-20140328/tools/pdf2txt.py -o %s %s/%s" %(output_name,directory_of_interest, file_to_analyze))
    
#gathers appropriate data from all .txt files and puts it into a dictionary (datadict), then deletes the .txt when done
datadict = {}

all_analyzed_files = os.listdir(os.curdir)
for filename in all_analyzed_files:
    ref_line = 0
    referece = 0
    avg = 0
    min = 0
    max = 0
    total_oligos = 0
    prcnt_no_reads = 0
    if ".txt" in filename:
        file = open("%s" %(filename), "rt")
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
            if line_count == reference:
                highest_key = filename + line[:-1]
                datadict[highest_key] = {} ##had trouble making csv play nice with defaultdict, sorry it's verbose
                datadict[highest_key]['Condition'] = filename.split('.txt')[0]
                datadict[highest_key]['Reference'] = line[:-1] #to get rid of \n 
            if line_count == avg:
                datadict[highest_key]['Average'] = line[:-1] #to get rid of \n
            if line_count == min:
                datadict[highest_key]['Min'] = line[:-1] #to get rid of \n 
            if line_count == max:
                datadict[highest_key]['Max'] = line[:-1] #to get rid of \n
            if line_count == total_oligos:
                datadict[highest_key]['Total Oligos'] = line[:-1] #to get rid of \n
            if line_count == prcnt_no_reads:
                datadict[highest_key]['Prcnt 0 Reads'] = line[:-1] #to get rid of \n
        # comment out the next line of code for debugging extraction
        os.remove(filename)

# makes the .csv with relevant data

### CHANGE .CSV DATA NAME HERE##
name_of_csv = "%s_analyzed.csv" %(directory_of_interest)
with open(name_of_csv, 'w') as csvfile:
    fieldnames = ['Condition', 'Reference', 'Average', 'Min', 'Max', 'Total Oligos', 'Prcnt 0 Reads']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for key in datadict:
        writer.writerow(datadict[key])
                            

