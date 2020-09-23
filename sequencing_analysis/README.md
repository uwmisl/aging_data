# aging_sequencing_analysis

All python files use Python 2.7

Each script has a more detailed description of its function, inputs and outputs at the top of the file.

..............

**_Common abbreviations_**

TP = time point (tp0 through tp5)

..............

**_Brief description of many files/file types:_**

**Seq_stats_compilation.xlsx** is an excel sheet with each timepoint's data. That includes each samples' average 
coverage ('Average'), the minimum and maximum coverage, the size of the file, and the 'Prcnt 0 reads' which is the 
percent of the file's strands that are missing.

**master_error_analysis.xlsx** an excel sheet that includes everything in Seq_stats_compilation.xlsx in addition to 
more detailed error averages, including averages for insertions, deletions, and substitutions. Contains each timepoint's
data with the run value as the leftmost column.

**pdf_to_csv_coverage_extraction.py** is a python script that helps scrape data from pdf documents (open the file 
for more info)

**generate_rates_subplots.py** generates the average ins/del/sub rate for each temperature and each file for a given
condition. Also shows the standard deviation between the data for the triplicates (though this is not necessarily 
visible given the scale)

**giant_generate_rates_subplots.py** same as generate_rates_subplots.py, except all conditions are plotted on
same 9 subplots

**pdf_to_csv_error_analysis_alternative.py** same as error_analysis (see above) but for files that were read by pdfminer
slightly differently than the majority. Scrapes same data as error_analysis.

**3mer_analysis/3mer_analysis_tables.py** Displays a heat map of p-values from t-test between trimer 
frequencies and GC:AT content of strands that stay missing versus strands with the top 'x'th percentile
coverage.

**missing_strands_analysis.py** Constructs a table for given file id and condition that displays
percentages of missing strands, missing strands that reappear, and missing strands that stay missing. 

**3mer_frequency_by_strand_f"x".pickle** A dictionary containing the trimer frequencies and GC:AT content
of each strand of a given file (f8 or f15)

**gc_content.py** calculates the GC : AT ratio of specified coverage levels for the specified files.

**3mer_analysis_without_primer.py** calculates the percentage of each possible trimer for specified low and high 
coverage levels and displays the results as a bar chart. NOTE: to see x labels optimally, it is best to open the chart
to full screen size.

**pdfminer-20140328** is the directory needed to run the pdf_to_csv python script. _ALWAYS keep in the same directory 
you run your pdf_to_csv python scripts!!_

**runxx** are directories with error analysis PDF files. 

**P2-20160229-1.6M-150-ID15-FP41-RP31.txt** is the text file where you can find each sequence of DNA that was supposed 
to be printed by the DNA synthesizer for file 15, file 8 sequences are not available.

**pickled_coverage_data** is the directory in which there's the pickled dictionary for each file and treatment where 
the sequence id is mapped to how many times it was seen by the sequencer (NOTE- this coverage number is an integer, 
NOT a float!)

**stay_missing_strands_fx.pickle** are both pickle files containing nested dictionaries. The structure is 
conditions: temperatures: timepoints: set of missing sequences that stayed missing from that time point onward. 
An example: `stay_missing_unpickled_data['Imagene']['65']['TP0'] = set([7680, 4100...4995])`
