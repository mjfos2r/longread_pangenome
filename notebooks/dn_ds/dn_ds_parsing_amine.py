import pandas
import re
import os
import glob
from collections import defaultdict

# ok we are gonna simplify the line parsing. feed it the list of lines,
# the current index 'i',
# and then where the second line needs to be pulled from 'y' (+1, +2, +3, etc)
def parse_lines(lines, i, y):
    line_x = lines[i+1].strip()
    log_l_x = line_x.split(',')[0].split('=')[1]
    aicc_x = line_x.split(',')[1].split('=')[1].split("(")[0]
    line_y = lines[i+y].strip()
    dnds_y = line_y.split('=')[1]
    return log_l_x.strip(), aicc_x.strip(), dnds_y.strip()

# Okay this is how we're gonna parse these files for Amine.
def process_file(file_path):
    # read in file path, crack it open, save lines to a var
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # pull gene name from file path
    gene = os.path.basename(file_path).split('.')[0]
    # set up new dict
    new_line = defaultdict()
    # set up gene key as the name we just parsed
    new_line["gene"] = gene
    # now lets set up all of our fields as None so we don't get any weird artifacts.
    a_gtr = None
    b_gtr = None
    c_gtr = None
    a_cod = None
    b_cod = None
    c_cod = None
    p_value_line = None
    # here are the specific strings that we need to be looking for.
    gtr_string = 'Obtaining the global omega estimate based on relative GTR branch lengths and nucleotide substitution biases'
    dn_ds_string = 'Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model'
    p_value_string = 'Branch-site unrestricted statistical test of episodic diversification'
    # init the counter to 0
    i = 0
    # start looping through all of the lines
    while i < len(lines):
        line = lines[i].strip()
        # if line starts with #, then start checking for those strings.
        if line.startswith('#'):
            # if it's a GTR line:
            if gtr_string in line:
                # parse lines +1 and +3, save to the corresponding vars and put em in the dict
                a_gtr,b_gtr,c_gtr = parse_lines(lines, i, 3)
                print(a_gtr,b_gtr,c_gtr)
                new_line["gtr_log_l"] = a_gtr
                new_line["gtr_aicc"] = b_gtr
                new_line["gtr_dnds"] = c_gtr
            elif dn_ds_string in line:
                # if the dn_ds for the full codon model, pull lines +1 and +2, save as corresponding vars and put em in the dict
                a_cod, b_cod, c_cod = parse_lines(lines, i, 2)
                print(a_cod,b_cod,c_cod)
                new_line["codon_log_l"] = a_cod
                new_line["codon_aicc"] = b_cod
                new_line["codon_dnds"] = c_cod
            elif p_value_string in line:
                # if this the p_value line, pull the next line, strip it down, get the p_value, and then save to var and put into dict
                p_value_line = lines[i+1].strip()
                p_val = p_value_line.split('=')[1].split('*')[0].strip()
                print(p_val)
                new_line["p_val"] = p_val
        # iterate counter
        i += 1
    # return the dict for this line.
    return new_line

######################################################################################################################################

# where the files?
dir_of_results = 'results_v2'

# gimme a list of all the results files
results_files = glob.glob(f'{dir_of_results}/*.txt')

# set up empty list of rows
rows = []

# parse through the results files
for infile in results_files:
    print(infile)
    # set up new_line as the output of our handy function
    new_line = process_file(infile)
    # add it to the list of rows
    rows.append(new_line)

# set the column names
col_names = ["gene","gtr_log_l","gtr_aicc","gtr_dnds","codon_log_l","codon_aicc","codon_dnds","p_val"]

# turn list of rows into dataframe
results_df = pandas.DataFrame(rows, columns=col_names)

# dump to csv
results_df.to_csv('results_busted_amine_v2.csv', sep=',')