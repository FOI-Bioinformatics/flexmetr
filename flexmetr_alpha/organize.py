#!/usr/bin/env python3

import os
import sys
import argparse
import shutil

import global_functions


### Parse input arguments
# setup
argparser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-i','--input',required=False,default='db:file_path',help='Path to directory of input files (default: use column "file_path" in database)')
argparser.add_argument('-o','--output',required=True,help='Path to output files')

# argparser.add_argument('-d','--database','-db','--db',required=False,default='flexmetr/db.sql',help='Path to database (default:flexmetr/db.sql)') # not implemented yet

argparser.add_argument('--select_column',required=True,help='Column in metadata to select')
argparser.add_argument('--select_values',required=True,help='Value in selected column to use. Multiple values may be specified, separated by comma')

argparser.add_argument('--group_by',required=True,help='Column in metadata to group output by')

argparser.add_argument('--metadata_file',required=False,default=None,help='Path to custom metdata file')
argparser.add_argument('--metadata_file_sep',required=False,default='\t',help='Separator to use in custom metadata file [default: tab]')
argparser.add_argument('--metadata_file_accession',required=False,default='\t',help='Column in custom metadata file that holds the accession number [default: first/0]')
#/
# parse input
args = argparser.parse_args()

input_dir = args.input
output_dir = args.output

metadata_db = args.database

select_column = args.select_column
select_values_raw = args.select_values

group_by = args.group_by

metadata_file = args.metadata_file
metadata_file_sep = args.metadata_file_sep
metadata_file_accession = args.metadata_file_accession
#/
###/

## Parse select_values
select_values = select_values_raw.split(',')
# remove preceding spaces
for enum,_ in enumerate(select_values):
    select_values[enum] = select_values[enum].replace(' ','')
#/
##/

## Parse metadata from DB
accession_metadata = {}
##/
## Parse metadata from custom file
if metadata_file:
    accession_metadata = global_functions.parse_metadata_file(metadata_file)
##/

## Traverse and identify + select accession numbers
input_accessions_paths = {} # accn->file_path
# Case1: Find files in DB
if input_dir == 'db:file_path':
    for accession,metadata in accession_metadata.items():
        # check if accession meets user-requirements
        add_accession = False
        if select_column in metadata:
            for select_value in select_values:
                if metadata[select_column] == select_value:
                    add_accession = True
                    break
        if add_accession:
            input_accessions_paths[accession] = metadata['file_path']
        #/
#/
# Case2: Find files in a local directory (user does not use the paths in the DB)
if input_dir != 'db:file_path':
    for path,dirs,files in os.walk(input_dir):
        for file_ in files:
            # parse accession number from file (files without accession number will return None)
            file_accession_number = global_functions.getAccessions(file_,return_first=True,suppress_warning=True)
            #/
            
            if file_accession_number:
                # check if accession has user-specified selection
                if file_accession_number in accession_metadata:
                    metadata = accession_metadata[file_accession_number]
                    add_accession = False
                    if select_column in metadata:
                        for select_value in select_values:
                            if metadata[select_column] == select_value:
                                add_accession = True
                                break
                    if add_accession:
                        file_path = path+'/'+file_
                        input_accessions_paths[file_accession_number] = file_path
                #/
##/

## Organize output
accessions_organized = {}
for accession in input_accessions_paths:
    metadata = accession_metadata[accession]
    group_by_val = metadata[group_by]
    if not group_by_val in accessions_organized:        accessions_organized[group_by_val] = []
    accessions_organized[group_by_val].append(accession)
##/

## Make directory for outputs and copy-in files
for group_by_val,accessions in accessions_organized.items():
    # define output dir for current group
    tmp_out = output_dir+'/'+group_by_val
    #/
    # check if previous dir exist, we do not expect this
    if os.path.exists(tmp_out):
        sys.exit('Warning: Output directory already exists! Please remove it before proceeding: '+tmp_out)
    #/
    # make the dir
    os.makedirs(tmp_out)
    #/
    # copy-in the files
    for accession in accessions:
        file_path = input_accessions_paths[accession]
        file_basename = os.path.basename(file_path)
        shutil.copy2(file_path,tmp_out+'/'+file_basename)
    #/
##/