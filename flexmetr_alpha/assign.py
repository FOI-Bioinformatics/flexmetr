#!/usr/bin/env python3

import os
import sys
import argparse
import shutil

import global_functions


### Parse input arguments
# setup
argparser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-i','--input',required=True,help='Path to input file or "-" to read from stream')
argparser.add_argument('-o','--output',required=False,default='-',help='Path to output file or "-" to print stream')

# argparser.add_argument('-d','--database','-db','--db',required=False,default='flexmetr/db.sql',help='Path to database (default:flexmetr/db.sql)') # not implemented yet

argparser.add_argument('-c','--column','--col','-col','--columns','--cols','-cols',required=True,help='Column in metadata to format output. Multiple columns may be specified by comma (e.g.: family,genus,species). May be specified as a formatter-string with # preceding the database column key (e.g.: #family_#genus_#species)')
argparser.add_argument('-s','--separator','--sep','-sep',required=False,default='_',help='Separator to use in output between multiple columns (default: underscore/_)')

argparser.add_argument('--clean_names',required=False,action='store_true',help='If specified, will attempt to clean input names from text additional to accession number')
argparser.add_argument('--clean_bvbrc',required=False,action='store_true',help='If specified, will attempt to clean various formats related to BVBRC, e.g. paranthesis and "accn|"')

argparser.add_argument('--scan_input',required=False,action='store_true',help='If specified, will assume the input is the path to a directory of file(s) and list them as input')
argparser.add_argument('--rename_files_dir',required=False,default=None,help='If specified with a path, will assume the input is a path of file(s) and attempt to rename them into the specified directory')

argparser.add_argument('--notify_missing',required=False,action='store_true',help='If specified, will assume each row in the input is an accession number and notify the user if there is no database entry')
argparser.add_argument('--skip_missing',required=False,action='store_true',help='If specified, will assume each row in the input is an accession number and skip the entry if there is no match for the entry in the database')
argparser.add_argument('--print_missing',required=False,action='store_true',help='If specified, only print the accession numbers that did not have a match in the database')

argparser.add_argument('--itol_names',required=False,action='store_true',help='If specified, output an ITOL-config-file to rename nodes')
argparser.add_argument('--itol_labels',required=False,action='store_true',help='If specified, output an ITOL-config-file to label nodes')
argparser.add_argument('--itol_colors',required=False,action='store_true',help='If specified, output a ITOL-config-file to color nodes')

argparser.add_argument('--keep_input_names',required=False,action='store_true',help='If specified, will keep the input formatting of names (as opposed to extracting the accession number)')
argparser.add_argument('--replace_spaces',required=False,default='',help='If specified with a character, will replace spaces in formatted text (spaces in metadata cells) (default:not set)')

argparser.add_argument('--metadata_file',required=False,default=None,help='Path to custom metdata file')
argparser.add_argument('--metadata_file_sep',required=False,default='\t',help='Separator to use in custom metadata file (default: tab)')
argparser.add_argument('--metadata_file_accession',required=False,type=int,default=0,help='Column in custom metadata file that holds the accession number (default: first/0)')
argparser.add_argument('--metadata_strip_quotes',required=False,action='store_true',default=False,help='If specified, will strip qutoes from metadata table cells')
argparser.add_argument('--metadata_replace_missing',required=False,default='',help='Replace missing entries in metadata with value (default:not set)')

argparser.add_argument('--id_list',required=False,default=None,help='[ALPHA] Path to custom ID-file to use in addition to NCBI accession numbers. One ID is expected per row in first element after using .split(<sep>)')
argparser.add_argument('--id_list_sep',required=False,default='\t',help='Separator to use in custom ID file (default: tab)')
argparser.add_argument('--id_list_additive',required=False,action='store_true',default=None,help='If specified, will load list and unite it with identifiers found in the default mode when not supplying a list, e.g. accession number. (default: only use identifiers in list)')
#/
# parse input
args = argparser.parse_args()

input_file = args.input
output_file = args.output

metadata_db = args.database
metadata_columns = args.column
out_separator = args.separator

clean_names = args.clean_names
clean_bvbrc = args.clean_bvbrc

scan_input = args.scan_input
rename_input_into_dir = args.rename_files_dir

notify_missing = args.notify_missing
skip_missing = args.skip_missing
print_missing = args.print_missing

itol_names_out = args.itol_names
itol_labels_out = args.itol_labels
itol_colors_out = args.itol_colors

keep_input_names = args.keep_input_names
replace_spaces = args.replace_spaces

metadata_file = args.metadata_file
metadata_file_sep = args.metadata_file_sep
metadata_file_accession = args.metadata_file_accession
metadata_strip_quotes = args.metadata_strip_quotes
metadata_replace_missing_with = args.metadata_replace_missing

custom_input_ID_list_path = args.id_list
custom_input_ID_list_sep = args.id_list_sep
id_list_additive = args.id_list_additive
#/
###/

## Read input (stdin or path[to scan] or file)
if input_file == '-':
    input_string = sys.stdin.read()
elif scan_input:
    # compile "input_string" as the equivalent of "ls <input_path>"
    input_string = '\n'.join(input_file+'/'+file_ for file_ in os.listdir(input_file))
else:
    with open(input_file,'r') as f:
        input_string = f.read()
##/

## Format default output (applicable for simple formatting. User-formatted output with # will replace this, below)
out_keys = metadata_columns.split(',')
# remove preceding spaces
for enum,_ in enumerate(out_keys):
    out_keys[enum] = out_keys[enum].replace(' ','')
#/
##/

## Format custom output
if metadata_columns.find('#') != -1:
    # parse split_key (will assume it is the first sign after the last '#'
    out_separator = ''
    hash_sign_passed = False
    for i in metadata_columns[::-1]:
        # toggle switch if a #-sign was traversed
        if i == '#':
            hash_sign_passed = True
            continue
        #/
        # when switch is toggled, parse first non-alphabetical letter or number
        if hash_sign_passed and not (i.isalpha() or i.isnumeric()):
            out_separator = i
            break
        #/
    #/
    # determine which keys to be out-formatted
    if out_separator:
        out_keys = metadata_columns.split(out_separator+'#') # when multiple entires exist (e.g. #family_#genus)
    else:
        out_keys = [metadata_columns] # when a single entry exist (e.g. #genus)
        
    for enum,_ in enumerate(out_keys):
        out_keys[enum] = out_keys[enum].replace('#','')
    #/
##/

## Check if parse custom IDs (in addition to NCBI accession numbers) from file
custom_input_ID_list = set()
if custom_input_ID_list_path:
    with open(custom_input_ID_list_path,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split(custom_input_ID_list_sep)
            custom_ID = line[0]
            custom_input_ID_list.add(custom_ID)
##/

## Run reg-ex over input, matching to GCx_NNNNNNNNN.V (x=F/A, N=1..9, V=1..9)
matches = global_functions.getAccessions(input_string,custom_input_ID_list=custom_input_ID_list,regex_and_list_ids_union=id_list_additive)
##/

## Determine accessions to import from DB
accessions_to_import = set()
for match in matches:
    accessions_to_import.add(match)
##/

## Parse metadata from DB
accession_metadata = {}
##/

## Check if parse metadata from custom file
if metadata_file:
    cols_to_import = set()#out_keys # parse only the keys specified by user input from DB
    accession_metadata = global_functions.parse_metadata_file(metadata_file,separator=metadata_file_sep,accessions_to_import=accessions_to_import,strip_quotes=metadata_strip_quotes,
                                                              cols_to_import=cols_to_import,accession_column=metadata_file_accession)
##/

## Check if user wants to replace missing values in metadata
if metadata_replace_missing_with:
    for accession,metadata in accession_metadata.items():
        for key,val in metadata.items():
            if val == '':
                metadata[key] = metadata_replace_missing_with
##/

## Check if user wants to check every input for database presence
if notify_missing or skip_missing or print_missing:
    # Assuming each row is an accession number in the database, test if it exists in the database
    missing = []
    found = []
    for row in input_string.split('\n'):
        if row:
            accession = global_functions.getAccessions(row,return_first=True,suppress_warning=False,custom_input_ID_list=custom_input_ID_list)
            if not accession or not accession in accession_metadata:
                missing.append(row)
            else:
                found.append(row)
    #/
    # Check if user wants to skip missing entires, then recompile the input string with only found entries
    if skip_missing:
        input_string = '\n'.join(found)+'\n'
    #/
    # Check if user wants to output missing entries, then recompile the input string with only missing entries
    if print_missing:
        input_string = '\n'.join(missing)+'\n'
    #/
##/

## Parse adjacent text to accesion numbers. Used to e.g. clean names, to keep original formatting in ITOL-out files.
matches_wAdj_text = global_functions.get_accessions_adjacent_text(input_string, matches)
##/

## Do tranformation of metadata via accession number in the input
output_string = input_string # make "copy" of input string. The "output_String" will hold the modifications
# Check if there is content adjacent to the accession number to remove (will remove everything upstream/downstream of the accession number until a special charater [ignoring underscore] is hit)
if clean_names:
    for match,match_wAdjText in matches_wAdj_text.items():
        output_string = output_string.replace(match_wAdjText,match)
#/
# Format output by replacing accession number with specified keys
for match in matches:
    # Only proceed if match exist in metadata
    if not match in accession_metadata: continue
    #/
    # Format output by the metadata according to user input
    output_formatted_arr = []
    for format_key in out_keys:
        format_value = accession_metadata[match][format_key]
        # check if we want to try to strip spaces (if they exist) from metadata cell values
        if replace_spaces:
            format_value = format_value.replace(' ',replace_spaces)
        #/
        # Check if cleanup bvbrc input
        if clean_bvbrc:
            format_value = format_value.replace('(',' ')
            format_value = format_value.replace(')',' ')
        #/
        output_formatted_arr.append(format_value)
    
    output_formatted = out_separator.join(map(str,output_formatted_arr))
    #/
    
    output_string = output_string.replace(match,output_formatted)
#/
##/

## Check if cleanup bvbrc input
if clean_bvbrc:
    output_string = output_string.replace('accn|','')
##/

## Check if user wanted to rename an input of files
if rename_input_into_dir:
    # get items in input/output. use the basename only in the output since it will be dumped in a new directory
    input_string_split = []
    input_string_dirnames = set() # keep track of dirnames. this is used for when keys in metadata include slashes.
    for entry in input_string.split('\n'):
        if entry:
            input_string_split.append(entry)
            input_string_dirnames.add(os.path.dirname(entry))
    output_string_split = []
    for entry in output_string.split('\n'):
        if entry:
            ## Remove double slashes if they exist
            entry = entry.replace('//','/')
            ##/
            ## Check if forbidden characters (if slashes are included in metadata this cause issues when getting paths)
            # find out which dirname this input has
            entry_dirname = None
            for input_dirname in sorted(input_string_dirnames,key=lambda x: len(x)): # sort so subdirectories are parsed last and will overwrite assignment
                if entry.find(input_dirname+'/') != -1:
                    entry_dirname = input_dirname+'/'
            #/
            # Determine basename (and clean forbidden characters if they exist)
            if entry_dirname == None:
                basename = os.path.basename(entry)
            else:
                basename = entry.replace(entry_dirname,'')
                
                for forbidden_char in ('/',' ','(',')','|',):
                    if basename.find(forbidden_char) != -1:
                        print('Replaced forbidden character "'+forbidden_char+'" that existed in metadata with "." before moving file: '+basename)
                        basename = basename.replace(forbidden_char,'.')
            #/
            ###/
            output_string_split.append(basename)
    #/
    # make output directory
    if os.path.exists(rename_input_into_dir):
        try:
            usr_inp = input('Warning: directory already exists! Proceed? (y/n) ')
            if not usr_inp.lower() in ('y','yes',):
                sys.exit('Terminated!')
        except EOFError: # this error occurs when running "assign"-module in a pipe. Handle it accordingly
            sys.exit('\rOutput directory already exists! Please remove it before proceeding')
        except Exception as e:
            sys.exit('Unknown error, terminated! ' +str(e))
    
    print('Making new directory: '+rename_input_into_dir)
    if not os.path.exists(rename_input_into_dir):       os.makedirs(rename_input_into_dir)
    #/
    # do copy-in with formatting in the new directory
    print('Copying and formatting input files into new directory')
    for i,_ in enumerate(input_string_split):
        source_file = input_string_split[i]
        target_file = rename_input_into_dir + '/' + output_string_split[i]
        shutil.copy(source_file,target_file)
    #/
    # reset output_string so it doesnt print
    output_string = ''
    #/
##/

## Check if user wanted to produce an ITOL-renaming file
if itol_names_out:
    # define base content of ITOL output
    itol_string = '\n'.join(['DATASET_TEXT','SEPARATOR COMMA','DATASET_LABEL,additional','DATA'])+'\n'
    #/
    # append ITOL output with "accession -> label" as "accession,label"
    for match in matches:
        # Format output by the metadata according to user input
        output_formatted_arr = []
        for format_key in out_keys:
            # set default value as the match
            value = match
            #/
            # update the value with metadata, if available
            had_value = False
            if match in accession_metadata:
                value = accession_metadata[match][format_key]
                had_value = True
            #/
            # save [unique] values: hinders accessions without metadata to be added multiple times
            if had_value or not value in output_formatted_arr:
                output_formatted_arr.append(value)
            #/
        output_formatted = out_separator.join(map(str,output_formatted_arr))
        #/
        # compile row
        map_from = match
        map_to = output_formatted
        
        #@ check if keep original input name (and not only the accession number)
        if keep_input_names:        map_from = matches_wAdj_text[match]
        #@/
        
        itol_string += ','.join([map_from,map_to,'-1','#000000','normal','1','0'])+'\n' # last 4 columns describe label and font formatting
        #/
    #/
    # set output string as the ITOL-content
    output_string = itol_string
    #/
##/

## Check if user wanted to produce an ITOL-LABEL-renaming file
if itol_labels_out:
    # define base content of ITOL output
    itol_string = '\n'.join(['LABELS','SEPARATOR COMMA','DATA'])+'\n'
    #/
    # append ITOL output with "accession -> label" as "accession,label"
    for match in matches:
        # Format output by the metadata according to user input
        output_formatted_arr = []
        for format_key in out_keys:
            output_formatted_arr.append(accession_metadata[match][format_key])
        
        output_formatted = out_separator.join(map(str,output_formatted_arr))
        #/
        
        itol_string += ','.join([match,output_formatted])+'\n'
    #/
    # set output string as the ITOL-content
    output_string = itol_string
    #/
##/

## Check if user wanted to produce an ITOL-coloring file
if itol_colors_out:
    # Import matplotlib colors
    try:
        import matplotlib.pyplot as plt
    except:
        sys.exit('Error: could not import matplotlib. To use this feature library "matplotlib" must be installed. Example: conda install matplotlib')
    #/
    # Determine what classifications exist amongst the datasets
    classifications_accessions = {} #classification -> accession
    accessions_classifications = {} #accession -> classification
    for match in matches:
        classis_raw = []
        for format_key in out_keys:
            # set default value as the match
            value = match
            #/
            # update the value with metadata, if available
            had_value = False
            if match in accession_metadata:
                value = accession_metadata[match][format_key]
                had_value = True
            #/
            # save [unique] values: hinders accessions without metadata to be added multiple times
            if had_value or not value in classis_raw:
                classis_raw.append(value)
            #/
        classi = '||'.join(map(str,classis_raw))
        
        if classi == '': continue # skip if empty AKA "no classification available"
        
        if not classi in classifications_accessions:      classifications_accessions[classi] = set()
        classifications_accessions[classi].add(match)
        
        accessions_classifications[match] = classi
    #/
    # Determine a color for each classification
    num_colors = len(classifications_accessions)
    classifications_colors = {}
    
    color_map = plt.get_cmap('rainbow')
    normalized_colors = plt.Normalize(vmin=0, vmax=num_colors-1)  # For some reason ITOL webviewer doesnt accept all colors. Need to normalize them first.
    color_arr = [plt.cm.colors.to_hex(color_map(normalized_colors(i))) for i in range(num_colors)]
    
    for enum,classi in enumerate(sorted(classifications_accessions)):
        classifications_colors[classi] = color_arr[enum]
    #/
    # define base content of ITOL output
    itol_string = '\n'.join(['DATASET_COLORSTRIP','SEPARATOR COMMA','DATASET_LABEL,NodeColors'])+'\n'
    #/
    # append ITOL output with the legend
    itol_string += '\n'.join(['LEGEND_TITLE,Nodes',
                              'LEGEND_SHAPES'+',1'*len(classifications_colors),
                              'LEGEND_LABELS,'+','.join([classi for classi,color in sorted(classifications_colors.items(),key=lambda x: x[0])]),
                              'LEGEND_COLORS,'+','.join([color for classi,color in sorted(classifications_colors.items(),key=lambda x: x[0])]) ])+'\n'
    #/
    # append ITOL output with "accession,color,classification"
    itol_string += 'DATA\n'
    for match in matches:
        if not match in accessions_classifications: continue # skip if it had no color assigned
        classi = accessions_classifications[match]
        color = classifications_colors[classi]
        
        # check if keep original input name (and not only the accession number)
        map_from = match
        if keep_input_names:        map_from = matches_wAdj_text[match]
        #/
        
        itol_string += ','.join([map_from,color,classi])+'\n'
    #/
    # set output string as the ITOL-content
    output_string = itol_string
    #/
##/

## Output (stdout or file)
if output_string:
    if output_file == '-':
        sys.stdout.write(output_string)
    else:
        with open(output_file,'w') as nf:
            nf.write(output_string)

if notify_missing:
    print('\nEntries in database (found,missing)=('+str(len(found))+','+str(len(missing))+')')
    if skip_missing:
        print('Skipped missing files!')
    else:
        print('Skipped files included and untouched!')
##/