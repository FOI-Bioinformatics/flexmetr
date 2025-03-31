
import os
import re

def parse_metadata_file(input_file,header_present=True,accession_column=0,separator='\t',
                        strip_quotes=False,accessions_to_import=set(),cols_to_import=set(),
                        discard_id_column=False):
    """
    Function to parse a user-specified metdata file.
    The first row defines column headers. If no header-row is given, columns will be enumerated.
    Assumes accession number is in the first column.
    """
    
    metadata = {}
    with open(input_file,'r') as f:
        
        # check if CSV-file, then use csv library to parse
        if separator == ',':
            import csv
            filehandle = csv.reader(f)
        else:
            filehandle = f
        #/
        
        header = None
        accession_key = None
        for enum,line in enumerate(filehandle):
            # Unless "comma" is used as separator, parse the line (csv library already parsed line for us)
            if not separator == ',':
                line = line.strip('\n')
                line = line.split(separator)
            #/
            
            # check if remove quotes (e.g. bvbrc metadatafile has quoted columns)
            if strip_quotes:
                for lineenum,_ in enumerate(line):
                    line[lineenum] = line[lineenum].replace('"','')
            #/
           
            # parse header
            if enum == 0:
                # Check if parse header from first row
                if header_present:
                    header = line
                #/
                # Else, make enumerated header
                else:
                    header = list(map(str,range(len(line))))
                #/
                # Parse accession number key
                accession_key = header[accession_column]
                #/
                continue
            #/
            
            # parse rows
            row_data = {}
            for colenum,entry in enumerate(line):
                column = header[colenum]
                # check if we know which columns to import, then skip current column if it is not part of that set
                if cols_to_import and (not column in cols_to_import and not column == accession_key): continue
                #/
                row_data[column] = entry
            #/
            
            # skip if no data parsed (i.e. if cols_to_import are specified and that column did not exist)
            if not row_data: continue
            #/
            
            # check if we know which accessions to import, then skip current accession if it is not part of that set
            if accessions_to_import and not row_data[accession_key] in accessions_to_import: continue
            #/
            
            # get ID to save at
            id_save = row_data[header[accession_column]]
            #/
            
            # check if skip id-column
            if discard_id_column:
               del row_data[header[accession_column]]
            #/
            
            # save row at accession number
            metadata[id_save] = row_data
            #/
    
    return metadata

def getAccessions(input_string,custom_input_ID_list=None,return_first=False,
                  suppress_warning=False,regex_and_list_ids_union=False):
    # run default matching using reg-ex against NCBI accession number
    regex_pattern = r"GC[A|F]_\d{9}\.\d*"
    matches = re.findall(regex_pattern, input_string)
    #/
    # check if user supplied a list and do not want to use regex matches
    if custom_input_ID_list != None and not regex_and_list_ids_union:
        matches = [] # reset matches as empty
    #/
    # Check if user had custom input list of IDs to scan
    if custom_input_ID_list != None:
        for custom_ID in custom_input_ID_list:
            if input_string.find(custom_ID) != -1:
                matches.append(custom_ID)
    #/
    # warning
    if not matches and not suppress_warning:
        print('Warning: No ID found at node: '+input_string)
    #/
    # Check if return all matches as a list
    if not return_first:
        return matches
    #/
    # Else, return first match (and warn if there was multiple occurrences of this match)
    else:
        if len(matches) > 1 and not suppress_warning:
            print('Warning: Had multiple regex matches at node')
        # Take first match as accession ID
        match = None
        if matches:
            match = matches[0]
        #/
        return match
    #/

def get_accessions_adjacent_text(input_string,matches):
    matches_wAdj_text = {} # match [e.g. accession number] in input -> match+adjacent text [e.g. <family>_<genus>_<species>_<accession>_<somemoretext>
    for match in matches:
        match_chunks_in_input_string = input_string.split(match)
        
        # Find locations of match in the input
        for enum1,_ in enumerate(match_chunks_in_input_string):
            enum2 = enum1+1
            if enum2 > len(match_chunks_in_input_string)-1: break # stop if we reached the end of splits
            
            chunk1 = match_chunks_in_input_string[enum1]
            chunk2 = match_chunks_in_input_string[enum2]
            
            # find leftside adjacent text
            left_adj = ''
            for character in chunk1[::-1]: # flip it from --left-->chunk--right--> to search <--left--chunk--right-->
                if character.isalpha() or character.isnumeric() or character in ('_','-','.'): # remove dots upwards of extension
                    left_adj += character
                else:
                    break
            left_adj = left_adj[::-1] # flip it back
            #/
            # find rightside adjacent text
            right_adj = ''
            for character in chunk2:
                if character.isalpha() or character.isnumeric() or character in ('_','-',): # do not scan dots: cant separate them from extension(s)
                    right_adj += character
                else:
                    break
            #/
            # compile match+adjacent text
            matches_wAdj_text[match] = left_adj + match + right_adj
            break # Break on first: Assume there is only one entry per accession number in the input.
            #/
        #/
    
    return matches_wAdj_text

