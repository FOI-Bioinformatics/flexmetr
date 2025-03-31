#!/usr/bin/env python3

import os
import sys
import re
import argparse
from random import randint
import matplotlib.pyplot as plt
from matplotlib import collections
import matplotlib.patches as patches

try:
    from Bio import Phylo
    from io import StringIO
except:
    sys.exit('Unable to import Biopython Phylo or StringIO package. Please make sure it has been installed.')

import global_functions

### Parse input arguments
# setup
argparser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-i','--input',required=True,help='Path to input file or "-" to read from stream')
argparser.add_argument('-o','--output',required=False,default='-',help='Path to output file or "-" to print stream')
argparser.add_argument('--plot',required=False,default=None,help='Path to output plot (default: do not output)')

# argparser.add_argument('-d','--database','-db','--db',required=False,default='flexmetr/db.sql',help='Path to database (default:flexmetr/db.sql)') # not implemented yet

argparser.add_argument('-c','--column','--col','-col','--columns','--cols','-cols',required=False,help='Column in metadata to format output. Multiple columns may be specified by comma (e.g.: family,genus,species)')
argparser.add_argument('--branch_columns',required=False,default=None,help='Columns to use to resolve metadata for branch-nodes (default: not set)')


argparser.add_argument('--outgroup',required=False,default=None,help='Set outgroup name (default: not set)')
argparser.add_argument('--remove_outgroup',required=False,action='store_true',default=None,help='If specified, will remove the outgroup from tree (default: not set)')

argparser.add_argument('--cansnper_column',required=False,default=None,help='If specified with a column name, will assume the column holds CanSNPer data (default: not set)')
argparser.add_argument('--cansnps',required=False,default=None,help='If specified with a comma-separated list without spaces, will import these CanSNPs (default: not set)')

argparser.add_argument('--flextaxd_outfiles_path',required=False,default=None,help='If specified with a path, will make a directory to output FlexTaxD-formatted files for arguments --taxonomy_file <treetaxonomy> --genomeid2taxid <genome_map> (default: not set)')
argparser.add_argument('--flextaxd_additional_genomes',required=False,default=None,help='If specified with a path, will load additional genomes from a FlexTaxD genome-map file. For example, with redundant genomes from RepGenR. Implies --flextaxd_outfiles_path (default: not set)')
argparser.add_argument('--flextaxd_root_name_linker',required=False,default=None,help='If specified, will add a "linker" with this name as root, connecting the target node in FlexTaxD to the imported tree (default: not set)')
argparser.add_argument('--flextaxd_apply_downstream',required=False,action='store_true',default=None,help='If specified, will apply branch metadata to all downstream branch nodes (default: not set)')
argparser.add_argument('--flextaxd_rank_columns',required=False,default=None,help='If specified, will resolve these columns as "rank" [third column in parent-child relations file] in FlexTaxD (default: not set)')
argparser.add_argument('--flextaxd_rank_columns_strip',required=False,action='store_true',default=None,help='If specified, will only use these columns to determine "rank" (default: not set)')
argparser.add_argument('--flextaxd_skip_missing_metadata',required=False,action='store_true',default=None,help='If specified, will not write e.g. "None" when a sample lacks a metadata column value (default: not set)')
argparser.add_argument('--flextaxd_keep_node_basenames',required=False,default=None,action='store_true',help='If specified, will keep branch and leaf node basename where metadata is available (default: only show metadata)')

argparser.add_argument('--branch_node_basename',required=False,default='branch',help='Specify basename used for branch labels. Expected format <branch_node_basename><number>')
argparser.add_argument('--leaf_node_basename',required=False,default='leaf',help='Specify basename used for leaf labels. Expected format <leaf_node_basename><number>')


argparser.add_argument('--branchify_leafs',required=False,action='store_true',default=None,help='If specified, will wedge in a branch-node upstream of leaf-nodes (default: not set)')

argparser.add_argument('--metadata_file',required=False,default=None,help='Path to custom metdata file')
argparser.add_argument('--metadata_file_sep',required=False,default='\t',help='Separator to use in custom metadata file [default: tab]')
argparser.add_argument('--metadata_file_accession',required=False,default='\t',help='Column in custom metadata file that holds the accession number [default: first/0]')
#/
# parse input
args = argparser.parse_args()

input_file = args.input
output_file = args.output
plot_output_file_path = args.plot

metadata_db = args.database
db_columns_raw = args.column
db_columns_branches_raw = args.branch_columns

outgroup_dataset = args.outgroup
remove_outgroup = args.remove_outgroup

cansnper_column = args.cansnper_column
cansnps_to_use = args.cansnps

flextaxd_outfiles_path = args.flextaxd_outfiles_path
flextaxd_additional_genomes = args.flextaxd_additional_genomes
flextaxd_root_name_linker = args.flextaxd_root_name_linker
flextaxd_apply_downstream = args.flextaxd_apply_downstream
flextaxd_rank_columns = args.flextaxd_rank_columns
flextaxd_rank_columns_strip = args.flextaxd_rank_columns_strip
flextaxd_skip_missing_metadata = args.flextaxd_skip_missing_metadata
flextaxd_keep_node_basenames = args.flextaxd_keep_node_basenames

branch_node_basename = args.branch_node_basename
leaf_node_basename = args.leaf_node_basename


branchify_leafs = args.branchify_leafs

metadata_file = args.metadata_file
metadata_file_sep = args.metadata_file_sep
metadata_file_accession = args.metadata_file_accession
#/
###/

##### FUNCTIONS
def metadata_parse_canSNPer_column(metadata_dict,db_columns_to_use,cansnps_to_use):
    ## See variable descriptions in argparse
    ## Metadata_dict is the imported metadata
    ## Returns input metadata_dict and db_columns_to_usem, now with expanded canSNPs
    
    # check if user want to parse specific canSNPs
    if cansnps_to_use != None:
        cansnps_to_use_parsed = []
        for cansnp_to_use in cansnps_to_use.split(','):
            if not cansnp_to_use in cansnps_to_use_parsed: # do not allow duplicates. Use array to maintain user input order of cansnps
                cansnps_to_use_parsed.append(cansnp_to_use)
        cansnps_to_use = cansnps_to_use_parsed
    #/
    # make new columns for each canSNP
    cansnps_parsed = set() # keep track of which cansnps were parsed and set these to "False" across samples without the cansnp
    for accession,metadata in metadata_dict.items():
        if cansnper_column in metadata:
            cansnp_assignments = metadata[cansnper_column].replace('"','')
            for cansnp_assignment in cansnp_assignments.split(';'):
                # ensure that there is no other metadata column with this value
                if cansnp_assignment in metadata:
                    print(f'FATAL: When attempting to parse CanSNP assignments from column {cansnper_column} there was already a column value imported for {cansnp_assignment}. Please check your input metadata and ensure that there are no conflicts.')
                    sys.exit()
                #/
                # check if user wants to use this canSNP
                if cansnps_to_use != None and not cansnp_assignment in cansnps_to_use: continue
                #/
                # save canSNP as new column
                metadata[cansnp_assignment] = True
                #/
                # save to tracker
                cansnps_parsed.add(cansnp_assignment)
                #/
    #/
    # exit if no cansnp was found in specified column
    if not cansnps_parsed:
        print(f'FATAL: No CanSNPs found in specified column {cansnper_column}. Please check the spelling and make sure that this column is imported (when using --columns parameter)')
        sys.exit()
    #/
    # set canSNPs to "False" where not present
    for cansnp_parsed in cansnps_parsed:
        for accession,metadata in metadata_dict.items():
            if not cansnp_parsed in metadata:
                metadata[cansnp_parsed] = False
    #/
    # remove cansnp assignment master column
    for accession,metadata in metadata_dict.items():
        del metadata[cansnper_column]
    #/
    # cleanup db_columns and add parsed cansnps while maintaining order
    db_columns_new = []
    for idx,val in enumerate(db_columns_to_use):
        if val == cansnper_column:
            for cansnp_assignment in cansnps_to_use:
                db_columns_new.append(cansnp_assignment)
        else:
            db_columns_new.append(val)
    db_columns_to_use = db_columns_new
    #/
    # return
    return metadata_dict,db_columns_to_use
    #/


def verify_imported_columns(expected_cols=None,imported_metadata=None,errormessage_description=None):
    db_columns_imported = set()
    for accn,data in imported_metadata.items():
        for db_col,value in data.items():
            db_columns_imported.add(db_col)
    
    if not set(expected_cols) == db_columns_imported:
        print('WARNING: Not all specified columns were found in the metadata input:')
        print('\n'.join(list(set(expected_cols).difference(db_columns_imported))))
        print(f'Please check your metadata and input arguments for {errormessage_description}')
        sys.exit()
#####/

### Import stuff
print('Begin import')

## Read input (stdin or file)
if input_file == '-':
    input_string = sys.stdin.read()
else:
    with open(input_file,'r') as f:
        input_string = f.read()
input_string = input_string.strip('\n')
##/

## Parse input into Biopython tree structure
tree_fileobject = StringIO(input_string)
tree = Phylo.read(tree_fileobject,'newick')
##/
## Check if user wants to branchify leaf nodes (e.g. to apply metadata at leaf-level)
if branchify_leafs:
    from Bio.Phylo.BaseTree import Clade
    print('Adding branch node upstream of leaf nodes now')
    
    # dummy-set a root to get path for highest nodes
    dummy_root = Clade(name='dummyroot')
    orig_root = tree.root
    dummy_root.clades.append(orig_root)
    tree.root = dummy_root
    #/
    # add branchnode upstream of leafnodes
    for leaf_node in tree.get_terminals():
        if leaf_node.name == 'dummyroot': continue
        leaf_path = tree.get_path(leaf_node.name)
        
        # bugcheck: all leaves should have a parent. If it does not, it means it is adjacent to undefined root
        if len(leaf_path) == 1:
            print('FATAL: Expected leaf-node to have a parent but it did not! Make fix for this!')
            sys.exit()
        #/
        # bugcheck: last node should be the leaf
        if not leaf_node.name == leaf_path[-1].name:
            print('FATAL: Expected last node in path to be the leaf. Make fix for this!')
            sys.exit()
        #/
        # get leaf parent
        leaf_parent = leaf_path[-2]
        #/
        # clear leaf from parent
        leaf_parent.clades.remove(leaf_node)
        #/
        # make new branch node
        new_branch_node = Clade(branch_length=0)
        #/
        # add leaf as child to new branchnode
        new_branch_node.clades.append(leaf_node)
        #/
        # add new branch node under parent
        leaf_parent.clades.append(new_branch_node)
        #/
    #/
    # when done, re-set the old root (i.e. a "no-root object". I cannot recreate this when setting tree.root = Clade(), even though orig_root is a "Clade()")
    print('Re-setting old root')
    tree.root = orig_root
    #/
##/
## Assign branch-node names (they are empty on import)
branch_node_enums = 0
for branch_node in tree.get_nonterminals():
    if branch_node.name == None or branch_node.name == '':
        branch_node.name = branch_node_basename+str(branch_node_enums)
        branch_node_enums += 1
##/
## Set outgroup (if supplied)
if outgroup_dataset:
    tree.root_with_outgroup(outgroup_dataset)
##/
## Remove outgroup (if specified)
if remove_outgroup:
    tree.prune(outgroup_dataset)
##/
## Get datasets (tree leaves)
datasets = set()
for leaf_node in tree.get_terminals():
    datasets.add(leaf_node.name)
##/
## Get branch node names of tree
branchnodes_names = set() # keep track of branch nodes
for branch_node in tree.get_nonterminals():
    branchnodes_names.add(branch_node.name)
##/
## Parse metadata from custom file
# parse columns
db_columns = []
for db_column in db_columns_raw.split(','):
    db_columns.append(db_column)
#/
# import metadata
accession_metadata = {}
if metadata_file:
    accession_metadata = global_functions.parse_metadata_file(metadata_file,accessions_to_import=datasets,discard_id_column=True,cols_to_import=db_columns)
print(f'Parsed leaf metadata, N={len(accession_metadata)}')
if len(accession_metadata) == 0:
    print('No metadata found. If you used custom columns, make sure that they can be exact-matched in the provided metadata. Terminating now!')
    sys.exit()
#/
# ensure that all columns were matched
verify_imported_columns(expected_cols=db_columns,imported_metadata=accession_metadata,errormessage_description='column metadata')
#/
##/
## Check if user supplied columns for branch-nodes to parse from metadata
# parse columns
db_columns_branches = []
if db_columns_branches_raw != None:
    for db_column in db_columns_branches_raw.split(','):
        db_columns_branches.append(db_column)
#/
# import metadata
branch_accession_metadata = {}
if db_columns_branches:
    if metadata_file:
        branch_accession_metadata = global_functions.parse_metadata_file(metadata_file,accessions_to_import=datasets,discard_id_column=True,cols_to_import=db_columns_branches)
    print(f'Parsed branch metadata, N={len(accession_metadata)}')
    if len(branch_accession_metadata) == 0:
        print('No metadata found for branches. If you used custom columns, make sure that they can be exact-matched in the provided metadata. Terminating now!')
        sys.exit()
#/
# ensure that all columns were matched
verify_imported_columns(expected_cols=db_columns_branches,imported_metadata=branch_accession_metadata,errormessage_description='branch metadata')
#/
##/

## Check if there is a cansnper path to expand
if cansnper_column != None:
    if cansnper_column in db_columns:
        print(f'Expanding CanSNPs in leaf metadata')
        accession_metadata,db_columns = metadata_parse_canSNPer_column(accession_metadata,db_columns,cansnps_to_use)
    if cansnper_column in db_columns_branches:
        print(f'Expanding CanSNPs in branch metadata')
        branch_accession_metadata,db_columns_branches = metadata_parse_canSNPer_column(branch_accession_metadata,db_columns_branches,cansnps_to_use)
##/
###/

### Assign node datasets
## Assign node datasets and keep track of parental nodes
datasets_nodes = {} # dataset -> nodes
branchNodes_parentNodes = {} # node -> parentNode
leafNodes_parentNodes = {} # leafNode -> parentNode
for dataset in datasets:
    dataset_ancestor_path = tree.get_path(dataset)
    for node_enum,_ in enumerate(dataset_ancestor_path):
        # get node
        node_name = dataset_ancestor_path[node_enum].name
        #/
        # get node parent (if possible)
        node_parent_name = None
        if node_enum > 0:
            node_parent_name = dataset_ancestor_path[node_enum-1].name
        #/
        # check if current node is a leaf (dataset). If so, then save its parent
        if node_name == dataset:
            leafNodes_parentNodes[node_name] = node_parent_name
        #/
        # save node datasets
        if not node_name in datasets_nodes:             datasets_nodes[node_name] = set()
        datasets_nodes[node_name].add(dataset)
        #/
        # save node parent
        if node_parent_name != None:
            branchNodes_parentNodes[node_name] = node_parent_name
        #/
##/

## Restructure: branchNode -> childNode
branchNodes_childNodes = {} # node -> childNode
for branchNode,parentNode in branchNodes_parentNodes.items():
    if not parentNode in branchNodes_childNodes:            branchNodes_childNodes[parentNode] = set()
    branchNodes_childNodes[parentNode].add(branchNode)
##/
###/

### Restructure: Metadata_column -> values -> datasets
metadata_cols_vals_datasets = {}
for dataset,metadata in accession_metadata.items():
    for column,value in metadata.items():
        if not column in metadata_cols_vals_datasets:               metadata_cols_vals_datasets[column] = {}
        if not value in metadata_cols_vals_datasets[column]:       metadata_cols_vals_datasets[column][value] = set()
        metadata_cols_vals_datasets[column][value].add(dataset)
print(f'There are N={len(metadata_cols_vals_datasets)} imported leaf metadata columns')
###/

### Use metadata to determine branch-node stuff
## Restructure: Metadata_column -> values -> datasets
branch_metadata_cols_vals_datasets = {}
for dataset,metadata in branch_accession_metadata.items():
    for column,value in metadata.items():
        if not column in branch_metadata_cols_vals_datasets:               branch_metadata_cols_vals_datasets[column] = {}
        if not value in branch_metadata_cols_vals_datasets[column]:       branch_metadata_cols_vals_datasets[column][value] = set()
        branch_metadata_cols_vals_datasets[column][value].add(dataset)
print(f'There are N={len(branch_metadata_cols_vals_datasets)} imported branch metadata columns')
##/

## Compute metadata at branchnodes
column_vals_branchnodes = {} # column -> values -> "branchnode where all branchnode_datasets have a metadata value"
for column in branch_metadata_cols_vals_datasets:
    # check if user supplied specific columns to use for branches only
    if db_columns_branches_raw != None:
        if not column in db_columns_branches: continue
    #/
    for value in branch_metadata_cols_vals_datasets[column]:
        datasets_with_col_val = branch_metadata_cols_vals_datasets[column][value]
        
        for branchnode,branchnode_datasets in datasets_nodes.items():
            if not branchnode in branchnodes_names: continue # skip if not a branch node
            if branchnode_datasets == datasets_with_col_val: # check if sets of datasets are identical between metadata and tree branch
                if not column in column_vals_branchnodes:               column_vals_branchnodes[column] = {}
                if not value in column_vals_branchnodes[column]:        column_vals_branchnodes[column][value] = set()
                column_vals_branchnodes[column][value].add(branchnode)
##/
###/

##### OUTPUTS
### Check if output FlexTaxD-formatted files (nodes renamed based on metadata)
if flextaxd_outfiles_path != None:
    print('Begin constructing FlexTaxD output')
    # Make (deep)copy of tree and modify node names for FTD-output
    import copy
    ftd_tree = copy.deepcopy(tree)
    #/
    #@ BRANCH NODES: Rename branch nodes
    print('Begin determine branch node names')
    # Determine new names for branch nodes based on metadata (in order, to get naming controllable. e.g.: species->subspecies->canSNP)
    branch_nodes_rename = {} # old name -> new name
    branch_nodes_ranks = {} # branchnode original name -> "rank"
    for db_column in db_columns_branches:
        if db_column in column_vals_branchnodes:
            for column_val in column_vals_branchnodes[db_column]:
                for branch_node in column_vals_branchnodes[db_column][column_val]:
                    ## Check if user wanted to use certain columns as rank for FTD
                    if flextaxd_rank_columns != None:
                        if db_column in flextaxd_rank_columns:
                            # save column name as rank at node
                            if not branch_node in branch_nodes_ranks:       branch_nodes_ranks[branch_node] = []
                            branch_nodes_ranks[branch_node].append(db_column)
                            #/
                            # check if user does not want to use rank columns as name annotation
                            if flextaxd_rank_columns_strip:
                                continue
                            #/
                    ##/
                    # init new branch node name, append with values
                    if not branch_node in branch_nodes_rename:          branch_nodes_rename[branch_node] = []
                    #/
                    # append value (if bool, then save the column name)
                    val_to_save = None
                    if type(column_val) == bool:
                        if column_val == True:
                            val_to_save = db_column
                    else:
                        if len(column_val) > 0: # only save if not an empty string
                            val_to_save = column_val
                    
                    #@ check if user wants to discard "None" assignments
                    if flextaxd_skip_missing_metadata == True:
                        if val_to_save == None:
                            continue
                    #@/
                    
                    branch_nodes_rename[branch_node].append(val_to_save)
                    #/
                    # check if user wants to apply rename to all downstream branch nodes
                    if flextaxd_apply_downstream:
                        # get the tree node
                        branch_node_treeObject = ftd_tree.find_any(branch_node)
                        #/
                        # get all downstreams (branch-nodes only)
                        downstream_nodes = branch_node_treeObject.get_nonterminals()
                        if downstream_nodes[0].name == branch_node:         downstream_nodes = downstream_nodes[1:] # remove current branchnode from "downstreams" (self is returned by biopython function and is first in list)
                        #/
                        print(f'Applying to downstream nodes: {val_to_save}')
                        # rename each downstream
                        for downstream_node in downstream_nodes:
                            downstream_node_name = downstream_node.name
                            if not downstream_node_name in branch_nodes_rename:      branch_nodes_rename[downstream_node_name] = []
                            branch_nodes_rename[downstream_node_name].append(val_to_save)
                        #/
                    #/
    
    branch_names_used = {} # keep track of which names were used to never put a duplicate name in case metadata is identical across nodes. Link to original tree node name
    branch_node_basename_sep = '_' # separator to use between branch basename and enumerate
    for branch_node in ftd_tree.get_nonterminals():
        branch_node_name = branch_node.name
        if branch_node_name in branch_nodes_rename:
            # determine new name at branch
            rename_to = branch_node_basename+branch_node_basename_sep+'_'.join(map(str,branch_nodes_rename[branch_node_name]))
            #/
            # set dummy name if empty
            if rename_to == branch_node_basename+branch_node_basename_sep:
                rename_to = branch_node_basename
            #/
            # remove "branch_node_basename" if there was metadata assigned
            if len(branch_nodes_rename[branch_node_name]) > 0 and not flextaxd_keep_node_basenames:
                rename_to = rename_to.replace(branch_node_basename+branch_node_basename_sep,'') # In e.g. "branchNode:metadata1_metadata2", replace "branchNode" with ""
            #/
            # if this rename already exists, append an enumerate until its unique
            enums_tested = 2 # start from 2 (1-based; the first node will be without _<enum>)
            try_rename_to = rename_to
            while try_rename_to in branch_names_used:
                try_rename_to = rename_to + f'_{enums_tested}'
                enums_tested += 1
            rename_to = try_rename_to
            #/
            # save
            branch_names_used[rename_to] = branch_node_name
            #/
    #/
    # Do branch-node renaming
    print('Renaming branches now')
    for branch_node in ftd_tree.get_nonterminals():
        branch_node_name = branch_node.name
        for rename_to,rename_from in branch_names_used.items():
            if rename_from == branch_node_name:
                # INFOprint
                print(f'Renamed branch {branch_node_name} -> {rename_to}')
                #/
                # set new name at node
                branch_node.name = rename_to
                #/
    #/
    #@/BRANCH NODES
    
    #@ LEAF NODES: rename leaf nodes
    print('Begin determine leaf node names')
    # Determine new names for leaf nodes
    leaf_names_used = {} # keep track of which names were used to never put a duplicate name in case metadata is identical across multiple nodes.
    leaf_node_basename_sep = '_'
    for leaf_node in ftd_tree.get_terminals():
        leaf_node_name = leaf_node.name
        if leaf_node_name in accession_metadata:
            metadata = accession_metadata[leaf_node_name]
            
            # get metadata to use in new name
            leaf_new_name = []
            for db_column in db_columns:
                if db_column in metadata:
                    column_val = metadata[db_column]
                    # skip value if empty
                    if column_val == '': continue
                    #/
                    # append value (if bool, then save the column name)
                    if type(column_val) == bool:
                        if column_val == True:
                            leaf_new_name.append(db_column)
                    else:
                        leaf_new_name.append(column_val)
                    #/
            #/
            # determine new name
            rename_to = leaf_node_basename+leaf_node_basename_sep+'_'.join(map(str,leaf_new_name))
            #/
            # set dummy name if empty
            if rename_to == leaf_node_basename+leaf_node_basename_sep:
                rename_to = leaf_node_basename
            #/
            # remove "branch_node_basename" if there was metadata assigned
            if len(leaf_new_name) > 0 and not flextaxd_keep_node_basenames:
                rename_to = rename_to.replace(leaf_node_basename+leaf_node_basename_sep,'') # In e.g. "branchNode:metadata1_metadata2", replace "branchNode" with ""
            #/
            # if this rename already exists, append an enumerate until its unique
            enums_tested = 2 # start from 2 (1-based; the first node will be without _<enum>)
            try_rename_to = rename_to
            while try_rename_to in leaf_names_used:
                try_rename_to = rename_to + f'_{enums_tested}'
                enums_tested += 1
            rename_to = try_rename_to
            #/
            # save
            leaf_names_used[rename_to] = leaf_node_name
            #/
    #/
    # Do leaf-node renaming
    print('Renaming leafs now')
    for leaf_node in ftd_tree.get_terminals():
        leaf_node_name = leaf_node.name
        for rename_to,rename_from in leaf_names_used.items():
            if rename_from == leaf_node_name:
                # INFOprint
                print(f'Renamed leaf {leaf_node_name} -> {rename_to}')
                #/
                # set new name at node
                leaf_node.name = rename_to
                #/
    #/
    #@/LEAF NODES
    ##/
    ## Check if user wants to load additional genomes (i.e. genomes_map with genomes that were redundant in RepGenR)
    additional_genomes = {} # genome_name -> accession
    if flextaxd_additional_genomes != None:
        print(f'Additional genomes supplied: {flextaxd_additional_genomes}')
        additional_genomes_including_dereplicated = 0
        with open(flextaxd_additional_genomes,'r') as f:
            for line in f:
                # parse line
                line = line.strip('\n')
                line = line.split('\t')
                #/
                # parse data
                genome_accn,genome_name = line
                #/
                # save data
                if not genome_name in additional_genomes:           additional_genomes[genome_name] = set()
                additional_genomes[genome_name].add(genome_accn)
                
                additional_genomes_including_dereplicated += 1
                #/
        print(f'Number of leafs imported N={len(additional_genomes)}, total number of genomes including dereplicated N={additional_genomes_including_dereplicated}')
    #/
    # rename imported leaf
    additional_genomes_renamed = {}
    for genome_name,accns in additional_genomes.items():
        # get new name, determined above
        genome_name_new = None
        for rename_to,rename_from in leaf_names_used.items():
            if rename_from == genome_name:
                genome_name_new = rename_to
        #/
        # save at new name
        additional_genomes_renamed[genome_name_new] = accns
        #/
    #/
    ##/
    ## Make FTD-output
    # init outdir (warn user if directory already exists)
    if os.path.exists(flextaxd_outfiles_path):
        raw_input = input(f'WARNING: Directory {flextaxd_outfiles_path} already exists. Do you want to write files here? (y/n): ')
        if not(raw_input and raw_input.lower() in ('yes','y',)):
            print(f'Answer given: {raw_input}, will terminate now!')
            sys.exit()
    if not os.path.exists(flextaxd_outfiles_path):        os.makedirs(flextaxd_outfiles_path)
    #/
    # Write tree structure (as parent-child, flextaxd format. AKA "tree2tax" argument "--mod_file")
    parent_childs_written = []
    print('Output tree relations now')
    with open(flextaxd_outfiles_path+'/'+'tree_parent_child_relations.tsv','w') as nf:
        # write header
        tmp_header = ['parent','child']
        
        if flextaxd_rank_columns != None:           tmp_header.append('rank') # if user wants to use metadata to determine rank, then add this column
        
        nf.write('\t'.join(tmp_header)+'\n')
        #/
        # write rows
        for leaf in ftd_tree.get_terminals():
            # get path with node names
            path = [] # root -> branchnode_N -> ... -> leaf
            if flextaxd_root_name_linker != None:      path.append(flextaxd_root_name_linker) # if user specified a root name to use
            for node in ftd_tree.get_path(leaf):
                path.append(node.name)
            #/
            # save parent-child relations
            for path_enum,parent in enumerate(path):
                if len(path)-1 == path_enum: continue # skip if we cannot get a child for current parent, i.e. out of boundary, we are finished
                child = path[path_enum+1]
                tmp_write = [parent,child]
                
                # check if user wanted to add a rank from metadata (applies to "child")
                if flextaxd_rank_columns != None:
                    rank_to_write = ''
                    
                    # check if this node should have a rank
                    if child in branch_nodes_ranks:
                        rank_to_write = branch_nodes_ranks[child] # check for non-renamed branch nodes
                    if child in branch_names_used:
                        old_parent_name = branch_names_used[child] # check for renamed branch nodes, get the original branch name
                        if old_parent_name in branch_nodes_ranks:
                            rank_to_write = branch_nodes_ranks[ old_parent_name ]
                    #/
                    # add rank to writeArr
                    if len(rank_to_write) > 1:
                        print(f'WARNING: Multiple ranks found for relation: {parent} {child}. This indicates that your metadata ambiguously describe your input (not expected for taxonomic ranks)')
                    
                    rank_to_write = ','.join(rank_to_write)
                    tmp_write.append(rank_to_write)
                    #/
                #/
                
                # only write this relationship if it was not already written in a previous path
                if not tmp_write in parent_childs_written:
                    # write row
                    nf.write('\t'.join(map(str,tmp_write))+'\n')
                    #/
                    # keep track of written relationships
                    parent_childs_written.append(tmp_write)
                    #/
                #/
            #/
        #/
    #/
    # Write nodes/genomes (AKA "genomeid2taxid" argument "--genomeid2taxid")
    print('Output genomes info now')
    with open(flextaxd_outfiles_path+'/'+'genome_id_map.tsv','w') as nf:
        # headerless file
        #/
        ## Write rows
        rows_written = [] # keep track of which rows were written. When user imported additional genomes, do not write leaf-nodes from tree twice
        # write rows for "tree leafs"
        for node_name,original_name in leaf_names_used.items():
            # parse accession id from original name (expected at <family>_<genus>_<species>_<GCx>_<number>.<v>)
            regex_pattern = r"GC[A|F]_\d{9}\.\d"
            match = re.search(regex_pattern,original_name)
            if match:
                matched_string = match.group()
                stripped_string = re.sub(f".*({regex_pattern}).*", r"\1", matched_string)
                accn = stripped_string # should be formatted as GCX_123456789.1
                
                writeArr = [accn,node_name]
                if not writeArr in rows_written:
                    nf.write('\t'.join(map(str,writeArr))+'\n')
                    rows_written.append(writeArr)
            #/
        #/
        # write rows loaded from "additional genomes"
        for node_name,accns in additional_genomes_renamed.items():
            for accn in accns:
                # format: col1=accession_number, col2=node_name
                writeArr = [accn,node_name]
                if not writeArr in rows_written:
                    nf.write('\t'.join(map(str,writeArr))+'\n')
                    rows_written.append(writeArr)
                #/
        #/
        ##/
    #/
    ##/
###/

### Do plotting
if plot_output_file_path:
    # init figure
    fig,(ax,ax_annotation) = plt.subplots(ncols=2, figsize=(30, 18), gridspec_kw={'width_ratios':[4,1]})
    #/
    # remove all stuff from ax_annotation
    if 1:
        for spine in ax_annotation.spines.values():
            spine.set_visible(False)
        ax_annotation.set_xticks([])
        ax_annotation.set_yticks([])
        ax_annotation.set_xticklabels([])
        ax_annotation.set_yticklabels([])
    #/
    # remove space between ax and ax_annotation
    plt.subplots_adjust(wspace=-1)
    #/
    # plot tree with biopython-Phylo
    Phylo.draw(tree, axes=ax,
               label_func=lambda leaf: leaf.name, # Replace Biopython label function that returns the full name (default, clip >40))
               do_show=False)
    #/
    
    #@@@@@@@ SECTION: apply metadata visualisations
    
    ## determine scale to use for plotting (text and elements)
    orig_font_size = plt.rcParams['font.size']
    plot_scaler = 1
    if 1:
        if len(datasets) > 200:
            plot_scaler = 0.3
        elif len(datasets) > 100:
            plot_scaler = 0.35
        elif len(datasets) > 50:
            plot_scaler = 0.4
        elif len(datasets) > 40:
            plot_scaler = 0.6
        elif len(datasets) > 30:
            plot_scaler = 0.8
        elif len(datasets) > 20:
            plot_scaler = 0.9
    ##/
    
    ## Get position of branch-nodes and leaf-nodes
    # get position of each dataset text label
    datasets_textlabel_pos = {} # dataset -> pos
    poses_vals_without_outgroup = []
    datasets_yposes = []
    for child in ax.get_children():
        if isinstance(child, plt.Text):
            label_text = child.get_text().lstrip(' ') # strip bullshit space infront of every text label
            # skip label if it is a branch-node or reference_sequence or is empty
            if label_text in branchnodes_names or label_text == '': continue
            #/
            # get pos and save
            position = child.get_position()
            datasets_textlabel_pos[label_text] = position
            #/
            # save coord of x offset (update: and Y-coordinate)
            if not (outgroup_dataset != None and label_text == outgroup_dataset):
                poses_vals_without_outgroup.append(position[0])
                datasets_yposes.append(position[1])
            #/
    #/
    # get position of branch node text labels
    branchnode_textlabel_pos = {} # branchnode -> pos
    for child in ax.get_children():
        if isinstance(child, plt.Text):
            label_text = child.get_text().lstrip(' ') # strip bullshit space infront of every text label
            if label_text in branchnodes_names:
                # get pos and save
                position = child.get_position()
                branchnode_textlabel_pos[label_text] = position
                #/
    #/
    # get all x_poses and sort (highest X-value last, for placement of dots)
    x_vals_sorted = sorted(poses_vals_without_outgroup) #sorted([x for x,y in datasets_textlabel_pos.values()])
    #/
    # re-position textlabels so that they (1) align and (2) position behind dot and (3) change font size
    datasets_textlabel_pos_repositioned = {}
    for child in ax.get_children():
        if isinstance(child, plt.Text):
            label_text = child.get_text().lstrip(' ') # strip bullshit space infront of every text label
            if label_text in datasets and (outgroup_dataset == None or label_text != outgroup_dataset):
                old_pos = child.get_position()
                child.set_position([x_vals_sorted[-1]+x_vals_sorted[-1]*0.04,old_pos[1]])
                child.set_fontsize(orig_font_size*(plot_scaler**8))
                
                # get pos and save
                position = child.get_position()
                datasets_textlabel_pos_repositioned[label_text] = position
                #/
    #/
    # compute average height between datasets (will plot ½ height above and below area-plots to make it look nicer)
    datasets_avg_plot_distance = (max(datasets_yposes)-min(datasets_yposes))/len(datasets_yposes)
    #/
    ##/
    
    ## Do plotting in annotation-axes
    # init-plot invisible line to create axes boundaries
    ax_annotation.plot([0,0],ax.get_ylim(),alpha=0)
    #/
    # init xoffset and set step size
    x_offset = 1*(plot_scaler**2)
    x_offset_init = x_offset
    x_offset_steps = 2*plot_scaler # increase offset by this step size
    #/
    #@ plot dots for metadata values (BRANCH-MATCHED DATASETS ONLY)
    if 0:
        for column in column_vals_branchnodes:
            for val in column_vals_branchnodes[column]:
                for branch_node in column_vals_branchnodes[column][val]:
                    branch_datasets = datasets_nodes[branch_node]
                    
                    for dataset in branch_datasets:
                        dataset_ycoord = datasets_textlabel_pos[dataset][1]
                        ax_annotation.scatter(x_offset,dataset_ycoord,color='red',s=100*(plot_scaler**2))
            x_offset += x_offset_steps # for each column, increase x_offset for next column/thing to plot
    #@/
    #@ plot dots for metadata values (ALL METADATA VALUES; NO GROUPING DONE)
    valtypes_colors = ['red', 'green', 'blue', 'gold', 'purple', 'crimson','grey']
    for column in metadata_cols_vals_datasets:
        # plot dot per value
        datasets_ycoords = []
        label_description_text = []
        for value_enum,value in enumerate(sorted(metadata_cols_vals_datasets[column])): # sort values so that e.g. values with TRUE/FALSE will have same color across different columns
            # plot scatters at datasets with this value
            for dataset in metadata_cols_vals_datasets[column][value]:
                dataset_ycoord = datasets_textlabel_pos[dataset][1]
                datasets_ycoords.append(dataset_ycoord)
                
                tmp_color = valtypes_colors[value_enum%len(valtypes_colors)]
                ax_annotation.scatter(x_offset,dataset_ycoord,color=tmp_color,s=100*(plot_scaler**2))
            #/
            # append description to label about this value
            label_description_text.append(f'{tmp_color}={value}')
            #/
        #/
        # write label
        tmp_label = column+' ('+', '.join(label_description_text)+')'
        plot_max_colum_label_length = 150
        if len(tmp_label) > plot_max_colum_label_length:
            print(f'Plot label for column {column} is very long. Will clip it now to {plot_max_colum_label_length} characters')
        ax_annotation.text(x_offset,min(datasets_ycoords)-1,tmp_label[:plot_max_colum_label_length],rotation=45,fontsize=orig_font_size*plot_scaler)
        #/
        x_offset += x_offset_steps # for each column, increase x_offset for next column/thing to plot
    #@/
    # plot dataset labels
    for dataset,textlabel_pos in datasets_textlabel_pos.items():
        ax_annotation.text(x_offset,textlabel_pos[1],dataset,fontsize=orig_font_size*plot_scaler)
    x_offset += x_offset_steps
    #/
    ##/
    
    ## Try to make visual guidelines for datasets<->leafnode connection
    # for each dataset, get position of text label and parent branch node. Use these coordinates to draw visual guideline
    for dataset in datasets:
        # get position of textlabel for dataset (look in "repositioned" first, then in "original". This makes sense for when the outgroup is not repositioned.)
        textlabel_pos = [0,0]
        if dataset in datasets_textlabel_pos_repositioned:
            textlabel_pos = datasets_textlabel_pos_repositioned[dataset]
        else:
            textlabel_pos = datasets_textlabel_pos[dataset]
        #/
        # get parent branchnode
        dataset_parent = leafNodes_parentNodes[dataset]
        #/
        # get parent branchnode position
        branchlabel_pos = [0,0] # the dataset with no parent will being in 0,0
        if dataset_parent != None:
            branchlabel_pos = branchnode_textlabel_pos[dataset_parent] # update position if dataset had a parent
        #/
        # draw guideline
        ax.plot([branchlabel_pos[0],ax.get_xlim()[1]],[textlabel_pos[1],textlabel_pos[1]],color='grey',linestyle=':',alpha=0.3,zorder=0) # from branchnode to max X value
        ax_annotation.plot([x_offset_init,x_offset-x_offset_steps],[textlabel_pos[1],textlabel_pos[1]],color='grey',linestyle=':',alpha=0.3,zorder=0) # from xlim start to last used x_offset, at branch_label y value
        #/
    #/
    ##/
    
    ## Do area-plotting in left axes
    area_colors = ['gray','green','red','orange','blue']
    num_area_draws = 0
    ax_texts_to_keep = []
    branch_nodes_numTexts = {} # keep track of texts at each branchnode to offset them when multiple texts are plotted
    for column in column_vals_branchnodes:
        for val in column_vals_branchnodes[column]:
            for branch_node in column_vals_branchnodes[column][val]:
                branch_datasets = datasets_nodes[branch_node]
                
                # get branch node x coord (use plot right border/xlim as xend)
                area_xstart = branchnode_textlabel_pos[branch_node][0]
                area_xend = ax.get_xlim()[1]
                #/
                # get dataset coords ycoords
                datasets_ycoords = []
                for dataset in branch_datasets:
                    dataset_pos = datasets_textlabel_pos[dataset]
                    datasets_ycoords.append(dataset_pos[1])
                area_ystart = min(datasets_ycoords)
                area_yend = max(datasets_ycoords)
                #/
                # plot area
                tmp_height = area_yend-area_ystart
                tmp_width = area_xend-area_xstart
                
                rect = patches.Rectangle((area_xstart,area_yend), tmp_width, -tmp_height, linewidth=0, facecolor=area_colors[num_area_draws%(len(area_colors)-1)],alpha=0.3) # x,y is bottom left of rectangle. next x2 numbers are width and height
                ax.add_patch(rect)
                
                num_area_draws += 1
                #/
                ## add text
                # get offset/init offset/iterate offset
                if not branch_node in branch_nodes_numTexts:            branch_nodes_numTexts[branch_node] = 0
                text_x_offset = branch_nodes_numTexts[branch_node]
                branch_nodes_numTexts[branch_node] += 1
                #/
                # determine coords for text
                text_x = area_xstart + (area_xend-area_xstart)*(0.1*(1+text_x_offset))
                text_y = (area_yend+area_ystart)/2
                #/
                # determine value to text
                text_to_plot = None
                if type(val) == bool and val == True: # if column was a bool, then plot the column name
                    text_to_plot = column
                else:
                    text_to_plot = val # else, plot the value
                #/
                # do text-plot
                txt = ax.text(text_x,text_y,text_to_plot,fontsize=orig_font_size*(plot_scaler*4),color='black')
                ax_texts_to_keep.append(txt)
                #/
                ##/
    ##/
    
    
    #@@@@@@@/
    
    ## cleanup plot
    # remove branch labels
    if 1 and 'cleanup branch labels':
        for text in ax.texts:
            #if text.get_text().lstrip(' ') in branchnodes_names:
            if not text in ax_texts_to_keep:
                text.remove()
    #/
    # cleanup display
    ax.spines['right'].set_visible(False) # remove right plot border
    ax.spines['top'].set_visible(False) # remove top plot border
    ax.set_xlabel(None) # remove ylabel
    ax.set_ylabel(None) # remove ylabel
    ax.yaxis.set_ticks([]) # remove yticks
    #/
    # set plot xand ylimits
    ax.set_xlim([ax.get_xlim()[0],max([xval for xval,yval in datasets_textlabel_pos.values()])]) # use start of xcoord for furthest right dataset label as end of plot
    #ax.set_ylim([ax.get_ylim()[0],ax.get_ylim()[1]-10])
    
    ax_annotation.set_xlim([0,10+x_offset])
    ax_annotation.set_ylim(ax.get_ylim())
    #/
    # apply layout
    plt.tight_layout()
    #/
    # Show the plot
    #plt.show()
    #/
    
    # Save as pdf
    if 1:
        png_output_path = plot_output_file_path
        print(f'Dumping plot as PDF: {png_output_path}')
        plt.savefig(png_output_path)
    #/
    
    # Convert the plot to interactive HTML (check if it is possible to get a richer interface)
    if 0 and 'needs improvments if it is going to be used':
        if plot_output_file_path != None and 'works':
            import mpld3
            html_str = mpld3.fig_to_html(fig)
            print(f'Dumping plot as html: {plot_output_file_path}')
            with open(plot_output_file_path,'w') as nf:
                nf.write(html_str)
    #/
    
    # output plot as html
    if 0 and 'does not work with lines plotted. possible to convert to supported format? (LineCollection to ... native line?)':
        import plotly.tools as tls
        plotly_fig = tls.mpl_to_plotly(fig)
        plotly_fig.write_html("interactive_plot_plotly.html")
    #/
###/
#####/OUTPUTS