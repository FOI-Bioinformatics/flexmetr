#!/usr/bin/env python3

import os
import sys
import argparse
try:            import ete3
except:         sys.exit('Unable to import ETE3 package. Please make sure it has been installed.')


import global_functions

### Parse input arguments
# setup
argparser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
argparser.add_argument('-i','--input',required=True,help='Path to input file or "-" to read from stream')
argparser.add_argument('-o','--output',required=False,default='-',help='Path to output file or "-" to print stream')

# argparser.add_argument('-d','--database','-db','--db',required=False,default='flexmetr/db.sql',help='Path to database (default:flexmetr/db.sql)') # not implemented yet

argparser.add_argument('-c','--column','--col','-col','--columns','--cols','-cols',required=False,help='Column in metadata to format output. Multiple columns may be specified by comma (e.g.: family,genus,species)')

argparser.add_argument('--debranch',required=False,type=float,default=None,help='Cuts branches below input distance (default: None)')

argparser.add_argument('--collapse',required=False,action='store_true',default=None,help='If specified, will collapse datasets at the same branch and select the dataset with shortest distance')

argparser.add_argument('--metadata_file',required=False,default=None,help='Path to custom metdata file')
argparser.add_argument('--metadata_file_sep',required=False,default='\t',help='Separator to use in custom metadata file [default: tab]')
argparser.add_argument('--metadata_file_accession',required=False,default='\t',help='Column in custom metadata file that holds the accession number [default: first/0]')

argparser.add_argument('--IDE_format_names',required=False,action='store_true',help='For developing purposes: when specified, will format node-names as "family_genus_species"')
#/
# parse input
args = argparser.parse_args()

input_file = args.input
output_file = args.output

metadata_db = args.database
db_columns = args.column

debranch_thresh = args.debranch
collapse_leafs = args.collapse

metadata_file = args.metadata_file
metadata_file_sep = args.metadata_file_sep
metadata_file_accession = args.metadata_file_accession

format_names = args.IDE_format_names
#/
###/

## Parse metadata from custom file
if metadata_file:
    accession_metadata = global_functions.parse_metadata_file(metadata_file)
##/

## Read input (stdin or file)
if input_file == '-':
    input_string = sys.stdin.read()
else:
    with open(input_file,'r') as f:
        input_string = f.read()
input_string = input_string.strip('\n')
##/

## Parse input into ete3 tree structure
tree = ete3.Tree(input_string)

if 0 and 'IDE, set root':
    try:
        root_name = 'GCF_003697165.2'
        print('[IDE] Attempting to set root: '+root_name)
        root_node = tree.search_nodes(name=root_name)[0]
        tree.set_outgroup(root_node)
    except:
        print('[IDE] Failed to set root')
##/

## Check if debranch branches below certain threshold
if debranch_thresh:
    # For every leaf-node, calculate distance to parent branch-nodes and re-position them at user-specified cutoff
    for node in tree.traverse():
        if node.is_leaf():
            # Get ancestors to node
            ancestors = node.get_ancestors()
            #/
            
            # IDE: check so ancestors are always sorted by distance
            ancestors_dist = [node.get_distance(ancestor) for ancestor in ancestors]
            if not ancestors_dist == sorted(ancestors_dist): print('IDE: Ancestors were not sorted by distance!')
            #/
            
            # Find the ancestor with acceptable distance: Check if distance to this ancestor is greater than the user-specified cutoff
            new_ancestor = None
            for ancestor in ancestors:
                leaf_ancestor_dist = node.get_distance(ancestor)
                
                
                if leaf_ancestor_dist > debranch_thresh:
                    new_ancestor = ancestor
                    break
            #/
            
            # Set node at new ancestor
            if new_ancestor != node._up:
                new_ancestor_dist = node.get_distance(new_ancestor)
                node.detach()
                new_ancestor.add_child(node)
                node.dist = new_ancestor_dist
            #/
    #/
    # Join branch-nodes that are single-connected
    for node in tree.traverse():
        # Get ancestors to node
        ancestors = node.get_ancestors()
        #/
        # Lift node to first ancestor that have at least two children
        new_ancestor = None
        for ancestor in ancestors:
            num_children = len(ancestor.get_children())
            if num_children >= 2:
                new_ancestor = ancestor
                break
        
        if new_ancestor != node._up:
            new_ancestor_dist = node.get_distance(new_ancestor)
            node.detach()
            new_ancestor.add_child(node)
            node.dist = new_ancestor_dist
        #/
    #/
    # Clean tree from dead branches (branch-nodes with no leaf-nodes of datasets)
    iterations = 0
    deletion_made = False
    while iterations < 1000000 and (deletion_made or iterations == 0):
        deletion_made = False
        for node in tree.traverse():
            if len(node.get_children()) == 0 and not node.name:
                node.delete()
                deletion_made = True
                break
        iterations += 1
    #/
##/
## Check if debranch based on database column
accession_metadata = None
if db_columns:
    ## Get accessions to import from DB
    accessions_to_import = set()
    for node in tree.traverse():
        name = node.name
        if name:
            accession = global_functions.getAccessions(name,return_first=True)
            
            if accession:
                accessions_to_import.add(accession)
    ##/
    
    ## Get columns to parse from database
    db_keys = db_columns.split(',')
    # remove preceding spaces
    for enum,_ in enumerate(db_keys):
        db_keys[enum] = db_keys[enum].replace(' ','')
    #/
    ##/
    
    ## Parse metadata from DB
    accession_metadata = {}
    ##/
    
    ## Check if parse metadata from custom file
    if metadata_file:
        cols_to_import = set()#db_keys # parse only the keys specified by user input from DB
        accession_metadata = global_functions.parse_metadata_file(metadata_file,accessions_to_import=accessions_to_import,cols_to_import=cols_to_import)
    ##/
    
    def get_node_classi(inp_node):
        name = inp_node.name
        accession = global_functions.getAccessions(name,return_first=True)
        
        # Save leaf-classifications at node
        vals = []
        for db_key in db_keys:
            value = accession_metadata[accession][db_key]
            if not value:
                value = 'NA'
            vals.append(value)
            
        vals_ID = '||'.join(map(str,vals))
        
        return vals_ID
    
    def get_node_children_classis_counts(inp_node):
        leafClassis_counts = {}
        childrens = inp_node.get_children()
        for child in childrens:
            if child.is_leaf() and child.name:
                name = child.name
                accession = global_functions.getAccessions(name,return_first=True)
                
                child_classi_ID = get_node_classi(child)
                if not child_classi_ID  in leafClassis_counts:     leafClassis_counts[child_classi_ID] = set()
                leafClassis_counts[child_classi_ID].add(accession)
                
        return leafClassis_counts
    
    # For each branch-node, raise its children if all leaf-nodes have the same classification
    iterations = 0
    deletion_made = False
    del_nodes = []
    while iterations < 1000000 and (deletion_made or iterations == 0):
        deletion_made = False
        for node in sorted(tree.traverse(),key=lambda x: x.get_distance(root_node)):
            if not node.is_leaf():
                # count the number of leaf-nodes per classification
                leafClassis_counts = get_node_children_classis_counts(node)
                #/
                # Sort children (datasets + branches) and leaf-children (datasets)
                childrens = node.get_children()
                childrens_leafs = []
                for child in childrens:
                    if child.is_leaf() and child.name:
                        childrens_leafs.append(child)
                #/
                # Check if only one classificaiton exist at node, if so then reposition leaf-nodes at parent
                if len(leafClassis_counts) == 1 and childrens == childrens_leafs:
                    delete_node = False
                    for child in childrens:
                        
                        if 0 :
                            ## IDE
                            name = child.name
                            accession = global_functions.getAccessions(name,return_first=True)
                            fam,gen,spe,accn = accession_metadata[accession]['family'],accession_metadata[accession]['genus'],accession_metadata[accession]['species'],accession
                            child_name = fam+'_'+gen+'_'+spe+'_'+accn
                            if child_name in ('Francisellaceae_Francisella_tularensis_GCA_000018925.1','Burkholderiaceae_Burkholderia_gladioli-A_GCA_009911875.1',):
                                sys.exit(child_name)
                            ##/
                            
                        parent = node._up
                        if parent:
                            # Check so there is not a conflict in classification at parent
                            parent_child_classis_counts = get_node_children_classis_counts(parent)
                            #/
                            # Require upstream node to have only same classifications as current node (leaf+branch), or no classifications at all (branch+branch)
                            if (parent_child_classis_counts.keys() == leafClassis_counts.keys()) or not parent_child_classis_counts:
                                
                                node_child_dist = node.get_distance(child)
                                node_parent_dist = node.get_distance(parent)
                                child_parent_dist = node_child_dist + node_parent_dist
                                
                                child.detach()
                                parent.add_child(child)
                                child.dist = child_parent_dist
                                delete_node = True # must delete node after repositioning all childs
                    
                    if delete_node:
                        node.delete()
                        deletion_made = True # toggle to restart while-loop
            #/
        iterations += 1
    #/
##/

## Check if reduce >2 leaf's to top2 (based on distance) leaf
if collapse_leafs:
    collapsed_childs = {} # "best_child" -> "collapsed childs"
    # Find branch-nodes with >1 leafs and select best leaf
    for node in tree.traverse():
        if not node.is_leaf():
            # get all childrens (leafs + branches)
            childrens = node.get_children()
            #/
            
            # sort leaf-childrens by classification
            childrens_classified = {}
            for child in childrens:
                if child.is_leaf() and len(child.get_children()) == 0:
                    # Check if we have metadata from DB imported. Else use default classification as "None"
                    if accession_metadata:
                        classi = get_node_classi(child)
                    else:
                        classi = None
                    #/
                    if not classi in childrens_classified:      childrens_classified[classi] = []
                    childrens_classified[classi].append(child)
            #/
            # for each classification of leafs-chindren, keep best child only
            for classi,leaf_childs in childrens_classified.items():
                # Check if >1 leaf-child, then select child with shortest distance to current branch-node and discard the rest
                if len(leaf_childs) > 1:
                    # determine which child to keep
                    best_child = sorted(leaf_childs,key=lambda x: x.dist)[0]
                    #/
                    # determine which childs to remove
                    discard_childs = []
                    for child in leaf_childs:
                        if not child == best_child:
                            discard_childs.append(child)
                    #/
                    # execute removal
                    for discard_child in discard_childs:
                        discard_child.detach()
                    #/
                    # save collapsed information
                    collapsed_childs[best_child.name] = [discard_child.name for discard_child in discard_childs]
                    #/
                #/
            #/
            # Check if branch now has a single orphan leaf, then move child to parent and readjust dist
            if len(node.get_children()) == 1:
                best_child = node.get_children()[0]
                parent = node._up

                best_child_parent_dist = parent.get_distance(best_child)
                
                best_child.detach()
                parent.add_child(best_child)
                best_child.dist = best_child_parent_dist
                
                node.delete()
            #/
    #/
##/


## IDE: Format names
if format_names and accession_metadata:
    for node in tree.traverse():
        if node.is_leaf():
            name = node.name
            accession = global_functions.getAccessions(name,return_first=True)
            fam,gen,spe,mycol1,accn = accession_metadata[accession]['family'],accession_metadata[accession]['genus'],accession_metadata[accession]['species'],accession_metadata[accession]['mycol1'],accession
            node.name = fam+'_'+gen+'_'+spe+'_'+mycol1+'_'+accn
##/


#print(tree)
print(tree.write(format=1))

