#!/usr/bin/env python3

## ___________________________________________________________________________
## Packages

import os





## ___________________________________________________________________________
## Functions, General

def argparse():
    """
    __________________________________________________________________________
    Description
        - Interface for collecting and organizing CLI input data
    """
    ## Packages
    import argparse

    ## Initialize parser
    parser = argparse.ArgumentParser(description="ATCCfinder")

    ## Define arguements
    parser.add_argument(
        "--api_key", help = "ATCC account-associated key used for accessing the ATCC servers.",
        required = False, default="" )

    parser.add_argument(
        "--download", help = "Specify which / any ATCC databases to download.",
        required = False, default = False, nargs = "*" )
    
    parser.add_argument(
        "--atcc_ids", help = "If provided, only downloads references associated with the specified ATCC ids. Input should be a newline-seperated .txt file.",
        required = False, default = "")

    parser.add_argument(
        "--overwrite", help="If downloading databases, specifies whether or not downloads in output folder should be kept or replaced.",
        required = False, default = False, action = argparse.BooleanOptionalAction)

    parser.add_argument(
        "--format", help="",
        required = False, default = False)

    parser.add_argument(
        "--out", help="",
        required = False, default = False)

    ## Store CLI input in vArgs dictionary {'option':'input'}
    args = vars(parser.parse_args())
    return(args)

def list_files_in_dir(dir, full_path=False, pattern=False, extension=False):
    """
    __________________________________________________________________________
    Description
        - List files in a directory (dir)
    __________________________________________________________________________
    Variables
        - dir       (string)    The target directory to be parsed
        - full_path (boolean)   If True, includes 'dir' path in output list items
        - pattern   (string)    Filters results by presence of this specified value
        - extension (boolean)   If True, includes file extensions in output list items
    """
    ## Packages
    import os
    ## Main
    if full_path:
        dirs = []
        for item in os.listdir(dir):
            dirs.append(os.path.join(dir, item))
    else:
        dirs = os.listdir(dir)

    ## Remove hidden files
    dirs_clean = []
    for file in dirs:
        if file.split("/")[-1][0] != ".":
            dirs_clean.append(file)

    ## Apply pattern matching to directories being returned
    if pattern:
        dirs_clean = [ x for x in dirs_clean if pattern in x ]

    ## Remove extension
    if not extension:
        dirs_clean = [ x.split(".")[0] for x in dirs_clean]

    return(dirs_clean)

def flatten_dictionary(dictionary, prefix='', sep="->"):
    """
    __________________________________________________________________________
    Description
        - Takes a nested dictionary an reduces nested keys into a single key.
    __________________________________________________________________________
    Variables
        - dictionary    Input nested dictionary
        - prefix        What the compressed key starts with
        - sep           Seperator for compressed key values
    """
    ## Intialize flattened dictionary
    flat = {}
    for rkey,val in dictionary.items():
        key = prefix + rkey
        if isinstance(val, dict):
            flat.update(flatten_dictionary(dictionary = val, prefix = key + sep))
        else:
            flat[key] = val
    return flat

def equalize_flat_dict(dict, sep="->", str="NA"):
    """
    __________________________________________________________________
    Description
        - Takes the output of function flatten_dictionary() and standardizes the lengths of flattened dictionary keys
    
    __________________________________________________________________
    Variables
        - dict  (dictionary)    The flattened dictionary to format
        - sep   (string)        The substring seperating flattened dictionary key values
        - str   (string)        The subtring to incorporate into equal-length formatted output dictionary keys
    """
    ## Identify longest key
    max_n = max( [len(x.split(sep)) for x in dict.keys()] )

    ## Initialize output dictionary
    dict_out = {}

    ## Iterate through keys and revise key lengths so that they are all the same split length
    for key in dict.keys():
        ## Identify the length of the key
        len_key = len(key.split(sep))
        ## Determine how much shorter key is than longest key
        blanks_to_add = max_n - len_key
        ## Define substring to equalize key length with maximum flattened key
        equal_filler = ( sep + "NA" ) * blanks_to_add
        ## Add equalizer substring to flattened to key
        key_standard = key + equal_filler
        ## Add to standardized output dictionary
        dict_out[key_standard] = dict[key]
    return(dict_out)

def list_unique_in_column(dir, col):
    """
    __________________________________________________________________
    Description
        - Returns a list of unique values present in a .tsv file column
        - Expects file to contain header

    __________________________________________________________________
    Variables
        - dir   (string)  
        - col   (string)  
    """
    pass

def export_gz(path, string):
    """
    Save file in gzip format
    """
    ## Packages
    import gzip
    ## Main
    with gzip.open(path, 'ab') as outfile:
        outfile.write(string.encode())

def test_if_fasta(file):

    """
    README: Tests if a file (x) is in fasta format by looking for a ">" first character
    """

    ## Packages
    import gzip

    ## Test if input is compressed (gzs)
    compressed = False
    if file.split(".")[-1] == "gz":
        compressed = True

    ## Main
    is_fasta=False 

    ## Scenario, Zipped file
    if compressed:
        with gzip.open(file) as infile:
            for n,line in enumerate(infile):
                if n==0:
                    if line.decode()[0]==">":
                        is_fasta = True

    ## Scenario, NOT zipped file
    else:
        with open(file) as infile:
            for n,line in enumerate(infile):
                if n==0:
                    if line[0]==">":
                        is_fasta = True

    return(is_fasta)

def load_atcc_ids(dir):
    """
    README: Load a newline-delimited textfile specifying ATCC genome ids to download.
    """

    ## Initialize list for storing ids
    ids = []

    ## Parse ids
    with open(dir) as infile:
        for id in infile:
            id_clean = id.strip()
            if id_clean != "":
                ids.append(id_clean)
            
    return(ids)





## ___________________________________________________________________________
## Functions, ATCC genome API wrapper

def ATCC_retrieve_ids(api_key):
    """
    __________________________________________________________________________
    Description
        - Uses the download_all_genomes() function to return a list of all ATCC genome IDs.
    __________________________________________________________________________
    Variables
        - api_key   (string)    ATCC-account affiliated API key used to access ATCC genome browser databases.
    """
    ## Packages
    import genome_portal_api
    import time

    ## Time, capture start time
    time_start = time.time()

    ## Initialize page counter
    page    = 1  ## Tracks the ATCC catalogue page (min = 1), each of which contains max 50 genomes to explore
    all_ids = [] ## List for storing all ATCC genome IDs, which will be returned by this function

    ## Iteratively browse pages until blank page is reached
    print("\nRetrieving ATCC Genome IDs...")
    while True:

        ## Report progress
        if page % 10 == 0:
            print("\t> ATCC pages parsed: ", page)

        ## Collect list of genomes on current page
        catalogue_page = genome_portal_api.download_all_genomes(api_key=api_key, page=page)

        ## Add genome IDs to export list or terminate retrival if last page has been reached
        if catalogue_page is not None:
            all_ids.extend([x["id"] for x in catalogue_page])
            page +=1
        else:
            break
    
    ## Time, capture end time and calculate elapsed
    time_end = time.time()
    delta_time = int(time_end - time_start)

    ## Report results and export
    print("\n\t>", str(len(all_ids)), " ATCC genome ids retrieved in ", delta_time, " seconds across ", page, " pages.")
    return(all_ids)

def ATCC_download_annotations(ids, api_key, dir_out, overwrite):
    """
    Description
        - Downloads ATCC genome annoptation files in the genbank .gbk format.

    Variables
        - api_key   (string)    ATCC-account affiliated API key used to access ATCC genome browser databases.
        - dir_out   (string)    Output directory where all ATCCfinder results will be stored (in auto-generated sub-folders)
        - ids       (list)      A list of ATCC internal IDs, which will be used to retrieve desired data for downloading
        - overwrite (bool)      Specifies whether or not pre-existing downloads should be kept. (False = keep previous downloads)
    """

    ## [Packages]
    import genome_portal_api
    import os

    ## [Paths] Define and create output directory
    dir_out_db = os.path.join(dir_out, "annotation")
    if not os.path.isdir(dir_out_db): os.makedirs(dir_out_db)

    ## [Overwrite] If downloads are NOT to be overwritten, then identify which files have already been downloaded and exclude these from the download list
    if not overwrite:
        ## Retrive downloaded items from folder
        downloaded = list_files_in_dir(dir=dir_out_db, full_path=False, pattern=False, extension=False)
        ## Identify which ids have been downloaded and exclude them from id download list
        ids = list( set(ids) - set(downloaded) )

    ## [Download] Iterate through ATCC genome ids and download associated annotations
    print("\nDownloading ATCC annotation database:")
    for id in ids:

        ## Download annotation from ATCC
        print(id)
        n_tries = 3
        for n in range(n_tries):

            ## Attempt download
            annotation = genome_portal_api.download_annotations(
                api_key              = api_key,
                id                   = id,
                download_link_only   = "False",
                download_annotations = "True" )

            ## Test if download was successful
            if "<Error>" not in annotation:
                ## Define output file
                file_out = os.path.join(dir_out_db, id + ".gbk")

                ## Export
                with open(file_out, "w+") as outfile:
                    outfile.write(annotation)

                ## Halt loop and move on to the next id
                break

def ATCC_download_catalogue(api_key, dir_out):
    """
    __________________________________________________________________
    Description
        - Function for downloading ATCC catalogue database.
    __________________________________________________________________
    Variables
        - ids       (list)      A list of ATCC internal IDs, which will be used to retrieve desired data for downloading
        - api_key   (string)    ATCC-account affiliated API key used to access ATCC genome browser databases.
        - dir_out   (string)    Output directory where all ATCCfinder results will be stored (in auto-generated sub-folders)
        - overwrite (bool)      Specifies whether or not pre-existing downloads should be kept. (False = keep previous downloads)
    """

    ## [Packages]
    import genome_portal_api
    import os

    print("\nDownloading ATCC catalogue database:")
    ## [Paths] Define and create output directory & pickle file
    dir_out_db = os.path.join(dir_out, "catalogue")
    if not os.path.isdir(dir_out_db): os.makedirs(dir_out_db)
    # file_out = os.path.join(dir_out_db, "catalogue.pkl")
    file_out = os.path.join(dir_out_db, "catalogue.tsv")

    ## [Download] Retrieve data
    catalogue = genome_portal_api.download_catalogue(api_key=api_key)

    ## [Parse] Initialize Variables
    category_levels = 0
    cat_a = ""
    cat_b = ""
    cat_c = ""

    with open(file_out, "w+") as outfile:

        ## Output, write header
        outfile.write("\t".join(["assembly_n", "atcc_id,", "product_id", "category_levels", "category_a", "category_b", "category_c", "value"]) + "\n")

        ## Iterate through genome meta data
        for n, data in enumerate(catalogue):

            ## Collect genome information
            n_str      = str(n+1) ## 0-indexed, so add 1 always
            atcc_id    = data['id']
            product_id = data['product_id']
            print(product_id)

            ## Iterate through highest level dictionary keys
            for key1 in data:

                ## Define current category
                cat_a = key1

                ## Identify items that are not nested at this level (level = 1)
                if type(data[key1]) is not dict:
                    category_levels = "1"
                    cat_b = "NA"
                    cat_c = "NA"
                    value = str(data[key1])
                    if value=="None": value = "NA"
                        
                    ## Export
                    outfile.write("\t".join([n_str, atcc_id, product_id, category_levels, cat_a, cat_b, cat_c, value]) + "\n")
                
                else:
                    for key2 in data[key1]:

                        ## Define current category
                        cat_b = key2

                        ## Identify items that are not nested at this level (level = 1)
                        if type(data[key1][key2]) is not dict:
                            category_levels = "2"
                            cat_c = "NA"
                            value = str(data[key1][key2])
                            if value=="None": value = "NA"
                                
                            ## Export
                            outfile.write("\t".join([n_str, atcc_id, product_id, category_levels, cat_a, cat_b, cat_c, value]) + "\n")

                        else:
                            for key3 in data[key1][key2]:

                                ## Define current category
                                cat_c = key3

                                ## No more items should be nested in dictionaries at this level
                                category_levels = "3"
                                value = str(data[key1][key2][key3])
                                if value=="None": value = "NA"
                                    
                                ## Export
                                outfile.write("\t".join([n_str, atcc_id, product_id, category_levels, cat_a, cat_b, cat_c, value]) + "\n")

def ATCC_download_metadata(ids, api_key, dir_out, overwrite):
    """
    __________________________________________________________________
    Description
        - Function for downloading ATCC metadata database.
    __________________________________________________________________
    Variables
        - ids       (list)      A list of ATCC internal IDs, which will be used to retrieve desired data for downloading
        - api_key   (string)    ATCC-account affiliated API key used to access ATCC genome browser databases.
        - dir_out   (string)    Output directory where all ATCCfinder results will be stored (in auto-generated sub-folders)
        - overwrite (bool)      Specifies whether or not pre-existing downloads should be kept. (False = keep previous downloads)
    """
    ## [Packages]
    import genome_portal_api
    import os

    ## [Paths] Define and create output directory & file
    print("\nDownloading ATCC metadata database:")
    dir_out_db = os.path.join(dir_out, "metadata")
    if not os.path.isdir(dir_out_db): os.makedirs(dir_out_db)

    ## If overwrite is disabled, then identify and exclude pre-downloaded ids from id list
    if not overwrite:
        ## Identify ids already present in download folder
        existing_ids = list_files_in_dir(dir=dir_out_db)
        ## Excluded downloaded files from id list
        ids = list( set(ids) - set(existing_ids) )
        print("\t> Overwrite disabled, ", str(len(ids)), " ids identified for downloading.")

    ## Iterate through IDs and retrieve metadata
    for n,id in enumerate(ids):

        ## Retrieve metadata dictionary from ATCC
        print( "\t> " + str(n) + " / ", str(len(ids)) )
        assembly_metadata_raw = genome_portal_api.download_metadata(api_key=api_key, id=id)

        ## Collapse nested dictionary pattern by combining nested keys into a single string
        assembly_metadata = flatten_dictionary(
            dictionary = assembly_metadata_raw,
            prefix     = '',
            sep        = "->" )

        ## Standardize length of flattened meta data dictionary keys
        assembly_metadata_std = equalize_flat_dict(
            dict = assembly_metadata,
            sep  = "->" )

        ## Define and initialize outfile
        file_out   = os.path.join(dir_out_db, id + ".tsv")
        with open(file_out, "w+") as outfile:
            
            ## Define header line
            outfile.write( "\t".join(["ID", "Category", "Value"]) + "\n" )

            ## Export meta data lines to file
            for key in assembly_metadata_std.keys():
                outfile.write( "\t".join([id, key, str(assembly_metadata_std[key])]) + "\n" )

def ATCC_download_assembly(ids, api_key, dir_out, overwrite):
    """
    __________________________________________________________________
    Description
        - Function for downloading ATCC reference sequences.
    __________________________________________________________________
    Variables
        - ids       (list)      A list of ATCC internal IDs, which will be used to retrieve desired data for downloading
        - api_key   (string)    ATCC-account affiliated API key used to access ATCC genome browser databases.
        - dir_out   (string)    Output directory where all ATCCfinder results will be stored (in auto-generated sub-folders)
        - overwrite (bool)      Specifies whether or not pre-existing downloads should be kept. (False = keep previous downloads)
    """

    ## [Packages]
    import genome_portal_api
    import glob
    import os

    ## [Paths] Define and create output directory & file
    print("\nDownloading ATCC assembly database:")
    dir_out_db = os.path.join(dir_out, "assembly")
    if not os.path.isdir(dir_out_db): os.makedirs(dir_out_db)

    ## Initialize variables
    page = 0 ## Tracks the ATCC catalogue page, each of which contains max 50 genomes to explore
    max_attempts = 3 ## Defines how many times to try and download a genome before giving up and moving on to the next item in the ATCC catalogue.

    ## If overwrite is disabled, then identify and exclude pre-downloaded ids from id list
    if not overwrite:
        ## Identify ids already present in download folder
        existing_ids = list_files_in_dir(dir=dir_out_db)
        ## Excluded downloaded files from id list
        ids = list( set(ids) - set(existing_ids) )
        print("\t> Overwrite disabled, ", str(len(ids)), " ids identified for downloading.")

    ## Iterate through genomes
    for id in ids:
        
        ## Initialize attempt counter
        downloading = True
        attempt     = 1

        ## Prior attempts to download genomes resulted in random occurences of genome ids not being recognized by ATCC.
        ## These same problematic ids would later be fine, necesitating this check to ensure the download succeeded.
        while downloading:

            ## Define genome outfile
            filename = id + ".fa.gz"
            dir_out = os.path.join(dir_out_db, filename)

            ## If genome has not been downloaded, then download it
            print(
                "____________________________________________________________"
                "\n"   + "Downloading ATCC assembly:"+
                "\n\t" + "> attempt "  + str(attempt) + " / " + str(max_attempts+1)+
                "\n\t" + "> ATCC ID: " + id +
                "\n"   )
            
            ## Download genome
            fasta = genome_portal_api.download_assembly(
                api_key            = api_key,
                id                 = id,
                download_link_only = "False",
                download_assembly  = "True" )

            ## Write genome to file
            if "<Error><Code>" not in "".join( fasta.keys() ):
                for key in fasta:
                    export_gz( path = dir_out, string = key + "\n" + fasta[key] + "\n" )
                downloading = False

def ATCC_download_manager(download, api_key, dir_out, overwrite):
    """
    __________________________________________________________________
    Description
        - Function for managing ATCC genome database downloads
    __________________________________________________________________
    Variables
        - download  (list)      Specifies which / any ATCC databases to download.
        - api_key   (string)    ATCC-account affiliated API key used to access ATCC genome browser databases.
        - dir_out   (string)    Output directory where all ATCCfinder results will be stored (in auto-generated sub-folders)
        - overwrite (bool)      Specifies whether or not pre-existing downloads should be kept. (False = keep previous downloads)
    """

    ##________________________________________________________________________
    ## Packages
    import genome_portal_api

    ##________________________________________________________________________
    ## Initialize

    if download:
        ## Identify which databases have been specified for downloading
        atcc_databases = {"annotation", "reference", "catalogue", "metadata"}
        dbs_to_dl      = set(download).intersection(atcc_databases)
        unknown_dbs    = set(download) - {"annotation", "reference", "catalogue", "metadata"}

        ## Report specified databases
        if dbs_to_dl:   print("\nDatabases specified for downloading:\n\t>",                    "\n\t> ".join(dbs_to_dl))
        if unknown_dbs: print("\nWARNING: Unknown databases specified for downloading:\n\t> " + "\n\t> ".join(unknown_dbs))

        ## Load ATCC genome IDs
        if set(download).intersection({"annotation", "reference", "metadata"}):

            if args['atcc_ids'] == "": ids = ATCC_retrieve_ids(api_key=api_key)
            if args['atcc_ids'] != "": ids = load_atcc_ids(dir=args['atcc_ids'])
        
        ##________________________________________________________________________
        ## [Download Database] Annotations
        if "annotation" in download:
            ATCC_download_annotations(
                ids       = ids       ,
                api_key   = api_key   ,
                dir_out   = dir_out   ,
                overwrite = overwrite )

        ##________________________________________________________________________
        ## [Download Database] reference
        if "reference" in download:
            ATCC_download_assembly(
                ids       = ids       ,
                api_key   = api_key   ,
                dir_out   = dir_out   ,
                overwrite = overwrite )
        
        ##________________________________________________________________________
        ## [Download Database] Catalogue
        if "catalogue" in download:
            ATCC_download_catalogue(
                api_key = api_key ,
                dir_out = dir_out )
        
        ##________________________________________________________________________
        ## [Download Database] Meta data
        if "metadata" in download:
            ATCC_download_metadata(
            ids       = ids       ,
            api_key   = api_key   ,
            dir_out   = dir_out   ,
            overwrite = overwrite )

def ATCC_format_manager(dir_in, dir_out):
    """
    __________________________________________________________________
    Description
        - Function for managing ATCC database formatting post-download
    __________________________________________________________________
    Variables
        - dir_out   (string)    The directory to export formatted databases
    """
    print("\nFormatting ATCC downloads:")
    format_references(dir_in=dir_in, dir_out=dir_out)

def list_refs_in_download(dir):
    """
    __________________________________________________________________________
    Description
        - Returns list of IDs downloaded from ATCC reference genome database
    __________________________________________________________________________
    Variables
        - dir   (string)    Path to reference database download directory
    """
    import os
    ids = []
    if os.path.exists(dir):
        ## collect names, which will include an underscore with species information
        ids = list_files_in_dir(dir=dir, full_path=False, pattern=".fa.gz", extension=False)
    return(ids)

def list_refs_in_compiled(dir):
    """
    __________________________________________________________________________
    Description
        - Returns a list of IDs found in fasta file of compiled ATCC reference genomes
    __________________________________________________________________________
    Variables
        - dir   (string)    Path to compiled reference fasta file
    """
    import os, re
    ids = []
    if os.path.exists(dir):
        with open(dir) as infile:
            for line in infile:
                if line.startswith(">"):
                    id_1 = re.split('genome_id="|";atcc_catalog_number='   , line)[1] ## ATCCGenome Id
                    id_2 = re.split('atcc_catalog_number="ATCC;|";species=', line)[1] ## ATCC Catalog number
                    id = id_1 + "_" + id_2
                    ids.append(id)
    return(ids)
    
def compile_refs(dir_download, dir_compiled, dir_headers):
    """
    __________________________________________________________________________
    Description
        - Takes multiple individual fasta files and combines them into a single query-able fasta file
    __________________________________________________________________________
    Variables
        - dir_download  Directory containing all downloaded aTCC reference assemblies
        - dir_compiled  File (.fasta) containing compiled ATCC refefence assemblies
        - dir_headers   Output path for file containing all fasta headers in compiled ATCC reference file
    """
    import gzip, os
    print("\nCompiling ATCC references..")

    ## Identify what IDs have been downloaded and compiled
    ids_download = list_refs_in_download(dir=dir_download) ## Downloaded IDs
    ids_compiled = list_refs_in_compiled(dir=dir_compiled) ## Compiled IDs
    print("\t- ATCC references downloaded:      ", len(ids_download))
    print("\t- ATCC references already compiled:", len(ids_compiled)  )

    ## Determine references not in compiled file
    ids_add = list( set(ids_download) - set(ids_compiled) )
    
    if len(ids_add) > 0:
        print("\t- References missing from compiled file, re-compiling..")

        ## Collect list of all downloaded reference files
        ref_names_path = list_files_in_dir(dir=dir_download, full_path=True,  pattern=".fa.gz", extension=True)

        ## Iterate through fa.gz files and merge
        with open(dir_compiled, "w+") as outfile, open(dir_headers, "w+") as out_headers:
            for n,file in enumerate(ref_names_path):

                ## Report progress
                if (n+1) % 200 == 0: print("\t- References merged: ", n)

                ## Iterate fasta file and download to files
                with gzip.open(file) as infile:
                    for line in infile:
                        
                        if line.decode().startswith(">") == True:
                            out_headers.write(line.decode())
                            outfile.write(line.decode().replace(' ', '_'))

                        if line.decode().startswith(">") == False:
                            outfile.write(line.decode())

                ## If last genome has been reached, do not include a newline at the end of the file.
                if n!=len(ref_names_path):
                    outfile.write("\n")
    
    print("\t- References successfully compiled")

def format_references(dir_in, dir_out):
    """
    __________________________________________________________________________
    Description
        - Combine ATCC reference genomes into a single query-able fasta file for use with minimap2
    """
    ## Import
    import os

    if os.path.exists(dir_in):

        ## [Paths] Define and create output directory
        out_compiled = os.path.join(dir_out, "atcc_references.fa")
        out_headers  = os.path.join(dir_out, "atcc_references_headers.txt")
        if not os.path.isdir(dir_out): os.makedirs(dir_out)

        ## Compile references, accounting for already-compiled IDs
        compile_refs(
            dir_download = dir_in,
            dir_compiled = out_compiled,
            dir_headers  = out_headers )





## ___________________________________________________________________________
## Main

## Load CLI
args = argparse()


## Download ATCC database(s), if specified
if args['download'] != False :

    ## Download any specified databases
    ATCC_download_manager(
        dir_out   = args['out'      ] ,
        download  = args['download' ] ,
        api_key   = args['api_key'  ] ,
        overwrite = args['overwrite'] )

    ## Format downloaded ATCC databases
    ATCC_format_manager(
        dir_in  = os.path.join(args['out'], "assembly"),
        dir_out = args['out'] )



if args['download'] == False and args['format'] != False:

    ## Format downloaded ATCC databases
    ATCC_format_manager(
        dir_in  = args['format'] ,
        dir_out = args['out'   ] )
