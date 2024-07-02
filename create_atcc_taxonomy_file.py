#!/usr/bin/env python3

## ___________________________________________________________________________
## Functions

def argparse():
    """
    __________________________________________________________________________
    Description
        - Interface for collecting and organizing CLI input data
    """
    ## Packages
    import argparse

    ## Initialize parser
    parser = argparse.ArgumentParser(description="A tool for generating a formatted ATCCfinder taxonomy table for ATCC genomes by parsing compiled ATCC multifasta reference genomes.")

    ## Define arguements
    parser.add_argument(
        "--infile", help = "Input multiline fasta file gernated from compoiling ATCC refernce genomes.",
        required = True, default="" )

    parser.add_argument(
        "--outfile", help = "Output file for parsed and formatted ATCC genome id and taxonomy associations.",
        required = True, default = False, nargs = "*" )

    ## Store CLI input in vArgs dictionary {'option':'input'}
    args = vars(parser.parse_args())
    return(args)

def create_atcc_taxonomy_file(infile, outfile):
    """
    """

    ## _______________________________________________________________________
    ## Packages

    import gzip
    import io
    import zipfile



    ## _______________________________________________________________________
    ## Extract headers

    ## Initialize disctionary for storing id : taxonomy pairs
    id_tax = {}

    ## Capture input fasta file format
    extension = infile.split(".")[-1]

    print("\n________________________________________")
    print("Parsing taxonomy information..")

    ## ____________________________________________________
    ## Scenario: Infile is uncompressed
    if extension in ["fasta", "fna", "fa"]:
        with open(infile) as inf:
            for n, line in enumerate(inf):
                ## Report progress
                if n%1000 == 0: print(" - " + str(n) + " unique IDs parsed")
                ## Parse and record unique data
                if line.startswith(">"):
                    id  = line.split()[0].replace(">","")
                    tax = line.split("species=")[1].split("contig_number=")[0].replace("\"","")
                    if id not in id_tax.keys():
                        id_tax[id] = tax

    ## ____________________________________________________
    ## Scenario: Infile has gzip compression
    if extension in ["tgz", "gz"]:
        with gzip.open(infile, "rt") as inf:
            for n, line in enumerate(inf):
                ## Report progress
                if n%1000 == 0: print(" - " + str(n) + " unique IDs parsed")
                ## Parse and record unique data
                if line.startswith(">"):
                    id  = line.split()[0].replace(">","")
                    tax = line.split("species=")[1].split("contig_number=")[0].replace("\"","")
                    if id not in id_tax.keys():
                        id_tax[id] = tax

    ## ____________________________________________________
    ## Scenario: Infile has zip compression
    if extension in ["zip"]:
        with zipfile.ZipFile(infile) as z:
            for n in z.namelist():
                with io.TextIOWrapper(z.open(n)) as inf:
                    for n, line in enumerate(inf):
                        ## Report progress
                        if n%1000 == 0: print(" - " + str(n) + " unique IDs parsed")
                        ## Parse and record unique data
                        id  = line.split()[0].replace(">","")
                        tax = line.split("species=")[1].split("contig_number=")[0].replace("\"","")
                        if id not in id_tax.keys():
                            id_tax[id] = tax



    ## _______________________________________________________________________
    ## Export taxonomy information
    
    with open(outfile, "w+") as outf:

        ## Initialize header
        outf.write( "\t".join( ["id", "taxonomy"] ) + "\n" )

        ## Write dictionary to output file
        for id in id_tax:
            outf.write( "\t".join( [id, id_tax[id]] ) + "\n" )





## ___________________________________________________________________________
## Main

if __name__ == "__main__":

    ## Load CLI inputs
    args = argparse()

    ## Execute
    create_atcc_taxonomy_file(
        infile  = args["infile"]  ,
        outfile = args["outfile"] )
