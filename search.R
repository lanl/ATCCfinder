#!/usr/bin/env Rscript

# ______________________________ ------------------------------------------
# Parameters --------------------------------------------------------------

args_interactive <- function() {
  "
  ______________________________________________________________________________
  Description
    - Define parameters here if running script from RStudio.
  "
  ## Define Arguments
  args <- list(
    target    = "/Users/skoehler/Desktop/project_2302_MO/Workflow_ATCC/atcc_downloads/assembly_compiled/atcc_references.fa",
    query     = "/Users/skoehler/Desktop/project_2302_MO/Workflow_ATCC/input/SK_Curated_Datasets/Amplicon_Sequences/Bacillus_subtilis_amplicon.fa",
    nhits     = "100",
    overwrite = "F",
    outdir    = "test_out" )
  
  ## Export list of arguments, for use when running script from RStudio
  return(args)
}





# ______________________________ ------------------------------------------
# Parameters, Commandline -------------------------------------------------

argparse <- function() {
  "
  ____________________________________________________________________________
  Description
    - Define and retrieve command-line input values
  "
  ## Packages
  suppressWarnings( library(argparse) )
  
  ## Define arguement parser
  parser <- ArgumentParser(description='Arguement Parser')
  
  ## Define Arguments
  parser$add_argument("--query"    , required=T, help="Query Sequence(s) fasta file")
  parser$add_argument("--target"   , required=T, help="Target database fasta file")
  parser$add_argument("--nhits"    , required=T, help="The maximum number of target hits to return per query", type="integer")
  parser$add_argument("--overwrite", required=T, help="(T/F) If alignment output already exists, should it be overwritten?", default="T")
  parser$add_argument("--outdir"   , required=T, help="Output directory")
  
  ## Load and format inputs
  args <- parser$parse_args()
  args["nhits"] <- as.numeric(args["nhits"])
  args["overwrite"] <- args["overwrite"] == "T"
  
  return(args)
}





# ______________________________ ------------------------------------------
# Functions ---------------------------------------------------------------

identify_run_environment_execute <- function() {
  "
  ______________________________________________________________________________
  Description
    - Manage argument parsing and main execution of software
  "
  ## _________________________________________________________________________
  ## Initialize Environment
  identify_and_set_wd_either()
  
  ## _________________________________________________________________________
  ## [- Interactive -] Use parameters defined in this file
  if ( interactive() == T ) { args <- args_interactive() }
  
  ## [- Command Line- ] Use parameters defined via command line argparse
  if ( interactive() == F ) {  args <- argparse() }
  
  ## Execute
}

make_path <- function(path) {
  dir.create(path=file.path(path), recursive=T, showWarnings=F)
}

ATCCfinder_align <- function(query, target, nhits, outdir, mode='minimap2') {
  "
  ____________________________________________________________________________
  Description
    - Wrapper function for running minimap2
  ____________________________________________________________________________
  Variables
    - query
    - target
    - outdir
    - mode    Defines what aligner should beused (minimap2, bwa, dynamic)
  ____________________________________________________________________________
  Notes
    - Requires an environment with minimap2 & bwa installed
    - Plan on implementing BWA:
        system( paste('bwa', 'index', target) )
        system( paste('bwa', 'aln', target, query, '>', 'reads.sai;', 'bwa', 'samse', target, 'reads.sai', query, '>', outfile) )
  "
  ## Align using minimap2
  if ( mode == "minimap2" ) {
    message("\nAligning query(s) to ATCC references via minimap2..\n")
    outfile <- file.path(outdir, "alignment.paf")
    # system( paste("minimap2", "-N", as.numeric(nhits)-1, "-x", "map-ont", "--paf-no-hit", target, query, ">", outfile) )
    system( paste("minimap2", "-c", "-N", as.numeric(nhits)-1, "-x", "map-ont", "--paf-no-hit", target, query, ">", outfile) )
    message("\nAlignment complete, paf file generated.\n")
  }
  
}

extract_fasta_headers <- function(fasta, outfile) {
  "
  ____________________________________________________________________________
  Description
    - Takes a fasta file and returns a file of all headers found in the file
  ____________________________________________________________________________
  Variables
    - fasta     The fasta file to extract headers from
    - outfile   Output file to store headers in
  "
  message("\nParsing target sequence meta data..")
  system(paste("grep", "-e", '">"', fasta, ">", outfile))
}

check_if_paf_exists <- function(file) {
  "
    __________________________________________________________________________
    Description
      - Tests T/F if a paf alignment file exists
    __________________________________________________________________________
    Variables
      - file    paf alignment file path
    "
  paf_exists <- file.exists(file)
  if (paf_exists) { message("\nNote: Alignment file already found in output directory. This result will be overwritten, if enabled.") }
  return( paf_exists )
}

load_parse_reference_headers <- function(file) {
  "
    __________________________________________________________________________
    Description
      - Takes ATCC reference header strings and parses them for ATCC-related meta data
    __________________________________________________________________________
    Variables
      - file    A tsv file of fasta headers (including '>' starting character) taken from ATCC references
  "
  ## Load file
  df <- read.delim(file=file, sep="\t", header=T)
  # df <- read.delim(file=file, sep="\t", header=F)
  # colnames(df) <- "raw"
  # 
  # ## Parse reference information from header
  # df[,"header_id"          ] <- gsub(x=df[,"raw"], pattern=paste0(c(">"                     , " .+"                    ), collapse="|"), replacement="")
  # df[,"assembly_id"        ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+assembly_id="        , " genome_id=.+"          ), collapse="|"), replacement="")
  # df[,"genome_id"          ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+genome_id="          , " atcc_catalog_number=.+"), collapse="|"), replacement="")
  # df[,"atcc_catalog_number"] <- gsub(x=df[,"raw"], pattern=paste0(c(".+atcc_catalog_number=", " species=.+"            ), collapse="|"), replacement="")
  # df[,"taxonomy_raw"       ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+species="            , " contig_number=.+"      ), collapse="|"), replacement="")
  # df[,"contig_number"      ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+contig_number="      , " topology=.+"           ), collapse="|"), replacement="")
  # df[,"topology"           ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+topology="           , ""                       ), collapse="|"), replacement="")
  # df[,"raw"] <- NULL
  # 
  # ## Parse taxonomy data
  # df[,"taxonomy_raw"] <- gsub(x=df[,"taxonomy_raw"], pattern="\\[|\\]", replacement="")
  # df[,"taxonomy_raw"] <- tolower(df[,"taxonomy_raw"])
  # 
  # df[,"n_taxname"]     <- unlist(lapply(X=df[,"taxonomy_raw"], FUN=function(x) length(strsplit(x, split=" ")[[1]])))
  # df[,"genus"] <- gsub(x=df[,"taxonomy_raw"], pattern=" .+", replacement="")
  # df[,"genus_species"] <- unlist(lapply(X=df[,"taxonomy_raw"], FUN=function(x) paste( strsplit(x, split=" ")[[1]][1:2], collapse=" ")))
  # 
  ## Remove duplicate rows
  df <- df[which(!duplicated(df)),]
  
  return(df)
}

load_format_paf <- function(file) {
  "
    __________________________________________________________________________
    Description
      - Load paf alignment file and apply standard formatting
    __________________________________________________________________________
    Variables
      - file    paf alignment file path
    "
  ## Headers
  cols_paf_mmap2 <- c("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "nmatch", "alen", "mapq")
  ## Load alignment results as single string to account for variable-length paf columns
  hits <- read.delim(file=file, header=F, sep="\n", fill=T, col.names = "raw")
  ## Parse first 12 columns of paf
  hits <- cbind(
    hits,
    data.frame(t(data.frame(lapply(
      X=hits[,"raw"], FUN = function(x) strsplit(x=x, split="\t")[[1]][1:12]
      ))))
    )
  row.names(hits) <- 1:nrow(hits)
  colnames(hits) <- c("raw", cols_paf_mmap2)
  
  ## Identify rows with 'AS' stat present & parse stat
  rows_as <- which( grepl(x=hits[,"raw"], pattern="AS:") )
  hits[rows_as,"AS"] <- as.numeric(unlist(lapply(
    X=hits[rows_as,"raw"], FUN = function(x) gsub(x=strsplit(x=x, split="\tAS:i:")[[1]][2], "", pattern="\t.+")
    )))
  
  ## Define column types
  hits$qlen   <- as.numeric(hits$qlen  )
  hits$qstart <- as.numeric(hits$qstart)
  hits$qend   <- as.numeric(hits$qend  )
  hits$tlen   <- as.numeric(hits$tlen  )
  hits$tstart <- as.numeric(hits$tstart)
  hits$tend   <- as.numeric(hits$tend  )
  hits$nmatch <- as.numeric(hits$nmatch)
  hits$alen   <- as.numeric(hits$alen  )
  hits$mapq   <- as.numeric(hits$mapq  )
  
  ## Identify hit type (primary, alternative, inversion)
  hits[,"tp"] <- ifelse(test = grepl(x=hits[,"raw"], pattern="tp:A:P"), yes="Primary", no="Secondary")
  
  ## Ensure no duplicate rows
  hits <- hits[which(!duplicated(hits)),]
  
  ## Define target (without contig information)
  hits[,"target"] <- gsub(x=hits[,"tname"], pattern="[_].+", "")
  
  return(hits)
}

ATCCfinder_summarize_paf <- function(target, query, paf, dir_out) {
  "
  ____________________________________________________________________________
  Description
    - Loads alignment results, scores hits, and identifies best reference match(es).
  ____________________________________________________________________________
  Variables
    - query       Fasta file that was queried to produce paf hit results
    - paf         Minimap2 alignment file results comparing user-specified query to ATCC reference target
    - outdir      Output directory
  "
  
  ## Load paf alignment results
  hits <- load_format_paf(file=paf)
  
  ## Load and parse full reference header data
  ref_info <- load_parse_reference_headers( file = file.path(dirname(target[[1]]), "atcc_reference_information.tsv") )
  
  ## Merge hits with header information
  hits_info <- merge(x=hits, y=ref_info, by.x="tname", by.y="header_id", all.x=T)
  
  ## Identify queries
  queries <- unique(hits[,"qname"])
  
  ## Iterate through queries and collect statistics
  hits_anot <- annotate_query_results(
    x       = hits_info[which(hits_info$tp=="Primary"),] ,
    queries = queries   )
  
  ## Summarize Query results
  result_summary <- ATCCfinder_summarize_query_results(
    x       = hits_anot ,
    queries = queries   )
  
  ## Export parsed results
  outfile_summary <- file.path(dir_out, "report_v2.tsv")
  write.table(x=result_summary, file=outfile_summary, sep="\t", row.names=F, quote=F)
  
  outfile_annot <- file.path(dir_out, "alignment_annotated_paf_P.csv")
  write.csv(x=hits_anot, file=outfile_annot, row.names=F, quote=F)
  
  ## Also export annotated hits with secodnary alignments
  ## Iterate through queries and collect statistics
  hits_anot_s <- annotate_query_results(
    x       = hits_info ,
    queries = queries   )
  outfile_annot_s <- file.path(dir_out, "alignment_annotated_paf_PS.csv")
  write.csv(x=hits_anot_s, file=outfile_annot_s, row.names=F, quote=F)
}

annotate_query_results <- function(x, queries) {
  "
    __________________________________________________________________________
    Description
      - Iterates through paf results by query, and analyzes results by various criteria
    __________________________________________________________________________
    Variables
      - x         Data frame of paf alignment results
      - queries   Vector of query strings to iterate through
    "
  n = 0
  for (q in queries) {
    n <- n + 1
    if (n%%100==0) {print(n)}
    ## subset hits to current query
    hits_query <- x[which(x[,"qname"]==q),]
    
    if ( nrow(hits_query) > 0 ) {
      ## Stats, Query
      query_name   <- q
      query_length <- max(hits_query$qlen)
      
      ## Stats, Target (general)
      n_hits <- nrow(hits_query)
      
      ## Stats, Target [ Longest Reference ]
      longest_ref <- max(hits_query$tlen)
      
      ## Stats, Target [ Longest Match ]
      max_nmatch   <- max(hits_query$nmatch)
      max_nmatch_p <- round( max_nmatch / query_length * 100, 2 )
      
      ## Stats, Target [ Minimum Contigs ]
      min_contigs <- as.numeric( min(hits_query$contig_number) )
      
      ## Stats, Target [ Maximum alignment score ]
      max_AS <- max(hits_query$AS)
      
      ## Stats, Target [ Maximum mapq ]
      max_q <- max(hits_query$mapq)
      
      ## Stats, Target [ Taxonomy Consensus ]
      taxa_hits <- data.frame(table(hits_query$genus_species))
      top_taxa_n <- max(taxa_hits$Freq)
      top_taxa <- taxa_hits[which(taxa_hits$Freq==top_taxa_n),"Var1"]
      
      ## Stats, Target [ Reference Consensus ]
      refs <- data.frame(table(hits_query$atcc_catalog_number))
      # refs$Var1 <- gsub(x=unique(refs$Var1), pattern="ATCC ", replacement="")
      top_ref_n <- max(refs$Freq)
      top_ref <- refs[which(refs$Freq==top_ref_n),"Var1"]
      
      ## Score hits based on metrics
      hits_query[,"stat_max_AS"    ] <- 0
      hits_query[,"stat_max_mapq"  ] <- 0
      hits_query[,"stat_top_taxa"  ] <- 0
      hits_query[,"stat_top_ref"   ] <- 0
      hits_query[,"stat_max_nmatch"] <- 0
      hits_query[,"stat_max_tlen"  ] <- 0
      hits_query[,"stat_min_cntigs"] <- 0
      hits_query[,"stat_algn_prim" ] <- 0
      
      hits_query[which(hits_query[,"AS"                 ] ==   max_AS     ), "stat_max_AS"    ] <- 1
      hits_query[which(hits_query[,"mapq"               ] ==   max_q      ), "stat_max_mapq"  ] <- 1
      hits_query[which(hits_query[,"genus_species"      ] %in% top_taxa   ), "stat_top_taxa"  ] <- 1
      hits_query[which(hits_query[,"atcc_catalog_number"] %in% top_ref    ), "stat_top_ref"   ] <- 1
      hits_query[which(hits_query[,"nmatch"             ] ==   max_nmatch ), "stat_max_nmatch"] <- 1
      hits_query[which(hits_query[,"tlen"               ] ==   longest_ref), "stat_max_tlen"  ] <- 1
      hits_query[which(hits_query[,"contig_number"      ] ==   min_contigs), "stat_min_cntigs"] <- 1
      hits_query[which(hits_query[,"tp"                 ] ==   "Primary"  ), "stat_algn_prim" ] <- 1
      
      ## Define cumulative score
      # hits_query[,"stat_cum_score"] <- rowSums(hits_query[,(ncol(hits_query)-5):ncol(hits_query)])
      # hits_query[,"stat_cum_score"   ] <- hits_query[,"stat_max_mapq"] + hits_query[,"stat_max_nmatch"]
      hits_query[,"stat_cum_score"] <- hits_query[,"stat_max_AS"  ] + hits_query[,"stat_max_nmatch"] + hits_query[,"stat_algn_prim"]
      
      ## Add information to compiled df
      ifelse(
        test = exists("paf_annot")                        ,
        yes  = paf_annot <- rbind( paf_annot, hits_query) ,
        no   = paf_annot <- hits_query                    )
    }
  }
  
  return(paf_annot)
}

ATCCfinder_summarize_query_results <- function(x, queries) {
  "
    __________________________________________________________________________
    Description
      - Iterates through paf results by query, and summarizes annotated results by various criteria
    __________________________________________________________________________
    Variables
      - x         Data frame of paf alignment results
      - queries   Vector of query strings to iterate through
    "
  ## Initialize output df
  cols_summary <- c("qname", "qlen", "n_total_hits", "max_mapq", "max_nmatch", "max_nmatch_%", "n_best_refs", "taxonomy_consensus", "best_refs", "best_refs_taxonomy")
  df_summary <- data.frame(matrix(ncol=length(cols_summary), nrow=0))
  colnames(df_summary) <- cols_summary
  
  for (q in queries) {
    
    ## subset hits to current query
    hits_query <- x[which(x[,"qname"]==q),]
    
    ## Remove rows with error asterisk ("*") target name
    hits_query <- hits_query[which(hits_query[,"tname"]!="*"),]
    
    if ( (nrow(hits_query) > 0) == T ) {
      ## Stats, Query
      query_name   <- q
      query_length <- max(hits_query$qlen)
      
      ## Stats, Target (general)
      n_hits <- nrow(hits_query)
      
      ## Stats, Target [ Longest Reference ]
      longest_ref <- max(hits_query$tlen)
      
      ## Stats, Target [ Longest Match ]
      max_nmatch   <- max(hits_query$nmatch)
      max_nmatch_p <- round( max_nmatch / query_length * 100, 2 )
      
      ## Stats, Target [ Minimum Contigs ]
      min_contigs <- as.numeric( min(hits_query$contig_number) )
      
      ## Stats, Target [ Maximum mapq ]
      max_q <- max(hits_query$mapq)
      
      ## Stats, Target [ Taxonomy Consensus ]
      taxa_hits <- data.frame(table(hits_query$genus_species))
      top_taxa_n <- max(taxa_hits$Freq)
      top_taxa <- taxa_hits[which(taxa_hits$Freq==top_taxa_n),"Var1"]
      
      ## Stats, Target [ Reference Consensus ]
      refs <- data.frame(table(hits_query$atcc_catalog_number))
      # refs$Var1 <- gsub(x=unique(refs$Var1), pattern="ATCC ", replacement="")
      top_ref_n <- max(refs$Freq)
      top_ref <- refs[which(refs$Freq==top_ref_n),"Var1"]
      
      
      
      ## Identify match(es) with highest cumulative score
      max_cumu <- max(hits_query$stat_cum_score)
      cumu <- hits_query[which(hits_query$stat_cum_score==max_cumu),]
      top_atcc <- sort(gsub(x=unique(cumu$atcc_catalog_number),"ATCC ", ""))
      top_hits_n <- length(top_atcc)
      
      cumu_tab <- data.frame(table(cumu$atcc_catalog_number, cumu$taxonomy_raw))
      cumu_tab <- cumu_tab[cumu_tab$Freq!=0,]
      colnames(cumu_tab) <- c("strain", "taxonomy", "n")
      cumu_tab$strain <- gsub(x=cumu_tab$strain, "ATCC ", "")
      cumu_tab <- cumu_tab[rev(order(cumu_tab$n)),]
      cumu_tab$strain_n <- paste0("[", cumu_tab$n, "]", cumu_tab$strain)
      cumu_tab$taxa_n   <- paste0("[", cumu_tab$n, "]", cumu_tab$taxonomy  )
      
      top_atcc_str <- paste0(cumu_tab$strain_n, collapse="; ")
      top_atcc_tax <- paste0(cumu_tab$taxa_n,   collapse="; ")
      
      ## Define unique top references taxonomy
      unique_top_taxa <- paste0(
        unique(unlist(lapply(
          X=as.character(cumu_tab$taxonomy),
          FUN=function(x) paste0(strsplit(x=x, split=" ")[[1]][1:2], collapse=" ") ))),
        collapse="; ")
      
      ## Add information to cumulative df
      df_summary[nrow(df_summary)+1,] <- c(query_name, query_length, n_hits, max_q, max_nmatch, max_nmatch_p, top_hits_n, unique_top_taxa, top_atcc_str, top_atcc_tax)
    }
    
    if ( (nrow(hits_query) > 0) == F ) {
      df_summary[nrow(df_summary)+1,] <- c(q, NA, 0, NA, NA, NA, NA, NA, NA, NA)
    }
    
  }
  return(df_summary)
}

summarize_top_hits <- funcion(x) {
  lapply(X=as.character(cumu_tab$taxonomy), FUN=function(x) paste0(strsplit(x=x, split=" ")[[1]][1:2], collapse=" ") )
}





# ______________________________ ------------------------------------------
# Globals -----------------------------------------------------------------

## Load CLI parameters
args <- argparse()

## Create output directory
make_path( path = args["outdir"] )

## Identify if alignment has already been performed
dir_paf <- file.path(args["outdir"], "alignment.paf")
paf_exists <- check_if_paf_exists( file = dir_paf )





# ______________________________ ------------------------------------------
# ATCCfinder, Search ------------------------------------------------------

## Search query sequence(s) in ATCC genomes
if ( args["overwrite"] == T | paf_exists == F ) {
  
  ATCCfinder_align(
    query  = args["query" ] ,
    target = args["target"] ,
    nhits  = args["nhits" ] ,
    outdir = args["outdir"] ,
    mode   = "minimap2"     )
  
}





# ______________________________ ------------------------------------------
# ATCCfinder, Report ------------------------------------------------------

## Process alignment results
ATCCfinder_summarize_paf(
  query   = args["query" ] ,
  target  = args["target"] ,
  paf     = dir_paf        ,
  dir_out = args["outdir"] )
