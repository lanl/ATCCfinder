#!/usr/bin/env Rscript

# ______________________________ ------------------------------------------
# Parameters --------------------------------------------------------------

args_rstudio <- function() {
  "
  ______________________________________________________________________________
  Description
    - Define parameters here if running script from RStudio.
  "
  
  set <- c("atcc", "staph")[2]
  
  ## Define Arguments
  if ( set == "atcc" ) {
    args <- list(
      query         = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/querying_a_custom_database/data/ncbi_Zymo_MCS_amplicons.fasta",
      target        = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/querying_a_custom_database/atcc_downloads/atcc_references.fa",
      target_meta   = "",
      nhits         = "100",
      note_nmatch_p = "95",
      overwrite     = "T",
      outdir        = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/querying_a_custom_database/data/output_whitelist" )
  }
  if ( set == "staph" ) {
    args <- list(
      query         = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/query_custom_database_Staphylococcus/MH798864.1.fasta",
      target        = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/query_custom_database_Staphylococcus/ncbi_genome_Staphylococcus.fa.gz",
      target_meta   = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/query_custom_database_Staphylococcus/ncbi_custom_target_database_Staphylococcus.tsv",
      nhits         = "100",
      note_nmatch_p = "95",
      overwrite     = "T",
      outdir        = "/Users/skoehler/Desktop/project_2302_MO/Writing/Publication_ATCCfinder/Journal_Revisions/figures/query_custom_database_Staphylococcus/output_Staphylococcus" )
  }
  
  ## Format parameters
  args <- format_parameters(args)
  args["note_nmatch_p"] <- as.numeric( args["note_nmatch_p"] )
  return(args)
}


args_cli <- function() {
  "
  ____________________________________________________________________________
  Description
    - Define and retrieve command-line input values
  "
  ## Packages
  suppressWarnings( library(argparse) )
  
  ## Define argument parser
  parser <- ArgumentParser(description='Arguement Parser')
  
  ## Define Arguments
  parser$add_argument("--query"        , required=T, help="Query Sequence(s) fasta file")
  parser$add_argument("--target"       , required=T, help="Target database fasta file")
  parser$add_argument("--target_meta"  , required=F, help="Required if using a custom target database. Provides meta data associated with sequences in a standardized format. Must contain columns 'header_id', 'reference_id', and 'taxonomy'. Additional columns may be included.", default="")
  parser$add_argument("--nhits"        , required=T, help="The maximum number of target hits to return per query")
  parser$add_argument("--note_nmatch_p", required=F, help="(0-100) Specify an alignment nmatch percentage at or above which results will be noted as potential matches despite their aggregate score not being the maximum.", default="95")
  parser$add_argument("--overwrite"    , required=F, help="(T/F) If alignment output already exists, should it be overwritten?", default="T")
  parser$add_argument("--outdir"       , required=T, help="Output directory")
  
  ## Load and format inputs
  args <- parser$parse_args()
  
  ## Format parameters
  args <- format_parameters(args)
  args["note_nmatch_p"] <- as.numeric( args["note_nmatch_p"] )
  return(args)
}

format_parameters <- function(args) {
  "
  ____________________________________________________________________________
  Description
    - Format input parameters to proper object types
  "
  
  args["nhits"    ] <- as.numeric(args["nhits"])
  args["target_meta"] <- ifelse( test = args["target_meta"] == "" , yes = F, no = args["target_meta"] )
  args["overwrite"  ] <- ifelse( test = args["overwrite"  ] == "T", yes = T, no = F)
  
  return(args)
}





# ______________________________ ------------------------------------------
# ATCCfinder --------------------------------------------------------------

ATCCfinder <- function() {
  "
  ______________________________________________________________________________
  Description
    - Manage argument parsing and main execution of software
  "
  
  ## _________________________________________________________________________
  ## [- Load parameters -]
  
  identify_and_set_wd_either()
  
  if ( interactive() == T ) { args <- args_rstudio() } ## RStudio
  if ( interactive() == F ) { args <-     args_cli() } ## Command line
  
  
  
  ## _________________________________________________________________________
  ## [- Define globals -]
  
  ## Files & directories
  dir_paf        <- file.path(args["outdir"], "alignment.paf")
  dir_db_headers <- file.path(args["outdir"], "database_headers.txt")
  dir_paf_anot   <- file.path(args["outdir"], "alignment_annotated.tsv")
  dir_summary    <- file.path(args["outdir"], "ATCCfinder_summary.tsv" )
  
  dir_query       = args["query"]
  dir_target      = args["target"]
  dir_target_meta = args["target_meta"]
  dir_out         = args["outdir"]
  nhits           = args["nhits"]
  note_nmatch_p   = args["note_nmatch_p"]
  
  
  
  ## _________________________________________________________________________
  ## Format output
  
  ## Create output directory
  make_path( path = args["outdir"] )
  
  
  
  ## _________________________________________________________________________
  ## ATCCfinder
  
  ## Search query sequence(s) against target database
  paf_exists <- check_if_paf_exists( file = dir_paf ) ## Identify if alignment has already been performed
  if ( args["overwrite"] == T | paf_exists == F ) { ## THIS LOGIC CHECK SHOULD BE PLACED WITHIN THE PROCEEDING FUNCTION
    
    ATCCfinder_align(
      query  = dir_query ,
      target = dir_target ,
      nhits  = nhits ,
      outdir = dir_out )
    
  }
  
  ## Annotate alignment results, identifying top hits
  ATCCfinder_annotate_alignment_stats(
    dir_paf      = dir_paf      ,
    dir_paf_anot = dir_paf_anot )
  
  ## Annotate alignment results, adding taxonomy information
  ATCCfinder_annotate_alignment_meta(
    dir_paf_anot    = dir_paf_anot    ,
    dir_target_meta = dir_target_meta ,
    dir_db_headers  = dir_db_headers  ,
    dir_target      = dir_target      )
  
  ## Summarize annotated alignment results for each query
  ATCCfinder_summarize_query_results(
    dir_paf_anot  = dir_paf_anot  ,
    dir_summary   = dir_summary   ,
    note_nmatch_p = note_nmatch_p )
  
  # ## Annotate alignments with target database meta data (nearly done)
  # ATCCfinder_annotate_paf(
  #   dir_paf         = dir_paf         ,
  #   dir_target      = dir_target      ,
  #   dir_target_meta = dir_target_meta ,
  #   dir_out         = dir_out         )
  
}





# ______________________________ ------------------------------------------
# Functions, General ------------------------------------------------------

identify_and_set_wd_either <- function(n_back=0) {
  "
  ____________________________________________________________________________
  Description
   - Identify and set the working directory when running code from command line
  ____________________________________________________________________________
  Parameters
   - nback  The number of directories to step back into when setting the working directory
  "
  ## Scenario: RStudio execution
  if ( interactive() == T ) {
    path_current <- strsplit(dirname(rstudioapi::getSourceEditorContext()$path), "/")[[1]] ## Identify the directory of the current R-Studio script
    path_base    <- Reduce(file.path, path_current[ 1:length(path_current) - n_back ])     ## Define the directory to be set
    setwd(path_base)                                                                       ## Set the working directory
  }
  
  ## Scenario: Command line execution
  if ( interactive() == F ) {
    initial.options <- commandArgs(trailingOnly=FALSE)
    file.arg.name   <- "--file="
    path_current    <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)]) ## Identify the directory of the current R-Studio script
    path_base       <- Reduce(file.path, path_current[1:length(path_current)-n_back])                ## Define the directory to be set
    setwd(system("pwd", intern=T) )
  }
  
  ## Parse data and return
  cat(
    "\n________________________________________\n\n\t",
    paste0('> Directory set: ', path_base))                                            ## Report the active directory
  
  return(path_base)                                                                    ## Return string of working directory path
  
}

make_path <- function(path) {
  dir.create(path=file.path(path), recursive=T, showWarnings=F)
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
  
  ## Remove raw column
  hits[,"raw"] <- NULL
  
  return(hits)
}





# ______________________________ ------------------------------------------
# Functions, ATCCfinder ---------------------------------------------------

ATCCfinder_align <- function(query, target, nhits, outdir) {
  "
  ____________________________________________________________________________
  Description
    - Wrapper function for running minimap2
  ____________________________________________________________________________
  Variables
    - query
    - target
    - outdir
  ____________________________________________________________________________
  Notes
    - Requires an environment with minimap2 installed
  "
  ## Align using minimap2
  message("\nAligning query(s) to ATCC references via minimap2..\n")
  outfile <- file.path(outdir, "alignment.paf")
  system( paste("minimap2", "-c", "-N", as.numeric(nhits)-1, "-x", "map-ont", "--paf-no-hit", target, query, ">", outfile) )
  message("\nAlignment complete, paf file generated.\n")
}

ATCCfinder_annotate_alignment_stats <- function(dir_paf, dir_paf_anot) {
  "
  ____________________________________________________________________________
  Description
    - Loads alignment results, scores hits, and identifies best reference match(es).
  ____________________________________________________________________________
  Variables
    - query       Fasta file that was queried to produce paf hit results
    - paf         Minimap2 alignment file results comparing user-specified query to reference target
    - outdir      Output directory
  "
  
  ## Load paf alignment results
  hits <- load_format_paf(file=dir_paf)
  
  ## Identify queries
  queries <- unique(hits[,"qname"])
  
  n = 0
  for (q in queries) {
    n <- n + 1
    if (n%%100==0) {print(n)}
    ## subset hits to current query
    hits_query <- hits[which(hits[,"qname"]==q),]
    
    if ( nrow(hits_query) > 0 ) {
      ## Stats, Query
      query_name   <- q
      query_length <- max(hits_query$qlen)
      
      ## Stats, Target (general)
      n_hits <- nrow(hits_query)
      
      ## Stats, Target [ Longest Match ]
      max_nmatch   <- max(hits_query$nmatch)
      max_nmatch_p <- round( max_nmatch / query_length * 100, 2 )
      
      ## Stats, Target [ Maximum alignment score ]
      max_AS <- max(hits_query$AS)
      
      ## Add column for percent match
      hits_query[,"nmatch_p"] <- round( hits_query$nmatch / hits_query$qlen * 100, 2 )
      
      ## Score hits based on metrics
      hits_query[,"stat_max_AS"    ] <- 0
      hits_query[,"stat_max_nmatch"] <- 0
      
      hits_query[which(hits_query[,"AS"     ] == max_AS     ), "stat_max_AS"    ] <- 1
      hits_query[which(hits_query[,"nmatch" ] == max_nmatch ), "stat_max_nmatch"] <- 1
      
      ## Define cumulative score
      hits_query[,"stat_cum_score"] <- hits_query[,"stat_max_AS"  ] + hits_query[,"stat_max_nmatch"]
      
      ## Add information to compiled df
      ifelse(
        test = exists("paf_anot")                        ,
        yes  = paf_anot <- rbind( paf_anot, hits_query) ,
        no   = paf_anot <- hits_query                    )
    }
  }
  
  ## Export annotated alignment results
  write.table(x=paf_anot, file=dir_paf_anot, sep="\t", row.names=F)
  
}

ATCCfinder_annotate_alignment_meta  <- function(dir_paf_anot, dir_target_meta, dir_db_headers, dir_target) {
  
  ## Load data
  paf <- read.delim(dir_paf_anot) ## Annotated alignment results
  
  ## Target database meta data
  if ( dir_target_meta[[1]] == F ) { meta <- load_target_metadata_atcc(dir_target, dir_db_headers) }
  if ( dir_target_meta[[1]] != F ) { meta <- read.delim(dir_target_meta[[1]])          }
  
  ## Ensure that target sequence names match meta data table entries
  paf_queries <- paf[which(paf$tname!="*"),"qname"]
  n_queries <- length(paf_queries)
  n_queries_in_meta <- sum(paf[which(paf$tname!="*"),"tname"]%in%meta$header_id, na.rm=T)
  dif <- n_queries - n_queries_in_meta
  if ( n_queries_in_meta == n_queries) { message("All target headers successfully found in target meta data table.") }
  if ( n_queries_in_meta != n_queries) { message(paste("[-WARNING-]", dif, "/", n_queries, "target headers missing from meta data table. Check formatting to ensure meta data headers match target database headers.")) }
  
  ## Merge & export
  paf_meta <- merge(x=paf, y=meta, by.x="tname", by.y="header_id", all.x=T)
  write.table(x=paf_meta, file=dir_paf_anot, sep="\t", row.names=F)
}

# ## Summarize Query results
# result_summary <- ATCCfinder_summarize_query_results(
#   x       = hits_anot ,
#   queries = queries   )
# 
# ## Export parsed results
# outfile_summary <- file.path(dir_out, "report_v2.tsv")
# write.table(x=result_summary, file=outfile_summary, sep="\t", row.names=F, quote=F)
# 
# 
# 

ATCCfinder_summarize_query_results <- function(dir_paf_anot, note_nmatch_p, dir_summary) {
  "
    __________________________________________________________________________
    Description
      - Iterates through paf results by query, and summarizes annotated results by various criteria
    __________________________________________________________________________
    Variables
      - x         Data frame of paf alignment results
      - queries   Vector of query strings to iterate through
    "
  
  ## Load annotated alignment results
  x <- read.delim(dir_paf_anot)
  
  ## Parse query sequence names
  queries <- unique(x$qname)
  
  ## Initialize output df
  cols_summary <- c("qname", "qlen", "n_total_hits", "max_mapq", "max_nmatch", "max_nmatch_%", "n_best_refs", "taxonomy_consensus", "best_refs", "best_refs_taxonomy", "note", "note_n_best_refs", "note_taxonomy_consensus", "note_best_refs", "note_best_refs_taxonomy")
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
      query_length <- unique(hits_query$qlen)
      
      ## Stats, Target (general)
      n_hits <- nrow(hits_query)
      
      ## Stats, Target [ Longest Reference ]
      longest_ref <- max(hits_query$tlen)
      
      ## Stats, Target [ Longest Match ]
      max_nmatch   <- max(hits_query$nmatch)
      max_nmatch_p <- round( max_nmatch / query_length * 100, 2 )
      
      ## Stats, Target [ Maximum mapq ]
      max_q <- max(hits_query$mapq)
      
      ## Stats, Target [ Taxonomy Consensus ]
      taxa_hits <- data.frame(table(hits_query$taxonomy))
      top_taxa_n <- max(taxa_hits$Freq)
      top_taxa <- taxa_hits[which(taxa_hits$Freq==top_taxa_n),"Var1"]
      
      ## Stats, Target [ Reference Consensus ]
      refs <- data.frame(table(hits_query$reference_id))
      top_ref_n <- max(refs$Freq)
      top_ref <- refs[which(refs$Freq==top_ref_n),"Var1"]
      
      
      
      ## Identify match(es) with highest cumulative score
      max_cumu <- max(hits_query$stat_cum_score)
      cumu <- hits_query[which(hits_query$stat_cum_score==max_cumu),]
      top_atcc <- sort(gsub(x=unique(cumu$reference_id),"ATCC ", ""))
      top_hits_n <- length(top_atcc)
      
      cumu_tab <- data.frame(table(cumu$reference_id, cumu$taxonomy))
      cumu_tab <- cumu_tab[cumu_tab$Freq!=0,]
      colnames(cumu_tab) <- c("strain", "taxonomy", "n")
      cumu_tab$strain <- gsub(x=cumu_tab$strain, "ATCC ", "")
      cumu_tab <- cumu_tab[rev(order(cumu_tab$n)),]
      cumu_tab$strain_n <- paste0("(", cumu_tab$n, ") ", cumu_tab$strain)
      cumu_tab$taxa_n   <- paste0("(", cumu_tab$n, ") ", cumu_tab$taxonomy  )
      
      top_atcc_str <- paste0(cumu_tab$strain_n, collapse="; ")
      top_atcc_tax <- paste0(cumu_tab$taxa_n,   collapse="; ")
      
      ## Define unique top references taxonomy
      unique_top_taxa <- paste0(
        unique(unlist(lapply(
          X=as.character(cumu_tab$taxonomy),
          FUN=function(x) paste0(strsplit(x=x, split=" ")[[1]][1:2], collapse=" ") ))),
        collapse="; ")
      
      
      
      ## Identify if any noteworthy hits exist without top score (No top score but have a large percentage of matching bases within alignment)
      not_cumu <- hits_query[which( hits_query$stat_cum_score != max_cumu & hits_query$nmatch_p >= note_nmatch_p ),]
      
      if ( (nrow(not_cumu) > 0) == T ) {
        
        note = paste0("NOTE: Alignments with non-maximum scores but high nmatch percentages (>=", note_nmatch_p, "%) exist. See following columns and output file 'alignment_annotated.tsv' for more details.")
        
        not_cumu <- not_cumu[which(not_cumu$nmatch_p>=note_nmatch_p),]
        not_cumu_top_atcc <- sort(gsub(x=unique(not_cumu$reference_id),"ATCC ", ""))
        not_cumu_top_hits_n <- length(not_cumu_top_atcc)
        
        not_cumu_tab <- data.frame(table(not_cumu$reference_id, not_cumu$taxonomy))
        not_cumu_tab <- not_cumu_tab[not_cumu_tab$Freq!=0,]
        colnames(not_cumu_tab) <- c("strain", "taxonomy", "n")
        not_cumu_tab$strain <- gsub(x=not_cumu_tab$strain, "ATCC ", "")
        not_cumu_tab <- not_cumu_tab[rev(order(not_cumu_tab$n)),]
        not_cumu_tab$strain_n <- paste0("(", not_cumu_tab$n, ") ", not_cumu_tab$strain)
        not_cumu_tab$taxa_n   <- paste0("(", not_cumu_tab$n, ") ", not_cumu_tab$taxonomy  )
        
        not_cumu_top_atcc_str <- paste0(not_cumu_tab$strain_n, collapse="; ")
        not_cumu_top_atcc_tax <- paste0(not_cumu_tab$taxa_n,   collapse="; ")
        
        ## Define unique top references taxonomy
        not_cumu_unique_top_taxa <- paste0(
          unique(unlist(lapply(
            X=as.character(not_cumu_tab$taxonomy),
            FUN=function(x) paste0(strsplit(x=x, split=" ")[[1]][1:2], collapse=" ") ))),
          collapse="; ")
      }
      
      if ( (nrow(not_cumu) > 0) == F ) {
        
        note = paste0("No alignments with non-maximum scores with noteworthy nmatch percentages (>=", note_nmatch_p, "%) to report.")
        
        not_cumu_top_hits_n      = "-"
        not_cumu_unique_top_taxa = "-"
        not_cumu_top_atcc_str    = "-"
        not_cumu_top_atcc_tax    = "-"
      }
      
      
      
      
      
      ## Add information to cumulative df
      df_summary[nrow(df_summary)+1,] <- c(query_name, query_length, n_hits, max_q, max_nmatch, max_nmatch_p, top_hits_n, unique_top_taxa, top_atcc_str, top_atcc_tax, note, not_cumu_top_hits_n, not_cumu_unique_top_taxa, not_cumu_top_atcc_str, not_cumu_top_atcc_tax)
    }
    
    if ( (nrow(hits_query) > 0) == F ) {
      df_summary[nrow(df_summary)+1,] <- c(q, NA, 0, NA, NA, NA, NA, NA, NA, NA)
    }
    
  }
  # return(df_summary)
  write.table(x=df_summary, file=dir_summary, sep="\t", quote=F, row.names=F)
}


# ATCCfinder_annotate_paf <- function(dir_paf, dir_target, dir_target_meta, dir_out) {
#   "
#   ____________________________________________________________________________
#   Description
#     - Annotate paf results with information about target sequences.
#   ____________________________________________________________________________
#   Variables
#     - dir_paf           Minimap2 alignment file results comparing user-specified query to reference target
#     - dir_target_meta   Provides meta data associated with sequences in a standardized format.
#     - dir_out           Output directory for all analysis results.
#   "
#   
#   ## Load paf alignment results
#   paf <- load_format_paf(file=dir_paf)
#   
#   ## Extract target database headers
#   parse_target_headers(dir_target=dir_target, dir_out=dir_out)
#   
#   ## Load target database meta data
#   #### If using ATCC database then automatically parse information from headers
#   if ( dir_target_meta == "" ) {
#     load_target_metadata_atcc(dir_db_headers=dir_db_headers)
#   }
#   
#   ## Load target database meta data
#   #### If using custom database, then load specified standardized metadata file
#   if ( dir_target_meta != "" ) {
#     load_target_metadata_custom(dir_target_meta=dir_target_meta)
#   }
#   
#   
# }

load_target_metadata_atcc <- function(dir_target, dir_db_headers) {
  "
    __________________________________________________________________________
    Description
      - Takes ATCC reference header strings and parses them for ATCC-related meta data
    __________________________________________________________________________
    Variables
      - file    A tsv file of fasta headers (including '>' starting character) taken from ATCC references
  "
  
  ## Identify infile directory
  # dir_target_headers <- file.path(dirname(dir_target[[1]]), "atcc_references_headers.txt")
  if ( grepl(x=dir_target[[1]], pattern="[.]gz") == TRUE  ) { system(paste( "gunzip -c", dir_target, "| grep '^>' >", dir_db_headers )) }
  if ( grepl(x=dir_target[[1]], pattern="[.]gz") == FALSE ) { system(paste( "cat"      , dir_target, "| grep '^>' >", dir_db_headers )) }
  
  ## Load file
  df <- read.delim(file=dir_db_headers, sep="\t", header=F)
  colnames(df) <- "raw"
  
  ## Parse reference information from header
  df[,"header_id"          ] <- gsub(x=df[,"raw"], pattern=paste0(c(">"                     , " .+"                    ), collapse="|"), replacement="")
  df[,"reference_id"       ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+assembly_id="        , " genome_id=.+"          ), collapse="|"), replacement="")
  df[,"taxonomy"           ] <- gsub(x=df[,"raw"], pattern=paste0(c(".+species="            , " assembly_id=.+"        ), collapse="|"), replacement="")
  df[,"raw"] <- NULL
  
  ## Parse taxonomy data
  df[,"taxonomy"] <- gsub(x=df[,"taxonomy"], pattern="\\[|\\]", replacement="")
  df[,"taxonomy"] <- tolower(df[,"taxonomy"])
  df[,"taxonomy"] <- paste0(toupper( substring(df[,"taxonomy"], 1, 1                    )) ,
                                     substring(df[,"taxonomy"], 2, nchar(df[,"taxonomy"])) )
  
  
  # 
  # df[,"n_taxname"]     <- unlist(lapply(X=df[,"taxonomy_raw"], FUN=function(x) length(strsplit(x, split=" ")[[1]])))
  # df[,"genus"] <- gsub(x=df[,"taxonomy_raw"], pattern=" .+", replacement="")
  # df[,"genus_species"] <- unlist(lapply(X=df[,"taxonomy_raw"], FUN=function(x) paste( strsplit(x, split=" ")[[1]][1:2], collapse=" ")))
  # 
  ## Remove duplicate rows
  df <- df[which(!duplicated(df)),]
  
  return(df)
}

load_target_metadata_custom <- function(dir_target_meta) {
  df <- read.delim(file=dir_target_meta[[1]], sep="\t", header=1)
}





# ______________________________ ------------------------------------------
# Main --------------------------------------------------------------------

ATCCfinder()




