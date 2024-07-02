# ATCCfinder

**Developed at Los Alamos National Laboratory (O# O4644)**

Download [ATCC Genome Portal](https://genomes.atcc.org/) microbial reference genomes and align query sequences.





<br /> <br />

# Table of Contents
1. [About](#About)
2. [Dependencies](#Dependencies)
    - [Downloading Databases](#Downloading)
    - [Searching Databases](#Search)
3. [Parameters](#Params)
4. [Example Usage](#Example)
5. [Output](#Output)





<br /> <br />

# About <a name="About"></a>

The American Type Culture Collection ([ATCC](https://www.atcc.org/)) sells a wide variety of microbes with strain-level taxonomy classification and associated sequenced reference genomes. ATCCfinder utilizes ATCC application interface software (API) to generate query-able databases from ATCC Genome resources. This tool provides the ability to generate databases of the four ATCC data types:
- Strain-specific genome assembly sequence data (reference)
- Information about how each strain was collected (meta, catalogue)
- Structural/functional information about genome assemblies (annotation).

ATCCfinder contains two core functionalities that may be used in conjunction or independently:
 1. Download ATCC references (with a valid API key)
 2. Query refernce sequences and report alignment results. The tool was built primarily for usage with ATCC refernce genomes, but custom sequence databases may also be searched against.

Once the ATCC reference genome database is retrieved by ATCCfinder, queries may be compared against ATCC reference genomes using the sequence alignment tool minimap2, whose results are then parsed to produce summary data describing what ATCC-available species and strain, if any, the query sequence matches.







<br /> <br />

# Dependencies <a name="Dependencies"></a>

**Downloading Databases** <a name="Downloading"></a>

If you plan to download databases from ATCC yourself, the following is required in an environment:
- [ATCC-Bioinformatics, Genome Portal API](https://github.com/ATCC-Bioinformatics/genome_portal_api)

The following is an example environment created for downloading databases from ATCC:

```
## Download genome_portal_api package
git clone https://github.com/ATCC-Bioinformatics/genome_portal_api.git

## Define environment
mamba create -n ATCCfinder_download
mamba activate ATCCfinder_download

## Install genome_portal_api package
mamba install git pip
pip install /path/to/genome_portal_api
```

Note that you will need to follow ATCC's instructions for generating your own account & API Key, see above Genome Portal API github link for instructions.



<br /> <br />
**Searching References** <a name="Search"></a>

The following software are required for performing alignment to ATCC references:
- [minimap2](https://github.com/lh3/minimap2)
- [R](https://www.r-project.org/)
- [samtools](https://github.com/samtools/samtools)

The following is an example environment created for searching databases from ATCC:

```
## Define environment
mamba create -n ATCCfinder_search
mamba activate ATCCfinder_search

## Install packages
mamba install -c bioconda minimap2
mamba install -c bioconda samtools
mamba install r-base r-argparse
```





<br /> <br />

# Parameters <a name="Params"></a>

**download.py**

| Parameter | Description |
| --- | --- |
| --help | Directory containing subread BAM file(s) |
| --download | Specify which / any ATCC databases to download |
| --overwrite, --no-overwrite | If downloading databases, specifies whether or not downloads in output folder should be kept or replaced |
| --format | Used to specify path to references (if, for example, references were downloaded from my database repository), which will be combined for searching |
| --out | Output folder |



<br /> <br />

**search.R**

| Parameter | Description |
| --- | --- |
| --query | Query Sequence(s) fasta file |
| --target | Target database fasta file |
| --nhits | The maximum number of target hits to return per query|
| --overwrite | (T/F) If alignment output already exists, should it be overwritten? |
| --outdir | Output directory |





<br /> <br />

# Example Usage <a name="Example"></a>

**Download ATCC Reference & Catalogue Databases:**

```
ATCCfinder/download.py \
--download reference catalogue \
--api_key <api_string> \
--overwrite
```

*Note that the combined ATCC reference genome file for use as alignment target is named `atcc_references.fa`*



<br /> <br />

**Align & Report a Query Sequence against ATCC References:**

```
ATCCfinder/search.R \
--target path/to/atcc_references.fa \
--query path/to/query.fa \
--nhits 100 \
--overwrite F \
--outdir path/for/output
```





# Output <a name="Output"></a>

`search.R` will return a file titled `report.tsv` contining the following columns summarizing alignment results:

| Column | Description |
| --- | --- |
| qname | Query sequence name|
| qlen | Query sequence length|
| n_total_hits | Total number of query alignments against reference |
| max_mapq | The maximum mapq score returned from any alignment hit |
| max_nmatch | The maximum number of matched basepairs from any alignment |
| max_nmatch_% | The maximum number of matched basepairs from any alignment, reported as percentage of `qlen` |
| n_best_refs | The number of alignment hits containing *both* `max_mapq` and `max_nmatch` scores |
| best_refs | The reference assembly(s) corresponding to `n_best_refs`, with bracketed values indicating the number of hits associated with this reference |
| best_refs_taxonomy | The taxonomy assignment corresponding to `n_best_refs`, with bracketed values indicating the number of hits associated with this classification |


