# ABOUT

Pipelines for Riboseq and RNA-seq profiling in prokaryotic systems.

# USAGE

NOTE: for all tools, run with no arguments to see usage information.

## Ribo-seq processing

```ribo-seq_sample_processing.sh```

### Outputs

This tool runs demuliplexing and data processing steps for ribosomal profiling of individual samples. It creates four outputs:

1. Tar archive of primary outputs (summary stats, gene counts, and bigWig profiles)
2. Manifest of files and formats in the tar archive
3. Tar archive of intermediate outputs, which are generally not needed in downstream analysis steps
4. Manifest of files and formats in the intermediate tar archive

Manifests list all of the files in each tar archive, along with which sample it's associated with and the file format. For example:

```
Mankin-01.results.tar.gz             All      tar.gz
Mankin-01sample1.stats.txt           sample1  Stats
Mankin-01sample1.genes.txt           sample1  Gene_counts
Mankin-01sample1.posRPM.txt          sample1  Gene_posRPM
Mankin-01sample1.negRPM.txt          sample1  Gene_negRPM
Mankin-01sample1.posRPM.bw           sample1  pos_bw
Mankin-01sample1.negRPM.bw           sample1  neg_bw
```

The first column are file names (e.g., Mankin-01ATCGT.stats.txt), the second column the sample name, and the third column the format/file type. The first row is the tar archive associated with this manifest (Mankin-01.results.tar.gz), all other rows are files within the tar archive.

### Required inputs

**Raw fastq data (-f)** Raw data (fastq or fastq.gz) for multiple samples from one sequencing run. This data will be demultiplexed based on the provided barcodes.

**Barcode list (-b)** List of sample barcodes and the corresponding names. Barcode sequences will be matched against the 3' end of the read.

The first column of this file should be the barcode sequences themselves, and the second column the name you want given to the file. If the second column is omitted, then the barcode sequence itself will be used as the name. Example:

```
ATCGT   Sample1
AGCTA   Sample2
CGTAA   Sample3
CTAGA   Sample4
```

Please use only alphanumeric characters, plus '-', '.', or '_' in the sample names (no spaces, slashes, or other special characters).

**Reference genome fasta file (-G)** Reference genome in fasta format.

**Reference annotation file (-A)** Gene annotations for your genome. This is a 5-column, tab-delimited file, formatted as:

```
[ID]  [start]  [end]  [strand]  [name/description]
```

For example:

```
ACT41903.1;thrL      190     255     +       thr operon leader peptide
ACT41904.1;thrA      336     2798    +       Bifunctional aspartokinase/homoserine dehydrogenase 1
ACT41905.1;thrB      2800    3732    +       homoserine kinase
ACT41906.1;thrC      3733    5019    +       L-threonine synthase
ACT41907.1;yaaX      5232    5528    +       DUF2502 family putative periplasmic protein
ACT41908.1;yaaA      5681    6457    -       peroxide resistance protein%2C lowers intracellular iron
```

Note that multi-chromosome genomes are not supported by this pipeline, as chromosome names are not included in this annotation file.

**Output prefix (-o)**

### Optional inputs

**Filter sequences fasta file (-F)** Reference sequences you want to exclude from downstream analysis, e.g., rRNAs, tRNAs, etc. Any reads mapping to these sequences will be excluded from downstream processing. This step can be omitted by not specifying a fasta file.

**Minimum number of reads per barcode (-min)** This parameter prevents any barcodes with very low coverage from being analyzed downstream, as any barcodes with fewer than this many reads will be ignored in downstream processing. Otherwise this tool will attempt to recover ALL of the samples in your barcode list, but if a sample fails or you include too many barcodes, you may still have a small number of reads for that sample, but not enough to consider analyzing.

**Run demultiplexing flexibly? (-flex)** If set to "Yes", will search for barcodes in the entire read, rather than just at the end; whichever barcode match is closest to the 3' of the read will be used. May be useful if reads are short and not enough adapter is included to be trimmed automatically. NOTE: use of this option will increase the compute time, as the barcode search will take longer.

Raw demultiplexed reads will appear in the secondary output, marked as "fastq" in the manifest.

**Bases to trim from 5' end of read (-hard5)** After demultiplexing, how many bases should be trimmed from the 5' end.

**Bases to trim from 3' end of read (-hard3)** After demultiplexing, how many bases should be trimmed from the 3' end.

Trimmed reads will appear in the secondary output, marked as "fastq-trimmed in the manifest.

**Bowtie2 seed length for filter sequences alignment (-Fseed)** Seed length to use in Bowtie2 alignment against the filter sequences (if provided).

Reads that pass this filtering step will appear in the secondary output, marked as "fastq-filtered" in the manifest.

**Bowtie2 seed length for genome sequences alignment (-Gseed)** Seed length to use in Bowtie2 alignment against the reference genome.

Mapped reads will appear in the secondary output, marked as "BAM" in the manifest.

In general, shorter seeds will provide more sensitive alignments, and will work better for short reads or short reference sequences, with some increase in execution time. The default setting has a shorter seed for the filter sequences in anticipation that many of those may be very short (e.g., tRNAs).

**Maximum mismathces allowed in alignment (-mis)** Any alignments with more than this many mismatches are discarded. They will be included in the .bam file, but not in the .assignCount.txt intermediate file (type "AssignCount" in the secondard output manifest).

**P-site offset (-Poff)** Offset to use from 3' end of read to find P-site.

**Minimum read length for P-site offset (-Pmin)** Reads shorter than this will be excluded from the analysis.

**Maximum read length for P-site offset (-Pmax)** Reads longer than this will be excluded from the analysis.

**Specify alternate path to look for dependencies in (-path)**



## Generate master table

```generate_master.sh```



## RNA-seq processing

```rna-seq_sample_processing.sh```



## Extract files from TAR outputs

```extract_from_tar.sh```



# SUPPORTING FILE FORMATS


# INSTALLATION


## DEPENDENCIES


## SUPPORTING SCRIPTS


