# ABOUT

Pipelines for Riboseq and RNA-seq profiling in prokaryotic systems.

# INSTALLATION

In addition to the scripts in this repository, the following third-party open source tools are required. These should be installed on your system, and you will need to set up your environment (PATH variable) to include them:

1) bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
2) samtools (http://www.htslib.org/download/)
3) wigToBigWig (Download binaries from http://hgdownload.soe.ucsc.edu/admin/exe/ – click on your system, then search for "wigToBigWig")

Otherwise this script will also use common linux utilities, such as echo, awk, cat, gzip, and sed.

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

**Output prefix (-o)**. Specific names for the four output files can also be given explicitly using flags -otar, -oman, -otmptar, and -tmpman (respectively, primary output tarball, primary output manifest, intermediate output tarball, intermedia output manifest)

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

### Outputs

This tool generates a master file based on specific treatment vs control contrasts, using results from riboseq processing with the 'Riboseq Sample Processing' tool.

Master file format:

```
GeneID_Treated	Position	codon_relative_rpm	codon_rpm	codon_count	codon_raw_density	GeneID_WT	Position[]	codon_relative_rpm	codon_rpm	codon_count	codon_raw_density	fd_codon_relative_rpm	fd_codon_rpm	fd_codon_count	fd_codon_raw_density	start_AA_cordinate	end_AA_cordinate	seq	rev_comp	AA	gene_start	gene_end
ACT43504.1;lpp	1702921	0.874999999997478	4337.49870908	42	42	ACT43504.1;lpp	1702921	0.666666666664663	1108.92065057	18	18	1.31250000000016	3.91145994697679	2.33333333333333	2.33333333333333	1702894	1702926	CTGGACAACATGGCTACTAAATACCGCAAGAAG		LDNMATKYRKK	1702690	1702926
ACT43504.1;lpp	1702744	0.0721153846151965	1549.10668182	15	15	ACT43504.1;lpp	1702744	0.0199999999998701	616.06702809	10	10	3.60576923078324	2.51450996594107	1.5	1.5	1702717	1702749	GCGGTAATCCTGGGTTCTACTCTGCTGGCAGGT		AVILGSTLLAG	1702690	1702926
ACT43504.1;lpp	1702843	0.0240384615383988	516.36889394	5	5	ACT43504.1;lpp	1702843	0.031999999999987	985.70724495	16	16	0.751201923075268	0.5238562429012	0.3125	0.3125	1702816	1702848	CAGCTGAGCAACGACGTGAACGCAATGCGTTCC		QLSNDVNAMRS	1702690	1702926
```

### Required inputs

**Control and Treatment sample manifests (-m1 and -m2)** These are the manifest output files for the primary outputs of the 'Riboseq Sample Processing' tool. They should have the follow file types in them (third column):

```
Gene_counts
Gene_posRPM
Gene_negRPM
```

**Control and Treatment sample names (-i1 and -i2)** The name of the control and treatment samples. These should match the sample names in the second column of the manifests.

The output master file will calculate fold-changes as treatment/control.

**Reference genome fasta file (-G)** Reference genome in fasta format.

**Reference annotation file (-A)** Gene annotations for your genome.

### Optional inputs

**Cutoff threshold for common genes (-c)** Minimum coverage threshold to include a gene in the master file, across both samples.

**Cutoff threshold for master file, control or treatment sample (-c1 and -c2)** Minimum codon coverage threshold (counts per codon) to include a codon in the master file, controlled separately for control or treatment samples.



## RNA-seq processing

```rna-seq_sample_processing.sh```

Inputs and outputs for this tool are identical to ribo-seq_sample_processing.sh, with the exception of flags -Poff, -Pmin, and -Pmax, which do not exist for this tool.

## Extract files from TAR outputs

```extract_from_tar.sh```

Extract specific files from a tar manifest. Basically a simple wrapper for tar -xf.

# INSTALLATION

To install, please download all scripts in this repository, e.g., with git clone. All supporting scripts noted below are expected to be in the same directory as the main analysis scripts described above.

Please also install the following software package and configure your PATH variable to include them.

## DEPENDENCIES

[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (and bowtie2-build)

[samtools](http://www.htslib.org/download/) (v 1.x)

[bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)

wigToBigWig and bedGraphToBigWig from [UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/) (browse to your Linux or Mac distribution to download these binaries).

## SUPPORTING SCRIPTS

The following scripts have been provided by collaborators, adapted from this reference:

Oh, E., Becker, A.H., Sandikci, A., Huber, D., Chaba, R., Gloge, F., Nichols, R.J., Typas, A., Gross, C.A., Kramer, G., et al. (2011). Selective ribosome profiling reveals the cotranslational chaperone action of trigger factor in vivo. Cell 147, 1295-1308.

Some have been minimally edited, but otherwise used as-is in the above pipelies. The pipelines above expect them to be in the same execution directoy as the pipeline.

```
assignCount3prime_flexible.py
BAM_To_AssignCount_input_format.pl
codon_rpm.pl
common_genes.pl
final.pl
fold_change_master_flexible.pl
GeneCount.py
rpm.pl
```

These scripts I have written, and are used internally in the pipeline. Some details:

```
hard_trim.sh
```

Very simple script to hard trim a given number of bases from the beginning or end of reads in fastq file.

```
split_fastq_by_barcode_in_read.pl
```

Demultiplexes a fastq file based on barcode sequences that are expected to be in the 3' position of each read. Can demultiplex a paired-end library, but demultiplexing is based on parsing R1.

Run without arguments to see usage.


