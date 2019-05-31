# Welcome to the mmannot Web page.

## What does mmannot do?

mmannot *annotates reads*, or *quantifies the features*.

For instance, suppose that you have sequenced your organism of interest with sRNA-Seq (RNA-Seq works too), and you want to know how many times you have sequenced miRNAs, rRNAs, tRNAs, etc.
This is what mmannot does.

## Why using this tool?

A huge proportion of the reads may actually map at several locations.
These multi-mapping reads are usually handled poorly by similar quantification tools.
In our methods, when a read maps at several locations, all these locations are inspected:

- If all these locations belong to the same feature (e.g. miRNAs, in case of a duplicated gene family), the read is still annotated as a miRNA.
- If the location belong to different features (e.g. 3'UTR and miRNA), the read is ambiguous, and is flagged as 3'UTR--miRNA.

In case 1, we say when have *rescued* a read.

## Why it matters?

* The method is unbiased, and extract all the information given by sequencing.
* A large proportion of the reads can be rescued.
* A read mapping two different features can help correcting the annotation, or finding a potential interaction between the features.


## Install

* Download the tool with git: `git clone https://github.com/mzytnicki/mmannot.git`
* The bundle contains:
  * the C++ file: `mmannot.cpp`,
  * a static build (if you cannot compile it): `mmannot_static`
  * the companion tool (see [*infra*](#add-the-nh-tag)) source code and static build: `addNH.cpp` and `addNH`,
  * several configuration files (see [*infra*](#configuration-file)), used for the publication.
* Compile everything with `make`.
* You will need a C++11 compiler, and zlib.


## How to use it?

Here are the options:

Compulsory options:

* `-a` *file*: annotation file in GTF format
* `-r` *file* [*file2*] ...: reads in BAM/SAM format

Main options:

* `-o` *output*: output file (see [*infra*](#output-file)). Default: standard output.
* `-n` *name1* [*name2*] ...: short name for each of the reads files
* `-c` *config*: configuration file (see [*infra*](#configuration-file)). Defaut: `config.txt`.
* `-s` *strand*: can be:
  * for paired-end reads: `U` (unknown), `FR` (forward-reverse), `RF` (reverse-forward), `FF` (forward-forward);
  * for single-end reads: `U` (unknown), `F` (forward), `R` (reverse);
  * Default: `F`.
* `-f` format (`SAM` or `BAM`): format of the read input files
* `-l` number overlap type (see [*infra*](#overlap)). Default: `-1`.

Configuration options:

* `-d` *integer*: upstream region size (default: 1000)
* `-D` *integer*: downstream region size (default: 1000)
* `-y` *strategy*: quantification strategy, mainly for comparison purposes.  Valid values are: `default`, `unique`, `random`, `ratio` (default: `default`)
* `-e` *integer*: attribute a read to a feature if at least *N*% of the hits map to the feature (default: 100%)

Output options:

* `-p`: print progress
* `-t` *integer*: number of threads (one thread per SAM/BAM file)
* `-m` *file*: mapping statistics of the reads (see [*infra*](#read-statistics)) 
* `-M` *file*: mapping statistics of the annotation (see [*infra*](#annotation-statistics))
* `-h`: help

### Annotation file

The annotation file should be in GTF. GFF might work too.  The feature should be provided in the second and/or the third field.

For instance, the following annotations are valid:

    chr1  miRAnnot        miRNA  100   120   .  +  .  ID "miRNA_1";
    chr1  geneAnnot       gene   1000  2000  .  +  .  ID "gene_1";
    chr1  protein_coding  gene   1000  2000  .  +  .  ID "gene_1";

However, this is not valid:

    chr1  annotationTool  gene   1000  2000  .  +  .  type "gene"; ID "gene_1";

Since the feature is provided in the last field.

### Configuration file

The annotation file describes which features you would like to quantify, and how.

#### Skeleton

The skeleton of a configuration is:

    Synonyms:
      (synonyms here)
    Introns:
      (introns here)
    Vicinity:
      (synonyms here)
    Order:
      (order here)

Replace text in parenthesis as desired.

#### Synonyms

The `synonym` section simply provides the list of synonyms.
It can be skipped.
For example, a transcript can be referred to as `mRNA` or `transcript`.
A 5' UTR can be `five_prime_UTR` or `5'UTR`.
You can state that both terms are equivalent, and should be handled similarly.

    Synonyms:
      mRNA:transcript
      five_prime_UTR:5'UTR
      three_prime_UTR:3'UTR
  
#### Introns

Introns are not always explicitly mentioned in an annotation file.
They are instead implicitly defined as the intervals that lie between consecutive exons of the same transcript.
If you want to quantify these introns this way.

  Introns:
    protein_coding:gene
    lincRNA:gene

Here, in `protein_coding:gene`, the first part (`protein_coding`) refers to an element of the second field of the GTF file, whereas the last part (`gene`) should be found in the third field of the GTF file.
You can then quantify these introns.
In this example, they can be referred to as `protein_coding:intron` and `lincRNA:intron` respectively.

#### Upstream and downstream regions

You can also include some upstream and downstream regions in you quantification.
The size of these regions should be provided as parameter of the program (`-d` and `-D` parameters).

    Vicinity:
      protein_coding:gene

This will create the following features: `protein_coding:upstream` and `protein_coding:downstream`.

#### Order

This is the main section.
Provide here the feature you want to quantify.
If you type:

    Order:
      protein_coding:CDS

You will count all the reads that overlap an annotation described this way:

    chr1  protein_coding  gene  1000  2000  .  +  .  ID=gene_1

If you want to count all the reads that overlap annotation, and that are located on the same strand, type something like:

    Order:
      protein_coding:CDS +

Notice that the order counts.  Suppose that you have a this configuration file:

    Order:
      protein_coding:CDS
      miRNA:primary_transcript
  
If a read maps on the CDS and on the miRNA (because the reads maps at two different locations, or two annotations overlap), it will be annotated as CDS only.
If you prefer to keep the ambiguity, put the two annotations on the same line (separated by a comma `,`):

    Order:
      protein_coding:CDS, miRNA:primary_transcript


#### Configuration file maker

If you have trouble creating the configuration, you can use a helper tool, `createConfigFile` (which requires Python >=2.5).

Type: 

    ./createConfigFile -i annotation.gtf -o config.txt

where `annotation.gtf` is your annotation file, and `config.txt` is the configuration file (created by `createConfigFile`).

Then, follow the instructions, which should be (hopefully) clear enough for you to create the configuration file, based on your annotation.


### Reads files

The reads should be given in SAM or BAM format.
The reads can be single end or paired-end.
mmannot must be informed about the number of times a read maps.
The usual read mapping tools can use two strategies to do so.

#### Reads mapped with BWA

BWA provides all the possible hits in the last field, with the `XA` tag.
mmannot works directly with the reads mapped by BWA.


#### Reads mapped with Bowtie or Bowtie2

Bowtie does not provides the number of times a read maps.
You should then add the `NH` flag, that provides the number of hits for each read.
There is a companion tool, `addNH`, that add this flag __on unsorted reads only__ (i.e. right after the mapping).
See [*infra*](#add-the-nh-tag) for more about this tool.

You should also check how your mapping tool handles multi-mapping reads (this can usually be tuned using the appropriate parameters).


### Output file

The output is a tab-separated file, and be used by EdgeR or DESeq, for instance. If the user provided *n* reads files, the output will contain *n*+1 columns:

    Gene                  sample_1  sample_2  ...
    feature_A             ...       ...       ...
    feature_B             ...       ...       ...
    feature_B--feature_C  ...       ...       ...

The first line is the list of the features. If a read maps several features (say, feature_B and feature_C), a new feature is added to the table, `feature_B--feature_C`. The reads that can be mapped to these features will be counted there (but not in the `feature_B` nor `feature_C` lines).


### Output stats

The output stats are given in standard error.

The general shape is

    Results for sample_A:
          # reads:                       N
          # uniquely mapped reads:       N (x%)
          # multi-mapping rescued reads: N (x%)
          # hits:                        N
          # ambiguous hits:              N (x%)
          # unassigned hits:             N (x%)

mmannot first provides the number of reads used, and the number of uniquely mapping reads.
The number of rescued reads is also provided.
Then, mmannot provides the number of hits (a read may have zero, one, or several hits).
Ambiguous hits and unassigned hits (i.e. not mapping any feature) is provided.


#### Read statistics

This file provides the list of features attributed to each reads.
The format is:

    read_1  N  feature_A: N_A    feature_B: N_B    ...

The description of each field follows:

* (`read_1`) the name of the read
* (`N`) the number of hits for this read
* (`feature_A: N_A`) the number of times the read overlaps an element of the feature A.

Note that some hits may be located away from any annotation, so *N_A + N_B + ... < = N*.


#### Annotation statistics

This file provides the number of reads per annotation.
The format is:

    ENSGXXXXX (feature_A)  N1
    ENSGXXXXX (feature_A)--ENSGYYYYY (feature_B)  N2

The first field provide the annotation id (`ENSGXXXXX`), its feature type (`feature_A`), and the number of reads matching this annotation (`N1`).
If some reads map several annotations, these annotation are merged (`ENSGXXXXX--ENSGYYYYY`) and the number of reads mapping these annotations is provided as well (`N2`).

#### Overlap

The way a read *R* is mapped to a feature *A* depends on the `-l` *n* value:

| if *n* is            | then *R* is mapped to *A* iff                        |
|----------------------|------------------------------------------------------|
| a negative value     | *R* is included in *A*                               |
| a positive integer   | they have at least *n* nucleotides in common         |
| a float value (0, 1) | *n*% of the nucleotides of *R* are shared with *A*   |


#### Add the `NH` tag

mmannot needs the `NH` tag in the SAM/BAM file to work.
This tag provide the number of hits for each read.
Unfortunately, some tools, such as Bowtie do not provide it, but it is not difficult to add it.
We provided a simple tool, `addNH`, to add it __in the unsorted SAM file__.

For instance, if you use bowtie1, you can use it this way:

    bowtie -S -a --best --strata database file.fastq | addNH > file.sam

`addNH` is provided as a statically compiled executable (thus usable on most Linux machines).  You can also compile the tool yourself: the code is provided.  Simply type `make addNH`.


## Troubleshooting

### Test data set

For your conveniance, test BAM and GTF files are available.

The GTF file is the annotation of the Y chromosome of the human chromosome, build 38.76.

The BAM file is the set of the reads available extracted from SRR5398624 (available from [GEO](https://www.ncbi.nlm.nih.gov/sra/SRR5398624/)), trimmed and mapped with [bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (with option `-a --best --strata -p 6 -m 20`).

You can test mmannot with

    ./mmannot -a test_dataset.gtf -r test_dataset.bam -c configHS38.txt
    

### Compilation error

My compiler says `cc1plus: error: unrecognized command line option ‘-std=c++11’`

The code is written in C++11, like many other bioinformatics tools. C++11 offers many new functions to the developper, and it seems that your compiler is too old for C++11. Starting from 4.7, released in 2012, GCC supports the C++11 features present in mmannot. You could consider upgrading to this version at least.

## Contact

Comment? Suggestion? Do not hesitate sending me an email: <matthias.zytnicki@inra.fr>.
