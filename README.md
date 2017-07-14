# hiplexdeplex
Tool suite to demultiplex hiplex data

## what is it?

it's a demultiplexer for sorting reads, such as from multiplexed samples from targeted-amplicon sequencing, based on their sequence tags.
It is capable of sorting out Chimeric reads and contains some base-quality and read-length filters, as well as some output statistics indicative of the success of the sequencing-process.
This should make it suitable for use as the only read-processing tool before aligning the reads to the reference genome.
Written in Python 2.7, depends on Biopython.

![what_does_it_do](./figures/detag_what_does_it_do.png)




## how does work?
The demultiplexer works by taking the read-structure into account, which is a result of the sequencing stragety employed to create the data.
Due to the two-step PCR, all reads are structured the same way:
They consist of the amplified region of interest, flanked on both sides with a primer-“heel” that is the same for all reads, and a sample-specific tag that is attached to all amplified regions of one sample.

![read-structure](./figures/Seq_strategy.png)

The demultiplexer exploits this by only aligning the sequences that are common to all reads, the heel-regions, to each read, inferring the tags of a read via their position relative to the heels, comparing the sequences found within the read to a dictionary of tag-sequences. This avoids aligning every single tag-sequences to the read, saving time in the process.


# Usage:

python detag.py [-h] --tags TAGS --heel HEEL [--pairedfastq]
                       [--forward FORWARD] [--reverse REVERSE]
                       [--unpairedfastq UNPAIREDFASTQ]
                       [--pairedparams PAIREDPARAMS] --outputprefix
                       OUTPUTPREFIX --h5score H5SCORE --h3score H3SCORE
                       --filtering FILTERING [--gsp GSP]

## arguments:
  -h, --help
                                      show this help message and exit

  --tags TAGS
                                      path/to/tagfile.txt.

  --heel HEEL
                                      path/to/heelfile.txt

  --pairedfastq
                                      input consists of paired fastq files, mutually exclusive with --unpairedfastq

  --forward FORWARD
                                      when --pairedfastq is specified, this flag is used to indicate the forward reads of the two paired fastq files.

  --reverse REVERSE
                                      when --pairedfastq is specified, this flag is used to indicate the reverse reads of the two paired fastq files.

  --unpairedfastq UNPAIREDFASTQ
                                      input consisting of one unpaired fastq file
  --pairedparams PAIREDPARAMS
                                      parameter for pairing the forward and reverse reads to make a consensus sequence.
                                      format is:
                                      overlap_kmer:overlap_hsp:overlap_min, default 5:7:10

                                      **overlap_kmer** specifies the kmer-size initially used to seed matches between the two paired reads.
                                      **overlap_hsp**  specifies the  minimum length of a High-Scoring Segment Pair (HSP) to be considered a match.
                                      **overlap_min**  specifies the minimum amount or minimum length of all HSP, in order for the read-pair to qualify for further processing.


  --outputprefix OUTPUTPREFIX         
                                      the prefix for the output. results in *outputprefix_tagnameXX.fastq* files.
                                      doesnt create folders. eg for *project/outputprefix_tagnameXX.fastq*, the folder project needs to be created beforehand.

  --h5score H5SCORE
                                      proportion of matching 5' heel_bases for acceptance.
                                      between 0 and 1. default 0.9

  --h3score H3SCORE
                                      proportion of matching 3' heel_bases for acceptance.
                                      between 0 and 1. default 0.9
  --filtering FILTERING
                                      parameters for filtering.
                                      **format:**
                                      strategy:min_len:max_len:mean_qual:min_qual
                                      **default:**
                                      0:70:150:30:20
                                      * **strategy** has a value of either  0, 1, 2 or 3.
                                        0  indicates full filtering, before the program searches for tags in the sequences.
                                        1 indicates no filtering of reads.
                                        2 indicates filtering down to high-quality regions
                                        3 indicates full filtering, but after the reads have been searched for the heel and tags.
                                      * **min_len** is the minimum length for a read to be not filtered out.
                                      * **max_len** is the maximum length for a read to be not filtered out.
                                      * **mean_qual** is the minimum average  per-base phred-score that a read needs to have in order be retained.
                                      * **min_qual** is the minimum phred-score a base needs to have in order to not be trimmed away in HQR-filtering.


  --gsp GSP             
                                      file with the gene specific primers, for additional
                                      diagnostics. same format as tag file, default = None
