# hiplexdeplex
Tool suite to demultiplex hiplex data

[generic description here]


 Written in Python 2.7, depends on Biopython

# Usage:

python detag.py [-h] --tags TAGS --heel HEEL [--pairedfastq]
                       [--forward FORWARD] [--reverse REVERSE]
                       [--unpairedfastq UNPAIREDFASTQ]
                       [--pairedparams PAIREDPARAMS] --outputprefix
                       OUTPUTPREFIX --h5score H5SCORE --h3score H3SCORE
                       --filtering FILTERING [--gsp GSP]

## optional arguments:
  -h, --help
                                      show this help message and exit

  --tags TAGS
                                      path/to/tagfile.txt.

  --heel HEEL
                                      path/to/heelfile.txt

  --pairedfastq
                                      input consists of --forward and --reverse

  --forward FORWARD
                                      forward reads of two paired fastq files

  --reverse REVERSE
                                      reverse reads of two paired fastq files

  --unpairedfastq UNPAIREDFASTQ
                                      input consisting of one unpaired fastq file
  --pairedparams PAIREDPARAMS
                                      overlap_kmer:overlap_hsp:overlap_min, default 5:7:10

  --outputprefix OUTPUTPREFIX         
                                      the prefix for the output. results in
                                      outputprefix_tagnameXX.fastq files.

  --h5score H5SCORE
                                      proportion of matching 5' heel_bases for acceptance.
                                      between 0 and 1. default 0.9

  --h3score H3SCORE
                                      proportion of matching 3' heel_bases for acceptance.
                                      between 0 and 1. default 0.9
  --filtering FILTERING
                                      parameters for filtering:
                                      status:min_len:max_len:mean_qual:min_qual default
                                      0:70:150:30:20
  --gsp GSP             
                                      file with the gene specific primers, for additional
                                      diagnostics. same format as tag file, default = None
