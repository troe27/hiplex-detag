# import modules
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from array import array
from Bio import SeqIO
import gzip
from collections import Counter

# commandlineparser


def cli_parser():

    '''
    - parses the command line interface.
    - outputs the specified parameters and filenames into a job-dictionary.
    '''
    parser_main = argparse.ArgumentParser(prog='python detag.py')
    parser_main.add_argument("--tags",
                             help="path/to/tagfile.txt.",
                             required=True)
    parser_main.add_argument("--heel",
                             help="path/to/heelfile.txt",
                             required=True)
    input_type = parser_main.add_mutually_exclusive_group(required=True)

    input_type.add_argument("--pairedfastq",
                            help="input consists of --forward and --reverse ",
                            action="store_true")
    parser_main.add_argument("--forward",
                             help="forward reads of two paired fastq files")
    parser_main.add_argument("--reverse",
                             help="reverse reads of two paired fastq files")
    input_type.add_argument("--unpairedfastq",
                            help="input consists of one  unpaired fastq file")
    parser_main.add_argument("--pairedparams",
                             help="overlap_kmer:overlap_hsp:overlap_min, default 5:7:10",
                             default="5:7:10")  # check sensible defaults!
    parser_main.add_argument("--outputprefix",
                             help="the prefix for the output. results in outputprefix_tagnameXX.fastq files.",
                             required=True)
    parser_main.add_argument("--h5score",
                             help="proportion of matching  5' heel_bases for acceptance. between 0 and 1. default 0.9",
                             required=True,
                             default=0.9)
    parser_main.add_argument("--h3score",
                             help="proportion of matching 3' heel_bases for acceptance. between 0 and 1. default 0.9",
                             required=True,
                             default=0.9)
    parser_main.add_argument("--filtering",
                             help="parameters for filtering: status:min_len:max_len:mean_qual:min_qual default 0:70:150:30:20",
                             required=True,
                             default="0:70:150:30:20")  # need good defaults
    parser_main.add_argument("--gsp",
                             help="file with the gene specific primers, for additional diagnostics. same format as tag file, default = None",
                             default=None)

    args = parser_main.parse_args()
    job_summary = ["starting Job with the following parameters:"]
    job = dict()
    job["heel_file"] = args.heel
    job["gsp"] = args.gsp
    job_summary.append("heel file is: "+str(job['heel_file']))
    job["tag_file"] = args.tags
    job_summary.append("tag file is: "+str(job['tag_file']))
    if args.h3score and args.h5score is None:
        job_summary.append("using default values for primerscores")
        job["h3score"] = float(0.9)
        job["h5score"] = float(0.9)
    else:
        job["h3score"] = float(args.h3score)
        job["h5score"] = float(args.h5score)

    job["output_prefix"] = args.outputprefix

    job_summary.append("primerscore (match/length ratio)  3'>{}  5'>{}".format(job["h3score"], job["h5score"]))
    if args.pairedfastq is True:
        job["filetype"] = "pairedfastq"
        job["forward_file"] = args.forward
        job["reverse_file"] = args.reverse
    if args.pairedfastq is not True:
        job["filetype"] = "unpairedfastq"
    job["unpaired_file"] = args.unpairedfastq
    if (args.unpairedfastq and args.forward) or (args.unpairedfastq and args.reverse )is not None:
        print " if the reads are not paired, you should not specify --forward and --reverse"
        exit()

    job_summary.append("output_prefix is {}".format(job["output_prefix"]))
    job_summary.append("file type is : {} ".format(job["filetype"]))
    job_summary.append("input files are : {} , {}, {}".format(job['forward_file'],
                                                              job['reverse_file'],
                                                              job['unpaired_file']))

    if args.pairedparams is not None and args.pairedfastq is True:
        param_list = args.pairedparams.split(":")
        job["kmer_overlap"] = int(param_list[0])
        job["hsp_overlap"] = int(param_list[1])
        job["min_overlap"] = int(param_list[2])
    elif args.pairedparams is None and args.pairedfastq is True:
        job_summary.append("no pairedparams found, using default parameters.")
        job["kmer_overlap"] = 7
        job["hsp_overlap"] = 5
        job["min_overlap"] = 10

    elif args.pairedparams is not None and args.pairedfastq is not True:
        job_summary.append("paired - parameters only relevant for paired fastq files, ignoring parameters.")
    if args.pairedfastq is True:
        job_summary.append('using paired-parameters kmer_overlap = {} , hsp_overlap = {} , min_overlap = {}'.format(job["kmer_overlap"],
                                                                                                                    job["hsp_overlap"],
                                                                                                                    job["min_overlap"]))
    filter_list = args.filtering.split(":")
    job["filter_strategy"] = int(filter_list[0])
    job["filter_minlen"] = int(filter_list[1])
    job["filter_maxlen"] = int(filter_list[2])
    job["filter_meanqual"] = int(filter_list[3])
    job["filter_minqual"] = int(filter_list[4])

    job_summary.append("filtering strategy {} , with minlen {} , maxlen {}, meanqual {} , minqual {}".format(job["filter_strategy"],
                                                                                                             job["filter_minlen"],
                                                                                                             job["filter_maxlen"],
                                                                                                             job["filter_meanqual"],
                                                                                                             job["filter_minqual"]))
    return job, job_summary

# file parser


def get_tagset(tag_input_file, end):
    '''
    tagfile in tagname;5'tag;3'tag\n order.
    End=5 makes a dict with the 5'tags as keys and matching 3'tags & tagname
    End=3 makes a dict with the 3'tags as keys and matching 5'tags & tagname
    '''
    tags = dict()
    tag_file = open(str(tag_input_file), "r")
    tag_length = 0
    bad_line_counter = 0
    ok = dict.fromkeys("GATC")
    min_len = None
    # Parse tag-file and perform sanity checks
    # checking for uniform length of tags:
    len_list = []
    tag_summary = ""
    for line in tag_file:
        line = line.rstrip()
        if line != "" and not line.startswith("#") and all(c in ok for c in line.split(";")[1] and line.split(";")[2]):
            len_list.append(len(line.split(";")[1]))
            len_list.append(len(line.split(";")[2]))
    if not all(length == len_list[0] for length in len_list):
        min_len = min(len_list)
        tag_summary += "difference in tag/primer length detected: truncating to {} bases \n".format(min_len)
    if end == "5":
        tag_summary += "tagname\t3'\t5'\n"
    else:
        tag_summary += "tagname\t5'\t3'\n"

    tag_file = open(str(tag_input_file), "r")
    # tag_file.seek(0)
    for line in tag_file:
        # print #line still good here
        line_list_a = list()
        line = line.rstrip()
        # print #line still good here
        if line == "":
            # print line nothing here
            continue  # skip empty lines
        if line.startswith("#"):
            # print line
            continue  # skip comments
        if all(c not in ok for c in line.split(";")[1] and line.split(";")[2]):
            # print line, "dummy"
            bad_line_counter += 1
            bad_line = line
            continue  # skip lines that are not in the right format
        # print line
        line_list_a = line.split(";")
            # now assumes definitive tag_combinations. e.g. tagname;5'tag;3'tag\n,
        if end == "3":  # if 3' , shuffle tags in order to enable gain of 3' and 5' dicts from the same file.
            order = [0, 2, 1]
            line_list = [line_list_a[i] for i in order]
        elif end == "5":
            line_list = line_list_a
        else:
            raise IOerror

        # print line_list
        if len(line_list) < 2:
            continue

        if min_len:

            if end == "5":
                tags[line_list[1][0:min_len]] = dict(tag=line_list[0],
                                                     matching=set())
                if tag_length == 0:
                    tag_length = int(min_len)

                if len(line_list) > 2:
                    for i in range(2, len(line_list)):
                        if line_list[i] != "":
                            rev_trunc_point = len(line_list[i]) - min_len
                            tags[line_list[1][0:min_len]]['matching'] |= set([line_list[i][rev_trunc_point:]])
            if end == "3":
                trunc_point = (len(line_list[1]) - min_len)
                tags[line_list[1][trunc_point:]] = dict(tag=line_list[0],
                                                        matching = set())
                if tag_length == 0:
                    tag_length = int(min_len)

                if len(line_list) > 2:
                    for i in range(2, len(line_list)):
                        if line_list[i] != "":
                            tags[line_list[1][trunc_point:]]['matching'] |= set([line_list[i][0:min_len]]) #:min_len
                            #print line_list[i]

        else:
            tags[line_list[1]] = dict(tag=line_list[0],
                                      matching = set() )
            if tag_length == 0:
                tag_length = len(line_list[1])

            if len(line_list) > 2:
                for i in range(2,len(line_list)):
                    if line_list[i] != "":
                        tags[line_list[1]]['matching'] |= set([line_list[i]])

        tags["_____length"] = tag_length

    for key in tags.keys() :
        if key is not "_____length":
            tagsum = str(tags[key]["tag"])+"\t"+str(key)+" "+str(" ".join(map(str,[x for x in list(tags[key]["matching"])]))+"\n")
            tag_summary+=tagsum
    if bad_line_counter > 3:
        tag_summary += "more than three lines in your tag-file dont conform to the required format, such as : \n"
        tag_summary += str(bad_line)
    #print tag_summary
    return tags, tag_summary

def get_heels(heel):
    heel5 = None
    heel3 = None
    if ";" in str(heel):
        heel5=heel.split(";")[0] #change here
        heel3=heel.split(";")[1] #change here
    elif ".txt" in str(heel):
        heel_file = open(str(heel), "r")
        counter=0
        for line in heel_file:
            line_list_a = list()
            line = line.rstrip()

            if line == "":
                continue
            if line.startswith("#"):
                continue
            if ";" not in line:
                continue
            line_list_a = line.split(";")
            counter+=1
            heel5 = line_list_a[0]
            heel3 = line_list_a[1]
            if counter >= 2:
                heel_summary= "please provide a valid heel file, number of valid lines too high"
            elif not heel3:
                heel_summary= "please provide a valid heel file, three_prime heel could not be found"
            elif  not heel5:
                heel_summary= "please provide a valid heel file, five_prime heel could not be found"
            else:
                heel_summary = "heels : 5' -> {} 3' -> {} ".format(heel5, heel3)
    return heel5, heel3, heel_summary

#formatting sequence_files

class Pair:

    def __init__(self, fastq_1, fastq_2,
                 kmer=7, hsp=5, min=10):
        self.kmer=kmer
        self.hsp = hsp
        self.min = min
        try:
            handle1 = gzip.open(fastq_1, "r")
            handle1.next() # This fails if not a gzip file, goes into the exception
            handle1 = gzip.open(fastq_1, "r")
            self.fastq1 = SeqIO.parse(handle1, "fastq")
            self.fastq2 = SeqIO.parse(gzip.open(fastq_2, "r"), "fastq")
        except IOError:
            self.fastq1 = SeqIO.parse(fastq_1, "fastq")
            self.fastq2 = SeqIO.parse(fastq_2, "fastq")
        self.qual_present = True


    def __iter__ (self):
        return self

    def next(self):
        s1 = self.fastq1.next()
        s2 = self.fastq2.next()

        return FastQPairQualSeq(s1, s2, self.kmer, self.hsp, self.min)

class QualSeq:
    def __init__(self, data):
        self.d = data

    def get(self):
        return self.d


class Qual:
    def __init__ (self, name, qual):
        self.name = name
        self.quals = qual


    def __getitem__ (self, item):
        return int(self.quals[item])

    def __repr__ (self):
        return "Qual(name = " + self.name + ", quals = [" + (", ".join([str(x) for x in self.quals])) + "]"

class FastQPairQualSeq(QualSeq):
    def __init__(self, s1, s2, kmer, hsp, min):
        self.s1=s1
        self.s2=s2
        self.kmer = kmer
        self.hsp = hsp
        self.min = min
        self.s = None
        self.quals = None

    def get(self):
        if not self.s:
            self._pair()
        return [self.s, self.quals]

    def _pair(self):
        s1 = str(self.s1.seq)
        s2_rec = self.s2.reverse_complement()
        s2 = str(s2_rec.seq)

        # Build kmer table of read 2
        kmers2 = dict()
        for i in range(len(s2) - self.kmer):
            k = s2[i:i + self.kmer]
            kmers2.update({k : (kmers2.get(k, []) + [i])})

        # Identify read1 kmers in read2
        kmer_pos = [[x[0], kmers2.get(x[1]), x[2]] for x in
                     [(s1[y], s1[y : y + self.kmer], y) for y in range(len(s1) - self.kmer)]]

        # Identify the longest kmer run
        runs = [ ]
        r = [ ]
        for i in range(len(kmer_pos)):
            if kmer_pos[i][1]:
                r.append(kmer_pos[i])
            else:
                if len(r):
                    if len(r) > self.hsp:
                        runs.append(r)
                r=[]
        runs.sort(key=lambda a: len(a), reverse=True)


        if not len(runs) or sum([len(r) for r in runs]) < self.min:
            return None
        run = runs[0]

        # Remove ambigous start/end
        while len(run):
            if len(run[0][1]) > 1:
                run = run[1:]
            else:
                break
        while len(run):
            if len(run[-1][1]) > 1:
                run = run[:-1]
            else:
                break
        if not len(run):
            return None
        # Splice together new sequence
        # TODO: This could be improved to keep the highest quality base
        # from each read.

        try:
            new_seq = s1[:run[-1][2]] + s2[run[-1][1][0]:]
            quals = self.s1.letter_annotations["phred_quality"][:run[-1][2]] + \
                s2_rec.letter_annotations["phred_quality"][run[-1][1][0]:]
        except IndexError:
            print "IndexError"
            print run, runs
            return None
        self.s=self.s2
        self.s.letter_annotations = dict()
        self.s.seq = Seq(new_seq, generic_dna)
        self.quals = Qual(self.s.id, quals)


# table for converting bases into integers suitable for bitwise &
A = 1
C = 2
G = 4
T = 8
trans_table = { 'A' : A,
                'C' : C,
                'G' : G,
                'T' : T,
                'R' : A | G,
                'Y' : C | T,
                'S' : G | C,
                'W' : A | T,
                'K' : G | T,
                'M' : A | C,
                'B' : C | G | T,
                'D' : A | G | T,
                'H' : A | C | T,
                'V' : A | C | G,
                'N' : A | C | T | G }

#detagging classes and functions

class DeTagSeq:
    # Translation of bases for primer identification.
    def __init__(self, p5, p5s, p3, p3s, t5, t3, gsp5=None, gsp3=None):
        self.p5s = p5s
        self.p3s = p3s
        self.t5 = t5
        self.t3 = t3
        #init gene specific primer part:
        self.gsp5 = gsp5
        self.gsp3 = gsp3
        #self.min_len = min_len
        # Translate primer sequence
        self.p3_rev = str(Seq(p3,generic_dna).reverse_complement()) ##change
        self.p5 = [trans_table[x] for x in p5.upper() if x in trans_table]
        self.p3 = [trans_table[x] for x in self.p3_rev.upper() if x in trans_table] ##change
        self.p5len = float(len(self.p5))
        self.p3len = float(len(self.p3))



    # Function to detag sequence and return a dict()
    # with sequence and metadata
    def detag_seq(self, seq_record, q=None):
        result = dict(rev=False)
        rev = False
        seq = seq_record.seq
        seq_str = str(seq)
        seq_list = [trans_table[x] if not x == 'N' else 0 for x in str(seq).upper() if x in trans_table]
        p5_pos = -1
        accepted_t3 = set()
        accepted_gsp3 = set()
        if len(self.p5):
            for x in range(0, min(len(seq_str), 2000 + len(self.p5) ) - len(self.p5)):
                m = 0
                for i in range(len(self.p5)):
                    if seq_list[x + i] & self.p5[i]:
                        m += 1
                if m / self.p5len > self.p5s:
                    p5_pos = x
                    break

            if p5_pos < 0:
                result["rev"] = True
                seq = seq.reverse_complement()
                seq_str = str(seq)
                seq_list = [trans_table[x] if not x == 'N' else 0 for x in str(seq).upper() if x in trans_table]

                for x in range(0, min(len(seq_str), 2000 + len(self.p5) ) - len(self.p5)):
                    m = 0
                    for i in range(len(self.p5)):
                        if seq_list[x + i] & self.p5[i]:
                            m += 1
                    if m / self.p5len > self.p5s:
                        p5_pos = x
                        break
            if p5_pos < 0:
                result["status"] = "no_primer5"
                if q:
                    return result, None # new change here
                else:
                    return result
        else:
            p5_pos=0
        #identify tag
        if len(self.t5):
            tag_len = self.t5["_____length"]
            tag_seq = seq_str[(p5_pos - tag_len):p5_pos]
            try:
                result["tag_name"] = self.t5[tag_seq]['tag']
                accepted_t3 = self.t5[tag_seq]['matching']
            except KeyError:
                result["status"] = "no_tag5"
                if q:
                    return result, None # new change here
                else:
                    return result
        else:
            result["tag_name"] = ""

        if self.gsp5:
            gsprim_len = self.gsp5["_____length"]
            gsprim_seq = seq_str[(p5_pos + len(self.p5)):(p5_pos + len(self.p5) + gsprim_len) ]
            try:
                result["gsprim_name"] = self.gsp5[gsprim_seq]['tag']
                accepted_gsp3 = self.gsp5[gsprim_seq]['matching']
            except KeyError:
                result["status_gsp"] = "no_gsprim5"
                #if q:
                #    return result, None
                #else:
                 #   return result
        else:
            result["gsprim_name"] = ""

        p3_pos = -1
        p3_matches = [ ]

        if len(self.p3):
            for x in range(len(seq) - len(self.p3) - 1, p5_pos, -1):
                m=0
                for i in range(len(self.p3)):
                    if seq_list[x + i] & self.p3[i]:
                        m += 1
                if m / self.p3len > self.p3s:
                    p3_pos = x
                    break
            if p3_pos < 0:
                result["status"] = "no_primer3"
                if q:
                    return result, None  # new change here
                else:
                    return result

            if len(self.t3):
                tag_len = self.t3["_____length"]
                tag_seq = seq_str[(p3_pos + len(self.p3)):(p3_pos + len(self.p3) + tag_len)]
                tag_seq = str(Seq(tag_seq,generic_dna).reverse_complement())
                try:
                    if len(accepted_t3) and tag_seq not in accepted_t3:     # changed here
                        result['status']='chimeric_tag'
                        result["tag_name"] += ";" + self.t3[tag_seq]['tag']  # added this here NUUCHANGE
                        return result, q  # new change here ## NUCHANGER
                    # result["tag_name"] += ("_" + self.t3[tag_seq]['tag']) # NUUCHANGE
                except KeyError:
                    result["status"] = "no_tag3"
                    if q:
                        return result, None  # new change here
            else:
                result["tag_name"] += ""

            if self.gsp3:
                gsprim_len = self.gsp3["_____length"]
                gsprim_seq = seq_str[(p3_pos - gsprim_len):p3_pos]
                gsprim_seq = str(Seq(gsprim_seq,generic_dna).reverse_complement())
                try:
                    if len(accepted_gsp3) and gsprim_seq not in accepted_gsp3:
                        result['status_gsp']='chimeric_gsprim'
                        result["gsprim_name"] += ";" + self.gsp3[gsprim_seq]['tag']
                        #return result, q

                except KeyError:
                    result["status_gsp"] = "no_gsprim3"
            else:
                result["gsprim_name"] += ""

            result["seq"] = seq_str[(p5_pos + len(self.p5)):p3_pos]
            seq_record.seq = seq[(p5_pos + len(self.p5)):p3_pos]
            result["seq_record"] = seq_record
        else:
            result["seq"] = seq_str[(p5_pos + len(self.p5)):]
            seq_record.seq = seq[(p5_pos + len(self.p5)):]
            result["seq_record"] = seq_record

        if not q:
            return result
        else:
            if p3_pos:
                q.quals = q.quals[(p5_pos + len(self.p5)):p3_pos]
            else:
                q.quals = q.quals[(p5_pos + len(self.p5)):]
            return result, q


#filter_functions

def filter_hqr(seq_record, qual, min_length, mean_min, min_qual):
    ''''''
    if len(seq_record.seq) < min_length:
        return [None, "too_short"]

    # Get regions with quality above minimum

    hqr = [ ]

    t_start=-1
    t_end=-1

    for i in range(len(qual.quals)):
        if t_start == -1: # Outside region
            if qual.quals[i] >= min_qual:
                t_start = i
        else: # Inside region
            if qual.quals[i] < min_qual:
                t_end = i
                hqr.append(dict(start=t_start,
                                end=t_end,
                                mean_quality=0,
                                length=0))
                t_start = -1
                t_end = -1
    if t_start != -1 and t_end == -1:
        hqr.append(dict(start=t_start,
                        end=len(qual.quals) - 1,
                        mean_quality=0,
                        length=0))
    # Check and modify regions to ensure mean quality is above limit
    # print "Before:", hqr
    for r in hqr:
        if r["end"] - r["start"] < min_length:
            continue
        while r["end"] > r["start"] and \
                (sum(qual.quals[r["start"]:r["end"]]) / float(r["end"] - r["start"])) \
                < mean_min:
            if qual.quals[r["start"]] < qual.quals[r["end"]]:
                r["start"] += 1
            else:
                r["end"] -= 1
        if r["end"] > r["start"]:
            r["mean_qual"] = (sum(qual.quals[r["start"]:r["end"]]) / float(r["end"] - r["start"]))
            r["length"] = r["end"] - r["start"]
    # print "After:", hqr

    hqr = filter(lambda a: a["length"] >= min_length and a["mean_qual"] >= mean_min, hqr)

    if not len(hqr):
        return [None, "low_mean_quality"]

    hqr.sort(key=lambda a: a["length"], reverse=True)

    if hqr[0]["end"] - hqr[0]["start"] < min_length:
        return [None, "low_mean_quality"]

    seq_record.seq = seq_record.seq[hqr[0]["start"]:hqr[0]["end"]]
    return [seq_record, None]


def filter_full(seq_record, qual, min_length, mean_min, min_qual):
    if len(seq_record.seq) < min_length:
        return [None, "too_short"]
    if float(sum(qual.quals)) / len(qual.quals) < mean_min:
        return [None, "low_mean_quality"]
    seq_list = list(seq_record.seq)
    seq_list = map(lambda b, q: b if q >= min_qual else 'N',
                   seq_list, qual.quals)
    if seq_list.count('N'):
        return [None, "low_min_quality"]
    seq_record.seq = Seq("".join(seq_list), seq_record.seq.alphabet)
    return [seq_record, None]


# execute and stats

class ResultProcessor(object):

    def __init__(self, tags, output_prefix=None, gsp=None):
        if gsp:
            self.gsp = [gsp[key]["tag"] for key in gsp.keys() if key is not "_____length"]
        else:
            self.gsp=None
        self.prefix = output_prefix
        self.total_reads = 0
        self.discarded_reads = 0
        self.kept_reads = 0
        self.tags = [tags[key]["tag"] for key in tags.keys() if key is not "_____length"]
        if self.gsp:
             self.stats_dict = {tag:dict(valid=0, gsps=Counter({gsp:0 for gsp in self.gsp}), failed=0, total=0, chimera_count=0) for tag in self.tags}
        else:
            self.stats_dict = {tag:dict(valid=0, failed = 0, total=0, chimera_count=0) for tag in self.tags}
        self.gsp_count = Counter()
        self.gsp_disc = Counter()
        self.status = {tag:Counter() for tag in self.tags}
        self.status["total"] = Counter()
        self.chimeras = Counter()
        self.chimeric_tags = Counter()
        #formatted stats:
        self.h_stats = list() #basic stats
        self.g_stats = list() #per tag per gsp  stats
        self.g_count = list()
    def open_tag_files(self):
        self.files = {a:open(str(self.prefix)+"_"+str(a)+".fq","w") for a in self.tags}
        for key in self.files.keys():
            self.files[key]

    def close_tag_files(self):
        for key in self.files.keys():
            self.files[key].close()
        #for tag in self.tags:
            #tag.close()
            #tag = close(str(self.prefix)+"_"+str(tag)+".fq")

    #def write_result(self, result, outfile):
        #rec = result["detagged_seq"]["seq_record"]
        #rec.letter_annotations["phred_quality"] = result["detag_qual"].quals
        #with open(outfile, "a") as handle:
            #SeqIO.write(rec, handle, "fastq")
        #SeqIO.write(rec, outfile, "fastq")


    def process_result(self, result):
        self.total_reads += 1

        if "detagged_seq" in result and "tag_name" in result["detagged_seq"]:
            if "chimeric_tag" in result["status"].keys():
                                                    self.stats_dict[result["detagged_seq"]["tag_name"].split(";")[0]]["total"] += 1
        # chimeric tags are counted for 5' tag in total
            else:
                self.stats_dict[result["detagged_seq"]["tag_name"]]["total"] += 1 # total reads per tag

        if "status" in result:
            self.status["total"].update(result["status"].keys()) #collect total status
                #print result["status"].keys()[0]
            if "detagged_seq" in result and "tag_name" in result["detagged_seq"]:
                if "chimeric_tag" in result["status"].keys():
                    self.status[result["detagged_seq"]["tag_name"].split(";")[0]].update(result["status"].keys())
                else:
                    self.status[result["detagged_seq"]["tag_name"]].update(result["status"].keys())


        if "keep" in result:
            rec = result["detagged_seq"]["seq_record"]
            rec.letter_annotations["phred_quality"] = result["detag_qual"].quals
            SeqIO.write(rec, self.files[result["detagged_seq"]["tag_name"]], "fastq")

            #self.write_result(result, str(self.prefix)+"_"+str(result["detagged_seq"]["tag_name"])+".fq" )
            #self.write_result(result, result["detagged_seq"]["tag_name"])
            self.kept_reads += 1
            self.stats_dict[result["detagged_seq"]["tag_name"]]["valid"] += 1
            if self.gsp:
                if "status_gsp" in result["detagged_seq"]:
                    self.gsp_count.update([result["detagged_seq"]["status_gsp"]])
                if "gsprim_name" in result["detagged_seq"] :
                    if "status_gsp" not in result["detagged_seq"]:
                        self.stats_dict[result["detagged_seq"]["tag_name"]]["gsps"].update([result["detagged_seq"]["gsprim_name"]])
                        self.gsp_count.update([result["detagged_seq"]["gsprim_name"]])
                    if "status_gsp" in result["detagged_seq"]:
                        if result["detagged_seq"]["status_gsp"] != "chimeric_gsprim":
                            self.stats_dict[result["detagged_seq"]["tag_name"]]["gsps"].update([result["detagged_seq"]["gsprim_name"]])
                            self.gsp_count.update([result["detagged_seq"]["gsprim_name"]])

        else:
            self.discarded_reads += 1
            # gsp
            if self.gsp:
                if self.gsp and "detagged_seq" in result:
                    if "status_gsp" in result["detagged_seq"]:
                        self.gsp_disc.update([result["detagged_seq"]["status_gsp"]])
            #chimeras
            if "status" in result:
                if "chimeric_tag" in  result["status"].keys():
                    self.chimeras.update([result["detagged_seq"]["tag_name"]])
                    self.chimeric_tags.update([result["detagged_seq"]["tag_name"].split(";")[0], result["detagged_seq"]["tag_name"].split(";")[1] ])

            #if result["detagged_seq"]["tag_name"] and not "chimeric_tag" in  result["status"].keys(): #wtf did you want to say here, tilman? this is why you shouldnt write code in the evenings.
                #self.stats_dict[result["detagged_seq"]["tag_name"]]["failed"] += 1

    def harvest_stats(self):
        self.h_stats.extend([("processed_reads_in_total:", self.total_reads)])
        self.h_stats.extend([("tag", "valid_reads_per_tag", "percentage_of_total_reads")])
        self.h_stats.extend([("kept:", self.kept_reads, (float(self.kept_reads)/float(self.total_reads))*float(100))])
        self.h_stats.extend([("discarded:", self.discarded_reads, (float(self.discarded_reads)/float(self.total_reads))*100)])
        self.h_stats.extend([("\t", "\t")])
        self.h_stats.extend([("tag", "valid_reads_per_tag")])
        self.h_stats.extend([(key,self.stats_dict[key]["valid"])for key in self.stats_dict.keys()])
        self.h_stats.extend([("\t", "\t")])
        self.h_stats.extend([("cause_for_failure", "amount_of_reads","percentage_of_total_reads" )])
        self.h_stats.extend([(key,self.status["total"][key],(float(self.status["total"][key])/float(self.total_reads))*float(100) ) for key in self.status["total"].keys()])
        self.h_stats.extend([("\t", "\t")])
        #self.h_stats.extend([("tag","most_common_failure_per_tag", "amount_of_reads")])
        #self.h_stats.extend([[key, self.status[key].most_common(1)[0][0],self.status[key].most_common(1)[0][1]] for key in self.status.keys() if key is not "total"])
        self.h_stats.extend([("\t", "\t")])
        self.h_stats.extend([("most_common chimeras", "amount_of_reads")])
        self.h_stats.extend([chimera for chimera in self.chimeras.most_common(10)])
        self.h_stats.extend([("\t", "\t")])
        self.h_stats.extend([("chimeric_tags", "amount_of_reads")])
        self.h_stats.extend([(key, self.chimeric_tags[key])for key in self.chimeric_tags.keys()])
        ##gsp_tab
        if self.gsp:
            #self.g_stats.append("primers")
            self.g_stats.append("tags;"+";".join(map(str,self.stats_dict[self.stats_dict.keys()[0]]["gsps"].keys())))
            for key in self.stats_dict.keys():
                self.g_stats.append(str(key)+";"+";".join(map(str,[self.stats_dict[key]["gsps"][primer] for primer in self.stats_dict[key]["gsps"].keys()])))
            self.g_count.append("gsp;read_count")
            for key in self.gsp_count.keys():
                self.g_count.append(str(key)+";"+str(self.gsp_count[key]))
            self.g_count.append(str(key))
        #return self.h_stats
    def write_stats(self):
        if self.gsp:
            with open(self.prefix+"_reads_per_gsp_per_tag.txt", "w") as handle:
                for l in self.g_stats:
                    handle.write(l+"\n")
            with open(self.prefix+"_stats.txt", "w") as handle:
                for line in self.h_stats:
                    handle.write("\t".join(map(str, line))+"\n")
            with open(self.prefix+"_gsp_counts.txt", "w") as handle:
                for line in self.g_count:
                    handle.write(line+"\n")
        else:
            with open(self.prefix+"_stats.txt", "w") as handle:
                for line in self.h_stats:
                    handle.write("\t".join(map(str, line))+"\n")

###need to fix write_logfile
def write_logfile(out_prefix,jobsum, heelsum, tagsum):
    with open(out_prefix+"_logfile.txt" ,"w") as handle:
        for i in jobsum:
            handle.write(i+"\n")
        handle.write(heelsum+"\n")
        handle.write(tagsum+"\n")

def main():
    job,jobsum = cli_parser()
    tags5,tagsum = get_tagset(tag_input_file=job["tag_file"],end='3')
    tags3 = get_tagset(tag_input_file=job["tag_file"],end='5')[0]
    if job["gsp"]:
    	gsp3,gspsum = get_tagset(tag_input_file=job["gsp"],end='3')
    	gsp5 = get_tagset(tag_input_file=job["gsp"],end='5')[0]

    heel5, heel3, heelsum = get_heels(job["heel_file"])
    seq_data = Pair(fastq_1=job["forward_file"],
                    fastq_2=job["reverse_file"],
                    hsp=job["hsp_overlap"],
                    kmer=job["kmer_overlap"],
                    min=job["min_overlap"])
    # init DeTagSeq
    if not job["gsp"]:
        dtseq = DeTagSeq(p3=heel3, p3s=job["h3score"], p5=heel5, p5s=job["h5score"], t3=tags3, t5=tags5)
    else:
        dtseq = DeTagSeq(p3=heel3, p3s=job["h3score"], p5=heel5, p5s=job["h5score"], t3=tags3, t5=tags5, gsp3=gsp3, gsp5=gsp5)

    # init ResultProcessor
    if not job["gsp"]:
        reproc = ResultProcessor(tags=tags3, output_prefix=job["output_prefix"])
    else:
        reproc = ResultProcessor(tags=tags3, output_prefix=job["output_prefix"], gsp=gsp3)

    # execute detag & reproc
    reproc.open_tag_files()

    for qseq in seq_data:
        r=qseq.get()
        res = dict(status = dict())
        seq = r[0]
        qual = r[1]
        detagged_seq = None
        if seq:
            res["id"] = seq.id
            if job["filter_strategy"] == 0: # Filter full ###
                if not qual:
                    raise Exception("Full sequence filtering requires quality data")

                seq, status = filter_full(seq_record=seq,
                                          qual=qual,
                                          min_length=job["filter_minlen"],
                                          mean_min=job["filter_meanqual"],
                                          min_qual=job["filter_minqual"])
            elif job["filter_strategy"] == 1: # Full sequence, no filtering###
                status = None
            elif job["filter_strategy"] == 2: # HQR###
                if not qual:
                    raise Exception("HQR filtering requires quality data")
                seq, status = filter_hqr(seq_record=seq,
                                         qual=qual,
                                         min_length=job["filter_minlen"],
                                         mean_min=job["filter_meanqual"],
                                         min_qual=job["filter_minqual"])
            elif job["filter_strategy"] == 3: # Primers first
                if not qual:
                    raise Exception("Amplicon quality filtering requires quality data")
                detagged_seq, qual = dtseq.detag_seq(seq,qual)
                if "status" in detagged_seq:
                    res["status"][detagged_seq["status"]] = 1
                    res["detagged_seq"] = detagged_seq
                    #result_list.append(res)
                    reproc.process_result(res) ###
                    continue

                seq, status = filter_full(seq_record=detagged_seq["seq_record"],
                                          qual=qual,
                                          min_length=job["filter_minlen"],
                                          mean_min=job["filter_meanqual"],
                                          min_qual=job["filter_minqual"])

            else:
                raise Exception("Unknown filter type")
        else:
            status = "failed_pair"
        if status:
            res["status"][status] = 1
            #result_list.append(res)
            reproc.process_result(res) ###
            continue
        if len(seq.seq) < job["filter_minlen"]:
            res["status"]["too_short"] = 1
            #result_list.append(res)
            reproc.process_result(res) ###
            continue

        if not detagged_seq:
            detagged_seq, qual = dtseq.detag_seq(seq, qual)
#	print detagged_seq

            res["detagged_seq"] = detagged_seq
        if "status" in detagged_seq:
#            print detagged_seq["status"]
            res["status"][detagged_seq["status"]] = 1
            #result_list.append(res)
            res["detagged_seq"] = detagged_seq
            reproc.process_result(res) ###
            continue
        if len(detagged_seq["seq"]) < job["filter_minlen"]:
            res["status"]["too_short"] = 1
            res["detagged_seq"] = detagged_seq
            #result_list.append(res)
            reproc.process_result(res) ###
            continue

        if job["filter_maxlen"] > 0 and len(detagged_seq["seq"]) > job["filter_maxlen"]:
            res["status"]["truncated"] = 1
            detagged_seq["seq"] = detagged_seq["seq"][:job["filter_maxlen"]]
            res["detagged_seq"] = detagged_seq
 #           continue

        if job["filter_maxlen"] < 0:
            res["status"]["truncated"] = 1
            detagged_seq["seq"] = detagged_seq["seq"][:job["filter_maxlen"]]
            res["detagged_seq"] = detagged_seq
#            continue

#        print res
        res["keep"] = True
        res["id"] = seq.id
        res["detag_qual"] = qual
        reproc.process_result(res)

    reproc.close_tag_files()
    reproc.harvest_stats()
    reproc.write_stats()
    write_logfile(job["output_prefix"],jobsum, heelsum, tagsum,)

#call main
main()
