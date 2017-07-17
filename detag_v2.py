#import modules
import argparse
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from array import array
from Bio import SeqIO
import gzip
from collections import Counter
import time
import datetime

#commandlineparser
def cli_parser():
    '''                                                                                            
    - parses the command line interface.                                                           
    - outputs the specified parameters and filenames into a job-dictionary.                        
    '''
    start_time_machine = time.time()
    start_time_human = datetime.datetime.utcnow()
    parser_main = argparse.ArgumentParser(prog='python detag.py') ### make sure the name is adequate when making the final script                      
    parser_main.add_argument("--tags",
                             help="path/to/tagfile.txt.",
                             required = True)
    parser_main.add_argument("--heel",
                             help="path/to/heelfile.txt",
                             required = True)
    input_type = parser_main.add_mutually_exclusive_group(required = True) # one of the input types need to be specified                                
                                                                                                                                                       
    input_type.add_argument("--pairedfastq",
                            help="input consists of --forward and --reverse ",
                            action="store_true")
    parser_main.add_argument("--forward",
                             help="the forward reads of two paired fastq files")
    parser_main.add_argument("--reverse",
                             help=" the reverse reads of two paired fastq files")
    input_type.add_argument("--unpairedfastq",
                            help="input consists of one  unpaired fastq file")
    parser_main.add_argument("--pairedparams",
                             help="overlap_kmer:overlap_hsp:overlap_min, default 5:7:10",
                             default="5:7:10")  ###check sensible defaults!                                                                               
    parser_main.add_argument("--outputprefix",
                            help="the prefix for the output. results in outputprefix_tagnameXX.fastq files.",
                            required = True)
    parser_main.add_argument("--h5score",
                             help="proportion of matching  5' heel_bases for acceptance. between 0 and 1. default 0.9",
                             required = True,
                             default=0.9) 
    parser_main.add_argument("--h3score",
                             help="proportion of matching 3' heel_bases for acceptance. between 0 and 1. default 0.9",
                             required = True,
                             default=0.9) 
    parser_main.add_argument("--filtering",
                             help="parameters for filtering: status:min_len:max_len:mean_qual:min_qual default 0:70:150:30:20",
                             required=True,
                             default="0:70:150:30:20")  # what are good defaults for filtering?
    args = parser_main.parse_args()

    job_summary = ["starting Job with the following parameters at {} :".format(start_time_human)]
    job = dict()
    job["heel_file"] = args.heel
    job_summary.append("heel file is: "+str(job['heel_file']))
    job["tag_file"] = args.tags
    job_summary.append("tag file is: "+str(job['tag_file']))
    if args.h3score and args.h5score is None:
        job_summary.append("using default values for primerscores")
        job["h3score"] =float(0.9) 
        job["h5score"] =float(0.9) 
    else:
        job["h3score"] = float(args.h3score)
        job["h5score"] = float(args.h5score)
    
    job["output_prefix"] = args.outputprefix 

    job_summary.append("primerscore (match/length ratio)  3'>{}  5'>{}".format(job["h3score"],job["h5score"]))
    if args.pairedfastq is True:
        job["filetype"] = "pairedfastq"
        job["forward_file"] = args.forward
        job["reverse_file"] = args.reverse
    if args.pairedfastq is not True:
        job["filetype"]= "unpairedfastq"
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
        param_list  = args.pairedparams.split(":")
        job["kmer_overlap"]= int(param_list[0])
        job["hsp_overlap"]= int(param_list[1])
        job["min_overlap"]= int(param_list[2])
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
    filter_list  = args.filtering.split(":")
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

#flprsr

def get_tagset(tag_input_file, end):
    '''                                                                                                                                                 
    tagfile in tagname;5'tag;3'tag\n order.
    End=3 makes a dict with the 3'tags as keys and matching 5'tags & tagname
    `'''
    tags = dict()
    tag_file = open(str(tag_input_file), "r")
    tag_length = 0

    # Parse tag-file and perform sanity checks
    for line in tag_file:
        line_list_a = list()
        line = line.rstrip()
        if line == "":
            continue
        if line.startswith("#"): # skip comments                                                                                                                                                               
            continue

        line_list_a = line.split(";")
        #now assumes definitive tag_combinations. e.g. tagname;5'tag;3'tag\n,                                                                                                                                  

        if end == "3": #if 3' , shuffle tags in order to enable gain of 3' and 5' dicts from the same file.                                                                                                    
            order = [0,2,1]
            line_list =[line_list_a[i] for i in order]
            tag_summary ="looking for the following tags: 3' 5'\n"
        elif end == "5":
            line_list = line_list_a
            tag_summary ="looking for the following tags: 5' 3' \n"

        #print line_list                                                                                                                                                                                       
        if len(line_list) < 2:
            continue
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
    def __init__(self, p5, p5s, p3, p3s, t5, t3):
        self.p5s = p5s
        self.p3s = p3s
        self.t5 = t5
        self.t3 = t3
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
                    return result, None 
                else:
                    return result
        else:
            p5_pos=0

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
                    return result, None 
                else:
                    return result 

            if len(self.t3):
                tag_len = self.t3["_____length"]
                tag_seq = seq_str[(p3_pos + len(self.p3)):(p3_pos + len(self.p3) + tag_len)]
                tag_seq = str(Seq(tag_seq,generic_dna).reverse_complement())
                try:
                    if len(accepted_t3) and tag_seq not in accepted_t3:     
                        result['status']='chimeric_tag'
                        result["tag_name"] += ";" + self.t3[tag_seq]['tag'] 
                        return result, q 
                except KeyError:
                    result["status"] = "no_tag3"
                    if q:
                        return result, None
            else:
                result["tag_name"] += ""

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
                hqr.append(dict(start = t_start,
                                end = t_end,
                                mean_quality = 0,
                                length = 0))
                t_start = -1
                t_end = -1
    if t_start != -1 and t_end == -1:
        hqr.append(dict(start = t_start,
                        end = len(qual.quals) - 1,
                        mean_quality = 0,
                        length = 0))
    # Check and modify regions to ensure mean quality is above limit                                                                                                                                           
    #print "Before:", hqr                                                                                                                                                                                      
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
    #print "After:", hqr                                                                                                                                                                                       

    hqr = filter(lambda a: a["length"] >= min_length and a["mean_qual"] >= mean_min,
                             hqr)

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


#execute and stats

def assess_results(result_list):
    processed_reads = len(result_list)
    keep = []
    nokeep = []
    for result in result_list:
        if "keep" in result:
            keep.append(result)
        else:
            nokeep.append(result)
    kept_reads = len(keep)
    tag_count = Counter([read["detagged_seq"]["tag_name"] for read in keep]).items()
    discarded_reads = len(nokeep)
    status_list = [i["status"].keys()[0] for i in nokeep]
    status_count = Counter(status_list).items()
    
    chim_list = [i["detagged_seq"]["tag_name"] for i in result_list if "chimeric_tag" in  i["status"].keys()]
    chim_combination_count = Counter(chim_list).most_common(10) 
    chim_tag_lists = [combination.split(";") for combination in chim_list]
    chim_tags = [tag for combination in chim_tag_lists for tag in combination]
    chim_tag_count = Counter(chim_tags).most_common(20)
    
    nokeepstats = [("processed_reads_in_total:", processed_reads)]
    nokeepstats.extend([("kept:", kept_reads)])
    nokeepstats.extend([("discarded:", discarded_reads)])
    nokeepstats.extend([("\t", "\t")])
    nokeepstats.extend([("tag:", "valid_reads_per_tag")])
    nokeepstats.extend(tag_count)
    nokeepstats.extend([("\t", "\t")])
    nokeepstats.extend([("cause_for_failure", "amount_of_reads")])
    nokeepstats.extend(status_count)
    nokeepstats.extend([("\t", "\t")])
    nokeepstats.extend([("most_common chimeras", "amount_of_reads")])
    nokeepstats.extend(chim_combination_count)
    nokeepstats.extend([("\t", "\t")])
    nokeepstats.extend([("most_common_chimeric_tags", "amount_of_reads")])
    nokeepstats.extend(chim_tag_count)
    return keep,  nokeepstats

def write_fastqs(keep, tags, out_prefix):
    key_list = [tags[key]["tag"] for key in tags.keys() if key is not "_____length"]
    reads_by_tag_dict={key:[] for key in key_list}
    for read in keep:
        rec = read["detagged_seq"]["seq_record"]
        rec.letter_annotations["phred_quality"] = read["detag_qual"].quals
        reads_by_tag_dict[read["detagged_seq"]["tag_name"]].append(rec)
    for key in reads_by_tag_dict.keys():
        with open(out_prefix+"_"+str(key)+".fastq", "w") as handle:
            SeqIO.write(reads_by_tag_dict[key], handle, "fastq")

def write_logfile(out_prefix,jobsum, heelsum, tagsum, stats ):
    with open(out_prefix+"_logfile.txt" ,"w") as handle:
        for i in jobsum:
            handle.write(i+"\n")
        handle.write(heelsum+"\n")
        handle.write(tagsum+"\n")
        for i in stats:
            handle.writelines(i[0]+"\t"+str(i[1])+"\n")

def main():
    job,jobsum = cli_parser()
    tags5,tagsum = get_tagset(tag_input_file=job["tag_file"],end='3')
    tags3 = get_tagset(tag_input_file=job["tag_file"],end='5')[0]
    heel5, heel3, heelsum = get_heels(job["heel_file"])
    seq_data = Pair(fastq_1=job["forward_file"],
                    fastq_2=job["reverse_file"],
                    hsp=job["hsp_overlap"], 
                    kmer=job["kmer_overlap"], 
                    min=job["min_overlap"])

    ##init DeTagSeq
    dtseq = DeTagSeq(p3=heel3,
                     p3s=job["h3score"],
                     p5=heel5,
                     p5s=job["h5score"],
                     t3=tags3,
                     t5=tags5)

    ## execute detag
    result_list = []

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
                    result_list.append(res)
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
            result_list.append(res)
            continue
        if len(seq.seq) < job["filter_minlen"]:
            res["status"]["too_short"] = 1
            result_list.append(res)
        
            continue
        if not detagged_seq:
            detagged_seq, qual = dtseq.detag_seq(seq, qual) 



    
        if "status" in detagged_seq:
            res["status"][detagged_seq["status"]] = 1
            result_list.append(res)
            res["detagged_seq"] = detagged_seq 
            continue
        if len(detagged_seq["seq"]) < job["filter_minlen"]: 
            res["status"]["too_short"] = 1
            res["detagged_seq"] = detagged_seq 
            result_list.append(res)
            continue

        if job["filter_maxlen"] > 0 and len(detagged_seq["seq"]) > job["filter_maxlen"]:
            res["status"]["truncated"] = 1
            detagged_seq["seq"] = detagged_seq["seq"][:job["filter_maxlen"]]
        if job["filter_maxlen"] < 0:
            res["status"]["truncated"] = 1
            detagged_seq["seq"] = detagged_seq["seq"][:job["filter_maxlen"]]

        res["detagged_seq"] = detagged_seq
        res["keep"] = True
        res["id"] = seq.id
        res["detag_qual"] = qual
        result_list.append(res)
    keep, stats = assess_results(result_list)
    write_fastqs(keep, tags3, job["output_prefix"])
    write_logfile(job["output_prefix"],jobsum, heelsum, tagsum, stats ) 

#call main
main()
