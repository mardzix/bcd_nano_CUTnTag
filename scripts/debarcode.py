import os
import argparse
from glob import glob
import sys
from re import split
import regex
import gzip
from contextlib import ExitStack

import yaml
from pysam import FastxFile
import Levenshtein                          # pip install python-Levenshtein

class bcdCT:
    def __init__(self,args):
        self.detect_input(args.input)
        self.detect_reads()
        self.single_cell=args.single_cell
        self.out_prefix=args.out_prefix

        if args.name:
            self.name = args.name
        else:
            self.autodetect_name()

        self.autodetect_barcodes(args)

    def detect_input(self,input):
        Error_message="*** Error: Wrong input files specified. The input must be either folder with _R1_*.fastq.gz _R2_*.fastq.gz _R3_*.fastq.gz files or paths to the files themselves ***\n" +\
        "The files should be placed in the same folder\n" +\
        "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/\n"

        input = [os.path.abspath(x) for x in input]
        if len(input) == 1 and os.path.isdir(input[0]):     # Case input is single directory
            self.input_dir = input[0]
            self.input_files = []
            self.input_files.extend(glob(self.input_dir + "/*.fastq.gz"))
            self.input_files.extend(glob(self.input_dir + "/*.fq.gz"))

        elif len(input) > 1:                                  # Case input are multiple files
            self.input_files = input
            self.input_dir = list(set([os.path.dirname(x) for x in self.input_files]))
            if not len(self.input_dir) == 1:
                sys.stderr.write(Error_message)
                sys.exit(1)
            if not sum([x.endswith('.fastq.gz') or x.endswith('.fq.gz') for x in self.input_files]) == len(self.input_files):
                sys.exit(1)
                sys.stderr.write(Error_message)
        else:
            sys.exit(1)
            sys.stderr.write(Error_message)

    def detect_reads(self):
        Error_message="*** Error: Please specify exactly one _R1_ _R2_ and _R3_ file or folder with exactly one of each files ***\n" + \
                      "e.g. /data/path_to_my_files/*L001*.fastq.gz or /data/path_to_my_files/\n"
        self.path_in = {}
        self.path_in['R1'] = [x for x in self.input_files if "_R1_" in x]
        self.path_in['R2'] = [x for x in self.input_files if "_R2_" in x]
        self.path_in['R3'] = [x for x in self.input_files if "_R3_" in x]

        if len(self.path_in['R1']) != 1 or len(self.path_in['R2']) != 1 or len(self.path_in['R3']) != 1:
            sys.stderr.write(Error_message)
            sys.exit(1)

        self.path_in = {key:self.path_in[key][0] for key in self.path_in.keys()}

    def autodetect_name(self):
        Error_message="*** Error: Prefix for R1 R2 R3 files not the same. Please use the same prefix for all the files or specify experiment name ***\n"

        self.name = [split("_R[0-9]_",str(x)) for x in self.path_in.values()]
        self.name = [x[0] for x in self.name]

        if len(list(set(self.name))) > 1:
            sys.stderr.write(Error_message)
            sys.exit(1)

        self.name=self.name[0].split("/")[-1]

    def in_handles(self,stack):
        in_stack = {x: stack.enter_context(FastxFile(self.path_in[x],'r')) for x in ['R1','R2','R3']}
        return in_stack

    def create_out_handles(self,stack):
        for bcd in self.picked_barcodes:
            os.makedirs(self.out_prefix + "/barcode_" + bcd, exist_ok=True)
        if self.single_cell:
            out_reads = ['R1','R2','R3']
        else:
            out_reads = ['R1','R3']

        self.path_out = {barcode: {} for barcode in self.picked_barcodes}
        self.path_out  = {barcode: {read: "{0}/barcode_{1}/{2}".format(self.out_prefix,barcode,os.path.basename(self.path_in[read])) for read in out_reads} for barcode in self.picked_barcodes}
        self.out_stack = {barcode: {read: stack.enter_context(gzip.open(self.path_out[barcode][read],'wt'))for read in out_reads} for barcode in self.picked_barcodes}


    def __iter__(self):
        with FastxFile(self.path_in['R1']) as f1, FastxFile(self.path_in['R2']) as f2, FastxFile(self.path_in['R3']) as f3:
            for r1,r2,r3 in zip(f1,f2,f3):
                yield r1, r2, r3

    def autodetect_barcodes(self,args):
        barcodes = {}
        n=0
        for read1,read2,read3 in self:
            hit = find_seq(args.pattern, read2.sequence, nmismatch=0)
            if not hit or hit == 'Multiple':
                continue
            hit = int(hit)
            read_barcode = get_read_barcode(read2, hit)
            try:
                barcodes[read_barcode] += 1
            except KeyError:
                barcodes[read_barcode] = 1
            n += 1
            if n == 50000:
                break

        top_barcodes = sorted(barcodes, key=barcodes.get, reverse=True)[:args.Nbarcodes]
        picked_barcodes = {key: barcodes[key] for key in top_barcodes}
        sys.stderr.write("\nDetected following most abundant barcodes out of first {} barcodes:\n{}\n".format(n, picked_barcodes))
        if args.barcode != "None":
            self.picked_barcodes = [args.barcode]
            sys.stderr.write("Barcode specified for demultiplexing [{barcode}] in top found barcodes: {bool} \n".format(bool = args.barcode in picked_barcodes.keys(), barcode = args.barcode))
            return

        self.picked_barcodes = picked_barcodes

def get_read_barcode(string,index):
    read_barcode = revcompl(string.sequence[index - 8:index])  # Get the barcode sequence
    return read_barcode

def extract_cell_barcode(read,index):
    cell_bcd_start = index + len(args.pattern)
    read.sequence = read.sequence[cell_bcd_start:cell_bcd_start + 16]  # Get the cell barcode
    read.quality  = read.quality[cell_bcd_start:cell_bcd_start + 16]   # Get corresponding Quality score
    return read

def revcompl(seq):
    revcomp_table = {
        "A": "T",
        "G": "C",
        "C": "G",
        "T": "A",
        "N": "N"
    }
    complement = "".join([revcomp_table[letter] for letter in seq.upper()])  # Complement
    return complement[::-1]  # Reverse

def rev(seq):
    return seq[::-1]

def find_seq(pattern, DNA_string, nmismatch=2):
    for n in range(0,nmismatch + 1):
        r = regex.compile('({0}){{e<={1}}}'.format(pattern, n))
        res = r.finditer(DNA_string)
        hit = [x.start() for x in res]
        if len(hit) == 0:
            continue
        if len(hit) > 1:
            return None
        if len(hit) == 1:
            return int(hit[0])
    return None

def main(args):
    exp = bcdCT(args)

    statistics = {
        "barcode_found": 0,
        "multiple_barcode_matches": 0,
        "no_barcode_match": 0,
        "no_spacer_found": 0,
        "too_short_read": 0
    }

    sys.stderr.write("Creating file output handles \n")
    with ExitStack() as stack:
        exp.create_out_handles(stack)
        n = 0
        sys.stderr.write("Starting demultiplexing \n")
        for read1,read2,read3 in exp:
            n+=1
            if n % 500000 == 0:
                sys.stderr.write("{} reads processed\n".format(n))
            assert (read1.name == read2.name == read3.name)                                                 # Make sure the fastq files are ok

            spacer_hit = find_seq(pattern=args.pattern,DNA_string=read2.sequence,nmismatch=2)
            if not spacer_hit:
                statistics["no_spacer_found"] +=1

            if spacer_hit:
                read_barcode = get_read_barcode(read2, spacer_hit)                                               # Returns only barcode e.g. ACTGACTG
                read_barcode_distance = {barcode: Levenshtein.distance(read_barcode,barcode) for barcode in exp.picked_barcodes}
                if sum([x <= int(args.mismatch) for x in read_barcode_distance.values()]) == 0:
                    statistics["no_barcode_match"] += 1
                    continue
                if sum([x <= args.mismatch for x in read_barcode_distance.values()]) > 1:
                    statistics["multiple_barcode_matches"] += 1
                    continue

                hit_barcode = min(read_barcode_distance,key=read_barcode_distance.get)

                if exp.single_cell:
                    read2 = extract_cell_barcode(read2, spacer_hit)  # Returns the whole read, only the cell barcode part
                    if len(read2.sequence) < 16:
                        statistics["too_short_read"] += 1
                        continue
                    exp.out_stack[hit_barcode]['R2'].write('{}\n'.format(str(read2)))

                # Write the outputs
                exp.out_stack[hit_barcode]['R1'].write('{}\n'.format(str(read1)))
                exp.out_stack[hit_barcode]['R3'].write('{}\n'.format(str(read3)))


                statistics["barcode_found"] += 1


    # Write the statistics file
    with open("{0}/{1}_statistics.yaml".format(exp.out_prefix,exp.name), 'w') as f:
        yaml.dump(statistics, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        nargs='+',
                        help='path to input R1,R2,R3 .fastq.gz files [3 files required]')

    parser.add_argument('-o', '--out_prefix',
                        type=str,
                        required=True,
                        help='Prefix to where to put the output files; Diretory will be created')

    parser.add_argument('-p', '--pattern',
                        type=str,
                        default="GCGTGGAGACGCTGCCGACGA",
                        help='Pattern that follows the antibody barcode \n \
                                  (Default: %(default)s)')

    parser.add_argument('--single_cell',
                        default=False,
                        action='store_true',
                        help='Data is single cell CUT&Tag (Default: %(default)s)')

    parser.add_argument('--name',
                        type=str,
                        default=None,
                        help='Custom name for the experiment (Default: Autodetect from filename)')

    parser.add_argument('--mismatch',
                        type=int,
                        default=1,
                        help='Maximum mismatches for sample barcode (Default: %(default)s)')

    parser.add_argument('--Nbarcodes',
                        type=int,
                        default=3,
                        help='Number of barcodes in experiment (Default: %(default)s)')

    parser.add_argument('--barcode',
                        type=str,
                        default='None',
                        help='Specific barcode to be extracted [e.g. ATAGAGGC] (Default: All barcodes [see --Nbarcodes])')

    args = parser.parse_args()
    main(args)