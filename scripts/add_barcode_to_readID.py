import gzip
import sys
import argparse
import pysam
import regex
import os
import yaml



def detect_files(files_list):
    reads = {'R1': [files_list[i] for i, x in enumerate(files_list) if '_R1_' in files_list[i]],
             'I2': [files_list[i] for i, x in enumerate(files_list) if '_R2_' in files_list[i]],
             'R2': [files_list[i] for i, x in enumerate(files_list) if '_R3_' in files_list[i]],
             }

    ab_barcode = [len(x) == 1 for x in reads.values()]      # Check if return exactly one file per R1,R2,R3
    if sum(ab_barcode) != 3:                                # Check if one file is returned for all 3 files R1,R2,R3
        return None

    reads = {key: reads[key][0] for key in reads}  # Unlist
    return reads


def find_seq(pattern, DNA_string, nmismatch):
    r = regex.compile('({0}){{e<={1}}}'.format(pattern, nmismatch))
    res = r.finditer(DNA_string)
    hit = [x.start() for x in res]
    if len(hit) > 1:
        return "Multiple"
    if len(hit) == 0:
        return None
    if len(hit) == 1:
        return hit[0]


def print_fastq(read):
    return '{}\n{}\n{}\n{}\n'.format(read[0], read[1], read[2], read[3])


def generate_output_folder(out_prefix, fastq_reads):
    out_reads = {key: out_prefix + os.path.basename(fastq_reads[key]).replace("_R3_", "_R2_") for key in
                 fastq_reads.keys() if key != "I2"}  # Treat R3 as R2, I2 is skipped
    out_reads = {key: out_reads[key] + '.gz' if not out_reads[key].endswith('.gz') else out_reads[key] for key in
                 out_reads.keys()}  # Add .gz if not there
    return out_reads


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

def get_barcodes(string,index,arguments):
    ab_barcode = revcompl(string.sequence[index - 8:index])  # Get the barcode sequence
    cell_barcode = None
    if arguments.single_cell:
        cell_bcd_start = index + len(arguments.pattern)
        cell_barcode = revcompl(string.sequence[cell_bcd_start:cell_bcd_start + 16])  # Get the cell barcode
    return(ab_barcode,cell_barcode)

def main(arguments):

    fastq_reads = detect_files(arguments.input)
    if not fastq_reads:
        sys.exit('*** ERROR: Provide single _R1_ _R2_ and _R3_ files ***')

    out_reads = generate_output_folder(arguments.out_prefix, fastq_reads)

    statistics = {
        "barcode_found": 0,
        "multiple_matches": 0,
        "no_barcode": 0
    }

    if not os.path.isdir(arguments.out_prefix):
        os.makedirs(arguments.out_prefix)

    with pysam.FastxFile(fastq_reads['R1']) as f1, \
            pysam.FastxFile(fastq_reads['I2']) as f2, \
            pysam.FastxFile(fastq_reads['R2']) as f3, \
            gzip.open(out_reads['R1'], 'wt') as out1, \
            gzip.open(out_reads['R2'], 'wt') as out2:

        for line1, line2, line3 in zip(f1, f2, f3):          # Iterate over 3 fastq files
            assert (line1.name == line2.name == line3.name)  # Make sure the fastq files are ok
            for n in range(0, 3):                            # Iterate through 0-1-2 mismatches, if hit found break
                hit = find_seq(arguments.pattern, line2.sequence, nmismatch=n)
                if not hit:
                    statistics["no_barcode"] += 1
                    continue
                if hit == 'Multiple':
                    statistics["multiple_matches"] += 1
                    sys.stderr.write("*** Multiple patterns found for read; skipping  ***")
                    break
                # Else number of hits == 1
                hit = int(hit)
                ab_barcode, cell_barcode = get_barcodes(line2,hit,arguments)

                line1.name = line1.name + "_" + ab_barcode  # Add barcode sequence to read name
                line3.name = line3.name + "_" + ab_barcode

                # Write the outputs
                out1.write('{}\n'.format(str(line1)))
                out2.write('{}\n'.format(str(line3)))
                statistics["barcode_found"] += 1
                break
    # Write the statistics file
    with open(arguments.out_prefix + '/statistics.yaml', 'w') as f:
        yaml.dump(statistics, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-i', '--input',
                        required=True,
                        type=str,
                        nargs='+',
                        help='path to input R1,R2,R3 .fastq files [3 files required]')


    parser.add_argument('-o', '--out_prefix',
                        type=str,
                        required=True,
                        help='Prefix to where to put the output files; Diretory will be created')

    parser.add_argument('-p', '--pattern',
                        type=str,
                        default="GCGTGGAGACGCTGCCGACGA",
                        help='Pattern that follows the antibody barcode \n \
                                  (Default: %(default)s)')

    parser.add_argument('--single-cell',
                        default=False,
                        action='store_true',
                        help='Data is single cell CUT&Tag')

    args = parser.parse_args()
    main(args)
