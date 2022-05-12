import sys
import pysam
from contextlib import ExitStack

bamFile             = sys.argv[1]      # Possorted bam file from 10x with barcodes           e.g. AAACGAAAGTTGCTTG-1
sampleID            = sys.argv[2]      # Sample name to prefix to bam barcode                e.g. e.g. bcdCT_MB21_02, making the final barcodes like: bcdCT_MB21_02_AAACGAAAGTTGCTTG-1
outFile             = sys.argv[3]      # bam file to write the output to                     e.g. /path/to/output/possorted_bam.bam

bamfile = pysam.AlignmentFile(bamFile, "rb")
header  = bamfile.header

n=0
with pysam.AlignmentFile(outFile,'wb',header=header) as f:
    for line in bamfile:
        n += 1
        if n % 1000000 == 0:
            sys.stderr.write("*** {} lines processed ***\n".format(n))
        try:
            barcode = line.get_tag("CB")
        except KeyError:
            continue
        line.set_tag("CB",sampleID + "_" + barcode)
        f.write(line)
