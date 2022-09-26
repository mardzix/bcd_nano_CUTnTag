import pysam
import sys

samfile = pysam.AlignmentFile(sys.argv[1], "rb")

reads_unique = {}
stats = {
    'not mapped': 0,
    'unique': 0,
    'LA duplicates': 0,
    'PCR duplicates': 0,
    'mate': 0
}

n = 0
# n_target = 10000 # test for 10 thousand reads
# n_target = 10000000 # test for 10 milion reads

read1 = False
read2 = False

for line in samfile:
    n += 1
    if n % 10000 == 0:
        sys.stderr.write('*** {} lines processed\n'.format(n))

    # Only take the first n_target reads # Debugging
    # if n == n_target:
    #     # print(reads_unique)
    #     print(stats)
    #     break

    read1 = line

    if read1.is_unmapped or read1.mate_is_unmapped:
        stats['not mapped'] += 1
        continue

    if read1.is_read1:
        stats['mate'] += 1
        continue

    # read2 = samfile.mate(line)
    cell_barcode = read1.get_tag('CR')

    read1_position = '{}_{}_{}'.format(read1.reference_name, + read1.reference_start, cell_barcode)
    read2_position = '{}_{}_{}'.format(read1.next_reference_name, + read1.next_reference_start, cell_barcode)

    try:  # Creates new fw read entry in the dictionary if not there
        reads_unique[read1_position]
    except KeyError:
        reads_unique[read1_position] = []

    # Case1 read2 list is empty == new read is observed             == new unique read
    if reads_unique[read1_position] == []:
        reads_unique[read1_position].append(read2_position)
        stats['unique'] += 1
        continue

    # Case read2 list is not empty and the read2 is not in list    == LA duplicate
    if read2_position not in reads_unique[read1_position]:
        reads_unique[read1_position].append(read2_position)
        stats['LA duplicates'] += 1
        continue

    # Case read2 list is not empty and read2 is in the list        == PCR duplicate
    if read2_position in reads_unique[read1_position]:
        stats['PCR duplicates'] += 1
        continue

    # Sanity check if something escapes
    sys.exit('found exception')

print(stats)