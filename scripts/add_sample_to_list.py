import sys

sample = sys.argv[1]

for line in sys.stdin:
    line     = line.rstrip()
    new_line = "{}_{}\n".format(sample,line)
    sys.stdout.write(new_line)
