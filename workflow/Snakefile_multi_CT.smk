gopalan_SRA = 'SRR14150925'

rule all:
    input:
        'foo'


rule gopalan_fastq_dump:
    output:
        "results/
    threads: 10
    resources:
        load = 50
    params:
        tmp = config['general']['tempdir'],
        out = "results/gopalan_multi_CT/fastq/",
    shell:
        "fasterq-dump -t {params.tmp} -f -e {threads} --split-files --include-technical -o {params.out} {gopalan_SRA}"
