include: 'Snakefile_prep.smk'
include: 'Snakefile_preprocess.smk'
include: 'Snakefile_pre_nbiotech.smk'

def get_fastq_basename_for_sample(prefix): #samples is str
    import glob,os
    files = glob.glob(prefix + '/**/*.fastq.gz',recursive=True)
    files = [os.path.basename(str(x)) for x in files]
    return(files)

def get_sample_id_from_fastq(prefix):
    import glob, os
    import re
    sample = glob.glob(prefix + '/**/*.fastq.gz',recursive=True)
    sample = [os.path.basename(x) for x in sample]
    sample = [re.split("_S[0-9]_",x)[0] for x in sample] # TODO fix as regexp
    sample = list(set(sample))
    return(" ".join(sample))

rule all:
    input:
        # Cellranger
        # ['results/downscale/cellranger/{antibody}/nbiotech/{sample}/{fastq}.log'.format(antibody='H3K27me3',sample=sample,fastq=f) for sample in ['H3K27me3_N1', 'H3K27me3_N2', 'H3K27me3_N3','H3K27me3_N4'] for f in get_fastq_basename_for_sample(sample, prefix = 'results/nbiotech_data/cellranger/fastq_final/')],
        expand('results/downscale/cellranger/H3K27me3/nbiotech/{sample}/outs/possorted_bam.bam',sample = ['H3K27me3_N1','H3K27me3_N2','H3K27me3_N3','H3K27me3_N4']),
        expand('results/downscale/cellranger/H3K27me3/nano_CT/{sample}/outs/possorted_bam.bam',sample = ['bcdCT_MB21_02','bcdCT_MB21_03','bcdCT_MB21_04','bcdCT_MB21_05'] ),
        # Mod fragments file
        expand('results/downscale/cellranger/{antibody}/{method}/{sample}/outs/fragments_mod.tsv.gz', antibody = 'H3K27me3', method = 'nano_CT',sample = ['bcdCT_MB21_02','bcdCT_MB21_03','bcdCT_MB21_04','bcdCT_MB21_05']),
        expand('results/downscale/cellranger/{antibody}/{method}/{sample}/outs/fragments_mod.tsv.gz', antibody = 'H3K27me3', method = 'nbiotech', sample = ['H3K27me3_N1','H3K27me3_N2','H3K27me3_N3','H3K27me3_N4']),


rule find_lowest_read_count: #
    input:
        glob.glob('results/nbiotech_data/cellranger/fastq_final/**/*.fastq.gz',recursive=True),
        glob.glob('results/multimodal_data/*/fastq_per_barcode/H3K27me3_*/**/*.fastq.gz',recursive=True),
    output:
        'results/downscale/summary.txt'
    shell:
        'ls {input} | while read line; do echo -n $line"; ";gunzip -c $line | wc -l ;  done > {output}'


rule downscale_fastq_nbiotech:
    input:
        'results/nbiotech_data/cellranger/fastq_final/{sample}/{fastq}',
    output:
        'results/downscale/fastq_nbiotech/{antibody}/{sample}/{fastq}',
    params:
        nreads = 30000000
    shell:
        'seqtk sample -s100 {input} {params.nreads} | gzip > {output} '

rule downscale_fastq_nano_CT:
    input:
        'results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_*/barcode_*/{fastq}'
        # lambda wildcards: glob.glob('results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_*/barcode_*/{fastq}'.format(antibody=wildcards.antibody, barcode =  barcodes_dict[wildcards.sample][wildcards.antibody], fastq=wildcards.fastq, sample = wildcards.sample )),
    output:
        'results/downscale/fastq_nano_CT/{antibody}/{sample}/{fastq}',
    params:
        nreads = 30000000
    shell:
        'seqtk sample -s100 {input} {params.nreads} | gzip > {output}'

rule cellranger_nbiotech:
    input:
        fastq = lambda wildcards: expand('results/downscale/fastq_nbiotech/{antibody}/{sample}/{{fastq}}'.format(antibody = wildcards.antibody, sample = wildcards.sample), fastq = get_fastq_basename_for_sample(prefix = 'results/nbiotech_data/cellranger/fastq_final/' + wildcards.sample))
    output:
        bam       = 'results/downscale/cellranger/{antibody}/nbiotech/{sample}/outs/possorted_bam.bam',
        fragments = 'results/downscale/cellranger/{antibody}/nbiotech/{sample}/outs/fragments.tsv.gz',
        # log         = 'results/downscale/cellranger/{antibody}/nbiotech/{sample}/log.txt'
    params:
        out_folder     = 'results/downscale/cellranger/{antibody}/nbiotech/',
        cellranger_ref = config['general']['cellranger_ref'],
        fastq_folder   = os.getcwd() + '/results/downscale/fastq_nbiotech/{antibody}/{sample}/',
        samples        = lambda wildcards: get_sample_id_from_fastq(os.getcwd() + '/results/nbiotech_data/cellranger/fastq_final/{sample}/'.format(sample = wildcards.sample))
    threads: 40
    shell:
        'rm -r {params.out_folder}/{wildcards.sample}; '
        'mkdir {params.out_folder}; '
        'cd {params.out_folder}; '
        '/data/bin/cellranger-atac count --id {wildcards.sample} --reference {params.cellranger_ref} --fastqs {params.fastq_folder} --sample {params.samples}'

rule cellranger_nano:
    input:
        fastq     = lambda wildcards: expand('results/downscale/fastq_nano_CT/{antibody}/{sample}/{{fastq}}'.format(antibody = wildcards.antibody,sample=wildcards.sample),fastq=get_fastq_basename_for_sample(prefix='results/nbiotech_data/cellranger/fastq_final/' + wildcards.sample))
    output:
        bam       = 'results/downscale/cellranger/{antibody}/nano_CT/{sample}/outs/possorted_bam.bam',
        fragments = 'results/downscale/cellranger/{antibody}/nano_CT/{sample}/outs/fragments.tsv.gz',
        # log = 'results/downscale/cellranger/{antibody}/nano_CT/{sample}/log.txt'
    params:
        out_folder     = 'results/downscale/cellranger/{antibody}/nano_CT/',
        cellranger_ref = config['general']['cellranger_ref'],
        fastq_folder   = os.getcwd() + '/results/downscale/fastq_nano_CT/{antibody}/{sample}/',
        samples        = lambda wildcards: get_sample_id_from_fastq(os.getcwd() + '/results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_*/barcode_*/'.format(antibody=wildcards.antibody,sample=wildcards.sample))
    threads: 40
    shell:
        'rm -rf {params.out_folder}/{wildcards.sample}; '
        'mkdir {params.out_folder}; '
        'cd {params.out_folder}; '
        '/data/bin/cellranger-atac count --id {wildcards.sample} --reference {params.cellranger_ref} --fastqs {params.fastq_folder} --sample {params.samples}'

rule add_sample_name_to_fragments:
    input:
        fragments = 'results/downscale/cellranger/{antibody}/{method}/{sample}/outs/fragments.tsv.gz'
    output:
        fragments = 'results/downscale/cellranger/{antibody}/{method}/{sample}/outs/fragments_mod.tsv.gz'
    params:
        script = workflow_dir + '/scripts/add_sample_to_fragments.py',
    shell:
        'python3 {params.script} {input.fragments} {wildcards.sample} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments}'