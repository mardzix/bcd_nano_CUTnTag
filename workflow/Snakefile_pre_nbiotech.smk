configfile: workflow.basedir + '/../config/config.yaml'

samples_list_nbiotech = [x for x in config['nbiotech_data']['samples'].keys()]
# ['H3K27ac_N1', 'H3K27ac_N2', 'H3K27me3_N1', 'H3K27me3_N2', 'H3K27me3_N3', 'H3K27me3_N4']

SRA_dict = config['nbiotech_data']['samples']
SRA_dict = {x: SRA_dict[x].split(" ") for x in SRA_dict.keys()}
# {'H3K27ac_N1': ['SRR12607042', 'SRR12607043'], 'H3K27ac_N2': ['SRR12607044', 'SRR12607045'], 'H3K27me3_N1': ['SRR12607034', 'SRR12607035'], 'H3K27me3_N2': ['SRR12607036', 'SRR12607037'], 'H3K27me3_N3': ['SRR12607038', 'SRR12607039'], 'H3K27me3_N4': ['SRR12607040', 'SRR12607041']}


shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; conda activate " + config['general']['conda_env']  + " ; ")

rename_fastq_dic = {"R1": "2", "R2": "3", "R3": "4", "R4": "1"}

localrules: nbiotech_fastq_dump

rule nbiotech_all:
    input:
        expand('results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam', sample= samples_list),
        expand('results/nbiotech_data/{sample}/seurat/bin_{binwidth}/Seurat_object.Rds', sample = samples_list,binwidth=5000),

        # Seurat direct download
        'results/nbiotech_data/data/seurat/H3K27me3_seurat_object.Rds'

rule nbiotech_fastq_dump:
    output:
        "results/nbiotech_data/cellranger/fastq/{sample}/{SRA}_1.fastq",
        "results/nbiotech_data/cellranger/fastq/{sample}/{SRA}_2.fastq",
        "results/nbiotech_data/cellranger/fastq/{sample}/{SRA}_3.fastq",
        "results/nbiotech_data/cellranger/fastq/{sample}/{SRA}_4.fastq",
    threads: 10
    params:
        tmp = config['general']['tempdir'],
        out = "results/nbiotech_data/cellranger/fastq/{sample}/{SRA}.fastq",
    shell:
        "fasterq-dump -t {params.tmp} -f -e {threads} --split-files --include-technical -o {params.out} {wildcards.SRA}"

rule nbiotech_gzip_fastq:
    input:
        "results/nbiotech_data/cellranger/fastq/{sample}/{file}.fastq",
    output:
        "results/nbiotech_data/cellranger/fastq/{sample}/{file}.fastq.gz"
    shell:
        "gzip {input}"

rule nbiotech_rename_fastq:
    input:
        lambda wildcards: expand("results/nbiotech_data/cellranger/fastq/{sample}/{file}_{n}.fastq.gz", n = rename_fastq_dic[wildcards.read], sample = wildcards.sample, file = wildcards.SRA)
    output:
        "results/nbiotech_data/cellranger/fastq_final/{sample}/{SRA}_S1_L001_{read}_001.fastq.gz"
    shell:
        "mv {input} {output}"

rule nbiotech_run_cellranger_nbiotech:
    input:
        lambda wildcards: expand("results/nbiotech_data/cellranger/fastq_final/{sample}/{SRA}_S1_L001_{read}_001.fastq.gz", sample=wildcards.sample, SRA = SRA_dict[wildcards.sample], read = ["R1","R2","R3"])
    output:
        bam  = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam',
        frag = 'results/nbiotech_data/cellranger/{sample}/outs/fragments.tsv.gz',
        meta = 'results/nbiotech_data/cellranger/{sample}/outs/singlecell.csv',
    params:
        cellranger_ref = '/data/ref/cellranger-atac/refdata-cellranger-atac-mm10-2020-A-2.0.0/',
        sample         = lambda wildcards: SRA_dict[wildcards.sample],
        fastq_dir      = lambda wildcards: os.getcwd() + "/results/nbiotech_data/cellranger/fastq_final/{sample}/".format(sample=wildcards.sample)
    threads: 40
    shell:
        'rm -r results/nbiotech_data/cellranger/{wildcards.sample}/; '
        'cd results/nbiotech_data/cellranger/; '
        '/data/bin/cellranger-atac count --id {wildcards.sample} --sample {params.sample} --reference {params.cellranger_ref} --fastqs {params.fastq_dir}'


rule nbiotech_bam_to_bw:
    input:
        cellranger_bam = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam'
    output:
        bigwig         = 'results/nbiotech_data/{sample}/bigwig/all_reads.bw'
    threads: 16
    shell:
        'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule nbiotech_run_macs_narrow:
    input:
        cellranger_bam = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam'
    output:
        narrow_peaks = 'results/nbiotech_data/{sample}/peaks/macs_narrow/{sample}_peaks.narrowPeak'
    params:
        macs_outdir  = 'results/nbiotech_data/{sample}/peaks/macs_narrow/'
    shell:
        'macs2 callpeak -t {input.cellranger_bam} -g mm -f BAMPE -n {wildcards.sample} '
        '--outdir {params.macs_outdir} --keep-dup=1 --llocal 100000 --cutoff-analysis --min-length 1000 --max-gap 1000  2>&1 '

rule nbiotech_run_macs_narrow_test:
    input:
        cellranger_bam = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam'
    output:
        narrow_peaks = 'results/nbiotech_data/{sample}/peaks/macs_narrow_llocal_{llocal}/{sample}_peaks.narrowPeak'
    params:
        macs_outdir  = 'results/nbiotech_data/{sample}/peaks/macs_narrow_llocal_{llocal}/'
    shell:
        'macs2 callpeak -t {input.cellranger_bam} -g mm -f BAMPE -n {wildcards.sample} '
        '--outdir {params.macs_outdir} --keep-dup=1 --llocal {wildcards.llocal} 2>&1 '

rule nbiotech_run_macs_broad:
    input:
        cellranger_bam = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam'
    output:
        broad_peaks = 'results/nbiotech_data/{sample}/peaks/macs_broad/{sample}_peaks.broadPeak'
    params:
        macs_outdir = 'results/nbiotech_data/{sample}/peaks/macs_broad/'
    shell:
        'macs2 callpeak -t {input} -g mm -f BAMPE -n {wildcards.sample} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup=1 --broad-cutoff=0.1 ' 
        '--min-length 1000 --max-gap 1000 --broad 2>&1 '

rule nbiotech_download_chromsizes:
    output:
        'results/nbiotech_data/mm10.chrom.sizes'
    params:
        url="http:/hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes"
    shell:
        "wget -O {output} {params.url}"

rule nbiotech_prep_SEACR_files:
    input:
        fragments = 'results/nbiotech_data/cellranger/{sample}/outs/fragments.tsv.gz',
        genome = 'results/nbiotech_data/mm10.chrom.sizes'
    output:
        'results/nbiotech_data/{sample}/peaks/SEACR/fragments.bg'
    shell:
        "bedtools genomecov -bg -g {input.genome} -i {input.fragments} > {output}"

rule nbiotech_run_SEACR:
    input:
        'results/nbiotech_data/{sample}/peaks/SEACR/fragments.bg',
    output:
        'results/nbiotech_data/{sample}/peaks/SEACR/peaks.relaxed.bed',
    params:
        out_prefix = 'results/nbiotech_data/{sample}/peaks/SEACR/peaks',
    shell:
        "~/bin/SEACR/SEACR_1.3.sh {input} 0.01 norm relaxed {params.out_prefix}"

rule nbiotech_add_barcode_fragments:
    input:
        fragments = 'results/nbiotech_data/cellranger/{sample}/outs/fragments.tsv.gz',
    output:
        fragments = 'results/nbiotech_data/{sample}/fragments/fragments.tsv.gz',
        index     = 'results/nbiotech_data/{sample}/fragments/fragments.tsv.gz.tbi',
    params:
        script    = os.path.dirname(workflow.basedir) + '/scripts/add_sample_to_fragments.py',
    shell:
        'python3 {params.script} {input.fragments} {wildcards.sample} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments}'

rule nbiotech_barcode_overlap_peaks:
    input:
        bam    = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam',
        peaks  = 'results/nbiotech_data/{sample}/peaks/macs_broad/{sample}_peaks.broadPeak',
    output:
        overlap = 'results/nbiotech_data/{sample}/barcode_metrics/peaks_barcodes.txt'
    params:
        get_cell_barcode     = os.path.dirname(workflow.basedir) + '/scripts/get_cell_barcode.awk',
        add_sample_to_list   = os.path.dirname(workflow.basedir) + '/scripts/add_sample_to_list.py',
        tmpdir               = config['general']['tempdir']
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" | python3 {params.add_sample_to_list} {wildcards.sample} | '
        'sort -T {params.tmpdir} | uniq -c > {output.overlap} && [[ -s {output.overlap} ]] ; '

rule nbiotech_barcode_metrics_all:
  input:
     bam       = 'results/nbiotech_data/cellranger/{sample}/outs/possorted_bam.bam',
  output:
    all_bcd    = 'results/nbiotech_data/{sample}/barcode_metrics/all_barcodes.txt'
  params:
      get_cell_barcode   = os.path.dirname(workflow.basedir) + '/scripts/get_cell_barcode.awk',
      add_sample_to_list = os.path.dirname(workflow.basedir) + '/scripts/add_sample_to_list.py',
      tmpdir             = config['general']['tempdir']
  shell:
    ' samtools view -f2 {input.bam}| '
    'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" | python3 {params.add_sample_to_list} {wildcards.sample} | '
    'sort -T {params.tmpdir} | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] ; '

####### CELLS SELECTION
rule nbiotech_cell_selection:
    input:
        bcd_all   = 'results/nbiotech_data/{sample}/barcode_metrics/all_barcodes.txt',
        bcd_peak  = 'results/nbiotech_data/{sample}/barcode_metrics/peaks_barcodes.txt',
        peaks     = 'results/nbiotech_data/{sample}/peaks/macs_broad/{sample}_peaks.broadPeak',
        metadata  = 'results/nbiotech_data/cellranger/{sample}/outs/singlecell.csv',
        fragments = 'results/nbiotech_data/{sample}/fragments/fragments.tsv.gz',
    output:
        'results/nbiotech_data/{sample}/cell_picking/cells_10x.png',
        'results/nbiotech_data/{sample}/cell_picking/cells_picked.png',
        'results/nbiotech_data/{sample}/cell_picking/cells_picked.bw',
        'results/nbiotech_data/{sample}/cell_picking/cells_not_picked.bw',
        'results/nbiotech_data/{sample}/cell_picking/metadata.csv',
    params:
        script      = os.path.dirname(workflow.basedir) + '/scripts/pick_cells.R',
        out_prefix  = 'results/nbiotech_data/{sample}/cell_picking/',
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --sample {wildcards.sample} --sample {wildcards.sample} --out_prefix {params.out_prefix}"

rule nbiotech_create_seurat_object:
    input:
        fragments = 'results/nbiotech_data/{sample}/fragments/fragments.tsv.gz',
        peaks     = 'results/nbiotech_data/{sample}/peaks/macs_broad/{sample}_peaks.broadPeak',
        metadata  = 'results/nbiotech_data/{sample}/cell_picking/metadata.csv',
        script    = os.path.dirname(workflow.basedir) + '/scripts/create_seurat_object.R',
    output:
        'results/nbiotech_data/{sample}/seurat/bin_{binwidth}/Seurat_object.Rds',
    params:
        out_prefix  = 'results/nbiotech_data/{sample}/seurat/bin_{binwidth}/',
        genome      = config['general']['genome'],
        antibody    = lambda wildcards: wildcards.sample.split("_")[0]
    shell:
        "Rscript {input.script} --sample {wildcards.sample} --antibody {params.antibody} --metadata {input.metadata}  --fragments {input.fragments} --peaks {input.peaks} --out_prefix {params.out_prefix} --window {wildcards.binwidth} --genome_version {params.genome}"

rule download_seurat_objects:
    output:
        'results/nbiotech_data/data/seurat/GSE163532_RAW.tar'
    params:
        url = config['nbiotech_data']['url']['seurat']
    shell:
        'wget -O {output} {params.url}'

rule untar_seurat:
    input:
        archive  = 'results/nbiotech_data/data/seurat/GSE157637_Seurat_v3_object.tar.gz'
    output:
        H3K27me3 = 'results/nbiotech_data/data/seurat/H3K27me3_seurat_object.Rds'
    params:
        outdir   = 'results/nbiotech_data/data/seurat/'
    shell:
        'tar -xvzf {input.archive}; '

