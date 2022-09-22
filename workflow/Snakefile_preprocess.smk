include: 'Snakefile_prep.smk'

rule all_preprocess:
    input:
        cellranger          = ['results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        bigwig_all          = ['results/multimodal_data/{sample}/{antibody}_{barcode}/bigwig/all_reads.bw'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        macs_narrow         = ['results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_narrow/{antibody}_peaks.narrowPeak'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        macs_broad          = ['results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_broad/{antibody}_peaks.broadPeak'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        chrom_sizes         = 'results/mm10.chrom.sizes',
        SEACR_peaks         = ['results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/SEACR/peaks.relaxed.bed'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        fragments           = ['results/multimodal_data/{sample}/{antibody}_{barcode}/fragments/fragments.tsv.gz'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        peaks_overlap       = ['results/multimodal_data/{sample}/{antibody}_{barcode}/barcode_metrics/peaks_barcodes.txt'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        barcodes_sum        = ['results/multimodal_data/{sample}/{antibody}_{barcode}/barcode_metrics/all_barcodes.txt'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        cell_pick           = ['results/multimodal_data/{sample}/{antibody}_{barcode}/cell_picking/metadata.csv'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        seurat              = ['results/multimodal_data/{sample}/{antibody}_{barcode}/seurat/bin_{binwidth}/Seurat_object.Rds'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody],binwidth = binwidth) for sample in samples_list for antibody in barcodes_dict[sample].keys() for binwidth in config['general']['binwidth']],
        seurat_peaks        = ['results/multimodal_data/{sample}/{antibody}_{barcode}/seurat/peaks/Seurat_object.Rds'.format(sample = sample, antibody = antibody, barcode = barcodes_dict[sample][antibody]) for sample in samples_list for antibody in barcodes_dict[sample].keys()],
        fragments_merged    = expand('results/multimodal_data/single_modality/{modality}/fragments/fragments.tsv.gz', modality= antibodies_list),# Merged fragments file
        fragments_biwgiw    = expand('results/multimodal_data/single_modality/{modality}/fragments/fragments.bw',modality=antibodies_list),
        peaks_merged        = expand('results/multimodal_data/single_modality/{modality}/peaks/macs_broad/{modality}_peaks.broadPeak',modality=antibodies_list),# Merged peaks file
        seurat_object       = expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object.Rds',modality=antibodies_list,feature=features),# Seurat object peaks
        seurat_clustered    = expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds',modality=antibodies_list,feature=features),

def get_peaks_file_from_modality(sample,modality,barcode):
    if modality == "ATAC":
        # Cellranger default peak calling
        # peaks = "results/multimodal_data/{sample}/cellranger/{sample}_{modality}_{barcode}/outs/peaks.bed".format(sample=sample, modality=modality, barcode=barcode)
        peaks = 'results/multimodal_data/{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample,modality=modality,barcode=barcode)
    elif modality == "H3K27ac":
        # Broad peaks
        peaks = 'results/multimodal_data/{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample, modality=modality, barcode=barcode)
    elif modality == "H3K27me3":
        # Broad peaks
        peaks = 'results/multimodal_data/{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample, modality=modality, barcode=barcode)
    else:
        peaks = 'results/multimodal_data/{sample}/{modality}_{barcode}/peaks/macs_broad/{modality}_peaks.broadPeak'.format(sample=sample, modality=modality, barcode=barcode)
    return peaks

def get_fragments_per_modality(modality, barcodes_dict):
    result = []
    for s in barcodes_dict:
        for m in barcodes_dict[s]:
            if modality == m:
               result.append('results/multimodal_data/{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz'.format(sample = s, modality = m, barcode=barcodes_dict[s][m]))
    return result
def get_seurat_per_modality(modality, barcodes_dict,feature):
    result = []
    for s in barcodes_dict:
        for m in barcodes_dict[s]:
            if modality == m:
               result.append('results/multimodal_data/{sample}/{modality}_{barcode}/seurat/{feature}/Seurat_object.Rds'.format(sample = s, modality = m, barcode=barcodes_dict[s][m], feature = feature))
    return result

rule demultiplex:
    input:
        script = workflow_dir + '/scripts/debarcode.py',
        fastq  = lambda wildcards: glob.glob(config['samples'][wildcards.sample]['fastq_path'] + '/**/*{lane}*R[123]*.fastq.gz'.format(lane=wildcards.lane),recursive=True)
    output:
        'results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_{barcode}/barcode_{barcode}/{id}_{number}_{lane}_R1_{suffix}',
        'results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_{barcode}/barcode_{barcode}/{id}_{number}_{lane}_R2_{suffix}',
        'results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_{barcode}/barcode_{barcode}/{id}_{number}_{lane}_R3_{suffix}',
    params:
        nbarcodes  = lambda wildcards: len(config['samples'][wildcards.sample]['barcodes']),
        out_folder = lambda wildcards: 'results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_{barcode}/'.format(sample=wildcards.sample, antibody=wildcards.antibody,barcode=wildcards.barcode),
    shell:
        "python3 {input.script} -i {input.fastq} -o {params.out_folder} --single_cell --barcode {wildcards.barcode} 2>&1"

rule run_cellranger:
    input:
        # Uses non-demultiplexed fastq files to figure out the seq id
        lambda wildcards: get_fastq_for_cellranger(config['samples'][wildcards.sample]['fastq_path'],sample=wildcards.sample,antibody=wildcards.antibody,barcode=wildcards.barcode)
    output:
        bam   = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam',
        frag  = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/fragments.tsv.gz',
        meta  = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/singlecell.csv',
        peaks = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/peaks.bed',
    params:
        cellranger_ref = config['general']['cellranger_ref'],
        fastq_folder   = lambda wildcards: os.getcwd() + '/results/multimodal_data/{sample}/fastq_per_barcode/{antibody}_{barcode}/barcode_{barcode}/'.format(sample=wildcards.sample, antibody=wildcards.antibody, barcode=wildcards.barcode)
    threads: 40
    shell:
        'rm -rf results/multimodal_data/{wildcards.sample}/cellranger/{wildcards.sample}_{wildcards.antibody}_{wildcards.barcode}/; '
        'cd results/multimodal_data/{wildcards.sample}/cellranger/; '
        '/data/bin/cellranger-atac count --id {wildcards.sample}_{wildcards.antibody}_{wildcards.barcode} --reference {params.cellranger_ref} --fastqs {params.fastq_folder}'

rule bam_to_bw:
    input:
        cellranger_bam = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam'
    output:
        bigwig         = 'results/multimodal_data/{sample}/{antibody}_{barcode}/bigwig/all_reads.bw'
    threads: 16
    shell:
        'bamCoverage -b {input.cellranger_bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule run_macs_narrow:
    input:
        cellranger_bam = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam'
    output:
        narrow_peaks = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_narrow/{antibody}_peaks.narrowPeak'
    params:
        macs_outdir  = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_narrow/'
    shell:
        'macs2 callpeak -t {input.cellranger_bam} -g mm -f BAMPE -n {wildcards.antibody} '
        '--outdir {params.macs_outdir} --keep-dup=1 --llocal 100000 --cutoff-analysis --min-length 1000 --max-gap 1000  2>&1 '

rule run_macs_narrow_test:
    input:
        cellranger_bam = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam'
    output:
        narrow_peaks = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_narrow_llocal_{llocal}/{antibody}_peaks.narrowPeak'
    params:
        macs_outdir  = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_narrow_llocal_{llocal}/'
    shell:
        'macs2 callpeak -t {input.cellranger_bam} -g mm -f BAMPE -n {wildcards.antibody} '
        '--outdir {params.macs_outdir} --keep-dup=1 --llocal {wildcards.llocal} 2>&1 '

rule run_macs_broad:
    input:
        cellranger_bam = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam'
    output:
        broad_peaks = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_broad/{antibody}_peaks.broadPeak'
    params:
        macs_outdir = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/macs_broad/'
    shell:
        'macs2 callpeak -t {input} -g mm -f BAMPE -n {wildcards.antibody} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup=1 --broad-cutoff=0.1 ' 
        '--min-length 1000 --max-gap 1000 --broad 2>&1 '

rule get_chromsizes:
    output:
        'results/mm10.chrom.sizes'
    params:
        faidx = config['general']['cellranger_ref'] + '/fasta/genome.fa.fai' # Assuming fasta index is in /fasta/genome.fa.fai, change if neccessary
    shell:
        'cut -f1,2 {params.faidx} > {output}'

rule prep_SEACR_files:
    input:
        fragments = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/fragments.tsv.gz',
        genome = 'results/mm10.chrom.sizes'
    output:
        'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/SEACR/fragments.bg'
    shell:
        "bedtools genomecov -bg -g {input.genome} -i {input.fragments} > {output}"

rule run_SEACR:
    input:
        'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/SEACR/fragments.bg',
    output:
        'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/SEACR/peaks.relaxed.bed',
    params:
        out_prefix = 'results/multimodal_data/{sample}/{antibody}_{barcode}/peaks/SEACR/peaks',
    shell:
        "~/bin/SEACR/SEACR_1.3.sh {input} 0.01 norm relaxed {params.out_prefix}"

rule add_barcode_fragments:
    input:
        fragments = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/fragments.tsv.gz',
    output:
        fragments = 'results/multimodal_data/{sample}/{antibody}_{barcode}/fragments/fragments.tsv.gz',
        index     = 'results/multimodal_data/{sample}/{antibody}_{barcode}/fragments/fragments.tsv.gz.tbi',
    params:
        script    = workflow_dir + '/scripts/add_sample_to_fragments.py',
    shell:
        'python3 {params.script} {input.fragments} {wildcards.sample} | bgzip > {output.fragments}; '
        'tabix -p bed {output.fragments}'

rule barcode_overlap_peaks:
    input:
        bam    = 'results/multimodal_data/{sample}/cellranger/{sample}_{modality}_{barcode}/outs/possorted_bam.bam',
        peaks  = lambda wildcards: get_peaks_file_from_modality(wildcards.sample,wildcards.modality,wildcards.barcode),
    output:
        overlap = 'results/multimodal_data/{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt'
    params:
        get_cell_barcode     = workflow_dir + '/scripts/get_cell_barcode.awk',
        add_sample_to_list   = workflow_dir + '/scripts/add_sample_to_list.py',
        tmpdir               = config['general']['tempdir']
    shell:
        'bedtools intersect -abam {input.bam} -b {input.peaks} -u | samtools view -f2 | '
        'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" | python3 {params.add_sample_to_list} {wildcards.sample} | '
        'sort -T {params.tmpdir} | uniq -c > {output.overlap} && [[ -s {output.overlap} ]] ; '

rule barcode_metrics_all:
  input:
     bam       = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam',
  output:
    all_bcd    = 'results/multimodal_data/{sample}/{antibody}_{barcode}/barcode_metrics/all_barcodes.txt'
  params:
      get_cell_barcode   = workflow_dir + '/scripts/get_cell_barcode.awk',
      add_sample_to_list = workflow_dir + '/scripts/add_sample_to_list.py',
      tmpdir             = config['general']['tempdir']
  shell:
    ' samtools view -f2 {input.bam}| '
    'awk -f {params.get_cell_barcode} | sed "s/CB:Z://g" | python3 {params.add_sample_to_list} {wildcards.sample} | '
    'sort -T {params.tmpdir} | uniq -c > {output.all_bcd} && [[ -s {output.all_bcd} ]] ; '

####### CELLS SELECTION
rule cell_selection:
    input:
        bcd_all   = 'results/multimodal_data/{sample}/{modality}_{barcode}/barcode_metrics/all_barcodes.txt',
        bcd_peak  = 'results/multimodal_data/{sample}/{modality}_{barcode}/barcode_metrics/peaks_barcodes.txt',
        peaks     = lambda wildcards: get_peaks_file_from_modality(wildcards.sample, wildcards.modality, wildcards.barcode),
        metadata  = 'results/multimodal_data/{sample}/cellranger/{sample}_{modality}_{barcode}/outs/singlecell.csv',
        fragments = 'results/multimodal_data/{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz',
    output:
        'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/cells_10x.png',
        'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/cells_picked.png',
        'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/cells_picked.bw',
        'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/cells_not_picked.bw',
        'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
    params:
        script      = workflow_dir + '/scripts/pick_cells.R',
        out_prefix  = 'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/',
    shell:
        "Rscript {params.script} --metadata {input.metadata} --fragments {input.fragments} --bcd_all {input.bcd_all} --bcd_peak {input.bcd_peak} --antibody {wildcards.modality} --sample {wildcards.sample} --out_prefix {params.out_prefix}"

rule create_seurat_object_bins:
    input:
        fragments = 'results/multimodal_data/{sample}/{antibody}_{barcode}/fragments/fragments.tsv.gz',
        metadata  = 'results/multimodal_data/{sample}/{antibody}_{barcode}/cell_picking/metadata.csv',
        script    = workflow_dir + '/scripts/create_seurat_object.R',
    output:
        'results/multimodal_data/{sample}/{antibody}_{barcode}/seurat/bin_{binwidth}/Seurat_object.Rds',
    params:
        out_prefix  = 'results/multimodal_data/{sample}/{antibody}_{barcode}/seurat/bin_{binwidth}/',
        genome      = config['general']['genome'],
    shell:
        "Rscript {input.script} --sample {wildcards.sample}   --antibody {wildcards.antibody} --metadata {input.metadata} --fragments {input.fragments} --out_prefix {params.out_prefix} --window {wildcards.binwidth} --genome_version {params.genome}"


rule create_seurat_object_peaks:
    input:
        fragments = 'results/multimodal_data/{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz',
        metadata  = 'results/multimodal_data/{sample}/{modality}_{barcode}/cell_picking/metadata.csv',
        peaks     = lambda wildcards: get_peaks_file_from_modality(wildcards.sample, wildcards.modality, wildcards.barcode),
        script    = workflow_dir + '/scripts/create_seurat_object.R',
    output:
        'results/multimodal_data/{sample}/{modality}_{barcode}/seurat/peaks/Seurat_object.Rds',
    params:
        out_prefix  = 'results/multimodal_data/{sample}/{modality}_{barcode}/seurat/peaks/',
        genome      = config['general']['genome'],
    shell:
        "Rscript {input.script} --sample {wildcards.sample}   --antibody {wildcards.modality} --metadata {input.metadata} --fragments {input.fragments} " \ 
        " --peaks {input.peaks} --out_prefix {params.out_prefix} --genome_version {params.genome}"


rule merge_fragments_file:
    input:
        fragments = lambda wildcards: get_fragments_per_modality(modality = wildcards.modality, barcodes_dict = barcodes_dict)
    output:
        fragments_merged = 'results/multimodal_data/single_modality/{modality}/fragments/fragments.tsv.gz',
        index            = 'results/multimodal_data/single_modality/{modality}/fragments/fragments.tsv.gz.tbi'
    params:
        tmpdir           = config['general']['tempdir']
    shell:
        'zcat {input.fragments} | sort -T {params.tmpdir} -k1,1 -k2,2n | bgzip > {output.fragments_merged} && tabix -p bed {output.fragments_merged}'

rule fragments_to_bw:
    input:
        fragments   = 'results/multimodal_data/single_modality/{modality}/fragments/fragments.tsv.gz',
        chrom_sizes = 'results/mm10.chrom.sizes'
    output:
        bam        = temp('results/multimodal_data/single_modality/{modality}/fragments/fragments.bam'),
        bam_sorted = temp('results/multimodal_data/single_modality/{modality}/fragments/fragments_sorted.bam'),
        index      = temp('results/multimodal_data/single_modality/{modality}/fragments/fragments_sorted.bam.bai'),
        bigwig     = 'results/multimodal_data/single_modality/{modality}/fragments/fragments.bw',
    threads: 8
    shell:
        "bedToBam -i {input.fragments} -g {input.chrom_sizes} > {output.bam} && "
        "samtools sort -@ {threads} -o {output.bam_sorted} {output.bam} &&"
        "samtools index {output.bam_sorted} && "
        "bamCoverage -b {output.bam_sorted} -o {output.bigwig} -p {threads} --minMappingQuality 5 --binSize 50 --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates"

rule run_macs_broad_merged:
    input:
        fragments = 'results/multimodal_data/single_modality/{modality}/fragments/fragments.tsv.gz'
    output:
        broad_peaks = 'results/multimodal_data/single_modality/{modality}/peaks/macs_broad/{modality}_peaks.broadPeak'
    params:
        macs_outdir = 'results/multimodal_data/single_modality/{modality}/peaks/macs_broad/'
    shell:
        'macs2 callpeak -t {input} -g mm -f BED -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup=1 --broad-cutoff=0.1 ' 
        '--min-length 1000 --max-gap 1000 --broad --nomodel 2>&1 '


rule merge_seurat_single_modality:
    input:
        seurat = lambda wildcards: get_seurat_per_modality(modality=wildcards.modality,barcodes_dict=barcodes_dict, feature=wildcards.feature),
        script = workflow_dir + '/scripts/merge_objects.R' # TODO uncomment
    output:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object.Rds'
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat}'


rule cluster:
    input:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object.Rds',
        # script = workflow_dir + '/scripts/UMAP_cluster.R' # TODO - temporarily commented out to prevent reruning big part of the pipeline
    output:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds'
    params:
        script = workflow_dir + '/scripts/UMAP_cluster.R',
        plot_group = 'sample'
    shell:
        'Rscript {params.script} -i {input.seurat} -o {output.seurat} -a {wildcards.feature} -d 40 -g {params.plot_group} '



