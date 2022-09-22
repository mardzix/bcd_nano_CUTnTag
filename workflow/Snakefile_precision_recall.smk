include: 'Snakefile_pre_nbiotech.smk'
include: 'Snakefile_single_modality.smk'
include: 'Snakefile_benchmarking.smk'


rule precision_recall_all:
    input:
        peaks        = expand('results/precision_recall/scCT/peaks/macs2_{q}/H3K27me3_peaks.broadPeak', q = [0.1,0.01,0.001]),
        peaks_nanoCT = expand('results/precision_recall/nanoCT/peaks/macs2_{q}/H3K27me3_peaks.broadPeak',q = [0.1,0.01,0.001]),
        bigwig       = expand('results/precision_recall/{s}/bigwig/merged_bam.bw', s =['nanoCT','scCT']),
        #
        meta3_nanoCT = 'results/precision_recall/3_techniques_metagene/nanoCT_metagene_matrix.txt.gz',
        meta3_scCT   = 'results/precision_recall/3_techniques_metagene/scCT_metagene_matrix.txt.gz',
        meta3_encode = 'results/precision_recall/3_techniques_metagene/encode_metagene_matrix.txt.gz',
        meta3_encode_plot = 'results/precision_recall/3_techniques_metagene/encode_metagene_matrix.png',


rule call_peaks_nbiotech:
    input:
        fragments = ['results/nbiotech_data/cellranger/H3K27me3_N1/outs/possorted_bam.bam',
                     'results/nbiotech_data/cellranger/H3K27me3_N3/outs/possorted_bam.bam',
                     'results/nbiotech_data/cellranger/H3K27me3_N2/outs/possorted_bam.bam',
                     'results/nbiotech_data/cellranger/H3K27me3_N4/outs/possorted_bam.bam']
    output:
        peaks = 'results/precision_recall/scCT/peaks/macs2_{q}/H3K27me3_peaks.broadPeak'
    params:
        macs_outdir = lambda wildcards: 'results/precision_recall/scCT/peaks/macs2_{q}/'.format(q = wildcards.q)
    shell:
        'macs2 callpeak -t {input.fragments} -g mm -f BAM -n H3K27me3 '
        '--outdir {params.macs_outdir} --llocal 1000000 --keep-dup=1 --broad-cutoff={wildcards.q} ' 
        '--qvalue {wildcards.q} --broad --nomodel --min-length 1000 --max-gap 1000  2>&1'

rule call_peaks_nanoCT:
    input:
        fragments = ['results/multimodal_data/bcdCT_MB21_02/cellranger/bcdCT_MB21_02_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
                     'results/multimodal_data/bcdCT_MB21_03/cellranger/bcdCT_MB21_03_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
                     'results/multimodal_data/bcdCT_MB21_04/cellranger/bcdCT_MB21_04_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
                     'results/multimodal_data/bcdCT_MB21_05/cellranger/bcdCT_MB21_05_H3K27me3_CCTATCCT/outs/possorted_bam.bam']
    output:
        peaks     = 'results/precision_recall/nanoCT/peaks/macs2_{q}/H3K27me3_peaks.broadPeak',
    params:
        macs_outdir = lambda wildcards: 'results/precision_recall/nanoCT/peaks/macs2_{q}/'.format(q = wildcards.q)
    shell:
        'macs2 callpeak -t {input.fragments} -g mm -f BAM -n H3K27me3 '
        '--outdir {params.macs_outdir} --llocal 1000000 --keep-dup=1 --broad-cutoff={wildcards.q} ' 
        '--qvalue {wildcards.q} --broad --nomodel --min-length 1000 --max-gap 1000  2>&1  '

rule merge_bam_nbiotech:
    input:
        fragments = ['results/nbiotech_data/cellranger/H3K27me3_N1/outs/possorted_bam.bam',
                     'results/nbiotech_data/cellranger/H3K27me3_N3/outs/possorted_bam.bam',
                     'results/nbiotech_data/cellranger/H3K27me3_N2/outs/possorted_bam.bam',
                     'results/nbiotech_data/cellranger/H3K27me3_N4/outs/possorted_bam.bam']
    output:
        bam = 'results/precision_recall/scCT/bigwig/merged_bam.bam'
    shell:
        'samtools merge {input.fragments} -o {output.bam}; '
        'samtools index {output.bam}'

rule merge_bam_nanoCT:
    input:
        fragments = ['results/multimodal_data/bcdCT_MB21_02/cellranger/bcdCT_MB21_02_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
                     'results/multimodal_data/bcdCT_MB21_03/cellranger/bcdCT_MB21_03_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
                     'results/multimodal_data/bcdCT_MB21_04/cellranger/bcdCT_MB21_04_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
                     'results/multimodal_data/bcdCT_MB21_05/cellranger/bcdCT_MB21_05_H3K27me3_CCTATCCT/outs/possorted_bam.bam']
    output:
        bam = 'results/precision_recall/nanoCT/bigwig/merged_bam.bam'
    shell:
        'samtools merge {input.fragments} -o {output.bam}; '
        'samtools index {output.bam}'

rule bam_to_bw_2:
    input:
        bam = 'results/precision_recall/{s}/bigwig/merged_bam.bam'
    output:
        bigwig = 'results/precision_recall/{s}/bigwig/merged_bam.bw'
    threads: 16
    shell:
        'bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --minMappingQuality 5 '
        ' --binSize 50 --centerReads --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates --extendReads'

rule download_blacklist:
    output:
        bed_gz = 'results/mm10-blacklist.v2.bed.gz',
        bed = 'results/mm10-blacklist.v2.bed',
    params:
        url = 'https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz?raw=true'
    shell:
        'wget -O {output.bed_gz} {params.url}; '
        'gunzip -c {output.bed_gz} > {output.bed}; '

rule metagene_compare_3_methods_a:
    input:
        nanoCT_bw = 'results/precision_recall/nanoCT/bigwig/merged_bam.bw',
        scCT_bw   = 'results/nbiotech_data/merged/bw/nbiotech_merged.bw',
        encode_bw = 'results/encode/bigwig/encode_merged.bw',
        nano_CT_peaks = 'results/precision_recall/nanoCT/peaks/macs2_0.01/H3K27me3_peaks.broadPeak',
        scCT_peaks    = 'results/precision_recall/scCT/peaks/macs2_0.01/H3K27me3_peaks.broadPeak',
        encode_peaks  = 'results/encode/peaks/peaks_all_merged.bed',
        blacklist     = 'results/mm10-blacklist.v2.bed.gz'
    output:
        matrix         = 'results/precision_recall/3_techniques_metagene/nanoCT_metagene_matrix.txt.gz',
        sorted_regions = 'results/precision_recall/3_techniques_metagene/nanoCT_sorted_regions.txt',
    threads: 16
    params:
        blacklist = 'results/precision_recall/3_techniques_metagene/mm10-blacklist.v2.bed.gz'
    shell:
        'computeMatrix reference-point '
        '-S {input.nanoCT_bw} '
        '-R {input.scCT_peaks} {input.encode_peaks} '
        '-o {output.matrix} '
        '-p {threads} '
        '--upstream 10000 --downstream 10000 --referencePoint center -bl {input.blacklist} '
        '--outFileSortedRegions {output.sorted_regions} --sortRegions descend'

rule metagene_compare_3_methods_b:
    input:
        nanoCT_bw = 'results/precision_recall/nanoCT/bigwig/merged_bam.bw',
        scCT_bw   = 'results/nbiotech_data/merged/bw/nbiotech_merged.bw',
        encode_bw = 'results/encode/bigwig/encode_merged.bw',
        nano_CT_peaks = 'results/precision_recall/nanoCT/peaks/macs2_0.01/H3K27me3_peaks.broadPeak',
        scCT_peaks    = 'results/precision_recall/scCT/peaks/macs2_0.01/H3K27me3_peaks.broadPeak',
        encode_peaks  = 'results/encode/peaks/peaks_all_merged.bed',
        blacklist= 'results/mm10-blacklist.v2.bed.gz'
    output:
        matrix         = 'results/precision_recall/3_techniques_metagene/scCT_metagene_matrix.txt.gz',
        sorted_regions = 'results/precision_recall/3_techniques_metagene/scCT_sorted_regions.txt',
    threads: 16
    params:
        blacklist = 'results/precision_recall/3_techniques_metagene/mm10-blacklist.v2.bed.gz'
    shell:
        'computeMatrix reference-point '
        '-S {input.scCT_bw} '
        '-R {input.scCT_peaks} {input.encode_peaks} '
        '-o {output.matrix} '
        '-p {threads} '
        '--upstream 10000 --downstream 10000 --referencePoint center -bl {input.blacklist} '
        '--outFileSortedRegions {output.sorted_regions} --sortRegions descend'

rule metagene_compare_3_methods_c:
    input:
        nanoCT_bw = 'results/precision_recall/nanoCT/bigwig/merged_bam.bw',
        scCT_bw   = 'results/nbiotech_data/merged/bw/nbiotech_merged.bw',
        encode_bw = 'results/encode/bigwig/encode_merged.bw',
        nano_CT_peaks = 'results/precision_recall/nanoCT/peaks/macs2_0.01/H3K27me3_peaks.broadPeak',
        scCT_peaks    = 'results/precision_recall/scCT/peaks/macs2_0.01/H3K27me3_peaks.broadPeak',
        encode_peaks  = 'results/encode/peaks/peaks_all_merged.bed',
        blacklist= 'results/mm10-blacklist.v2.bed.gz'
    output:
        matrix         = 'results/precision_recall/3_techniques_metagene/encode_metagene_matrix.txt.gz',
        sorted_regions = 'results/precision_recall/3_techniques_metagene/encode_sorted_regions.txt',
    threads: 16
    params:
        blacklist = 'results/precision_recall/3_techniques_metagene/mm10-blacklist.v2.bed.gz'
    shell:
        'computeMatrix reference-point '
        '-S {input.encode_bw} '
        '-R {input.scCT_peaks} {input.encode_peaks} '
        '-o {output.matrix} '
        '-p {threads} '
        '--upstream 10000 --downstream 10000 --referencePoint center -bl {input.blacklist} '
        '--outFileSortedRegions {output.sorted_regions} --sortRegions descend'

rule plot_3_metagene:
    input:
        encode = 'results/precision_recall/3_techniques_metagene/encode_metagene_matrix.txt.gz',
        scCT   = 'results/precision_recall/3_techniques_metagene/scCT_metagene_matrix.txt.gz',
        nanoCT = 'results/precision_recall/3_techniques_metagene/nanoCT_metagene_matrix.txt.gz',
    output:
        encode = 'results/precision_recall/3_techniques_metagene/encode_metagene_matrix.png',
        scCT   = 'results/precision_recall/3_techniques_metagene/scCT_metagene_matrix.png',
        nanoCT = 'results/precision_recall/3_techniques_metagene/nanoCT_metagene_matrix.png',
        encode_pdf = 'results/precision_recall/3_techniques_metagene/encode_metagene_matrix.pdf',
        scCT_pdf   = 'results/precision_recall/3_techniques_metagene/scCT_metagene_matrix.pdf',
        nanoCT_pdf = 'results/precision_recall/3_techniques_metagene/nanoCT_metagene_matrix.pdf',
    shell:
        'plotHeatmap -m {input.nanoCT} -o {output.nanoCT} --zMax 10 ;  '
        'plotHeatmap -m {input.scCT} -o {output.scCT} --zMax 25 ; '
        'plotHeatmap -m {input.encode} -o {output.encode} --zMax 40; '
        'plotHeatmap -m {input.nanoCT} -o {output.nanoCT_pdf} --zMax 10 ;  '
        'plotHeatmap -m {input.scCT} -o {output.scCT_pdf} --zMax 25 ; '
        'plotHeatmap -m {input.encode} -o {output.encode_pdf} --zMax 40; '



