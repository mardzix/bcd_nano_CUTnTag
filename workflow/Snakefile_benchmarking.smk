include: "Snakefile_preprocess.smk"
include: "Snakefile_pre_nbiotech.smk"
include: "Snakefile_single_modality.smk"


rule all_benchmarks:
    input:
        # PCA together with nbiotech data
        expand('results/benchmarks/{feature}/PCA_on_markers/{idents}/{m}_matrix.npz',feature='peaks',idents='idents_short',m=['markers', 'markers_positive']),

        # metagene H3K27me3/H3K27ac comparison with nbiotech
        'results/benchmarks/peaks/specificity_benchmark/peaks_merged.bed',
        'results/benchmarks/peaks/specificity_benchmark/peaks_matrix.npz',
        'results/benchmarks/peaks/specificity_benchmark/peaks_specific/',
        'results/benchmarks/peaks/specificity_benchmark/metaplot.mtx',
        'results/benchmarks/peaks/specificity_benchmark/metaplot.pdf',

rule markers_to_bed:
    input:
        csv = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/{m}.csv',
        script = workflow_dir + '/scripts/markers_to_bed.R'
    output:
        bed = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/{m}.bed',
    params:
        nmarkers = 50
    shell:
        'Rscript {input.script} --input {input.csv} --nmarkers {params.nmarkers} --output {output.bed}'

rule merge_markers_bed:
    input:
        markers = lambda wildcards: expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/{m}.bed',modality =
            ['ATAC', 'H3K27ac', 'H3K27me3'],feature=wildcards.feature,idents=wildcards.idents,m=wildcards.m),
    output:
        'results/benchmarks/{feature}/PCA_on_markers/{idents}/{m}_merged.bed'
    shell:
        'cat {input.markers} | sort -k1,1 -k2,2n | bedtools merge -i - > {output}'

rule matrix_PCA:
    input:
        markers='results/benchmarks/{feature}/PCA_on_markers/{idents}/{m}_merged.bed',
        bigwig_ATAC='results/multimodal_data/single_modality/ATAC/seurat/{feature}/bam_per_cluster/{idents}/bigwig/',
        bigwig_K27ac='results/multimodal_data/single_modality/H3K27ac/seurat/{feature}/bam_per_cluster/{idents}/bigwig/',
        bigwig_K27me3='results/multimodal_data/single_modality/H3K27me3/seurat/{feature}/bam_per_cluster/{idents}/bigwig/',
        bigwig_nbiot='results/nbiotech_data/data/bigwig/',
    output:
        'results/benchmarks/{feature}/PCA_on_markers/{idents}/{m}_matrix.npz'
    threads: 8
    shell:
        'multiBigwigSummary BED-file -b {input.bigwig_ATAC}/*.bw {input.bigwig_K27ac}/*.bw {input.bigwig_K27me3}/*.bw '
        '{input.bigwig_nbiot}/H3K27ac*.bw {input.bigwig_nbiot}/H3K27me3*.bw -o {output} --BED {input.markers} -p {threads}'


rule merge_all_peaks:
    input:
        'results/multimodal_data/single_modality/H3K27ac/{feature}/macs_broad/H3K27ac_peaks.broadPeak',
        'results/multimodal_data/single_modality/H3K27me3/{feature}/macs_broad/H3K27me3_peaks.broadPeak'
    output:
        'results/benchmarks/{feature}/specificity_benchmark/peaks_merged.bed'
    shell:
        'cat {input} | sort -k1,1 -k2,2n | bedtools merge -i - > {output}'

rule create_matrix_H3K27ac_H3K27me3_metagene:
    input:
        bed          = 'results/benchmarks/{feature}/specificity_benchmark/peaks_merged.bed',
        nbiot_K27ac  = 'results/nbiotech_data/data/bigwig/H3K27ac_Astrocytes.bw',
        nbiot_K27me3 = 'results/nbiotech_data/data/bigwig/H3K27me3_Astrocytes.bw',
        bcd_K27ac    = 'results/multimodal_data/single_modality/H3K27ac/seurat/peaks/bam_per_cluster/idents_short/bigwig/',
        bcd_K27me3   = 'results/multimodal_data/single_modality/H3K27me3/seurat/peaks/bam_per_cluster/idents_short/bigwig/'
    params:
        pattern = 'AST_[TN][ET].bw'
    output:
        np_mtx  = 'results/benchmarks/{feature}/specificity_benchmark/peaks_matrix.npz',
        txt_mtx = 'results/benchmarks/{feature}/specificity_benchmark/peaks_matrix.txt',
    threads: 16
    shell:
        'multiBigwigSummary BED-file -b {input.nbiot_K27ac} {input.nbiot_K27me3} {input.bcd_K27ac}/{params.pattern} {input.bcd_K27me3}/{params.pattern} '
        '--BED {input.bed} -p {threads} -o {output.np_mtx} --outRawCounts	{output.txt_mtx}'

rule filter_peaks_by_matrix:
    input:
        matrix_txt = 'results/benchmarks/{feature}/specificity_benchmark/peaks_matrix.txt',
        matrix_npz= 'results/benchmarks/{feature}/specificity_benchmark/peaks_matrix.npz',
        script = workflow_dir + '/scripts/filter_peaks_by_matrix.R'
    output:
        directory('results/benchmarks/{feature}/specificity_benchmark/peaks_specific/')
    shell:
        'Rscript {input.script} --matrix_txt {input.matrix_txt} --matrix_npz {input.matrix_npz} --output {output}'

rule create_meta_matrix:
    input:
        peaks     = 'results/benchmarks/{feature}/specificity_benchmark/peaks_specific/',
        bw_K27ac  = 'results/multimodal_data/single_modality/H3K27ac/seurat/{feature}/bam_per_cluster/idents_short/bigwig/',
        bw_K27me3 = 'results/multimodal_data/single_modality/H3K27me3/seurat/{feature}/bam_per_cluster/idents_short/bigwig/'
    output:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot.mtx',
    params:
        pattern = 'AST_[TN][ET].bw'
    threads: 8
    shell:
        'computeMatrix reference-point -S {input.bw_K27ac}/{params.pattern} {input.bw_K27me3}/{params.pattern} '
        ' -R {input.peaks}/*.bed -o {output} -p {threads} --referencePoint center '
        '--upstream 20000 --downstream 20000'

rule plot_heatmap:
    input:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot.mtx',
    output:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot.pdf',
    shell:
       'plotHeatmap -m {input} -o {output}'
