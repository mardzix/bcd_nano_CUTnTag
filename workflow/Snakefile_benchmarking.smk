include: "Snakefile_preprocess.smk"
include: "Snakefile_pre_nbiotech.smk"
include: "Snakefile_single_modality.smk"


rule all_benchmarks:
    input:
        # PCA together with nbiotech data
        expand('results/benchmarks/{feature}/PCA_on_markers/{idents}/{m}_matrix.npz',feature='peaks',idents='idents_short',m=['markers', 'markers_positive']),

        # metagene H3K27me3/H3K27ac comparison with nbiotech
#         'results/benchmarks/peaks/specificity_benchmark/peaks_merged.bed',
        'results/benchmarks/peaks/specificity_benchmark/peaks_matrix.npz',
        'results/benchmarks/peaks/specificity_benchmark/peaks_specific/',
        'results/benchmarks/peaks/specificity_benchmark/metaplot.pdf',
        'results/benchmarks/peaks/specificity_benchmark/metaplot_2.pdf',
        'results/benchmarks/peaks/specificity_benchmark/fingerpring.pdf',
        # Encode
        'results/encode/download/ENCFF508DLX.bam',
        'results/encode/peaks/peaks_all_merged.bed',
        # Devel
        'results/encode/bigwig/encode_merged.bw'

# rule markers_to_bed_benchmark:
#     input:
#         csv = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/{m}.csv',
#         script = workflow_dir + '/scripts/markers_to_bed.R'
#     output:
#         bed = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/{m}.bed',
#     params:
#         nmarkers = 50
#     shell:
#         'Rscript {input.script} --input {input.csv} --nmarkers {params.nmarkers} --output {output.bed}'

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
        bw_K27me3 = 'results/multimodal_data/single_modality/H3K27me3/seurat/{feature}/bam_per_cluster/idents_short/bigwig/',
        nbiot_K27ac= 'results/nbiotech_data/data/bigwig/H3K27ac_Astrocytes.bw',
        nbiot_K27me3='results/nbiotech_data/data/bigwig/H3K27me3_Astrocytes.bw',
    output:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot.mtx',
    params:
        pattern = 'AST_TE.bw'
    threads: 8
    shell:
        'computeMatrix reference-point -S {input.bw_K27ac}/{params.pattern} {input.bw_K27me3}/{params.pattern} {input.nbiot_K27ac} {input.nbiot_K27me3}'
        ' -R {input.peaks}/*.bed -o {output} -p {threads} --referencePoint center '
        '--upstream 20000 --downstream 20000'

rule plot_heatmap:
    input:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot.mtx',
    output:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot.pdf',
    shell:
       'plotHeatmap -m {input} -o {output} --colorMap Purples Greens Purples Greens --zMax 5.5 6.5 3.5 4.5'

rule create_meta_matrix_2:
    input:
        peaks     = '/data/ref/cellranger-atac/refdata-cellranger-atac-mm10-1.2.0/genes/genes.gtf',
        bw_K27ac  = 'results/multimodal_data/single_modality/H3K27ac/seurat/{feature}/bam_per_cluster/idents_short/bigwig/',
        bw_K27me3 = 'results/multimodal_data/single_modality/H3K27me3/seurat/{feature}/bam_per_cluster/idents_short/bigwig/',
        nbiot_K27ac= 'results/nbiotech_data/data/bigwig/H3K27ac_Astrocytes.bw',
        nbiot_K27me3='results/nbiotech_data/data/bigwig/H3K27me3_Astrocytes.bw',
    output:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot_2.mtx',
    params:
        pattern = 'AST_TE.bw'
    threads: 8
    shell:
        'computeMatrix reference-point -S {input.bw_K27ac}/{params.pattern}  {input.nbiot_K27ac} '
        ' -R {input.peaks} -o {output} -p {threads} '
        ' --upstream 100000 --downstream 200000'

rule plot_heatmap_2:
    input:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot_2.mtx',
    output:
        'results/benchmarks/{feature}/specificity_benchmark/metaplot_2.pdf',
    shell:
       'plotHeatmap -m {input} -o {output} --colorMap Purples Greens'


rule download_encode:
    input:
        script =  os.path.dirname(workflow.basedir) + '/scripts/download_encode.sh'
    output:
        # P0 forebrain
        'results/encode/download/ENCFF508DLX.bam',
        'results/encode/download/ENCFF338UAL.bam',
        'results/encode/download/ENCFF290TRM.bam',
        'results/encode/download/ENCFF852BHI.bam',
        # P0 hindbrain
        'results/encode/download/ENCFF014BZV.bam',
        'results/encode/download/ENCFF298PYM.bam',
        'results/encode/download/ENCFF837IMP.bam',
        'results/encode/download/ENCFF685NLF.bam',
        # P0 midbrain
        'results/encode/download/ENCFF134BLZ.bam',
        'results/encode/download/ENCFF687PXP.bam',
        'results/encode/download/ENCFF072WBY.bam',
        'results/encode/download/ENCFF981PDE.bam',
        # P0 Hindbrain peaks
        'results/encode/download/ENCFF953WTE.bed.gz',
        'results/encode/download/ENCFF301ZAK.bed.gz',
        # P0 Midbrain peaks
        'results/encode/download/ENCFF911MII.bed.gz',
        'results/encode/download/ENCFF973VFO.bed.gz',
        # P0 Forebrain peaks
        'results/encode/download/ENCFF434VST.bed.gz',
        'results/encode/download/ENCFF769WNP.bed.gz',
    params:
        out_folder = 'results/encode/download/'
    shell:
        'cd {params.out_folder}; '
        'sh  {input.script}; '
        'ls *.bam | while read line; do samtools index $line; done'

rule plot_fingerprint:
    input:
        'results/encode/download/ENCFF508DLX.bam',
        'results/encode/download/ENCFF338UAL.bam',
        'results/encode/download/ENCFF290TRM.bam',
        'results/encode/download/ENCFF852BHI.bam',
        'results/nbiotech_data/cellranger/H3K27me3_N1/outs/possorted_bam.bam',
        'results/nbiotech_data/cellranger/H3K27me3_N2/outs/possorted_bam.bam',
        'results/nbiotech_data/cellranger/H3K27me3_N3/outs/possorted_bam.bam',
        'results/nbiotech_data/cellranger/H3K27me3_N4/outs/possorted_bam.bam',
        'results/multimodal_data/bcdCT_MB21_02/cellranger/bcdCT_MB21_02_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
        'results/multimodal_data/bcdCT_MB21_03/cellranger/bcdCT_MB21_03_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
        'results/multimodal_data/bcdCT_MB21_04/cellranger/bcdCT_MB21_04_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
        'results/multimodal_data/bcdCT_MB21_05/cellranger/bcdCT_MB21_05_H3K27me3_CCTATCCT/outs/possorted_bam.bam',
    output:
        counts = 'results/benchmarks/{feature}/specificity_benchmark/fingerpring.txt',
        pdf    = 'results/benchmarks/{feature}/specificity_benchmark/fingerpring.pdf'
    threads: 16
    shell:
        'plotFingerprint --bamfiles {input} --outRawCounts {output.counts} --minMappingQuality 30 -plot {output.pdf}'

rule merge_encode_peaks:
    input:
        # P0 Hindbrain peaks
        'results/encode/download/ENCFF953WTE.bed.gz',
        'results/encode/download/ENCFF301ZAK.bed.gz',
        # P0 Midbrain peaks
        'results/encode/download/ENCFF911MII.bed.gz',
        'results/encode/download/ENCFF973VFO.bed.gz',
        # P0 Forebrain peaks
        'results/encode/download/ENCFF434VST.bed.gz',
        'results/encode/download/ENCFF769WNP.bed.gz',
    output:
        peaks_all       = temp('results/encode/peaks/peaks_all.bed'),
        peaks_all_merge = 'results/encode/peaks/peaks_all_merged.bed'
    shell:
        'ls results/encode/download/*.bed.gz | while read line; do zcat $line >> {output.peaks_all}; done; '
        'cat {output.peaks_all} | sort -k1,1 -k2,2n | bedtools merge -i stdin > {output.peaks_all_merge}; '

rule encode_bam_to_bw:
    input:
        'results/encode/download/'
    output:
        'results/encode/bigwig/encode_merged.bw'
    threads: 16
    shell:
        'samtools merge -@ {threads} -o results/encode/download/encode_merged.bam `ls *.bam`; '
        'samtools index results/encode/download/encode_merged.bam; '
        '  bamCoverage -b results/encode/download/encode_merged.bam -o output -p {threads}; '
        '  done'


