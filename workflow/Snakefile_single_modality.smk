include: "Snakefile_preprocess.smk"


def get_fragments_per_modality(modality, barcodes_dict):
    result = []
    for s in barcodes_dict:
        for m in barcodes_dict[s]:
            if modality == m:
               result.append('results/{sample}/{modality}_{barcode}/fragments/fragments.tsv.gz'.format(sample = s, modality = m, barcode=barcodes_dict[s][m]))
    return result

def get_seurat_per_modality(modality, barcodes_dict,feature):
    result = []
    for s in barcodes_dict:
        for m in barcodes_dict[s]:
            if modality == m:
               result.append('results/{sample}/{modality}_{barcode}/seurat/{feature}/Seurat_object.Rds'.format(sample = s, modality = m, barcode=barcodes_dict[s][m], feature = feature))
    return result


rule all_single_modality:
    input:
        expand('results/single_modality/{modality}/fragments/fragments.tsv.gz',modality = antibodies_list),                                         # Merged fragments file
        expand('results/single_modality/{modality}/fragments/fragments.bw', modality = antibodies_list),
        expand('results/single_modality/{modality}/peaks/macs_broad/{modality}_peaks.broadPeak', modality = antibodies_list),                       # Merged peaks file
        expand('results/single_modality/{modality}/seurat/{feature}/Seurat_object.Rds', modality = antibodies_list, feature = features),                           # Seurat object peaks
        expand('results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds', modality= antibodies_list, feature = features),
        expand('results/single_modality/{modality}/seurat/{feature}/integration/integration_RNA.Rds', modality = antibodies_list, feature = features),
        expand('results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',modality = antibodies_list, feature = 'peaks'), # TODO: use 'peaks' as variable
        expand('results/single_modality/{modality}/seurat/{feature}/bigwig/{idents}/', modality = antibodies_list, feature = 'peaks', idents = ['idents_L1','idents_L2','idents_L3','seurat_clusters']), # TODO: use 'peaks' and idents as variable
        expand('results/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers.csv', modality = antibodies_list, feature = 'peaks', idents = ['idents_L1','idents_L2','idents_L3','seurat_clusters']),

rule merge_fragments_file:
    input:
        fragments = lambda wildcards: get_fragments_per_modality(modality = wildcards.modality, barcodes_dict = barcodes_dict)
    output:
        fragments_merged = 'results/single_modality/{modality}/fragments/fragments.tsv.gz',
        index            = 'results/single_modality/{modality}/fragments/fragments.tsv.gz.tbi'
    params:
        tmpdir           = config['general']['tempdir']
    shell:
        'zcat {input.fragments} | sort -T {params.tmpdir} -k1,1 -k2,2n | bgzip > {output.fragments_merged} && tabix -p bed {output.fragments_merged}'

rule fragments_to_bw:
    input:
        fragments   = 'results/single_modality/{modality}/fragments/fragments.tsv.gz',
        chrom_sizes = 'results/mm10.chrom.sizes'
    output:
        bam        = temp('results/single_modality/{modality}/fragments/fragments.bam'),
        bam_sorted = temp('results/single_modality/{modality}/fragments/fragments_sorted.bam'),
        index      = temp('results/single_modality/{modality}/fragments/fragments_sorted.bam.bai'),
        bigwig     = 'results/single_modality/{modality}/fragments/fragments.bw',
    threads: 8
    shell:
        "bedToBam -i {input.fragments} -g {input.chrom_sizes} > {output.bam} && "
        "samtools sort -@ {threads} -o {output.bam_sorted} {output.bam} &&"
        "samtools index {output.bam_sorted} && "
        "bamCoverage -b {output.bam_sorted} -o {output.bigwig} -p {threads} --minMappingQuality 5 --binSize 50 --smoothLength 250 --normalizeUsing RPKM --ignoreDuplicates"

rule run_macs_broad_merged:
    input:
        fragments = 'results/single_modality/{modality}/fragments/fragments.tsv.gz'
    output:
        broad_peaks = 'results/single_modality/{modality}/peaks/macs_broad/{modality}_peaks.broadPeak'
    params:
        macs_outdir = 'results/single_modality/{modality}/peaks/macs_broad/'
    shell:
        'macs2 callpeak -t {input} -g mm -f BED -n {wildcards.modality} '
        '--outdir {params.macs_outdir} --llocal 100000 --keep-dup=1 --broad-cutoff=0.1 ' 
        '--min-length 1000 --max-gap 1000 --broad --nomodel 2>&1 '


rule merge_seurat_single_modality:
    input:
        seurat = lambda wildcards: get_seurat_per_modality(modality=wildcards.modality,barcodes_dict=barcodes_dict, feature=wildcards.feature),
        script = workflow_dir + '/scripts/merge_objects.R'
    output:
        seurat = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object.Rds'
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat}'


rule cluster:
    input:
        seurat = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object.Rds',
        script = workflow_dir + '/scripts/UMAP_cluster.R'
    output:
        seurat = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds'
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat} -a {wildcards.feature} -d 40 '


rule integrate_with_scRNA:
    input:
        seurat = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds',
        rna    = '/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/Sten_RNA/clustering/01.clustering_20000cells.Rds', # TODO fix path here
        script = workflow_dir + '/scripts/integrate_CT_RNAseq.R'
    output:
        'results/single_modality/{modality}/seurat/{feature}/integration/integration_RNA.Rds'
    shell:
        "Rscript {input.script} -i {input.seurat} -r {input.rna} -o {output}"


rule cluster_final_and_rename:
    input:
        seurat   = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds',
        notebook =  workflow_dir + '/notebooks/single_modality/{modality}_rename_clusters.Rmd',
    output:
        'results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds'
    params:
        report     = os.getcwd() + '/results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.html',
        out_prefix = os.getcwd() + '/results/'
    shell:
        "Rscript -e \"rmarkdown::render(input='{input.notebook}',\
                                        output_file = '{params.report}', \
                                        params=list(out_prefix = '{params.out_prefix}',modality = '{wildcards.modality}', feature = '{wildcards.feature}'))\" "


rule export_bw:
    input:
        seurat    = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',
        fragments = 'results/single_modality/{modality}/fragments/fragments.tsv.gz',
        script    = workflow_dir + '/scripts/export_bw.R'
    output:
        bigwig    = directory('results/single_modality/{modality}/seurat/{feature}/bigwig/{idents}/')
    shell:
        "Rscript {input.script}  --input {input.seurat} --fragments {input.fragments} --output_folder {output.bigwig} --idents {wildcards.idents} "

rule find_markers:
    input:
        seurat = 'results/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',
        script = workflow_dir + '/scripts/find_markers.R'
    output:
        markers = 'results/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers.csv',
    shell:
        "Rscript {input.script} -i {input.seurat} -o {output.markers} --idents {wildcards.idents}"