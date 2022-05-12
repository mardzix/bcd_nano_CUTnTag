include: "Snakefile_preprocess.smk"
include: "Snakefile_pre_nbiotech.smk"

idents = ['idents_L1','idents_L2','idents_L3','seurat_clusters','idents_short']

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


rule all_single_modality:
    input:
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/integration/integration_RNA.Rds', modality = antibodies_list, feature = 'peaks'),
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',modality = antibodies_list, feature = 'peaks'), # TODO: use 'peaks' as variable
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/bigwig/{idents}/', modality = antibodies_list, feature = 'peaks', idents = idents), # TODO: use 'peaks' and idents as variable
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers.csv', modality = antibodies_list, feature = 'peaks', idents = idents),
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers.bed', modality= antibodies_list, feature = 'peaks', idents = idents),
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers_positive.bed', modality= antibodies_list, feature = 'peaks', idents = idents),
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/L3_niche_markers2/L3_markers.csv', modality = antibodies_list, feature = 'peaks'),
        # Bam files
        # ['results/{sample}/{antibody}_{barcode}/bam/possorted_bam_sampleID.bam'.format(sample=sample,antibody=antibody,barcode=barcodes_dict[sample][antibody]) for sample in samples_list  for antibody in barcodes_dict[sample].keys()],
        expand('results/multimodal_data/single_modality/{modality}/bam/possorted_bam_sampleID.bam',modality = antibodies_list),
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/bam/',modality=antibodies_list,feature='peaks',ident = idents),
        expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/bigwig/',modality=antibodies_list,feature='peaks',ident = idents),



rule integrate_with_scRNA:
    input:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds',
        rna    = '/data/proj/GCB_MB/single-cell-CUT-Tag/nbiotech_paper/analysis/results/Sten_RNA/clustering/01.clustering_20000cells.Rds', # TODO fix path here -> Config
        script = workflow_dir + '/scripts/integrate_CT_RNAseq.R'
    output:
        'results/multimodal_data/single_modality/{modality}/seurat/{feature}/integration/integration_RNA.Rds'
    shell:
        "Rscript {input.script} -i {input.seurat} -r {input.rna} -o {output}"


rule cluster_final_and_rename:
    input:
        seurat     = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds',
        notebook   =  workflow_dir + '/notebooks/single_modality/{modality}_rename_clusters.Rmd',
        integrated =  'results/multimodal_data/single_modality/{modality}/seurat/{feature}/integration/integration_RNA.Rds',
    output:
        seurat   = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds'
    params:
        report     = os.getcwd() + '/results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.html',
        out_prefix = os.getcwd() + '/',
    shell:
        "Rscript -e \"rmarkdown::render(input='{input.notebook}', "
        "                                output_file = '{params.report}', "
        "                                params=list(out_prefix = '{params.out_prefix}', "
        "                                           modality = '{wildcards.modality}', "
        "                                           feature = '{wildcards.feature}', "
        "                                           input = '{params.out_prefix}{input.seurat}', "
        "                                           integrated = '{params.out_prefix}{input.integrated}', "                # TODO - fix absolute paths integration                         
        "                                           output = '{params.out_prefix}{output.seurat}'))\" "


rule export_bw:
    input:
        seurat    = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',
        fragments = 'results/multimodal_data/single_modality/{modality}/fragments/fragments.tsv.gz',
        script    = workflow_dir + '/scripts/export_bw.R'
    output:
        bigwig    = directory('results/multimodal_data/single_modality/{modality}/seurat/{feature}/bigwig/{idents}/')
    shell:
        "Rscript {input.script}  --input {input.seurat} --fragments {input.fragments} --output_folder {output.bigwig} --idents {wildcards.idents} "

rule find_markers:
    input:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',
        script = workflow_dir + '/scripts/find_markers.R'
    output:
        markers     = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers.csv',
        markers_pos = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/markers/{idents}/markers_positive.csv',
    shell:
        "Rscript {input.script} -i {input.seurat} -o {output.markers} --idents {wildcards.idents}"

rule find_L3_markers:
    input:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',
        script = workflow_dir + '/scripts/find_L3_markers.R'
    output:
        'results/multimodal_data/single_modality/{modality}/seurat/{feature}/L3_niche_markers2/L3_markers.csv'
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output}'

rule export_cluster_barcode_table:
    input:
        seurat = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds',
        script = workflow_dir + '/scripts/export_cluster_barcode_table.R'
    output:
        csv = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/cluster_barcode_table.csv',
    shell:
        "Rscript {input.script} -i {input.seurat} -o {output.csv} -d {wildcards.ident}"

rule add_sampleID_to_bam:
    input:
        bam = 'results/multimodal_data/{sample}/cellranger/{sample}_{antibody}_{barcode}/outs/possorted_bam.bam', \
        script = workflow_dir + '/scripts/add_sampleID_to_bam.py'
    output:
        bam = temp('results/multimodal_data/{sample}/{antibody}_{barcode}/bam/possorted_bam_sampleID.bam'),
    shell:
        "python3 {input.script} {input.bam} {wildcards.sample} {output.bam}"


def get_bamfiles_per_modality(modality, barcodes_dict):
    bam_files = []
    for sample in barcodes_dict.keys():
        if modality in barcodes_dict[sample].keys():
            bam = 'results/multimodal_data/{sample}/{modality}_{barcode}/bam/possorted_bam_sampleID.bam'.format(sample=sample,modality=modality,barcode=
            barcodes_dict[sample][modality])
            bam_files.append(bam)
    return bam_files


rule merge_bam_accross_samples:
    input:
        bam=lambda wildcards: get_bamfiles_per_modality(wildcards.modality,barcodes_dict=barcodes_dict)
    output:
        bam='results/multimodal_data/single_modality/{modality}/bam/possorted_bam_sampleID.bam',
    threads: 8
    shell:
        "samtools merge -@ {threads} {output.bam} {input.bam}"

rule export_bam_per_cluster:
    input:
        bam    = 'results/multimodal_data/single_modality/{modality}/bam/possorted_bam_sampleID.bam',
        table  = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/cluster_barcode_table.csv',
        script = workflow_dir + '/scripts/filter_bam_by_barcode.py'
    output:
        bam_files = directory('results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/bam/'),
    shell:
        "python3 {input.script} {input.bam} {input.table} NA {output.bam_files}"

rule bam_to_bw_per_cluster:
    input:
        bam = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/bam/',
        script = workflow_dir + '/scripts/all_bam_to_bw.sh'
    output:
        bw = directory('results/multimodal_data/single_modality/{modality}/seurat/{feature}/bam_per_cluster/{ident}/bigwig/')
    threads: 8
    shell:
        'sh {input.script} {input.bam} {output.bw} {threads}'
