include: "Snakefile_prep.smk"
include: "Snakefile_single_modality.smk"

samples = ['nonN','GABA','GLUT']

def get_filenames_from_url_txt(path):
    import os
    result = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            result.append(os.path.basename(line))
    return result

rule scATAC_all:
    input:
#        expand('results/scATAC_bingren/download/snap/{file}',file = get_filenames_from_url_txt(workflow_dir + '/config/bingren_snap_files.txt')),              # SNAP files
        "results/scATAC_bingren/seurat/Seurat_merged_clustered.Rds",                                                                                          # Seurat object
        # Integration bingren data
        expand("results/multimodal_data/single_modality/{modality}/seurat/peaks/integration/integrated_with_ATAC_{reference}.Rds", modality= ['ATAC','H3K27ac'],reference = ['GABA','GLUT','nonN','merged_clustered']) # TODO move to the snakefile_single_modality


############################# Bing ren mouse brain
# Li, Y.E., Preissl, S., Hou, X. et al. An atlas of gene regulatory elements in adult mouse cerebrum. Nature 598, 129–136 (2021). https://doi.org/10.1038/s41586-021-03604-1

rule download_tsv_matrix:
    output:
        mat_nonN  = "results/scATAC_bingren/download/matrix_nonN.tsv.gz",
        meta_nonN = "results/scATAC_bingren/download/meta_nonN.tsv",
        mat_GABA  = "results/scATAC_bingren/download/matrix_GABA.tsv.gz",
        meta_GABA = "results/scATAC_bingren/download/meta_GABA.tsv",
        mat_GLUT  = "results/scATAC_bingren/download/matrix_GLUT.tsv.gz",
        meta_GLUT = "results/scATAC_bingren/download/meta_GLUT.tsv",
    params:
        mat_nonN  = "http://catlas.org/mousebrain/cellbrowser/NonN_all/exprMatrix.tsv.gz",
        meta_nonN = "http://catlas.org/mousebrain/cellbrowser/NonN_all/meta.tsv",
        mat_GABA  = "http://catlas.org/mousebrain/cellbrowser/GABA_all/exprMatrix.tsv.gz",
        meta_GABA = "http://catlas.org/mousebrain/cellbrowser/GABA_all/meta.tsv",
        mat_GLUT  = "http://catlas.org/mousebrain/cellbrowser/Glutamate_all/exprMatrix.tsv.gz",
        meta_GLUT = "http://catlas.org/mousebrain/cellbrowser/Glutamate_all/meta.tsv",
    shell:
        "wget -O {output.mat_nonN} {params.mat_nonN}; "
        "wget -O {output.mat_GABA} {params.mat_GABA}; "
        "wget -O {output.mat_GLUT} {params.mat_GLUT}; "
        "wget -O {output.meta_nonN} {params.meta_nonN}; "
        "wget -O {output.meta_GABA} {params.meta_GABA}; "
        "wget -O {output.meta_GLUT} {params.meta_GLUT}; "

rule download_snap_matrix:
    input:
        urls = workflow_dir + '/config/bingren_snap_files.txt'
    output:
        expand('results/scATAC_bingren/download/snap/{file}',file = get_filenames_from_url_txt(workflow_dir + '/config/bingren_snap_files.txt'))
    params:
        outdir = 'results/scATAC_bingren/download/snap/'
    shell:
        "cd {params.outdir}; "
        "wget -i {input.urls}; "


rule gunzip_data:
    input:
        "results/scATAC_bingren/download/matrix_{sample}.tsv.gz"
    output:
        "results/scATAC_bingren/download/matrix_{sample}.tsv"
    shell:
        "gunzip -c {input} > {output}"

rule subset_matrix:
    input:
        matrix = "results/scATAC_bingren/download/matrix_{sample}.tsv",
        script = workflow_dir + '/scripts/subset_bingren_scATAC.py'
    output:
        "results/scATAC_bingren/download/matrix_{sample}_subseted.tsv"
    shell:
        "python3 {input.script} {input.matrix} 0.1 {output}"

rule create_seurat_object_bingren:
    input:
        matrix = "results/scATAC_bingren/download/matrix_{sample}_subseted.tsv",
        metadata = "results/scATAC_bingren/download/meta_{sample}.tsv",
        script = workflow_dir + '/scripts/create_seurat_bingren.R'
    output:
        seurat = "results/scATAC_bingren/seurat/Seurat_{sample}.Rds"
    params:
        modality = 'ATAC'
    shell:
        "Rscript {input.script} --matrix {input.matrix} --metadata {input.metadata} --out {output.seurat} --modality {params.modality}"

rule merge_seurat_bingren:
    input:
        seurat = expand("results/scATAC_bingren/seurat/Seurat_{sample}.Rds",sample = ['nonN','GABA','GLUT']),
        script = workflow_dir + '/scripts/merge_objects.R'
    output:
        seurat = "results/scATAC_bingren/seurat/Seurat_merged.Rds"
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat}'

rule cluster_merged_bingren:
    input:
        seurat = "results/scATAC_bingren/seurat/Seurat_merged.Rds",
        script = workflow_dir + '/scripts/UMAP_cluster.R'
    output:
        seurat = "results/scATAC_bingren/seurat/Seurat_merged_clustered.Rds"
    params:
        assay = "GA",
        ndim  = 40
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat} -a {params.assay} -d {params.ndim} '

rule integrate_bingren_CT:
    input:
        script    = workflow_dir + '/scripts/integrate_CT_ATAC.R',
        reference = "results/scATAC_bingren/seurat/Seurat_{reference}.Rds",
        query     = "results/multimodal_data/single_modality/{modality}/seurat/peaks/Seurat_object_clustered_renamed.Rds"
    output:
        seurat = "results/multimodal_data/single_modality/{modality}/seurat/peaks/integration/integrated_with_ATAC_{reference}.Rds"
    params:
        reference_assay = 'GA',
        query_assay     = 'GA'
    shell:
        "Rscript {input.script} --reference {input.reference} --reference_assay {params.reference_assay} --query {input.query} --query_assay {params.query_assay} --out {output.seurat}"
