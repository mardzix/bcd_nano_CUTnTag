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
        'results/scATAC_bingren/seurat/Seurat_ATAC_clustered.Rds',
        'results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed.gz',
        # expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/integration/integrated_with_ATAC_{nfeatures}.Rds', modality = ['ATAC','H3K27ac'], feature = 'peaks', nfeatures = [1000,5000,10000,20000,40000,50000,75000,100000,125000]),


############################# Bing ren mouse brain
# Li, Y.E., Preissl, S., Hou, X. et al. An atlas of gene regulatory elements in adult mouse cerebrum. Nature 598, 129–136 (2021). https://doi.org/10.1038/s41586-021-03604-1

rule download_metadata:
    output:
        'results/scATAC_bingren/download/snap/metadata.tsv'
    params:
        out_zip = 'results/scATAC_bingren/download/snap/metadata.tsv.zip',
        url = 'http://renlab.sdsc.edu/yangli/downloads/mousebrain/Supplementary_Tables/Supplementary%20Table%202%20-%20Metatable%20of%20nuclei.tsv.zip'
    shell:
        'wget -O {params.out_zip} {params.url} && '
        'gunzip -c {params.out_zip} > {output}'

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

rule create_seurat_from_snap:
    input:
        snap = expand('results/scATAC_bingren/download/snap/{file}',file = get_filenames_from_url_txt(workflow_dir + '/config/bingren_snap_files.txt')),
        metadata = 'results/scATAC_bingren/download/snap/metadata.tsv',
        script = workflow_dir + '/scripts/snap_to_seurat.R'
    output:
        seurat = 'results/scATAC_bingren/seurat/Seurat_ATAC.Rds',
        snap = 'results/scATAC_bingren/seurat/Seurat_ATAC_snap.Rds',
        fragments = temp('results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed'),
    params:
        ncells = 50000,
        log= 'results/scATAC_bingren/seurat/seurat.log'
    shell:
        "Rscript {input.script} --snap {input.snap} --metadata {input.metadata} --ncells {params.ncells} --output {output.seurat} > {params.log}"

rule process_fragments_bed:
    input:
        'results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed'
    output:
        bed = 'results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed.gz',
        index = 'results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed.gz.tbi',
    shell:
        'cut -f1-5 {input} | bgzip > {output.bed} && '
        'tabix -p bed {output.bed}'

rule cluster_bingren:
    input:
        seurat = 'results/scATAC_bingren/seurat/Seurat_ATAC.Rds',
        script = workflow_dir + '/scripts/UMAP_cluster.R'
    output:
        seurat = 'results/scATAC_bingren/seurat/Seurat_ATAC_clustered.Rds',
    params:
        assay = 'peaks',
        ndim   = 40
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat} -a {params.assay} -d {params.ndim} '

rule integrate_with_CT:
    input:
        atac = 'results/scATAC_bingren/seurat/Seurat_ATAC_clustered.Rds',
        ct = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered.Rds',
        script = workflow_dir + '/scripts/integrate_CT_ATAC.R',
        ref_frag = '/data/proj/GCB_MB/bcd_CT/single-cell/results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed.gz'
    output:
        integrated = 'results/multimodal_data/single_modality/{modality}/seurat/{feature}/integration/integrated_with_ATAC_{nfeatures}.Rds'
    threads: 16
    shell:
        'Rscript {input.script} --reference {input.atac} --query {input.ct} '
        '--reference_assay peaks --query_assay peaks '
        '--reference_fragments {input.ref_frag} --out {output.integrated} '
        '--reference_group MajorType --query_group seurat_clusters --downsample_features {wildcards.nfeatures}'
