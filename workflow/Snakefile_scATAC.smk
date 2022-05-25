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
        # Bing Ren's scATAC-seq
        'results/scATAC_bingren/seurat/Seurat_ATAC_clustered.Rds',
        'results/scATAC_bingren/seurat/Seurat_ATAC_fragments.bed.gz',
        # expand('results/multimodal_data/single_modality/{modality}/seurat/{feature}/integration/integrated_with_ATAC_{nfeatures}.Rds', modality = ['ATAC','H3K27ac'], feature = 'peaks', nfeatures = [1000,5000,10000,20000,40000,50000,75000,100000,125000]),

        # 10x scATAC-seq
        'results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_filtered_peak_bc_matrix.h5',
        'results/scATAC_10x/seurat/Seurat_object.Rds',
        # expand('results/scATAC_10x/seurat/bigwig/{idents}/', idents = ['predicted.id']),
        expand('results/scATAC_10x/seurat/bam_per_cluster/{idents}/bam/', idents = ['predicted.id']),
        expand('results/scATAC_10x/seurat/bam_per_cluster/{idents}/bigwig/', idents = ['predicted.id']),
        expand('results/scATAC_10x/seurat/bam_per_cluster/{idents}/peaks/', idents =['predicted.id']),

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

rule download_10x_ATAC:
    input:
        script = workflow_dir + '/scripts/download_10x_scATAC.sh'
    output:
        matrix    = 'results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_filtered_peak_bc_matrix.h5',
        metadata  = 'results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_singlecell.csv',
        fragments = 'results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_fragments.tsv.gz',
        allen_RNA = 'results/scATAC_10x/download/allen_brain.rds',
    params:
        folder   = 'results/scATAC_10x/download/'
    shell:
        'cd {params.folder}; '
        'sh {input.script}'

rule run_10x_ATAC:
    input:
        matrix    = os.getcwd() + '/results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_filtered_peak_bc_matrix.h5',
        metadata  = os.getcwd() + '/results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_singlecell.csv',
        fragments = os.getcwd() + '/results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_fragments.tsv.gz',
        allen_RNA = os.getcwd() + '/results/scATAC_10x/download/allen_brain.rds',
        script    = workflow_dir + '/notebooks/Adult_ATAC_10x.ipynb',
    output:
        seurat = 'results/scATAC_10x/seurat/Seurat_object.Rds',
    params:
        out_prefix = os.getcwd() + '/',
        report     = os.getcwd() + '/results/scATAC_10x/seurat/seurat_report.ipynb',
        rmd        = workflow_dir + '/notebooks/Adult_ATAC_10x.Rmd',
    shell:
        "papermill {input.script} {params.report} "
        "                         -p out_prefix {params.out_prefix} "
        "                         -p modality ATAC "
        "                         -p feature peaks "
        "                         -p data_10x {input.matrix} "
        "                         -p metadata {input.metadata} "
        "                         -p fragments {input.fragments} "
        "                         -p allen_RNA {input.allen_RNA} "
        "                         -p out_seurat {output.seurat}"

rule export_cluster_barcode_table_scATAC_10x:
    input:
        seurat = 'results/scATAC_10x/seurat/Seurat_object.Rds',
        script = workflow_dir + '/scripts/export_cluster_barcode_table.R'
    output:
        csv = 'results/scATAC_10x/seurat/{idents}/cluster_barcode_table.csv',
    shell:
        "Rscript {input.script} -i {input.seurat} -o {output.csv} -d {wildcards.idents}"

# rule export_bw_scATAC_10x:
#     input:
#         seurat    = 'results/scATAC_10x/seurat/Seurat_object.Rds',
#         fragments = 'results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_fragments.tsv.gz',
#         script    = workflow_dir + '/scripts/export_bw.R'
#     output:
#         bigwig    = directory('results/scATAC_10x/seurat/bigwig/{idents}/')
#     shell:
#         "Rscript {input.script}  --input {input.seurat} --fragments {input.fragments} --output_folder {output.bigwig} --idents {wildcards.idents} "


rule export_bam_per_cluster_scATAC_10x:
    input:
        bam    = 'results/scATAC_10x/download/8k_mouse_cortex_ATACv1p1_nextgem_Chromium_X_possorted_bam.bam',
        table  = 'results/scATAC_10x/seurat/{idents}/cluster_barcode_table.csv',
        script = workflow_dir + '/scripts/filter_bam_by_barcode.py'
    output:
        bam_files = directory('results/scATAC_10x/seurat/bam_per_cluster/{idents}/bam/'),
    shell:
        "python3 {input.script} {input.bam} {input.table} NA {output.bam_files}"

rule bam_to_bw_per_cluster_scATAC_10x:
    input:
        bam = 'results/scATAC_10x/seurat/bam_per_cluster/{idents}/bam/',
        script = workflow_dir + '/scripts/all_bam_to_bw.sh'
    output:
        bw = directory('results/scATAC_10x/seurat/bam_per_cluster/{idents}/bigwig/')
    threads: 8
    shell:
        'sh {input.script} {input.bam} {output.bw} {threads}'

rule bam_to_peaks_per_cluster_scATAC_10x:
    input:
        bam    = 'results/scATAC_10x/seurat/bam_per_cluster/{idents}/bam/',
        script = workflow_dir + '/scripts/all_bam_to_peaks.sh'
    output:
        peaks  = directory('results/scATAC_10x/seurat/bam_per_cluster/{idents}/peaks/'),
    shell:
        'sh {input.script}  {input.bam} {output.peaks}'
