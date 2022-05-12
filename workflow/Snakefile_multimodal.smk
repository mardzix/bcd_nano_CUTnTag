include: "Snakefile_preprocess.smk"
include: "Snakefile_single_modality.smk"

def get_input_for_multiple_modalities_merge(combination, feature):
    combination = combination.split("_")
    files = ["results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds".format( modality = modality, feature = feature) for modality in combination]
    return files

rule all_multimodal:
    input:
       # [expand('results/multiple_modalities/{combination}/seurat_multimodal/bin_{feature}/Seurat_object.Rds', combination="_".join(combination),feature=config['general']['binwidth']) for combination in modalities_combinations],
        [expand('results/multiple_modalities/{combination}/seurat_multimodal/{feature}/Seurat_object.Rds', combination="_".join(combination), feature = ['peaks']) for combination in modalities_combinations],
        "results/multimodal_data/WNN/seurat/Seurat_object_WNN.Rds",
        "results/multimodal_data/WNN/seurat/pseudotime/Seurat_object_WNN_pseudotime.Rds",
        expand('results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/scvelo/{modalities}/Seurat_object_{clusters}.h5Seurat',modalities=['ATAC_H3K27ac','ATAC_H3K27me3','H3K27me3_H3K27ac'], clusters = ['all_clusters','OPC_MOL'])

rule merged_multiple_modalities:
    input:
        seurat = lambda wildcards: get_input_for_multiple_modalities_merge(combination=wildcards.combination,
                                                                           feature = wildcards.feature),
        script = os.path.dirname(workflow.basedir) + '/scripts/merge_modalities2.R'
    output:
        'results/multiple_modalities/{combination}/seurat_multimodal/{feature}/Seurat_object.Rds'

    shell:
        'Rscript {input.script} -i {input.seurat} -m {wildcards.combination} -o {output}'

rule run_WNN:
    input:
        seurat = "results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/Seurat_object.Rds",
        notebook = workflow_dir + '/notebooks/multimodal/WNN_Seurat.Rmd'
    output:
        seurat_out = "results/multimodal_data/WNN/seurat/Seurat_object_WNN.Rds",
    params:
        report     = os.getcwd() + '/results/multimodal_data/WNN/seurat/Seurat_object_WNN.html',
        seurat_in  = os.getcwd() + "/results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/Seurat_object.Rds",
        seurat_out = os.getcwd() + "/results/multimodal_data/WNN/seurat/Seurat_object_WNN.Rds",
        notebook   = workflow_dir + '/code/notebooks/multimodal/WNN_Seurat.Rmd'
    shell:
        "Rscript -e \"rmarkdown::render(input='{input.notebook}', "
        "                                output_file = '{params.report}', "
        "                                params=list(seurat = '{params.seurat_in}', "
        "                                           output = '{params.seurat_out}'))\" "


rule pseudotime_WNN:
    input:
        seurat = "results/multimodal_data/WNN/seurat/Seurat_object_WNN.Rds",
        script = workflow_dir + '/scripts/run_slingshot.R'
    output:
        seurat = "results/multimodal_data/WNN/seurat/pseudotime/Seurat_object_WNN_pseudotime.Rds",
    shell:
        "Rscript {input.script} --seurat {input.seurat} --clusters OPC MOL --reduction wnn.umap --out {output.seurat}"

rule convert_seurat_to_anndata:
    input:
        seurat =    'results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/Seurat_object.Rds',
        script =    workflow_dir + '/scripts/seurat_to_scvelo.R',
        seurat_pt = 'results/multimodal_data/WNN/seurat/pseudotime/Seurat_object_WNN_pseudotime.Rds'
    output:
        seurat = 'results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/scvelo/{modalities}/Seurat_object_{clusters}.h5Seurat',
        h5ann  = 'results/multiple_modalities/ATAC_H3K27ac_H3K27me3/seurat_multimodal/peaks/scvelo/{modalities}/Seurat_object_{clusters}.h5ad'
    params:
        clusters   = lambda wildcards: '' if wildcards.clusters == 'all_clusters' else '-c ' + ' '.join(wildcards.clusters.split('_')) ,
        idents     = 'idents_short',
        modalities = lambda wildcards: ' '.join(wildcards.modalities.split('_')),
        assay      = 'GA'
    shell:
        'Rscript {input.script} -i {input.seurat} -o {output.seurat} {params.clusters} -d {params.idents} -m {params.modalities} -a {params.assay} -t {input.seurat_pt}'