include: "Snakefile_preprocess.smk"
include: "Snakefile_single_modality.smk"

def get_input_for_multiple_modalities_merge(combination, feature):
    combination = combination.split("_")
    files = ["results/multimodal_data/single_modality/{modality}/seurat/{feature}/Seurat_object_clustered_renamed.Rds".format( modality = modality, feature = feature) for modality in combination]
    return files

rule all_multimodal:
    input:
       # [expand('results/multiple_modalities/{combination}/seurat_multimodal/bin_{feature}/Seurat_object.Rds', combination="_".join(combination),feature=config['general']['binwidth']) for combination in modalities_combinations],
        [expand('results/multiple_modalities/{combination}/seurat_multimodal/{feature}/Seurat_object.Rds', combination="_".join(combination), feature = ['peaks']) for combination in modalities_combinations]

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
