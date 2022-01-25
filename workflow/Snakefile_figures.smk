include: "Snakefile_prep.smk"
include: "Snakefile_preprocess.smk"
include: "Snakefile_single_modality.smk"



rule figures_all:
    input:
        'results/figures/figure_1/UMAP_integrated.pdf'


rule figures_figure1:
    input:
        sequential = '/data/proj/GCB_MB/single-cell-CUT-Tag/sequential/mouse_brain/results/H3K27me3/clustering/01.clustering.Rds', # TODO integrate better into the pipelina
        nbiotech   = '/data/proj/GCB_MB/bcd_CT/single-cell/results/nbiotech_data/data/seurat/H3K27me3_seurat_object.Rds',
        notebook   = '/data/proj/GCB_MB/bcd_CT/single-cell/code/notebooks/figures/figure1.Rmd',
    output:
        'results/figures/figure_1/UMAP_integrated.pdf'
    params:
        report     = os.getcwd() + '/results/figures/figure_1/figure_1.html',
        out_folder = os.getcwd() + '/results/'
    shell:
        "Rscript -e \"rmarkdown::render(input='{input.notebook}',\
                                        output_file = '{params.report}', \
                                        params=list('nbiotech'='{input.nbiotech}', 'sequential'='{input.sequential}','out_folder'='{params.out_folder}'))\" "
