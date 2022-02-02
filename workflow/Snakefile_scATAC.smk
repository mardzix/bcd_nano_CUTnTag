rule scATAC_all:
    input:
        "results/scATAC_bingren/meta_nonN.tsv"

# rule download:
#     output:
#         'results/scATAC/download/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5',
#         'results/scATAC/download/atac_v1_adult_brain_fresh_5k_singlecell.csv',
#         'results/scATAC/download/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz',
#         'results/scATAC/download/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz.tbi',
#     params:
#         outdir = 'results/scATAC/download/'
#     shell:
#         "cd {params.outdir}; "
#         "wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5; "
#         "wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_singlecell.csv; "
#         "wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz; "
#         "wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz.tbi; "


############################# Bing ren mouse brain
# Li, Y.E., Preissl, S., Hou, X. et al. An atlas of gene regulatory elements in adult mouse cerebrum. Nature 598, 129–136 (2021). https://doi.org/10.1038/s41586-021-03604-1

rule download:
    output:
        mat_nonN  = "results/scATAC_bingren/matrix_nonN.tsv.gz",
        meta_nonN = "results/scATAC_bingren/meta_nonN.tsv",
        mat_GABA  = "results/scATAC_bingren/matrix_GABA.tsv.gz",
        meta_GABA = "results/scATAC_bingren/meta_GABA.tsv",
        mat_GLUT  = "results/scATAC_bingren/matrix_GLUT.tsv.gz",
        meta_GLUT = "results/scATAC_bingren/meta_GLUT.tsv",
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
