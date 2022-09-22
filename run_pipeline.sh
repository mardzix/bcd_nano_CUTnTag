#!/bin/bash
snakemake --snakefile ../code/workflow/Snakefile_single_modality.smk  --cores 16 --profile htcondor -p
snakemake --snakefile ../code/workflow/Snakefile_multimodal.smk  --cores 16 --profile htcondor -p
snakemake --snakefile ../code/workflow/Snakefile_nbiotech.smk  --cores 4 --profile htcondor -p --resources load=50 --rerun-incomplete
snakemake --snakefile ../code/workflow/Snakefile_scATAC.smk  --cores 16 --profile htcondor -p
