#!/bin/bash
# run_pipeline.sh

# activate conda env
source /home/hkaufm49/anaconda3/etc/profile.d/conda.sh

conda activate regnoise_env

Rscript pipeline/scripts/bp_cells_time_measure.R

