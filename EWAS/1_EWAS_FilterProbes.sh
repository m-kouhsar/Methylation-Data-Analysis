#!/bin/bash

beta_file="$1"
OutputPrefix="$2"

ScriptDir=/lustre/projects/Research_Project-191391/Morteza/github/Methylation-Data-Analysis/EWAS
#################################################################################################

Rscript ${ScriptDir}/1_EWAS_FilterProbes.R "$beta_file" "$OutputPrefix" "$ScriptDir"

echo "All done!"