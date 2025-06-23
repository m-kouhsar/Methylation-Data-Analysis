#!/bin/bash

beta_file="/lustre/projects/Research_Project-191391/Morteza/KCL_FACS/Raw/KCL_FACS.betas.rds"
OutputPrefix=/lustre/projects/Research_Project-191391/Morteza/KCL_FACS/Results/EWAS/KCL_FACS

ScriptDir=/lustre/projects/Research_Project-191391/Morteza/github/Methylation-Data-Analysis/EWAS

Rscript ${ScriptDir}/1_EWAS_FilterProbes.R "$beta_file" "$OutputPrefix" "$ScriptDir"
