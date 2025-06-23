#!/bin/bash

beta_file="/lustre/projects/Research_Project-191391/Morteza/KCL_FACS/Results/QC/KCL_FACS_EPICV1.Normalized.rdat"
pheno_file="/lustre/projects/Research_Project-191391/Morteza/KCL_FACS/Raw/FACS_samples_all_plates_combined.csv"
EPIC_manifest="/lustre/projects/Research_Project-191391/Morteza/Ref/MethylationEPIC_v-1-0_B4.csv"
OutputPrefix=/lustre/projects/Research_Project-191391/Morteza/KCL_FACS/Results/EWAS/KCL_FACS

ScriptDir=/lustre/projects/Research_Project-191391/Morteza/github/Methylation-Data-Analysis

Rscript ${ScriptDir}/EWAS_PrepData.R "$normalized_file" "$EPIC_manifest" "$OutputPrefix" "$ScriptDir"
