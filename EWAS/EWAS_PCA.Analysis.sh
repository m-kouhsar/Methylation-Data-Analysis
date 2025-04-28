#!/bin/bash

##beta_file: Name of beta values file in rds format
##pheno_file: name of phenotype information file in csv format
##variables_numeric: names of numeric variables (name of pheno_file columns) that you want to see in the correlation plot
##variables_factor: names of non-numeric variables (name of pheno_file columns) that you want to see in the correlation plot
##out_prefix: Prefix string for naming the outputs

beta_file=10396_Normalised.NoControl.0.2.prepared.rds
pheno_file=Pitts.Pheno.NoControl.0.2.csv
variables_numeric=Age,Prop
variables_factor=Plate,Sex,Tissue_Type,SentrixID,Psychosis
out_prefix=10396_Normalised.NoControl.0.2

ScriptDir=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Scripts

Rscript ${ScriptDir}/EWAS_PCA.Analysis.R ${beta_file} ${pheno_file} ${variables_numeric} ${variables_factor} ${out_prefix}


