#!/bin/bash


in_pref=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/
out_pref=Feb2023/EWAS/Pitts.Pheno.NoControl.0.2.CDR

beta_file=${in_pref}10396_Normalised.NoControl.0.2.prepared.rds
pheno_file=${in_pref}Pitts.Pheno.NoControl.0.2.csv
trait=Psychosis
covars_fact=Sex,Tissue_Type,SentrixID,Plate,Basename,last_CDR
covars_num=Age,Prop
batchs=0
model_combat=0
model_sva=~Psychosis+Age+Sex+Prop+Plate+Tissue_Type+last_CDR
max_sor_var=10
model_lm=~Psychosis+Age+Sex+Prop+Plate+Tissue_Type+last_CDR
calc_sv=1
sv_file=0

script_dir=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Scripts

Rscript ${script_dir}/EWAS_CheckSurr.R $beta_file $pheno_file $trait $covars_fact $covars_num $batchs $model_combat $model_sva $max_sor_var $model_lm $out_pref $calc_sv $sv_file

