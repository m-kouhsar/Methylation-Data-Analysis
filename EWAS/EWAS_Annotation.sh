#!/bin/bash


in_pref=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Feb2023/EWAS/CDR/Pitts.Pheno.NoControl.0.2.CDR
out_pref=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Feb2023/Annotation/CDR/Pitts.Pheno.NoControl.0.2.CDR
run_lm=0
model_lm=~Psychosis+Age+Sex+Prop+Plate+Tissue_Type
sva_file=${in_pref}_sva.Rdata
num_sur=2
ewas_file=${in_pref}_EWAS_RESULTS_SV2.Rdata
EPIC_ref=/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv

ScriptDir=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Scripts

for i in {0..10}
do
  num_sur=$i
  ewas_file=${in_pref}_EWAS_RESULTS_SV${i}.Rdata
  
  echo $ewas_file
  Rscript ${ScriptDir}/EWAS_Annotation.R $out_pref $run_lm $model_lm $num_sur $sva_file $ewas_file $EPIC_ref
done

#Rscript ${ScriptDir}/EWAS_Annotation.R  $wd $out_pref $run_lm $model_lm $num_sur $sva_file $ewas_file $EPIC_ref

