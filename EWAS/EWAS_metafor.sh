#!/bin/bash


wd=/mnt/data1/Morteza/Methylation_Data/ewas
ScriptDir=/mnt/data1/Morteza/Methylation_Data/ewas/Scripts

for i in {0..10}
do
  study1=Pitts/Quantile/Pitts8_EWAS_RESULTS_Annot_${i}.Rdata
  for j in {0..10}
  do
    study2=BDR/Dec2022/Quantile/BDR10.quantile_EWAS_RESULTS_Annot_${j}.Rdata
    out_pref=metafor/Quantile/result_quantile_Pitts8.${i}.BDR10.${j}
    echo $study1 , $study2
    Rscript ${ScriptDir}/EWAS_metafor.R  $wd $study1 $study2 $out_pref
  done
done

#study1=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/final_Nov10_2022/Pitts8_EWAS_RESULTS_Annot_1.Rdata
#study2=/mnt/data1/Morteza/Methylation_Data/ewas/BDR/final_Nov10_2022/BDR10.1_EWAS_RESULTS_Annot_2.Rdata
#out_pref=/mnt/data1/Morteza/Methylation_Data/ewas/BDR/final_Nov10_2022/BDR10.1.sv0.Pitts8.sv1.
#  echo $study1 , $study2
#Rscript /mnt/data1/Morteza/Methylation_Data/ewas/Scripts/EWAS_metafor.R  $wd $study1 $study2 $out_pref
