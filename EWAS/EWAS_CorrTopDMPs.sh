#!/bin/bash

in_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f
out_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f

script_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f/Scripts

EWAS_file1=${in_dir}/Pitts.Psycho0.2.SV1.CohensF_EWAS_RESULTS_SV1CohensF.Rdata
EWAS_file2=${in_dir}/BDR.Selected.Psycho0.2.SV1.CohensF_EWAS_RESULTS_SV2CohensF.Rdata

Num_cpg=100
out_prefix=${out_dir}/Pitts.Psycho.0.2.BDR.Selected.CohensF


#for i in ${in_dir}/*.Rdata
#do
  #echo $i
#  file_name=${i%".Rdata"}
#  file_name=${file_name#"${in_dir}/"}
#  echo $file_name
  
#	for j in {100,200,300,400,500}
#	do
 #    out_prefix=${out_dir}/${file_name}
     #echo $out_prefix
#		 Rscript ${script_dir}/EWAS_CorrTopDMPs.R $i $EWAS_file2 $j $out_prefix
#	done
#done

Rscript ${script_dir}/EWAS_CorrTopDMPs.R $EWAS_file2  $EWAS_file2 $Num_cpg  $out_prefix

echo "All done!"

