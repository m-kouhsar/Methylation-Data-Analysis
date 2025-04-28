#!/bin/bash

in_dir=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Feb2023/Annotation/Braak
out_dir=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Feb2023/comb-p/Braak

cd $wd

for i in ${in_dir}/*_dmr_*.txt
do
  file_name=${i%".txt"}
  file_name=${file_name#"${in_dir}/"}
  out_dir1=${out_dir}/${file_name}
  mkdir -p $out_dir1
  echo $file_name
  #echo $out_dir1
  comb-p pipeline  -c 4 --seed 1e-3 --dist 1000 -p ${out_dir1}/${file_name} --anno hg19 $i
done

echo "All done!"
