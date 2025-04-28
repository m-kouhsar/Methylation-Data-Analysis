#!/bin/bash

normalized_file=10396_Normalised.NoControl.0.2.rdat
cross_hyd_ref=/mnt/data1/EPIC_reference/CrossHydridisingProbes_McCartney.txt
SNPPprob_ref=/mnt/data1/EPIC_reference/SNPProbes_McCartney.txt
EPIC_ref=/mnt/data1/EPIC_reference/MethylationEPIC_v-1-0_B4.csv
out_file=10396_Normalised.NoControl.0.2

ScriptDir=/mnt/data1/Morteza/Methylation_Data/ewas/Pitts/Pitts.All/Scripts

Rscript ${ScriptDir}/EWAS_PrepData.R ${normalized_file} ${cross_hyd_ref} ${SNPPprob_ref} ${EPIC_ref} ${out_file}
