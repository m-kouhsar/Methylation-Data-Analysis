#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=4:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=QC.EPIC.V1.%j.out

source $1

OutDir="$(dirname $OutPrefix)"
OutName="$(basename $OutPrefix)"
mkdir -p $OutDir

Rscript -e "rmarkdown::render('$ScriptDir/DNAm_QC_EPIC_V1.Rmd', output_dir='${OutDir}',output_file='${OutName}_DNAm_QC_EPIC_V1.html',\
            output_options = list(pandoc_args = c(paste0('--include-before-body=', '$ScriptDir/References/DNAm_QC_Header.html'),paste0('--include-after-body=','$ScriptDir/References/DNAm_QC_Header.html'))), \
            params=list(Title='${Title}',Study_Name='${Study_Name}',User_Array='${User_Array}',User_QC='${User_QC}'))" \
            --args "${ScriptDir}" "${OutPrefix}" "${idatPath}" "${SampleSheet}" "${IntensityThreshold}" 


