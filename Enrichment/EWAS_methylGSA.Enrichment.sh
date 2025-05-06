#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address

in_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f
EWAS_Result_file=${in_dir}/Pitts.Psycho0.2.SV1.CohensF_allcpg.csv
pval_col=4
sig_cpg_file=${in_dir}/Pitts.Psycho0.2.Sig.CohensF0.25.Pval.e-3.csv
array_type=EPIC
method=ORA
GS_type=KEGG
minsize=10
maxsize=500
out_pref=${in_dir}/Pitts.Psycho0.2

Rscript Scripts/EWAS_methylGSA.Enrichment.R $EWAS_Result_file $pval_col $sig_cpg_file $array_type $method $GS_type $minsize $maxsize $out_pref 
