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
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

module load R

out_pref=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f/Blood_EWAS/Blood.EWAS.CohensF.SV3

beta_file=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f/Blood_EWAS/10399_betas.rds
pheno_file=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f/Blood_EWAS/10399_pheno.csv
trait=Phenotype
covars_fact=Sex,Plate
covars_num=Age,CD8T,CD4T,Bcell,Gran,Mono,NK
batchs=0
model_combat=0
model_sva=~Phenotype+Age+Sex+Plate+CD8T+CD4T+Bcell+Gran+Mono+NK
num_sor_var=3
model_lm=~Phenotype+Age+Sex+Plate+CD8T+CD4T+Bcell+Gran+Mono+NK
calc_sv=0
sv_file=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f/Blood_EWAS/10399_sva.Rdata
num_threads=16

script_dir=/lustre/projects/Research_Project-191391/Morteza/Genotyping/Pitts.All/cohens_f/Scripts

Rscript ${script_dir}/EWAS_CohensF.R $beta_file $pheno_file $trait $covars_fact $covars_num $batchs $model_combat $model_sva $num_sor_var $model_lm $out_pref $calc_sv $sv_file  $num_threads
