#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=escobosaolmo@gmail.com
#SBATCH --mem=8G
#SBATCH --job-name=crispr-a-web

#SBATCH -e stderr_nf_%j.err
#SBATCH -o stdout_nf_%j.out
module load  Singularity/3.7.1-nopaths-foss-2016b
module load nextflow/20.07.1-Java-11.0.4

SINGULARITY_IMAGES=/scratch/lab_mguell/shared_data/crispr-ga_nextflow/Singularity



# Run nextflow
nextflow run /scratch/lab_mguell/shared_data/crispr-ga_nextflow/crispr-ga_umis_subs.nf -c nextflow_input.config > NF_output.txt;

awk -v date="$(date)" -v pwd="" '$0 ~ "^Completed" && !found {print "Time: "date"\t"pwd"\tCompleted"; found++ } END { if (!found) print "Time: "date"\t"pwd"\tNotCompleted";}' NF_output.txt >> /scratch/lab_mguell/shared_data/crispr-ga_nextflow/runs.log
cp /scratch/lab_mguell/shared_data/crispr-ga_nextflow/upload_second-part.php index.php
touch end/

