#!/bin/bash --login

#SBATCH --job-name=EDTA
#SBATCH --time=167:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem=300G
#SBATCH --output=genome_EDTA.SLURMout
#SBATCH -o "stdout.%x.%j.%N"
#SBATCH -e "stderr.%x.%j.%N"

# You can use a conda env, module load it, or whatever works. This is just my preferred method.
# You can also use EDTA1.
# You can nohup this too, but it can take awhile to run.

# Pull EDTA2 container:
#singularity pull edta_2.2.2.sif docker://quay.io/biocontainers/edta:2.2.2--hdfd78af_1

cd /path/to/working_directory

singularity exec edta_2.2.2.sif EDTA.pl --genome genome.fasta --anno 1 -t 80
