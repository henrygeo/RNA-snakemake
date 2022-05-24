#!/bin/bash
#SBATCH --time=00:60:00
#SBATCH --mem=8000
#SBATCH --account=def-jstinchc
#SBATCH --mail-user=georgia.henry@mail.utoronto.ca
#SBATCH --mail-type=FAIL

## Loading modules required for script commands
module load python
module load scipy-stack

source ENV/bin/activate #activate the virtual python environment
pip install multiqc --no-index #load in from the wheel

## Running MultiQC on FastQC files
multiqc /scratch/henrygeo/Results --module fastqc --module star