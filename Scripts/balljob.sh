#!/bin/bash
#SBATCH --account=def-jstinchc   
#SBATCH --mem=16000      
#SBATCH --time=00:30:00
module load gcc/9.3.0 r/4.1.2              # Adjust version and add the gcc module used for installing packages.
Rscript /scratch/henrygeo/Scripts/ballg.R