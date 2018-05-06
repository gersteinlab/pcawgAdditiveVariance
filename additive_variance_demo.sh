#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=additive_variance_pipeline
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=10000 
#SBATCH --time=100:00:00

matlab -nodisplay -nosplash -nojvm < additive_variance_demo.m > additive_variance_demo.log

