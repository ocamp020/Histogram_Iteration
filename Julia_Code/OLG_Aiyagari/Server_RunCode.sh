#!/bin/bash

#SBATCH --job-name=HistIterOLG
#SBATCH --account=def-robin370
#SBATCH --mail-type=ALL                  
#SBATCH --mail-user=brobin63@uwo.ca
#SBATCH --time=0-23:55:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=SIM-Main.out 
#SBATCH --error=SIM-Error.err

module load julia/1.7.0

julia Main_Solve_OLG_Model.jl


module unload julia/1.7.0