#!/bin/bash
#SBATCH -J PELE
#SBATCH --output=${system}.out
#SBATCH --error=${system}.err
#SBATCH --ntasks=$ncpus

module purge
export PELE="${pele_path}"
export SCHRODINGER="${schrodinger_path}"
export PATH=${conda_path}:$$PATH
module load intel mkl impi gcc # 2> /dev/null
module load boost/1.64.0
${conda_path}/python3.8 -m pele_platform.main $system_yml
