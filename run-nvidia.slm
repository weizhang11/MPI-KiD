#!/bin/sh
#SBATCH --account=NWP501
#SBATCH --cluster-constraint=green
#SBATCH --output=outnvidia.rpt
#SBATCH --clusters=miller
#SBATCH --partition=ampere
#SBATCH --time=00:40:00
#SBATCH --nodes=1
##SBATCH --ntasks=8
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
##SBATCH --gpus-per-node=4
#SBATCH --mail-type=ALL

date

srun ./bin/KiD_CU_2D.exe namelists/CU_2D.nml output/CU_2D_nvidia_org.nc

ncdump -v mean_cloud_mass_path ./output/CU_2D_nvidia_org.nc

