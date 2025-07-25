#!/bin/bash
#SBATCH --output=/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/correlation_function_logs/displacement_shear_correlation_%j.out
#SBATCH --error=/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/correlation_function_logs/displacement_shear_correlation_%j.err
#SBATCH --partition=comp,pscomp
#SBATCH --job-name=displacement_shear_correlation
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00

# Change to script directory
cd /home/murray/intrinsic_alignments/scripts

# Activate conda environment
source ~/.bashrc
conda activate calum_conda
echo "Conda activated"
echo "Working directory: $(pwd)"
echo "Python path: $(which python)"

# Specify config file as a variable for easy change
CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_observed_displacement_shear.ini
#CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_observed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_reconstructed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_rsd_removed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_reconstructed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_rsd_removed.ini

python run_displacement_shear.py --config_file $CONFIG_FILE

exit 0