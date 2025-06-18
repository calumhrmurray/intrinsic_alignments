#!/bin/zsh
#SBATCH --output=/n17data/murray/desi_data/DESI/correlation_function_logs/lrg_run_correlation_function.out
#SBATCH --error=/n17data/murray/desi_data/DESI/correlation_function_logs/lrg_run_correlation_function.err
#SBATCH --partition=comp,pscomp
#SBATCH --job-name=run_correlation_function
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=5:00:00

# Activate conda environment
source ~/.bashrc
conda activate calum_conda

# Specify config file as a variable for easy change
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_observed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_reconstructed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_rsd_removed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_observed.ini
# CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_reconstructed.ini
CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_rsd_removed.ini

python run_correlation_functions.py --config_file $CONFIG_FILE

exit 0