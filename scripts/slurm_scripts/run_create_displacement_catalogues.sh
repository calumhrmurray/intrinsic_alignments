#!/bin/zsh
#SBATCH --output=/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/correlation_function_logs/create_displacement_catalogues.out
#SBATCH --error=/n17../results/correlation_functions/murray/desi_../results/correlation_functions/DESI/correlation_function_logs/create_displacement_catalogues.err
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
CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/lrg_displacement_fields.ini

python create_displacement_fields.py --config_file $CONFIG_FILE

exit 0