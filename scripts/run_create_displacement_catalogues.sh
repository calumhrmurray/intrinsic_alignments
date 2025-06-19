#!/bin/zsh
#SBATCH --output=/n17data/murray/desi_data/DESI/correlation_function_logs/create_displacement_catalogues.out
#SBATCH --error=/n17data/murray/desi_data/DESI/correlation_function_logs/create_displacement_catalogues.err
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
CONFIG_FILE=/home/murray/intrinsic_alignments/scripts/config_files/bgs_displacement_fields.ini

python create_displacement_fields.py --config_file $CONFIG_FILE

exit 0