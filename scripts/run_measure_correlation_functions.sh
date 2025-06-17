#!/bin/zsh
#SBATCH --output=/n17data/murray/desi_data/DESI/correlation_function_logs/run_correlation_function.out
#SBATCH --error=/n17data/murray/desi_data/DESI/correlation_function_logs/run_correlation_function.err
#SBATCH --partition=comp,pscomp
#SBATCH --job-name=run_correlation_function
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=4:00:00

python run_correlation_functions.py --config_file /home/murray/intrinsic_alignments/scripts/config_files/bgs_recon.ini

exit 0