#!/bin/zsh
#SBATCH --output=/n09data/guerrini/glass_mock/run_glass_mock_10.out
#SBATCH --error=/n09data/guerrini/glass_mock/run_glass_mock_10.err
#SBATCH --partition=comp,pscomp
#SBATCH --job-name=run_glass_mock
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=48:00:00

module purge
module load intelpython

module load openmpi

source activate glass

cd ~/sp_validation/glass_mock

python make_unions_glass_sim.py -p /n09data/guerrini/glass_mock/ -n 4096 -N 10

#cosmosis cosmosis_config/cosmosis_pipeline_SP_v1.3_LFmask_8k.ini


exit 0