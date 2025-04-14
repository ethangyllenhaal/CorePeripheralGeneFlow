#!/bin/bash

#SBATCH --ntasks=20
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --job-name=cont_msp
#SBATCH --output=_trees_out_%j
#SBATCH --error=_trees_error_%j
#SBATCH --partition=normal
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=address@email.edu

# Load miniconda and parallel modules
module load miniconda3
module load parallel

# Activate the environment and load env_parallel
source activate msptree-env
source $(which env_parallel.bash)

# Set a directory path to working directory, then make node list for parallelizing across nodes
dir=$SLURM_SUBMIT_DIR
scontrol show hostname > $dir/node_list_${SLURM_JOB_ID}

export LANG=C

# runs simulations 8 per node
# first runs python script, then makes a variant only fasta, and finally deletes original fasta
env_parallel -j 8 --delay 1 --joblog $dir/msp_log --sshloginfile $dir/node_list_${SLURM_JOB_ID} \
             'python $dir/core_phylo_noanc.py -r {6} -m {1} -n {2} -e {3} -i {4} -s {5} \
                     -o $dir/output/sim/coresim_m{1}_n{2}_e{3}_i{4}_{5}_r{6}
              sh $dir/variant_fasta.sh $dir/output/sim/coresim_m{1}_n{2}_e{3}_i{4}_{5}_r{6}.fasta \
                 $dir/output/sim/coresim_m{1}_n{2}_e{3}_i{4}_{5}_r{6}_var.fasta
              rm $dir/output/sim/coresim_m{1}_n{2}_e{3}_i{4}_{5}_r{6}.fasta' \
              ::: 0 0.1 0.2 0.4 0.8 ::: 50000 100000 ::: 200000 400000 ::: 50000 100000 200000 ::: CAB ACB ABC ::: {1..20}

# runs trees 2 per node with 4 threads
# --redo just in case run is interupted
env_parallel -j 2 --delay 1 --joblog $dir/iqtree_log --sshloginfile $dir/node_list_${SLURM_JOB_ID} \
	     'iqtree -s $dir/output/sim/coresim_m{1}_n{2}_e{3}_i{4}_{5}_r{6}_var.fasta \
                     --prefix $dir/output/tree/coresim_m{1}_n{2}_e{3}_i{4}_{5}_r{6} \
                     -m JC -B 1000 -T 4 --redo' \
              ::: 0 0.1 0.2 0.4 0.8 ::: 50000 100000 ::: 200000 400000 ::: 50000 100000 200000 ::: CAB ACB ABC ::: {1..20}
