#!/bin/bash
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=1
#SBATCH -t 95:59:00
#SBATCH --mem=10GB
# sends mail when process begins, and 
# when it ends. Make sure you define your email 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=yarmola@princeton.edu

/u/yarmola/momsearch/bin/simple_tree_cat -r data_5.24 '' | ./verify_cusp '' > test.log
