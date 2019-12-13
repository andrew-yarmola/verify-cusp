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


base_dir="/u/yarmola/verify-cusp"
bin_dir="$base_dir/bin"
data_dir="/scratch/network/yarmola/run_cleanup/tree_tar"

cd $bin_dir

./rootcat "$data_dir" | ./verify_cusp '' > "$base_dir/test.log"
