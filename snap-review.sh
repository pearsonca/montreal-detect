cat <<EOF
#!/bin/bash
#SBATCH -r n
#SBATCH --job-name=$1
#SBATCH -o $1.o
#SBATCH -e $1.err
#SBATCH -m a
#SBATCH --mail-user=cap10@ufl.edu
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=1
#SBATCH -N1
#SBATCH --mem-per-cpu=2gb

module load gcc/5.2.0 R/3.2.2
cd /ufrc/singer/cap10/montreal-detect
make input/detection/$2/snapFTPR.rds
EOF
