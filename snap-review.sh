cat <<EOF
#!/bin/bash
#PBS -r n
#PBS -N $1
#PBS -o $1.o
#PBS -e $1.err
#PBS -m a
#PBS -M cap10@ufl.edu
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2gb

module load gcc/5.2.0 R/3.2.2
cd /scratch/lfs/cap10/montreal-detect
make input/detection/$2/snapFTPR.rds
EOF
