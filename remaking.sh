#!/bin/bash
#SBATCH -r n
#SBATCH --job-name=remaking
#SBATCH -o remaking.o
#SBATCH -e remaking.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=cap10@ufl.edu
#SBATCH -t=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH -N 1
#SBATCH --mem-per-cpu=2gb
#SBATCH --array=1-17

module load gcc/5.2.0 R/3.2.2
cd /lfs/scratch/cap10/montreal-detect

tar=`cat notes_remaking | tail -n +$SLURM_ARRAY_TASK_ID | head -1`

rm $tar

func() {
  resval="$(make $tar 2>&1 > /dev/null)"
  testval=$?
}

func

while [ $testval != 0 ]
do
rm $tar
basetar=`echo $resval | grep -Po ": \./input/.+/base\.rds"`
basetar=${basetar:4}
rm $basetar
make $basetar
func
done