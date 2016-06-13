#!/bin/bash
#PBS -r n
#PBS -N remaking
#PBS -o remaking.o
#PBS -e remaking.err
#PBS -m a
#PBS -M cap10@ufl.edu
#PBS -l walltime=2:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=2gb
#PBS -t 1-17

cd /lfs/scratch/cap10/montreal-detect

tar=`cat notes_remaking | tail -n +$PBS_ARRAYID | head -1`

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