#!/bin/sh

#PBS -A GEO111
#PBS -N CG_direction
#PBS -j oe
#PBS -o job_cg_direction.$PBS_JOBID.o
#PBS -l walltime=1:00:00
#PBS -l nodes=12

NPROC=384

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

iter_old=M16
iter_new=M17

xoutput_tag=XTAG_$iter_new

# dir for last iteration, cg_direction.bp
direction_0="/lustre/atlas/proj-shared/geo111/rawdata/$iter_old/kernels/post_processing/smooth/direction_M16"
# dir for this iteration, cg_direction.bp
direction_1="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_$iter_new/kernel_post/direction_$iter_new/"
# dir for last iteration, kernels_precond.bp
gradient_0="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_$iter_new/kernel_post"
# dir for this iteration, kernels_precond.bp
gradient_1="/lustre/atlas/proj-shared/geo111/Wenjie/kernels_$iter_new/kernel_post"

if [ ! -d $direction_0 ]; then
  echo Error: set $direction_0
  exit
fi
if [ ! -d $gradient_0 ]; then
  echo Error: set $gradient_0
  exit
fi
if [ ! -d $gradient_1 ]; then
  echo Error: set $gradient_1
  exit
fi

if [ ! -d $direction_1 ]; then
  echo MKDIR $direction_1
  mkdir $direction_1
fi

echo submit compute cg direction
aprun -n $NPROC ./xcompute_cg_direction $direction_0 $direction_1 $gradient_0 $gradient_1  > $xoutput_tag
echo done successfully
