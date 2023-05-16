#!/bin/bash
#SBATCH --job-name=run
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --output=/gpfs/projects/MirandaGroup/jonathan/cola_projects/LCDM/a/slurm_files/out/run_%a.out
#SBATCH --error=/gpfs/projects/MirandaGroup/jonathan/cola_projects/LCDM/a/slurm_files/err/run_%a.err
#SBATCH --partition=long-24core
#SBATCH -t 4:00:00

colasolver=/gpfs/projects/MirandaGroup/jonathan/FML/FML/COLASolver/nbody
param_file=/gpfs/projects/MirandaGroup/jonathan/cola_projects/LCDM/a/lua_files/parameter_file${SLURM_ARRAY_TASK_ID}.lua
module purge > /dev/null 2>&1
module load slurm
source /gpfs/home/jsgordon/miniconda/etc/profile.d/conda.sh
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Number of task is $SLURM_NTASKS
echo Number of cpus per task is $SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR
conda activate cola
source start_cola

export OMP_PROC_BIND=close
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
  export OMP_NUM_THREADS=1
fi


mpirun -n ${SLURM_NTASKS} --report-bindings --mca btl tcp,self --bind-to core --map-by numa:pe=${OMP_NUM_THREADS} $colasolver $param_file