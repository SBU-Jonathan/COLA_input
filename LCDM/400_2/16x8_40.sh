#!/bin/bash
#SBATCH --job-name=5D_emu2
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=8
#SBATCH --output=/gpfs/projects/MirandaGroup/jonathan/cola_projects/5D_emulator2/slurm_files/out/run_%a.out
#SBATCH --error=/gpfs/projects/MirandaGroup/jonathan/cola_projects/5D_emulator2/slurm_files/err/run_%a.err
#SBATCH --partition=medium-40core
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=1

colasolver=/gpfs/projects/MirandaGroup/jonathan/FML/FML/COLASolver/nbody
param_file=/gpfs/projects/MirandaGroup/jonathan/cola_projects/5D_emulator2/lua_files/parameter_file${SLURM_ARRAY_TASK_ID}.lua
module purge > /dev/null 2>&1
module load slurm/seawulf3/21.08.8
source /gpfs/home/jsgordon/miniconda/etc/profile.d/conda.sh
echo Running on host `hostname`
echo Job started at `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
echo Number of task is $SLURM_NTASKS
echo Number of cpus per task is $SLURM_CPUS_PER_TASK

cd $SLURM_SUBMIT_DIR
conda activate cola
source start_cola

export OMP_PROC_BIND=close
export OMP_NUM_THREADS=1

mpirun -n ${SLURM_NTASKS} --report-bindings --mca btl tcp,self --bind-to core --map-by numa:pe=${OMP_NUM_THREADS} $colasolver $param_file

echo Job ended at `date`