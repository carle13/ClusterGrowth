#!/bin/bash
#SBATCH --ntasks=40                # Nombre total de processus MPI
#SBATCH --ntasks-per-node=40       # Nombre de processus MPI par noeud
#SBATCH --hint=nomultithread       # 1 processus MPI par coeur physique (pas d'hyperthreading)
#SBATCH -A nyu@cpu
#SBATCH --time=20:00:00            # Temps d’exécution maximum demande (HH:MM:SS)


module purge
module load lammps/20210929-mpi

echo $tempG
echo $randomSeed

srun /gpfswork/rech/nyu/uvm82kt/.Software/1_Lammps_withLassoLars_BIS/src/lmp_mpi -i input.lmp -log seed${randomSeed}.log -var seed ${randomSeed} -var T ${tempG}
