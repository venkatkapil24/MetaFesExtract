#!//bin/bash

#$ -N META-FES-EXTRACT

#$ -S /bin/bash

#$ -cwd

#$ -e sge-OUTPUTFILES/$JOB_NAME.err

#$ -o sge-OUTPUTFILES/$JOB_NAME.out

#$ -q new.q

#$ -pe mpi 8

unset SGE_ROOT
export F_UFMTENDIAN=big

/opt/mpi/openmpi/1.3.3/intel/bin/mpirun  -np $NSLOTS ./META-FES-EXTRACT.x
