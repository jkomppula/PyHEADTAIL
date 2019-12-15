#!/bin/bash
# A shell which is used to launch PyHEADTAIL simulations in a MPI environment 
# from a Jupyter Notebook

nnodes=$1
job_id=$2
algorithm=$3
n_particles=$4
n_slices=$5
n_turns=$6
# Path to the installed environment
# (same as ENVPATH in the file install_mpi_environment.sh)
# ENVIRONMENT=/home/jani/environments/mpi_dev

# source $ENVIRONMENT/virtualenvs/py2.7/bin/activate
# which python
mpirun -np $nnodes python main.py $job_id $algorithm $n_particles $n_slices $n_turns
