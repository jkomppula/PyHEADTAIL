import subprocess
print "start"
# The first argument is a number of processors given to the mpirun and the second one is the case index


algoritmhs = [
        'linear_mpi_optimized_fft',
        'memory_optimized',
#        'linear_mpi_full_ring_fft'
        ]

n_particles = 10000
n_turns = 150
n_slices = 100

job_counter = 0

for algoritmh in algoritmhs:
            subprocess.call("./run_local_job.sh 4 " + str(job_counter) +  " " + algoritmh +  " " + str(n_particles) +  " "  + str(n_slices)+  " " + str(n_turns), shell=True)
            job_counter += 1
print "end"
