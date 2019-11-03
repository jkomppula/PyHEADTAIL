import subprocess
print "start"
# The first argument is a number of processors given to the mpirun and the second one is the case index


n_slices = [
#        25,
#        50,
        100
        ]

filling_schemes = [
        'train',
        'random'
]

algoritmhs = [
        'memory_optimized',
        'linear_mpi_optimized_fft',
        'linear_mpi_full_ring_fft'
        ]

n_bunches = [
        4,
        6,
        8,
        12,
        16,
        24,
        32,
        46,
        64,
        90,
        128,
        180,
        256,
        364,
        512,
        724,
        1024,
        1448,
        2048
        ]

for n_s in n_slices:
    job_counter = 0 
    for algoritmh in algoritmhs:
        for filling_scheme in filling_schemes:
            for n_b in n_bunches:
                subprocess.call("./run_local_job.sh 4 " + str(job_counter) +  " " + algoritmh +  " " + filling_scheme +  " " + str(n_b) +  " " + str(n_s), shell=True)
                job_counter += 1
                if (algoritmh == 'memory_optimized') and (n_b == 724):
                    break
                
                if ((n_s==100) and (algoritmh == 'memory_optimized')) and (n_b == 256):
                    break

print "end"