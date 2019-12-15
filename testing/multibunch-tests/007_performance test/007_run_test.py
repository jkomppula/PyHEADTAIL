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
        'linear_mpi_optimized_fft',
        'memory_optimized',
#        'linear_mpi_full_ring_fft'
        ]

n_bunches = [
        4,
        6,
        8,
        12,
        16,
        24,
        32,
        44,
        64,
        92,
        128,
        180,
        256,
        364,
        512,
        724,
        1024,
        1448,
        2048,
        2896
        ]

n_bunches = n_bunches[::-1]

for n_s in n_slices:
    job_counter = 0 
    for algoritmh in algoritmhs:
        for filling_scheme in filling_schemes:
            for n_b in n_bunches:
                if (algoritmh == 'linear_mpi_optimized_fft') or ((n_s==100) and (algoritmh == 'memory_optimized') and (n_b <= 512)):
                    n_turns = 100
                else:
                    n_turns = 12
                subprocess.call("./run_local_job.sh 4 " + str(job_counter) +  " " + algoritmh +  " " + filling_scheme +  " " + str(n_b) +  " " + str(n_s)+  " " + str(n_turns), shell=True)
                job_counter += 1
print "end"