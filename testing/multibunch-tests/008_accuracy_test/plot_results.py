import subprocess
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
print "start"
# The first argument is a number of processors given to the mpirun and the second one is the case index


circumference = 26658.883
h_RF = 3564

sigma_z = 0.09
n_slices = 100
slices = np.linspace(-4.*sigma_z, 4.*sigma_z, n_slices)

filling_scheme = []
for i in range(2):
    for j in range(72):
        filling_scheme.append(i*80 + j)


n_s = 100

d1_1 = np.loadtxt('data/linear_mpi_optimized_fft_rank_0__slice_by_slice.dat')
d1_2 = np.loadtxt('data/linear_mpi_optimized_fft_rank_1__slice_by_slice.dat')
d1_3 = np.loadtxt('data/linear_mpi_optimized_fft_rank_2__slice_by_slice.dat')
d1_4 = np.loadtxt('data/linear_mpi_optimized_fft_rank_3__slice_by_slice.dat')

d2_1 = np.loadtxt('data/memory_optimized_rank_0__slice_by_slice.dat')
d2_2 = np.loadtxt('data/memory_optimized_rank_1__slice_by_slice.dat')
d2_3 = np.loadtxt('data/memory_optimized_rank_2__slice_by_slice.dat')
d2_4 = np.loadtxt('data/memory_optimized_rank_3__slice_by_slice.dat')

d1 = d1_1[:,1]*d1_1[:,4]
d1 = np.append(d1, d1_2[:,1]*d1_2[:,4])
d1 = np.append(d1, d1_3[:,1]*d1_3[:,4])
d1 = np.append(d1, d1_4[:,1]*d1_4[:,4])

d2 = d2_1[:,1]*d2_1[:,4]
d2 = np.append(d2, d2_2[:,1]*d2_2[:,4])
d2 = np.append(d2, d2_3[:,1]*d2_3[:,4])
d2 = np.append(d2, d2_4[:,1]*d2_4[:,4])


#
#d2_1 = np.loadtxt('data/linear_mpi_full_ring_fft_rank_0__slice_by_slice.dat')
#d2_2 = np.loadtxt('data/linear_mpi_full_ring_fft_rank_1__slice_by_slice.dat')
#d2_3 = np.loadtxt('data/linear_mpi_full_ring_fft_rank_2__slice_by_slice.dat')
#d2_4 = np.loadtxt('data/linear_mpi_full_ring_fft_rank_3__slice_by_slice.dat')
fig = plt.figure(constrained_layout=True)
gs = GridSpec(2, 3, figure=fig)
# identical to ax1 = plt.subplot(gs.new_subplotspec((0, 0), colspan=3))
ax1 = fig.add_subplot(gs[0, :2])
ax2 = fig.add_subplot(gs[1, :2])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[1, 2])

for i, bunch in enumerate(filling_scheme):
    i_from = (len(filling_scheme)-i-1)*n_slices 
    i_to = ((len(filling_scheme)-i-1)+1)*n_slices
    z = slices+bunch*circumference/float(h_RF)
    if i == 0:
        ax1.plot(z,d1[i_from:i_to], 'C0-', label="Old")
        ax1.plot(z,d2[i_from:i_to], 'C1:', label="New")
    else:
        ax1.plot(z,d1[i_from:i_to], 'C0-')
        ax1.plot(z,d2[i_from:i_to], 'C1:')
    
    ax2.plot(z,d1[i_from:i_to]-d2[i_from:i_to], 'C2-')
    
    
    if i == 80:
        ax3.plot(z,d1[i_from:i_to], 'C0-')
        ax3.plot(z,d2[i_from:i_to], 'C1:')
        ax4.plot(z,d1[i_from:i_to]-d2[i_from:i_to], 'C2-')


#ax.plot(d1_1[:,1]*d1_1[:,4])
#ax.plot(d2_1[:,1]*d2_1[:,4])

#ax.plot(d1_1[:,4])
#ax.plot(d2_1[:,4])
#
#ax.loglog(d_normal_train[:,0],d_normal_train[:,1],'C0o-',label="Original algorithm, bunch tain")
#ax.loglog(d_normal_random[:,0],d_normal_random[:,1],'C0d--',label="Original algorithm, random filling")
#ax.loglog(d_optimized_train[:,0],d_optimized_train[:,1],'C1o-',label="New algorithm, bunch tain")
#ax.loglog(d_optimized_random[:,0],d_optimized_random[:,1],'C1d--',label="New algorithm, random filling")
#ax.loglog(d_normal_train[:,0],n2_line,'k:',label="$\propto x^2$")
#ax.loglog(d_optimized_train[:,0],n_line,'k--',label="$\propto x$")

ax2.set_xlabel('Distance from the first bunch [m]')
#ax4.set_xlabel('Distance from the first bunch [m]')
ax1.set_ylabel('Slice oscillation [m]')
ax2.set_ylabel('Algorithm difference')
#ax.set_ylim(0.1e-2,1000)
#ax.set_xlim(4,2100)
ax1.legend(frameon=False)

plt.show()
#fig.savefig("PyHT_on_laptop__100_slices_per_bunch.png", dpi=300)
print "end"
