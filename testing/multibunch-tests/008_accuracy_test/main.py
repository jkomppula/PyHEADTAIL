from __future__ import division

import sys, os, time
import numpy as np
from mpi4py import MPI
from scipy.constants import c, e, m_p


BIN = os.path.expanduser("../../../../")
sys.path.append(BIN)


def run(argv):
    
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    np.random.seed(0+rank) # different macroparticle coordinates on each rank
    
    job_id = int(argv[0])
    algorithm = argv[1]
    n_particles = int(argv[2])
    n_slices = int(argv[3])
    n_turns = int(argv[4])
    
    if algorithm == 'memory_optimized':
        mpi_settings = algorithm
    elif algorithm == 'linear_mpi_full_ring_fft':
        mpi_settings = algorithm
    elif algorithm == 'linear_mpi_optimized_fft':
        mpi_settings = algorithm
    else:
        raise ValueError('Unknown algorithm.')
    
    
    
    BIN = os.path.expanduser("../")
    sys.path.append(BIN)
    
    from PyHEADTAIL.particles.slicing import UniformBinSlicer, UniformChargeSlicer
    from PyHEADTAIL.impedances.wakes import CircularResonator, WakeField
    from PyHEADTAIL.monitors.monitors import BunchMonitor, SliceMonitor
    from PyHEADTAIL.machines.synchrotron import Synchrotron
    
    chroma = 3


    # Random number seed must be defined for the bunch-by-bunch trace comparison
    # but different seed must be used in each rank in order to avoid weird coherent
    # oscillations at the beginning.

#    np.random.seed(0)
    
    print('I am rank ' + str(rank) + ' out of ' + str(size))
    print('  ')

    # BEAM AND MACHNINE PARAMETERS
    # ============================

    n_macroparticles = n_particles
    intensity = 2.3e11

    charge = e
    mass = m_p
    alpha = 53.86**-2

    p0 = 7000e9 * e / c

    accQ_x = 62.31
    accQ_y = 60.32
    Q_s = 2.1e-9
    circumference = 26658.883
    s = None
    alpha_x = None
    alpha_y = None
    beta_x = circumference / (2.*np.pi*accQ_x)
    beta_y = circumference / (2.*np.pi*accQ_y)
    D_x = 0
    D_y = 0
    optics_mode = 'smooth'
    name = None
    n_segments = 1

    # detunings
    Qp_x = chroma
    Qp_y = chroma


    app_x = 0
    app_y = 0
    app_xy = 0

    longitudinal_mode = 'linear'

# Number of bunches simulated
    h_RF = 3564


    wrap_z = False

    machine = Synchrotron(
            optics_mode=optics_mode, circumference=circumference,
            n_segments=n_segments, s=s, name=name,
            alpha_x=alpha_x, beta_x=beta_x, D_x=D_x,
            alpha_y=alpha_y, beta_y=beta_y, D_y=D_y,
            accQ_x=accQ_x, accQ_y=accQ_y, Qp_x=Qp_x, Qp_y=Qp_y,
            app_x=app_x, app_y=app_y, app_xy=app_xy,
            alpha_mom_compaction=alpha, longitudinal_mode=longitudinal_mode,
            h_RF=np.atleast_1d(h_RF), p0=p0,
            charge=charge, mass=mass, wrap_z=wrap_z, Q_s=Q_s)


    
    filling_scheme = []
    for i in range(2):
        for j in range(72):
            filling_scheme.append(i*80 + j)

    # BEAM
    # ====
    epsn_x = 2e-6
    epsn_y = 2e-6
    sigma_z = 0.09
    allbunches = machine.generate_6D_Gaussian_bunch(
        n_macroparticles, intensity, epsn_x, epsn_y, sigma_z=sigma_z,
        filling_scheme=filling_scheme, matched=False)
    

    # CREATE BEAM SLICERS
    # ===================
    slicer_for_wakefields = UniformBinSlicer(n_slices, z_cuts=(-4.*sigma_z, 4.*sigma_z),
                               circumference=machine.circumference, h_bunch=h_RF)

    # WAKE PARAMETERS
    # ============================
    
    # Remember: n_turns_wake wake must be several times the decay time of the 
    # the wake function if slice-by-slice locations will be comparared on 
    # the floating point rounding error level 
    
    
    n_turns_wake = 4
    R_shunt = 135e6 * 40
    frequency = 1.97e9
    Q = 31000/100.
    wakes = CircularResonator(R_shunt, frequency, Q, n_turns_wake=n_turns_wake)

    wake_field = WakeField(slicer_for_wakefields, wakes, mpi=mpi_settings,
                           Q_x=accQ_x, Q_y=accQ_y, beta_x=beta_x, beta_y=beta_y)


#    if case == 1:
    machine.one_turn_map.append(wake_field)

    # SIMULATION PARAMETERS
    # ============================
#    n_turns = 100



    tracking_data = np.zeros((n_turns,4))
    

    # TRACKING LOOP
    # =============

    header_string = str(n_slices) + " slices, job " + str(job_id) + ": " + mpi_settings + " algorithm, " 
    for i in range(n_turns):
        t0 = time.clock()

        machine.track(allbunches)
        tracking_data[i,0] = i
        tracking_data[i,1] = time.clock() - t0
        tracking_data[i,2] = allbunches.mean_x()
        tracking_data[i,3] = allbunches.mean_y()
        

        if (rank == 0) and (i%10 == 0):
            if i == 0:
                print(header_string)
    
            bunch_list = allbunches.split_to_views()
            bunch = bunch_list[0]
            print 'Bunch: {:4d} \t {:+3e} \t {:+3e} \t {:+3e} \t {:3e} \t {:3e} \t {:3f} \t {:3f} \t {:3f} \t {:3s}'.format(i, bunch.mean_x(), bunch.mean_y(), bunch.mean_z(), bunch.epsn_x(), bunch.epsn_y(), bunch.epsn_z(), bunch.sigma_z(), bunch.sigma_dp(), str(time.clock() - t0))


    bunch_list = allbunches.split_to_views()
    output_data = np.zeros((len(bunch_list)*n_slices, 5))
    
    for i, b in enumerate(bunch_list):
        s = b.get_slices(slicer_for_wakefields)
        
        i_from = i*n_slices
        i_to = (i+1)*n_slices
        
        output_data[i_from:i_to,0] = np.ones(n_slices) * i
        output_data[i_from:i_to,1] = s.mean_x
        output_data[i_from:i_to,2] = s.mean_y
#        output_data[i_from:i_to,3] = s.z
        output_data[i_from:i_to,4] = s.n_macroparticles_per_slice
        
        
    
    

    file_prefix = "data/" + str(mpi_settings) + "_rank_" + str(rank) + '__'
    
    np.savetxt(file_prefix + 'slice_by_slice.dat',output_data, header=header_string)
    
    if rank == 0:
        print '\n*** Successfully completed!'


if __name__=="__main__":
	run(sys.argv[1:])
