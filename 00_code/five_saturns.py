import rebound as rb
import numpy as np
from celmech.nbody_simulation_utilities import align_simulation

inc_scale = 1e-3

def lammers_inst_law(t,m):
    A = 11.9
    C = 5.2
    log_delta = (np.log10(t * (m/3e-6)) - C) / A
    da_by_a = 10**log_delta * m**(0.25)
    return da_by_a

def get_sim(t,m):
    sim = rb.Simulation()
    sim.units = ("Msun","AU","yr")
    sim.integrator = 'ias15'
    sim.add(m=1)
    da_by_a = lammers_inst_law(t,m)
    ai = 1
    for i in range(5):
        inc = np.random.rayleigh(scale=inc_scale)
        sim.add(m=m,a=ai,l='uniform',inc=inc,Omega='uniform')
        ai = ai * (1+da_by_a)/(1-da_by_a)
    sim.move_to_com()
    align_simulation(sim)
    return sim

#if __name__=="__main__":
#    import sys
#    i = int(sys.argv[1])
#    print(rb.__version__)
#    T_inst = 10**np.random.uniform(4,5)
#    sim = get_sim(T_inst,5.15e-4)
#    savedir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
#    sim.save_to_file(savedir + "five_saturn_sim_{}.sa".format(i),walltime=10)
#    sim.integrate(1e8)

if __name__=="__main__":
    import sys
    import os
    i = int(sys.argv[1])
    srcdir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
    file_name = srcdir + "five_saturn_sim_{}.sa".format(i)
    if os.path.isfile(file_name):
        print("Opening file {}".format(file_name))
        sim = rb.Simulation(file_name)
    else:
        print("No file named {}".format(file_name))
        print("Creating new simulation...")
        T_inst = 10**np.random.uniform(3.,4.)
        mass = 5.15e-4
        sim = get_sim(T_inst,mass)
    sim.save_to_file(file_name, walltime=30.,delete_file=False)
    print("Integrating...")
    sim.integrate(1e8)
