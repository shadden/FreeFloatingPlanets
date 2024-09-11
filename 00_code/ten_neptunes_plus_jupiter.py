import rebound as rb
import numpy as np
from celmech.nbody_simulation_utilities import align_simulation

# test 1,2,3
inc_scale = 1e-3

def lammers_inst_law(t,m):
    A = 11.9
    C = 5.2
    log_delta = (np.log10(t * (m/3e-6)) - C) / A
    da_by_a = 10**log_delta * m**(0.25)
    return da_by_a
mNep = 5.15e-5
def get_sim(t):
    sim = rb.Simulation()
    sim.units = ("Msun","AU","yr")
    sim.integrator = 'ias15'
    sim.add(m=1)
    sim.add(m=1e-3,a=1) # add a jupiter
    da_by_a = lammers_inst_law(t,mNep)
    da_by_a_jup = lammers_inst_law(t,1e-3)
    ai = (1+da_by_a_jup)/(1-da_by_a_jup)
    for i in range(10):
        inc = np.random.rayleigh(scale=inc_scale)
        sim.add(m=mNep,a=ai,l='uniform',inc=inc,Omega='uniform')
        ai = ai * (1+da_by_a)/(1-da_by_a)
    sim.move_to_com()
    align_simulation(sim)
    return sim

if __name__=="__main__":
    import sys
    i = int(sys.argv[1])
    srcdir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
    file_name = srcdir + "ten_neptunes_plus_jupiter_sim_{}.sa".format(i)
    print("Opening file {}".format(file_name))
    sim = rb.Simulation(file_name)
    sim.save_to_file(file_name, walltime=25.,delete_file=False)
    print("Integrating...")
    sim.integrate(1e8)
