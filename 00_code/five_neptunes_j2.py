import rebound as rb
import reboundx as rbx
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

def get_sim(t,m,J2,Req):
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
    extras = rbx.Extras(sim)
    gh = extras.load_force("gravitational_harmonics")
    extras.add_force(gh)
    sim.particles[0].params["J2"] = J2
    sim.particles[0].params["R_eq"] = Req

    return sim,extras

if __name__=="__main__":
    import sys
    i = int(sys.argv[1])
    #T_inst = 10**np.random.uniform(4,5)
    srcdir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
    file_name = srcdir + "j2_five_neptune_sim_{}.sa".format(i)
    extras_file_name = srcdir + "j2_five_neptune_sim_{}_extras.bin".format(i)

    try:
        print("Opening file {}".format(file_name))
        sim = rb.Simulation(file_name)
        extras = rbx.Extras(sim,extras_file_name)
    except:
        print("No file named '{}'...".format(file_name))
        print("creating new simulation")
        ain = 0.1 / 10. # planet at 0.1 AU relative to 10 AU Neptune
        mIn = 3e-4
        Req = ain
        J2 = 0.5 * mIn * (ain/Req)**2
        sim, extras = get_sim(3e4,5.15e-5,J2,Req)
        extras.save(extras_file_name)
    sim.save_to_file(file_name, walltime=25.,delete_file=False)
    print("Integrating...")
    sim.integrate(1e8)
