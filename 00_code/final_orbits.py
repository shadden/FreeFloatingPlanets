import numpy as np
import rebound as rb
savedir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
def orbits_of_sim(sim):
    orbits = sim.orbits(sim.particles[0])
    orbit_data = np.zeros((len(orbits),6))
    for i,o in enumerate(orbits):
        orbit_data[i] = o.a,o.e,o.inc,o.M,o.omega,o.Omega
    return orbit_data
import sys
planet_type = sys.argv[1]
I = int(sys.argv[2])

archive_file = savedir + "five_{}_sim_{}.sa".format(planet_type,I)
sim = rb.Simulation(archive_file)
datadir = savedir
orbits_arr = orbits_of_sim(sim)
np.save(datadir+"{}_{}_final_orbits".format(planet_type,I),orbits_arr)
