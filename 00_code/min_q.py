import numpy as np
import rebound as rb
savedir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
def minq_of_sim(sim):
    orbits = sim.orbits(sim.particles[0])
    minq = np.inf
    for o in orbits:
        if o.a > 0:
            q = o.a * (1 - o.e)
            minq = np.min((minq,q))
    return minq

import sys
planet_type = sys.argv[1]
I = int(sys.argv[2])

archive_file = savedir + "five_{}_sim_{}.sa".format(planet_type,I)
sa = rb.Simulationarchive(archive_file)
times,minq,energy = np.transpose([(sim.t,minq_of_sim(sim),sim.energy()) for sim in sa])
abs_dE = np.abs(energy[1:]-energy[:-1])
max_dEbyE = np.max(abs_dE/np.abs(energy[0]))
min_q = np.min(minq)
datadir = "/cita/h/home-2/hadden/Projects/15_FreeFloatingPlanetProduction/03_data/"
np.save(datadir+"{}_{}_tfin_max_dE_min_q".format(planet_type,I),np.array((times[-1],max_dEbyE,min_q)))