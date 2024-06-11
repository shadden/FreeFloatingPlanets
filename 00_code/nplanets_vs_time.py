import numpy as np
import rebound as rb
import sys
savedir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
def count_bound(sim):
    orbits = sim.orbits(sim.particles[0])
    return np.sum([o.a>0 for o in orbits])

planet_type = sys.argv[1]
I = int(sys.argv[2])
archive_file = savedir + "five_{}_sim_{}.sa".format(planet_type,I)
sa = rb.Simulationarchive(archive_file)
Nlast = 5
E0 = sa[0].energy()
times,counts,dE = [],[],[]
for sim in sa:
    N = count_bound(sim)
    if N != Nlast:
        times.append(sim.t)
        counts.append(N)
        dE.append(np.abs(sim.energy()/E0 - 1))
    Nlast = N
counts_and_dE_vs_time = np.transpose((times,counts,dE))
# datadir = "/cita/h/home-2/hadden/Projects/15_FreeFloatingPlanetProduction/03_data/"
# np.save(datadir+"{}_{}_final_orbits".format(planet_type,I),orbits_arr)