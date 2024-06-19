import numpy as np
import rebound as rb
import sys
from celmech.miscellaneous import linking_l
savedir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
planet_type = sys.argv[1]
def count_links(sim):
    orbits = sim.orbits(primary=sim.particles[0])
    links = 0
    for i in range(sim.N-1):
        oi = orbits[i]
        for j in range(i+1,sim.N-1):
            oj = orbits[j]
            links += linking_l(oi,oj)<0
    return links

def crossing_time(sa):
    for sim in sa:
        if count_links(sim)>0:
            return sim.t
    else:
        return np.nan
archive_fi_str=savedir+"five_{}_sim_{}.sa"
crossing_times = np.zeros(100)
for i in range(100):
    print(i)
    sa = rb.Simulationarchive(archive_fi_str.format(planet_type,i))
    crossing_times[i] = crossing_time(sa)
np.save(savedir+"{}_crossing_times".format(planet_type),crossing_times)