import numpy as np
from celmech.nbody_simulation_utilities import get_simarchive_integration_results
from glob import glob
#from matplotlib import pyplot as plt
import rebound as rb
from celmech.miscellaneous import linking_l
def process_sa(sa,Nout):
    target_times = np.linspace(sa.tmin,sa.tmax,Nout)
    times = np.zeros(Nout)
    energy = np.zeros(Nout)
    a,e,inc,links=np.zeros((4,Nout,5))
    for i,ttarget in enumerate(target_times):
        sim = sa.getSimulation(ttarget)
        orbits = sim.orbits(primary=sim.particles[0])
        a[i]=[o.a for o in orbits]
        e[i]=[o.e for o in orbits]
        inc[i]=[o.inc for o in orbits]
        times[i] = sim.t
        energy[i]= sim.energy()
        for j in range(5):
            if orbits[j].a<0:
                continue
            for k in range(j+1,5):
                if orbits[k].a<0:
                    continue
                if linking_l(orbits[j],orbits[k])<0:
                    links[i,j]+=1
                    links[i,k]+=1
    return {'times':times,'a':a,'e':e,'inc':inc,'links':links,'energy':energy}

if __name__=="__main__":
    import sys
    i = int(sys.argv[1])
    srcdir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
    for planet_type in ("neptune","earth","saturn","3neptune"):
        fi = srcdir + "five_{}_sim_{}.sa".format(planet_type,i)
        sa = rb.Simulationarchive(fi)
        results = process_sa(sa,2048)
        outfile = srcdir + "five_{}_sim_{}_result.npz".format(planet_type,i)
        np.savez_compressed(outfile,**results)
