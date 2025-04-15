import numpy as np
import rebound as rb
from celmech import Poincare
from celmech.miscellaneous import critical_relative_AMD
from celmech.nbody_simulation_utilities import align_simulation
from celmech.nbody_simulation_utilities import get_canonical_heliocentric_orbits

def bound_only_sim_copy(sim):
    copy_sim = rb.Simulation()
    copy_sim.units=sim.units
    orbits = get_canonical_heliocentric_orbits(sim)
    copy_sim.add(sim.particles[0].copy())
    for o,p in zip(orbits,sim.particles[1:]):
        if o.a>0:
            copy_sim.add(p.copy())
    copy_sim.move_to_com()
    align_simulation(copy_sim)
    return copy_sim,np.array([o.a>0 for o in orbits])

def calculate_AMDs_and_AMD_crits(sim):
    pvars = Poincare.from_Simulation(sim)
    amds = np.array([p.Gamma + p.Q for p in pvars.particles[1:]])
    i_inner = np.argmin([p.a for p in pvars.particles[1:]])
    inner_planet = pvars.particles[i_inner+1]    
    AMD_crits = np.array([
        (p.Lambda * critical_relative_AMD(inner_planet.a/p.a,inner_planet.m/p.m) if i != i_inner else 0)
        for i,p in enumerate(pvars.particles[1:]) 
    ])
    return amds, AMD_crits
if __name__=="__main__":
    savedir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
    import sys
    SIMTYPE=int(sys.argv[1])
    sim_family = ("five_neptune_sim","five_3neptune_sim","five_saturn_sim","five_30neptune_sim","ten_neptunes_plus_jupiter_sim")[SIMTYPE]
    for I in range(100):
        print(I)
        archive_file = savedir + "{}_{}.sa".format(sim_family,I)
        sa = rb.Simulationarchive(archive_file)
        sim0 = sa[0]
        Npl = sim0.N-1
        times = np.geomspace(1e5,sa.tmax,600)    
        AMD_crit = np.zeros((times.size,Npl))
        AMD_pp = np.zeros((times.size,Npl))
        for i, sim in enumerate(sa.getSimulations(times)):
            sim_bound,bound_mask = bound_only_sim_copy(sim)
            AMD_pp[i,bound_mask],AMD_crit[i,bound_mask] = calculate_AMDs_and_AMD_crits(sim_bound)

        np.savez_compressed(savedir+f"amd_data_{sim_family}_{I}",AMD_pp = AMD_pp, AMD_crit = AMD_crit)