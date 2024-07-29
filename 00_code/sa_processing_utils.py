import numpy as np
import rebound as rb
from celmech.miscellaneous import AMD_stability_coefficients, AMD_stable_Q
from celmech.nbody_simulation_utilities import align_simulation
from celmech.miscellaneous import linking_l

def orbits_of_sim(sim):
    """
    get the orbital elements of planets as an array
    """
    orbits = sim.orbits(sim.particles[0])
    orbit_data = np.zeros((len(orbits),6))
    for i,o in enumerate(orbits):
        orbit_data[i] = o.a,o.e,o.inc,o.M,o.omega,o.Omega
    return orbit_data

def amd_stability_data(sim):
    """
    Given a simulation, return the number of bound planets, whether the system is AMD-stable, and the pairwise AMD stability coefficients.
    """
    copy_sim = rb.Simulation()
    copy_sim.units=sim.units
    orbits = sim.orbits(primary=sim.particles[0])
    copy_sim.add(sim.particles[0].copy())
    for o,p in zip(orbits,sim.particles[1:]):
        if o.a>0:
            copy_sim.add(p.copy())
    copy_sim.move_to_com()
    align_simulation(copy_sim)
    Nbound = copy_sim.N-1
    stableQ = AMD_stable_Q(copy_sim)
    coeffs = AMD_stability_coefficients(copy_sim)
    return (Nbound,stableQ,coeffs)

def count_links(sim):
    orbits = sim.orbits(primary=sim.particles[0])
    links = 0
    for i in range(sim.N-1):
        oi = orbits[i]
        for j in range(i+1,sim.N-1):
            oj = orbits[j]
            links += linking_l(oi,oj)<0
    return links

def count_bound(sim):
    orbits = sim.orbits(sim.particles[0])
    return np.sum([o.a>0 for o in orbits])

def get_min_q(sim):
    orbits = sim.orbits(sim.particles[0])
    return np.min([o.a*(1-o.e) for o in orbits if o.a>0])


def process_sa(sa,outfile):
    # - 1
    # crossing times
    # count data
    # final orbit
    # max dE
    # min q
    sim0 = sa[0]
    E0 = sim0.energy()
    tcross=np.inf
    
    Nbound = 5
    eject_times=[0.]
    counts=[Nbound]

    
    dE_by_E_max = 0
    dE_times = []
    dE_values = []

    min_q = np.inf    
    min_q_times = []
    min_q_values = []

    for sim in sa:
        
        # get first orbit crossing
        if  tcross==np.inf and count_links(sim)>0:
            tcross=sim.t
        
        # track number of bound planets
        N1 = count_bound(sim)
        if N1 < Nbound:
            Nbound = N1
            eject_times.append(sim.t)
            counts.append(Nbound)
        
        # track max energy error
        dE_by_E = np.abs(sim.energy()/E0 - 1)
        if dE_by_E > dE_by_E_max:
            dE_by_E_max = dE_by_E
            dE_times.append(sim.t)
            dE_values.append(dE_by_E_max)


        # 
        q = get_min_q(sim)
        if q < min_q:
            min_q = q
            min_q_times.append(sim.t)
            min_q_values.append(min_q)

    final_sim = sa[-1]
    _,amd_stableQ,amd_coeffs = amd_stability_data(final_sim)
    final_time = final_sim.t

    min_q_times.append(final_time)
    min_q_values.append(min_q)

    dE_times.append(final_time)
    dE_values.append(dE_by_E_max)

    eject_times.append(final_sim.t)    
    counts.append(Nbound)

    np.savez_compressed(outfile,
                        min_q_times=min_q_times,
                        min_q_values=min_q_values,
                        dE_times=dE_times,
                        dE_values=dE_values,
                        eject_times=eject_times,
                        Nbound=counts,
                        final_amd_coeffs=amd_coeffs,
                        final_orbits=orbits_of_sim(final_sim),
                        tcross=tcross,
                        amd_stableQ=amd_stableQ
    )
