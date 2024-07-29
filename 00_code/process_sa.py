import rebound as rb
import numpy as np
from sa_processing_utils import process_sa
import sys
pl_type = sys.argv[1]
sim_id = int(sys.argv[2])

in_dir = "/fs/lustre/cita/hadden/06_free_floating_planets/03_simulations/"
out_dir = "/fs/lustre/cita/hadden/06_free_floating_planets/05_summary_data/"
in_file = in_dir+"five_{}_sim_{}.sa".format(pl_type,sim_id)
out_file = out_dir+"{}_{}_summary_data".format(pl_type,sim_id)
sa = rb.Simulationarchive(in_file)
process_sa(sa,out_file)