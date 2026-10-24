#!/bin/bash
# Thin wrapper. The mesh-independence study now lives in run_domain_sweep.sh
# (SWEEP=mesh), so it shares all the safety features: error logging, skipping
# finished cases, mesh backup, time-folder pruning, angle-of-attack loop.
#
#   UP_FIXED=-10 DOWN_FIXED=27 HALF_H=8 AOA="0 15 25" bash run_mesh_sweep.sh
#
# UP_FIXED / DOWN_FIXED / HALF_H = the converged domain from the earlier studies.
SWEEP=mesh exec bash "$(dirname "$0")/run_domain_sweep.sh"
