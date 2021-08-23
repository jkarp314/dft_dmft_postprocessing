#example script to calculate the density matrix after a converged DMFT calculation
#uses triqs 3.0
# this example is for a 3 orbital calculation for Sr2RuO4


import numpy as np
from h5 import HDFArchive
from triqs_dft_tools.sumk_dft import *
from triqs.gf import *
from triqs.operators.util import *
from triqs_cthyb import *
import triqs.utility.mpi as mpi

U = 2.3
J = .4
beta = 50

dc_type = 1                      # DC type: 0 FLL, 1 Held, 2 AMF
use_blocks = True                # use block structure from DFT input
prec_mu = 1e-5

outfile = 'J0.4_out' #output file containing converged DMFT solution

SK = SumkDFT(hdf_file = '../SRO400k.h5',use_dft_blocks=use_blocks) #load input from Wannier Hamiltonian


p = {}
# solver
p["random_seed"] = 123 * mpi.rank + 567
p["length_cycle"] = 200
p["n_warmup_cycles"] = int(5e4)
p["n_cycles"] = int(2e8/mpi.size)
p["move_double"]= True
p["perform_tail_fit"] = True
p["fit_max_moment"] = 4
p["fit_min_w"] = 5
p["fit_max_w"] = 15
p["measure_density_matrix"] = True
p["use_norm_as_weight"] = True

n_orb = SK.corr_shells[0]['dim']
spin_names = ["up","down"]
orb_names = [i for i in range(n_orb)]
gf_struct = SK.gf_struct_solver_list[0]
Umat, Upmat = U_matrix_kanamori(n_orb=n_orb, U_int=U, J_hund=J)
h_int = h_int_kanamori(spin_names, orb_names, map_operator_structure=SK.sumk_to_solver[0], U=Umat, Uprime=Upmat, J_hund = J)
S = Solver(beta=beta, gf_struct=gf_struct)

previous_runs = 0
if mpi.is_master_node():
    ar = HDFArchive(outfile+'.h5','a')
    previous_runs = ar['iterations']
    S.Sigma_iw = ar['Sigma_iw']
    SK.chemical_potential = ar['mu-%d'%previous_runs]
    SK.dc_imp = ar['dc_imp']
    del ar
previous_runs    = mpi.bcast(previous_runs)
S.Sigma_iw = mpi.bcast(S.Sigma_iw)
SK.chemical_potential = mpi.bcast(SK.chemical_potential)
SK.dc_imp = mpi.bcast(SK.dc_imp)

SK.symm_deg_gf(S.Sigma_iw,ish=0)                        # symmetrizing Sigma
SK.set_Sigma([S.Sigma_iw])
chemical_potential = SK.calc_mu( precision = prec_mu )  # find the chemical potential for given density
Gloc = SK.extract_G_loc()[0]                       # calc the local Green function

nlat = Gloc.total_density().real

S.G0_iw << inverse(S.Sigma_iw + inverse(Gloc))

S.solve(h_int=h_int, **p)
for name, s_iw in S.Sigma_iw:
    S.Sigma_iw[name] = make_hermitian(s_iw)

nimp = S.G_iw.total_density().real

if mpi.is_master_node():
    ar = HDFArchive('density_matrix.h5','a')
    ar['G_0'] = S.G0_iw
    ar['G_tau'] = S.G_tau
    ar['G_iw'] = S.G_iw
    ar['Sigma_iw'] = S.Sigma_iw
    ar['nimp'] = nimp
    ar['nlat'] = nlat
    ar['mu'] = SK.chemical_potential.real
    ar['rho'] = S.density_matrix
    ar['h_loc_diag'] = S.h_loc_diagonalization
    del ar
