import json
import math
import numpy as np


gam_a = 1.4

D = 0.1
L = 10 * D

M = 1.4
Re = 200.0

P = 101325.0
rho = 1.225

v1 = M * np.sqrt(gam_a * P / rho)
mu = rho * v1 * D / Re

# print('mu: ', mu)
# print('v1: ', v1)
# print('rho: ', rho)
# print('Kn = ' + str( np.sqrt(np.pi*gam_a/2)*(M/Re) )) # Kn < 0.01 = continuum flow

dt = 0.8E-06
Nt = int(L * 4 / v1 / dt)
t_save = 25 #int(Nt/4)
t_step_start_stats = int(Nt/2)

# print(Nt, t_step_start_stats)

Nx = 127
Ny = Nx
Nz = Ny

# load initial sphere locations
sphere_loc = np.loadtxt('sphere_array_locations.txt')
N_sphere = len(sphere_loc)

particle_vf = N_sphere * (4.0/3.0*np.pi*(D/2.0)**3) / (L**3)
fluid_vf = 1.0 - particle_vf

# immersed boundary dictionary
ib_dict = {}
for i in range(N_sphere):
    ib_dict.update({
        f"patch_ib({i+1})%geometry": 8,
        f"patch_ib({i+1})%x_centroid": sphere_loc[i, 0],
        f"patch_ib({i+1})%y_centroid": sphere_loc[i, 1],
        f"patch_ib({i+1})%z_centroid": sphere_loc[i, 2],
        f"patch_ib({i+1})%radius": D / 2,
        f"patch_ib({i+1})%slip": "F",
    })

# Configuring case dictionary
case_dict = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    # x direction
    "x_domain%beg": -L/2,
    "x_domain%end": L/2,
    # y direction
    "y_domain%beg": -L/2,
    "y_domain%end": L/2,
    # z direction
    "z_domain%beg": -L/2,
    "z_domain%end": L/2,
    "cyl_coord": "F",
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,  
    "t_step_save": t_save,  
    "t_step_start_stats": t_step_start_stats,
    # Simulation Algorithm Parameters
    # Only one patches are necessary, the air tube
    "num_patches": 1,
    # Use the 5 equation model
    "model_eqns": 2,
    # 6 equations model does not need the K \div(u) term
    "alt_soundspeed": "F",
    # One fluids: air
    "num_fluids": 1,
    # Correct errors when computing speed of sound
    "mixture_err": "T",
    # Use TVD RK3 for time marching
    "time_stepper": 3,
    # Reconstruct the primitive variables to minimize spurious
    # Use WENO5
    "weno_order": 5,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "F",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": 2,
    "low_Mach": 0,
    "wave_speeds": 1,
    # periodic bc
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    # Set IB to True and add 1 patch
    "ib": "T",
    "num_ibs": N_sphere,
    "viscous": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    #"q_filtered_wrt": "T",
    "parallel_io": "T",
    "forcing_wrt": "T",
    # Patch: Constant Tube filled with air
    # Specify the cylindrical air tube grid geometry
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    # Uniform medium density, centroid is at the center of the domain
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": 10 * D,
    "patch_icpp(1)%length_y": 10 * D,
    "patch_icpp(1)%length_z": 10 * D,
    # Specify the patch primitive variables
    "patch_icpp(1)%vel(1)": v1,
    "patch_icpp(1)%vel(2)": 0.0e00,
    "patch_icpp(1)%vel(3)": 0.0e00,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0e00,
    # Patch: Sphere Immersed Boundary
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,

    # new case additions
    "periodic_forcing": "T",
    "u_inf_ref": v1,
    "rho_inf_ref": rho,
    "P_inf_ref": P,
    "mom_f_idx": 1,
    "forcing_window": 1,
    "forcing_dt": 1.0 / (2.5 * dt),
    "fluid_volume_fraction": fluid_vf,  # 1 - particle volume fraction
    # compute unclosed terms in volume filtered momentum equation

    #"volume_filter_momentum_eqn": "T",
    
    "filter_width": 3.0 * D / 2 * np.sqrt(2 / (9 * np.pi)),
    "compute_particle_drag": "T",
    # MPI domain decomposition into slabs instead of blocks
    "slab_domain_decomposition": "T",
}

case_dict.update(ib_dict)

print(json.dumps(case_dict))
