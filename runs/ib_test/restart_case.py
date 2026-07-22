import json
import math
import numpy as np

gam_a = 1.4

# particle diameter
D = 0.1
R = D/2.0

# domain length
L = 10.0 * D

# particle params
rho_s = 10.0 
vol_s = 4.0/3.0 * np.pi * R**3
mass_s = rho_s * vol_s
N_s = 2
particle_vf = (N_s * vol_s) / (L**3) 
fluid_vf = 1.0 - particle_vf

# fluid params
M = 2.0
Re = 500.0

P = 101325
rho = 1.225

v1 = M * np.sqrt(gam_a * P / rho) 
mu = rho * v1 * D / Re

v12 = v1#/2

#print('mu: ', mu)
#print('v1: ', v1)
#print('rho: ', rho)
#print('Kn = ' + str( np.sqrt(np.pi*gam_a/2)*(M/Re) )) # Kn < 0.01 = continuum flow

dt = 5.0E-06
Nt1 = 200 
Nt2 = 500
t_save = 10

Nx = 127
Ny = Nx
Nz = Ny

# print(f'CFL: {v1*dt/(L/(Nx+1))}')

# immersed boundary dictionary
ib_dict = {}
ib_dict.update({
    f"patch_ib({1})%geometry": 8,
    f"patch_ib({1})%x_centroid": 0.0,
    f"patch_ib({1})%y_centroid": -0.43,
    f"patch_ib({1})%z_centroid": 0.0,
    f"patch_ib({1})%vel(2)": -10.0,
    f"patch_ib({1})%radius": D / 2,
    f"patch_ib({1})%slip": "F",
    f"patch_ib({1})%moving_ibm": 2,
    f"patch_ib({1})%mass": mass_s,

    f"patch_ib({2})%geometry": 8,
    f"patch_ib({2})%x_centroid": 0.0,
    f"patch_ib({2})%y_centroid": +0.43,
    f"patch_ib({2})%z_centroid": 0.0,
    f"patch_ib({2})%vel(2)": +10.0,
    f"patch_ib({2})%radius": D / 2,
    f"patch_ib({2})%slip": "F",
    f"patch_ib({2})%moving_ibm": 2,
    f"patch_ib({2})%mass": mass_s,
    })

# Configuring case dictionary
case_dict = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    'old_ic': 'T',
    'old_grid': 'T',
    't_step_old': 0,
    # x direction
    "x_domain%beg": -5.0 * D,
    "x_domain%end": 5.0 * D,
    # y direction
    "y_domain%beg": -5.0 * D,
    "y_domain%end": 5.0 * D,
    # z direction
    "z_domain%beg": -5.0 * D,
    "z_domain%end": 5.0 * D,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": Nt1,
    "t_step_stop": Nt2,  
    "t_step_save": t_save,  
    # Simulation Algorithm Parameters
    "num_patches": 0,
    # Use the 5 equation model
    "model_eqns": 2,
    # 6 equations model does not need the K \div(u) term
    "alt_soundspeed": "F",
    # One fluids: air
    "num_fluids": 1,
    # time step
    "mpp_lim": "F",
    # Correct errors when computing speed of sound
    "mixture_err": "T",
    # Use TVD RK3 for time marching
    "time_stepper": 3,
    # Reconstruct the primitive variables to minimize spurious
    # Use WENO5
    "weno_order": 5,
    "weno_eps": 1.0e-14,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 2,
    "mapped_weno": "T",
    "null_weights": "F",
    "mp_weno": "T",
    "riemann_solver": 2,
    "wave_speeds": 1,
    # periodic bc
    "bc_x%beg": -3,
    "bc_x%end": -3,
    "bc_y%beg": -7, # -7, -11
    "bc_y%end": -8, # -8, -12
    "bc_z%beg": -3,
    "bc_z%end": -3,

    'bc_y%grcbc_in': "T", 
    'bc_y%grcbc_out': "F", 
    "bc_y%vel_in(1)": 0.0,
    "bc_y%vel_in(2)": v12,
    "bc_y%vel_in(3)": 0.0,
    "bc_y%pres_in": P,
    "bc_y%alpha_rho_in(1)": rho,
    "bc_y%alpha_in(1)": 1.0,

    # Set IB to True and add 1 patch
    "ib": "T",
    "num_ibs": N_s,
    "viscous": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "parallel_io": "T",
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,
    }

case_dict.update(ib_dict)

print(json.dumps(case_dict))
