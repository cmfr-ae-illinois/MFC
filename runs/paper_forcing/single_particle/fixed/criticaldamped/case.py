import json
import numpy as np

# load initial sphere locations
N_s = 1

gam_a = 1.4

# particle diameter
D = 0.1
R = D/2.0

# domain length
Lx = 10. * D
Ly = 10. * D 
Lz = 10. * D

# particle params
rho_s = 10.0 
vol_s = 4.0/3.0 * np.pi * R**3
mass_s = rho_s * vol_s
particle_vf = (N_s * vol_s) / (Lx * Ly * Lz) 
fluid_vf = 1.0 - particle_vf

# fluid params
M_tgt = 1.4
Re_tgt = 500.0

P_tgt = 101325.0
rho_tgt = 1.225

v_tgt   = M_tgt * np.sqrt(gam_a * P_tgt / rho_tgt) 
mu = rho_tgt * v_tgt * D / Re_tgt

factor = 0.2
rho_start = factor * rho_tgt
v_start = v_tgt
P_start = factor * P_tgt

# print('mu: ', mu)
# print('v_tgt: ', v_tgt)
# print('rho: ', rho)
# print('Kn = ' + str( np.sqrt(np.pi*gam_a/2)*(M_tgt/Re) )) # Kn < 0.01 = continuum flow

dt = 1.0e-6
t_final = 4 * Ly / v_tgt
Nt = int(t_final / dt)
t_step_save = Nt // 200

Nx = 127
Ny = Nx
Nz = Ny

alpha_forcing = 0.627
forcing_start = -1

# immersed boundary dictionary
ib_dict = {}
for i in range(N_s):
  ib_dict.update({
      f"patch_ib({i+1})%geometry": 8,
      f"patch_ib({i+1})%x_centroid": 0.0,
      f"patch_ib({i+1})%y_centroid": 0.0,
      f"patch_ib({i+1})%z_centroid": 0.0,
      f"patch_ib({i+1})%radius": D / 2,
      f"patch_ib({i+1})%slip": "F",
      })

# Configuring case dictionary
case_dict = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    # x direction
    "x_domain%beg": -Lx/2,
    "x_domain%end": +Lx/2,
    # y direction
    "y_domain%beg": -Ly/2,
    "y_domain%end": +Ly/2,
    # z direction
    "z_domain%beg": -Lz/2,
    "z_domain%end": +Lz/2,
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": t_step_save,  
    # Simulation Algorithm Parameters
    "num_patches": 1,
    # Use the 5 equation model
    "model_eqns": 2,
    # 6 equations model does not need the K \div(u) term
    "alt_soundspeed": "F",
    # One fluid
    "num_fluids": 1,
    # Correct errors when computing speed of sound
    "mixture_err": "T",
    # Use TVD RK3 for time marching
    "time_stepper": 3,
    # Reconstruct the primitive variables to minimize spurious

    # Use MUSCL
    "recon_type": 2,
    "muscl_order": 2,
    "muscl_lim": 1,

    # Use WENO5
    # "weno_order": 5,
    # "weno_eps": 1.0e-16,
    # "weno_Re_flux": "T",
    # "weno_avg": "T",
    # "mp_weno": "T",

    "riemann_solver": 2,
    "wave_speeds": 1,
    "avg_state": 1,
    # periodic bc
    "bc_x%beg": -1,
    "bc_x%end": -1,
    "bc_y%beg": -1,
    "bc_y%end": -1,
    "bc_z%beg": -1,
    "bc_z%end": -1,
    # Set IB to True and add 1 patch
    "ib": "T",
    "num_ibs": N_s,
    "viscous": "T",
    "fd_order": 4,
    "ib_neighborhood_radius": 20,
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "cons_vars_wrt": "T", 
    "E_wrt": "T",
    "parallel_io": "T",
    "ib_state_wrt": "T",
    # fluid patch parameters
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": Lx,
    "patch_icpp(1)%length_y": Ly,
    "patch_icpp(1)%length_z": Lz,
    # Specify the patch primitive variables
    "patch_icpp(1)%vel(1)": 0.0e00,
    "patch_icpp(1)%vel(2)": v_start,
    "patch_icpp(1)%vel(3)": 0.0e00,
    "patch_icpp(1)%pres": P_start,
    "patch_icpp(1)%alpha_rho(1)": rho_start,
    "patch_icpp(1)%alpha(1)": 1.0e00,
    # Patch: Sphere Immersed Boundary
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,

    # periodic forcing
    "periodic_forcing": "T",
    "u_inf_ref": v_tgt,
    "rho_inf_ref": rho_tgt,
    "P_inf_ref": P_tgt,
    "mom_f_idx": 2,
    "forcing_window": 1,
    "forcing_dt": 1.0/(alpha_forcing*dt),
    "fluid_volume_fraction": fluid_vf,
    "forcing_wrt": "T",
    "forcing_start": forcing_start,#Nt//10,

    }

case_dict.update(ib_dict)

print(json.dumps(case_dict))
