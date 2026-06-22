import json
import math
import numpy as np

# load initial sphere locations
sphere_loc = np.loadtxt('sphere_array_locations.txt')
N_s = len(sphere_loc)

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
particle_vf = (N_s * vol_s) / (L**3) 
fluid_vf = 1.0 - particle_vf

# fluid params
M_start = 2.0
M_tgt = 2.0
Re = 500.0

P = 101325
rho = 1.225

v_start = M_start * np.sqrt(gam_a * P / rho) 
v_tgt   = M_tgt   * np.sqrt(gam_a * P / rho) 
mu = rho * v_tgt * D / Re

#print('mu: ', mu)
#print('v1: ', v1)
#print('rho: ', rho)
#print('Kn = ' + str( np.sqrt(np.pi*gam_a/2)*(M/Re) )) # Kn < 0.01 = continuum flow

dt = 1.0e-6
t_final = 6 * L / v_tgt
Nt = int(t_final / dt)
t_step_save = 1 #Nt // 300

Nx = 127
Ny = Nx
Nz = Ny

# immersed boundary dictionary
ib_dict = {}
for i in range(N_s):
  ib_dict.update({
      f"patch_ib({i+1})%geometry": 8,
      f"patch_ib({i+1})%x_centroid": sphere_loc[i, 0],
      f"patch_ib({i+1})%y_centroid": sphere_loc[i, 1],
      f"patch_ib({i+1})%z_centroid": sphere_loc[i, 2],
      f"patch_ib({i+1})%vel(2)": 0.0,
      f"patch_ib({i+1})%radius": D / 2,
      f"patch_ib({i+1})%slip": "F",
      f"patch_ib({i+1})%moving_ibm": 0, #2
      f"patch_ib({i+1})%mass": mass_s,
      })

# Configuring case dictionary
case_dict = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
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
    "t_step_start": 0,
    "t_step_stop": Nt,
    "t_step_save": t_step_save,  
    # Simulation Algorithm Parameters
    "num_patches": 1,
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
    "weno_order": 1,
    "weno_eps": 1.0e-16,
    "weno_Re_flux": "T",
    "weno_avg": "T",
    "avg_state": 1,
    # "mp_weno": "T",
    "riemann_solver": 2,
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
    "num_ibs": N_s,
    "viscous": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "parallel_io": "T",
    # fluid patch parameters
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": 0.0,
    # Uniform medium density, centroid is at the center of the domain
    "patch_icpp(1)%y_centroid": 0.0,
    "patch_icpp(1)%z_centroid": 0.0,
    "patch_icpp(1)%length_x": 10 * D,
    "patch_icpp(1)%length_y": 10 * D,
    "patch_icpp(1)%length_z": 10 * D,
    # Specify the patch primitive variables
    "patch_icpp(1)%vel(1)": 0.0e00,
    "patch_icpp(1)%vel(2)": v_start,
    "patch_icpp(1)%vel(3)": 0.0e00,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho,
    "patch_icpp(1)%alpha(1)": 1.0e00,
    # Patch: Sphere Immersed Boundary
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,

    # periodic forcing
    "periodic_forcing": "T",
    "u_inf_ref": v_tgt,
    "rho_inf_ref": rho,
    "P_inf_ref": P,
    "mom_f_idx": 2,
    "forcing_window": 1,
    "forcing_dt": 1.0/(8.0*dt),
    "fluid_volume_fraction": fluid_vf,
    "forcing_wrt": "T",
    "forcing_start": 0,
    }

case_dict.update(ib_dict)

print(json.dumps(case_dict))
