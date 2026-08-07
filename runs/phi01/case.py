import json
import numpy as np
from MFC_particle_forces import Osnes_CD

# load initial sphere locations
sphere_loc = np.loadtxt('sphere_array_locations.txt')
N_s = len(sphere_loc)

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
M = 0.8
Re = 1000.0

P = 101325
rho = 1.225

v1 = M * np.sqrt(gam_a * P / rho) 
mu = rho * v1 * D / Re

# body force
CD = Osnes_CD(particle_vf, Re, M, gam_a)
drag = 0.5 * rho * v1**2 * np.pi * R**2 * CD
g0 = drag / mass_s
fRe = CD * Re / 24.0
tau_p = 2.0/9.0 * rho_s * R**2 / (mu * fRe)

#print('mu: ', mu)
#print('v1: ', v1)
#print('rho: ', rho)
#print('Kn = ' + str( np.sqrt(np.pi*gam_a/2)*(M/Re) )) # Kn < 0.01 = continuum flow

dt = 2.0E-06
Nt = 20 #int(1 * L / v1 / dt)
t_save = 2
t_step_start_stats = 2

Nx = 127
Ny = Nx
Nz = Ny

Cg = 1.2; Cp = 100.0
K_Pg = -1.0/(Cg*tau_p)
K_Dg = -0.25
K_Pp = -2.0*P/(Cp*M)

collision_time = 5.0 * dt

# immersed boundary dictionary
ib_dict = {}
for i in range(N_s):
  ib_dict.update({
      f"patch_ib({i+1})%geometry": 8,
      f"patch_ib({i+1})%x_centroid": sphere_loc[i, 0],
      f"patch_ib({i+1})%y_centroid": sphere_loc[i, 1],
      f"patch_ib({i+1})%z_centroid": sphere_loc[i, 2],
      f"patch_ib({i+1})%radius": D / 2,
      f"patch_ib({i+1})%slip": "F",
      f"patch_ib({i+1})%moving_ibm": 2, 
      f"patch_ib({i+1})%mass": mass_s,
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
    # time step
    "mpp_lim": "F",
    # Correct errors when computing speed of sound
    "mixture_err": "T",
    # Use TVD RK3 for time marching
    "time_stepper": 3,
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
    "fd_order": 4,
    "viscous": "T",
    # Formatted Database Files Structure Parameters
    "format": 1,
    "precision": 2,
    "prim_vars_wrt": "T",
    "E_wrt": "T",
    "q_filtered_wrt": "T",
    "ib_state_wrt": "T", 
    "parallel_io": "T",
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
    "forcing_dt": 1.0/(2.0*dt),
    "fluid_volume_fraction": 0.9,

    "volume_filter_momentum_eqn": "T",
    "volume_filter_width": 3.0*D/2 * np.sqrt(2/(9*np.pi)),
    "slab_domain_decomposition": "T", 

    # controls
    "particle_control": "T",
    "particle_control_start": 5, #Nt//10,
    "particle_bf": -g0,

    "cntrl_p%Re_tgt": Re,
    "cntrl_p%M_tgt": M,
    "cntrl_p%K_Pg": K_Pg,
    "cntrl_p%K_Dg": K_Dg,
    "cntrl_p%K_Pp": K_Pp,
    "cntrl_p%window_size": 1,

    # collisions
    "collision_model": 1,  # soft-sphere collision model
    "ib_coefficient_of_friction": 0.1,
    "collision_time": collision_time,
    "coefficient_of_restitution": 0.95,  # almost perfectly elastic
    }

case_dict.update(ib_dict)

print(json.dumps(case_dict))
