import json
import math
import numpy as np
from MFC_particle_forces import Osnes_CD

gam_a = 1.4


# Liao parameters
mu = 0.001 
rho_f = 1000
rho_s = 1.14 * rho_f
D = 0.01/6

x0 = 0; x1 = 6*D
y0 = -6*D; y1 = 30*D
z0 = 0; z1 = 6*D

sphere_loc1 = [3.03*D, 21*D, 3.03*D]
sphere_loc2 = [2.97*D, 18.96*D, 2.97*D]

Re = 120

N_s = 2

# calculations from Liao params
vel = Re * mu / (rho_f * D)

R = D / 2.0

vol_s = 4.0/3.0 * np.pi * R**3
mass_s = rho_s * vol_s
particle_vf = (N_s * vol_s) / ((x1-x0)*(y1-y0)*(z1-z0)) 
fluid_vf = 1.0 - particle_vf

# fluid params
M = 0.1
cs = vel / M
P = cs**2 * rho_f / gam_a

# control params
CD = Osnes_CD(particle_vf, Re, M, gam_a)
drag = 0.5 * rho_f * vel**2 * np.pi * R**2 * CD
g0 = drag / mass_s
fRe = CD * Re / 24.0
tau_p = 2.0/9.0 * rho_s * R**2 / (mu * fRe)

Cg = 2.0; Cp = 1000.0

K_Pg = -1.0/(Cg*tau_p)
K_Dg = -0.2
K_Pp = -2.0*P/(Cp*M)


Nx = 63
Ny = 511
Nz = 63

CFL = 0.45
dx = (x1 - x0)/Nx
dy = (y1 - y0)/Ny
dz = (z1 - z0)/Nz

t_final = 0.8
dt = CFL * min(dx/cs, dy/(vel+cs), dz/cs)
Nt = int(t_final/dt)
t_save = 2# Nt//300

collision_time = 10.0 * dt

# print('vel', vel)
# print('cs', cs)
# print(dt)
# print(D/dx, D/dy, D/dz)

# immersed boundary dictionary
ib_dict = {}
ib_dict.update({
    f"patch_ib({1})%geometry": 8,
    f"patch_ib({1})%x_centroid": sphere_loc1[0],
    f"patch_ib({1})%y_centroid": sphere_loc1[1],
    f"patch_ib({1})%z_centroid": sphere_loc1[2],
    f"patch_ib({1})%radius": D / 2,
    f"patch_ib({1})%slip": "F",
    f"patch_ib({1})%moving_ibm": 2,
    f"patch_ib({1})%mass": mass_s,

    f"patch_ib({2})%geometry": 8,
    f"patch_ib({2})%x_centroid": sphere_loc2[0],
    f"patch_ib({2})%y_centroid": sphere_loc2[1],
    f"patch_ib({2})%z_centroid": sphere_loc2[2],
    f"patch_ib({2})%radius": D / 2,
    f"patch_ib({2})%slip": "F",
    f"patch_ib({2})%moving_ibm": 2,
    f"patch_ib({2})%mass": mass_s,
    })
# start two spheres vertically stacked

# Configuring case dictionary
case_dict = {
    # Logistics
    "run_time_info": "T",
    # Computational Domain Parameters
    # x direction
    "x_domain%beg": x0,
    "x_domain%end": x1,
    # y direction
    "y_domain%beg": y0,
    "y_domain%end": y1,
    # z direction
    "z_domain%beg": z0,
    "z_domain%end": z1,
    "cyl_coord": "F",
    "m": Nx,
    "n": Ny,
    "p": Nz,
    "dt": dt,
    "t_step_start": 0,
    "t_step_stop": Nt,  
    "t_step_save": t_save,  
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
    "wave_speeds": 1,
    "low_Mach": 2,
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
    # Patch: Constant Tube filled with air
    # Specify the cylindrical air tube grid geometry
    "patch_icpp(1)%geometry": 9,
    "patch_icpp(1)%x_centroid": (x1+x0)/2,
    # Uniform medium density, centroid is at the center of the domain
    "patch_icpp(1)%y_centroid": (y1+y0)/2,
    "patch_icpp(1)%z_centroid": (z1+z0)/2,
    "patch_icpp(1)%length_x": x1-x0,
    "patch_icpp(1)%length_y": y1-y0,
    "patch_icpp(1)%length_z": z1-z0,
    # Specify the patch primitive variables
    "patch_icpp(1)%vel(1)": 0.0e00,
    "patch_icpp(1)%vel(2)": vel,
    "patch_icpp(1)%vel(3)": 0.0e00,
    "patch_icpp(1)%pres": P,
    "patch_icpp(1)%alpha_rho(1)": rho_f,
    "patch_icpp(1)%alpha(1)": 1.0e00,
    # Patch: Sphere Immersed Boundary
    # Fluids Physical Parameters
    "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),  # 2.50(Not 1.40)
    "fluid_pp(1)%pi_inf": 0,
    "fluid_pp(1)%Re(1)": 1.0 / mu,

    # periodic forcing
    "periodic_forcing": "T",
    "u_inf_ref": vel,
    "rho_inf_ref": rho_f,
    "P_inf_ref": P,
    "mom_f_idx": 2,
    "forcing_window": 1,
    "forcing_dt": 1.0/(0.5*dt),
    "fluid_volume_fraction": fluid_vf,
    "forcing_wrt": "T",

    # controls
    "particle_control": "T",
    "particle_bf": -g0,

    "cntrl_p%Re_tgt": Re,
    "cntrl_p%M_tgt": M,
    "cntrl_p%K_Pg": K_Pg,
    "cntrl_p%K_Dg": K_Dg,
    "cntrl_p%K_Pp": K_Pp,
    "cntrl_p%window_size": 1,

    # collisions
    "collision_model": 1,  # soft-sphere collision model
    "ib_coefficient_of_friction": 0.00,
    "collision_time": collision_time,
    "coefficient_of_restitution": 1.0,  
    "ib_state_wrt": "T",

    }

case_dict.update(ib_dict)

print(json.dumps(case_dict))