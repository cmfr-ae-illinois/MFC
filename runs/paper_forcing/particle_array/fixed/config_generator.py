import copy
import json
import sys
from pathlib import Path
import numpy as np

# freestream conditions
M = 1.2
Re = 1500.0

P = 101325
rho = 1.225

# particle diameter
D = 0.1

# physical params
gamma = 1.4
R = 287.0
Pr = 0.72

cp = gamma * R / (gamma - 1.0)
v = M * np.sqrt(gamma * P / rho) 
mu = rho * v * D / Re
k = mu * cp / Pr

# filter width 
filter_width = 1.5*D * np.sqrt(2.0/(9.0*np.pi))

# saved file indices
FIRST_STEP = 100
LAST_STEP = 1000
STEP_INTERVAL = 100

DEFAULT_MAP_KIND = "box"

BOX_MAP = {
    "kind": "box",
    "length": 0.4,  
    "points": 50,  
}

AXISYMMETRIC_MAP = {
    "kind": "axisymmetric",
    "axis": 1, # 0 = x, 1 = y, 2 = z
    "axial_length": 0.2,  
    "axial_points": 101,
    "radius": 0.1,  
    "radial_points": 51,
    "azimuthal_points": 64,  
}

CONFIG = {
    "periodic": [True, True, True],
    "grid": {
        "shape": [128, 128, 128],  
        "lengths": [1.0, 1.0, 1.0],
        "origin": [0.0, 0.0, 0.0],
    },
    "physics": {
        "eos": "ideal_gas",
        "transport": "constant",
        "gamma": gamma,  
        "gas_constant": R,  
        "mu": mu,  
        "conductivity": k,
    },
    "io": {
        "format": "shared",
        "dtype": "<f8",
        "particle_dtype": "<f8",
        "nvars": 6,  
        "indices": {"rho": 1, "momentum": [2, 3, 4], "energy": 5},
        "flow": "restart_data/lustre_{step}.dat",
        "particles": "restart_data/ib_state_{step}.dat",
    },
    "sigmas": [filter_width],  
    "alpha_min": 1e-16,
    "fft_workers": 4,
    "surface": {
        "polar_points": 32,
        "offset_cells": 2.0,
        "velocity_bc": "no_slip",
        "thermal_bc": "adiabatic",
    },
    "fields": [], # all outputs
    "output": "closure_statistics.npz",
}


def build_config(map_kind):
    maps = {"box": BOX_MAP, "axisymmetric": AXISYMMETRIC_MAP}
    config = copy.deepcopy(CONFIG)
    config["snapshots"] = list(range(FIRST_STEP, LAST_STEP + 1, STEP_INTERVAL))
    config["map"] = copy.deepcopy(maps[map_kind])
    return config


def main():

    map_kind = "box"
    
    path = Path("post_process_config.json")

    config = build_config(map_kind)
    text = json.dumps(config, indent=2, allow_nan=False) + "\n"
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    print(text, end="")
    print(f"Saved configuration to {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
