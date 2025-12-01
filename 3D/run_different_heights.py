import os
import numpy as np
import yaml
import shutil
import pandas as pd
import re
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
from geo_3D import geometry
from setup import simulation_pyelmer
import T_boundary_from2Dto3D



def extract_length(sim_name):
    # Looks for 'length=...' and returns float value
    m = re.search(r"length=([0-9.]+)", sim_name)
    if m:
        return float(m.group(1))
    else:
        raise ValueError(f"Could not extract length from '{sim_name}'")
    


# Load base configurations
with open("config_geometry.yml") as f:
    config_geo_base = yaml.safe_load(f)
with open("config_sim.yml") as f:
    config_sim = yaml.safe_load(f)
with open("config_mat.yml") as f:
    config_mat = yaml.safe_load(f)


# simulations = [
#     "06_reference_case_isotropic_45mm"
# ] # reference


simulations = [
    # "length=0.030_vpull=0.2",
    "length=0.045_vpull=0.2",
    # "length=0.060_vpull=0.2",
    # "length=0.075_vpull=0.2",
    # "length=0.090_vpull=0.2",
    # "length=0.117_vpull=0.2",
] # CsI

# simulations = [
#     "length=0.030_vpull=0.033",
#     "length=0.045_vpull=0.033",
#     "length=0.060_vpull=0.033",
#     "length=0.075_vpull=0.033",
#     "length=0.090_vpull=0.033",
#     "length=0.117_vpull=0.033",
# ] # oxides



        #crystal_length = 0.045
        crystal_length = extract_length(simulaion)


        model = geometry(
            config_geo_base,
            sim_dir=sim_dir,
            dir_2D=dir_2D,
            #name=f"sim_{idx}",
            visualize=False,
            include_atmosphere=True,
            crystal_height_cut=crystal_length 
        )

        T_boundary_from2Dto3D.make_T_files(dir_2D + "/results",sim_dir)

        simulation_pyelmer(
            model,
            config_sim,
            sim_dir=sim_dir,
            config_mat=config_mat
        )


        run_elmer_grid(sim_dir, "mesh.msh")
        run_elmer_solver(sim_dir)
        err, warn, stats = scan_logfile(sim_dir)
        print("Errors:", err)
        print("Warnings:", warn)
        print("Statistics:", stats)
        print("-" * 30)
    
    except Exception as e:
        print(f"Simulation {simulaion} FAILED with error:\n{e}")
        print("-" * 30)