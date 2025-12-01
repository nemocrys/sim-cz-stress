import os
import numpy as np
import yaml
import shutil
import pandas as pd
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
from geo_3D import geometry
from setup import simulation_pyelmer
import T_boundary_from2Dto3D


# Load base configurations
with open("config_geometry.yml") as f:
    config_geo_base = yaml.safe_load(f)
with open("config_sim.yml") as f:
    config_sim = yaml.safe_load(f)
with open("config_mat.yml") as f:
    config_mat = yaml.safe_load(f)

def run_simulations(simulations):
    for simulation in simulations:
        print(f"Starting simulation: {simulation}")

        dir_2D = "../2D/"+simulation +"/simdata/01_case" # 2D simulation path
        sim_dir = "./simdata/"+simulation + "/"  # where 3D simulation will be stored
        os.makedirs(sim_dir, exist_ok=True)

        model = geometry(
            config_geo_base,
            sim_dir=sim_dir,
            dir_2D=dir_2D,
            #name=f"sim_{idx}",
            visualize=False,
            include_atmosphere=True
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


if __name__ == "__main__":

    simulations = [
        "Csi_reference_case",
        # 'Csi_optimum'
    ]
    run_simulations(simulations)