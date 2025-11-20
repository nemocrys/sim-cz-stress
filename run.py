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

simulations = [
    "2D_Csi_reference_case"
]


# simulations = [
#     "2D_shield_h_0.008_r_0.006",
#     "2D_shield_h_0.012_r_0.006",
#     "2D_shield_h_0.012_r_0.010",
#     "2D_shield_h_0.018_r_0.006",
# ]

# simulations = [
#     "2D_top_susceptor_h_0.047",
#     "2D_top_susceptor_h_0.067",
#     "2D_top_susceptor_h_0.087",
# ]

# simulations = [
#     "2D_side_insulation_h_0.025",
#     "2D_side_insulation_h_0.050",
# ]

# simulations = [
#     "2D_top_coil_1",
#     "2D_top_coil_2",
#     "2D_top_coil_3",
# ]

# simulations = [
#     "2D_top_susceptor_h_0.087_iridium",
#     "2D_top_susceptor_h_0.087_tungsten",
# ]

# simulations = [
#     '2D_iridium_thinnerCrucible',
#     # "2D_iridium",
#     "2D_tungsten",
# ]

# simulations = [
#     #"2D_optimum_CsI",
#     #"2D_optimum_Sapphire_vol2",
#     "2D_optimum_Ga2O3",
# ]


# simulations = [
#     "2D_top_coil_0_susceptor",
#     # "2D_top_coil_1_susceptor",
#     # "2D_top_coil_2_susceptor",
#     # "2D_top_coil_3_susceptor",
# ]

# simulations = [
#     #"2D_optimum_Sapphire_thin_crc_update_properties_topCoil_shield",
#     #"2D_optimum_Ga2O3_thin_crc_update_properties_shield"
#     "2D_optimum_CsI_shieldAdjust"
# ]

# simulations = [
#     #"2D_optimum_Sapphire_thin_crc_update_properties_topCoil_shield",
#     #"2D_optimum_Ga2O3_thin_crc_update_properties_shield"
#     "2D_optimum_CsI_shieldAdjust"
# ]

for simulaion in simulations:
    print(f"Starting simulation: {simulaion}")

    sim_dir = "./simdata/11_" +simulaion + "/"  # where 3D simulation will be stored
    dir_2D = "./2D/"+simulaion +"/simdata/01_case" # 2D simulation path
    # dir_2D = "./2D/"+simulaion +"/simdata/length=0.058_vpull=0.2" # 2D simulation path
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

    