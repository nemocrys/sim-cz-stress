import os
import numpy as np
import yaml
import shutil
import pandas as pd
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
from geometry import geometry
from setup import simulation_pyelmer
import post

# Load base configurations
with open("config_geometry.yml") as f:
    config_geo_base = yaml.safe_load(f)
with open("config_sim.yml") as f:
    config_sim = yaml.safe_load(f)
with open("config_mat.yml") as f:
    config_mat = yaml.safe_load(f)

total_height = sum(config_geo_base["crystal"]["lengths"])
#cut_heights = np.arange(0.010, 0.1167, 0.01)
# cut_heights = np.arange(0.010, 0.1167, 0.001)


cut_heights = np.arange(0.028, total_height, 0.001)
if not np.any(np.isclose(cut_heights, total_height)):
    cut_heights = np.append(cut_heights, total_height)


v_pull = 0.2
sim_dirs = []
lengths = []
vpulls = []
times = []

global_csv_dir = "./simdata/macplasData/macplasCSVs"
os.makedirs(global_csv_dir, exist_ok=True)

for idx, h_cut in enumerate(cut_heights, start=1):
    try:
        # Deep copy geometry config and override cut height
        config_geo = yaml.safe_load(yaml.dump(config_geo_base))
        config_geo["crystal"]["height_cut"] = float(h_cut)

        # Set up a unique simulation directory for each run
        sim_dir = os.path.join("./simdata/macplasData", f"length={h_cut:.3f}_vpull={v_pull}")
        os.makedirs(sim_dir, exist_ok=True)

        # sim_dirs.append(f"length={h_cut:.3f}_vpull={v_pull}")
        # lengths.append(h_cut * 10)
        # vpulls.append(v_pull)
        # times.append(h_cut * 10 / v_pull)

        print(f"=== Running simulation {idx}/{len(cut_heights)}: height_cut = {h_cut:.3f} --> {sim_dir} ===")

        # Generate geometry
        model = geometry(
            config_geo,
            sim_dir=sim_dir,
            name=f"sim_{idx}",
            visualize=False,
            include_atmosphere=True
        )

        # Launch Elmer simulation
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

        post.run_post(sim_dir)
        post.csv_for_macplas(sim_dir, name=f"length={h_cut:.3f}_vpull={v_pull}")

        csv_name = f"length={h_cut:.3f}_vpull={v_pull}.csv"
        src_csv = os.path.join(sim_dir, csv_name)
        if os.path.exists(src_csv):
            shutil.copy(src_csv, os.path.join(global_csv_dir, csv_name))
        else:
            print(f"Warning: expected CSV '{src_csv}' not found.")

        post.evaluate_heat_flux(sim_dir, config_mat["csi_solid"]["Heat Conductivity"])




        sim_dirs.append(f"length={h_cut:.3f}_vpull={v_pull}")
        lengths.append(h_cut * 1e3)
        vpulls.append(v_pull)
        times.append((h_cut * 1e3) / v_pull)
        print(f"--- Completed simulation {idx} ---\n")


    except Exception as e:
        print(f"!!! Simulation {idx} FAILED: height_cut = {h_cut:.3f} --> {e}\n")
        continue  # Skip to the next cut height

# Save summary CSV of all runs
df = pd.DataFrame({
    "filename": sim_dirs,
    "crystal length in mm": lengths,
    "v_pull in mm/min": vpulls,
    "time in min": times
})
output_csv = os.path.join(global_csv_dir, "time.csv")
df.to_csv(output_csv, index=False)
print(f"Summary CSV written to: {output_csv}")
