import os
import numpy as np
from pyelmer import elmerkw as elmer
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
import yaml
from copy import deepcopy
import T_boundary_from2Dto3D

from geo_3D import geometry


def simulation_pyelmer(
    model, config, sim_dir="./simdata", config_mat={}, elmer_config_file="config_elmer.yml"
):

    ### SOLVERS SETUP

    sim = elmer.load_simulation("3D_steady", elmer_config_file)
    solver_heat = elmer.load_solver("HeatSolver", sim, elmer_config_file)
    solver_stress = elmer.load_solver("StressSolver", sim, elmer_config_file)
    flux_solver = elmer.load_solver("FluxSolver", sim, elmer_config_file)
    elmer.load_solver("ResultOutputSolver", sim, elmer_config_file)#
    elmer.load_solver("save_scalars", sim, elmer_config_file) 
    elmer.load_solver("save_line_solver_axial", sim, elmer_config_file)  
    elmer.load_solver("save_line_solver_radial", sim, elmer_config_file)  


    axial_line = elmer.load_solver("save_line_solver_axial",sim,elmer_config_file)
    axial_line.data.update(
        {
            "Equation": '"SaveLineAxial"',
            "Polyline Coordinates(2,3) ":"0.0 0.0 0.02 0.0 0.0 0.142",
            "Filename": '"axial.dat"',
        }
    )
    radial_line = elmer.load_solver("save_line_solver_radial",sim,elmer_config_file)
    radial_line.data.update(
        {
            "Equation": '"SaveLineRadial"',
            "Polyline Coordinates(2,3) ":"0.0 -0.01 0.07 0.0 0.01 0.07",
            "Filename": '"radial.dat"',
        }
    )

    ### EQUATIONS SETUP

    eqn_main = elmer.Equation(sim, "eqn_main", [solver_heat,solver_stress, flux_solver])



    # add crystal
    crystal = elmer.Body(sim,'crystal', [model["crystal"].ph_id],{ "crystal" : "Logical True"})
    material_name = model["crystal"].params.material
    mat = elmer.Material(sim, material_name, config_mat[material_name])
    melting_point = mat.data["Melting Point"] 
    ic = elmer.InitialCondition(sim, "T_crystal", {"Temperature": model["crystal"].params.T_init})

    crystal.equation = eqn_main
    crystal.material = mat
    crystal.initial_condition = ic


    boundries_list =[] #save all boundries for post-processing


    #Crystal bnd

    include_file = ["phase-if.dat","crystal_T_boundary_side.dat","crystal_T_boundary_top.dat"]
    interpolation_coordinate = [1,2,1] # T(x),T(y),T(x)
    #interpolation_coordinate_at_sif = [1,3,1] # T(x),T(z),T(x)
    interpolation_coordinate_at_sif = [3,3,1] # T(x),T(z),T(x)

    for index, bnd in enumerate ([
        "bnd_crystal_btm",
        "bnd_crystal_side",
        "bnd_crystal_top"
    ],start=0):
        boundries_list.append(bnd) # for post-pro

        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.data = {
            "Temperature " : f"""Variable Coordinate  {interpolation_coordinate_at_sif[index]} 
    Real
        include {include_file[index]} 
    End"""     
        }

        if index ==2:
            bnd.data.update(
                {
                    #"Displacement 1":0.0,
                    #"Displacement 2":0.0,
                    "Displacement 3":0.0,
                }
            )




    with open( os.path.join(sim_dir, "boundaries.txt") , 'w') as file:
        for line in boundries_list:
            file.write(line + "\n")

    sim.write_sif(sim_dir)


if __name__ == "__main__":
    # sim_dir = "./simdata/12_reference_diff_heights_reference/"
    # dir_2D = "./2D/2D_Csi_reference_case/simdata/06_reference_case_isotropic_60mm"

    # dir_2D = "./2D/Andrejs_plot/CsI_optimised/macplasData/length=0.060_vpull=0.2"  
    # sim_dir = "./simdata/12_reference_diff_heights_CsI_opt/"
    # dir_2D = "./2D/Andrejs_plot/Ga2O3_optimised/macplasData/length=0.060_vpull=0.033"  
    # sim_dir = "./simdata/12_reference_diff_heights_Ga2O3/"
    # dir_2D = "./2D/Andrejs_plot/Al2O3_optimised/macplasData/length=0.060_vpull=0.033"  
    # sim_dir = "./simdata/12_reference_diff_heights_Al2O3/"

    dir_2D = "./2D/2D_Csi_reference_case/simdata/05_reference_case_isotropic"  
    sim_dir = "./simdata/13_CsI_reference_mesh_0.8/"

    if os.path.exists(sim_dir):
        raise ValueError("Please remove the old simulation directory.")

    with open("config_geometry.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir,dir_2D=dir_2D,crystal_height_cut=None)

    with open("config_sim.yml") as f:
        config_sim = yaml.safe_load(f)
    with open("config_mat.yml") as f:
        config_mat = yaml.safe_load(f)


    

    T_boundary_from2Dto3D.make_T_files(dir_2D + "/results",sim_dir)

    simulation_pyelmer(model, config_sim, sim_dir, config_mat)

    run_elmer_grid(sim_dir, "mesh.msh")
    run_elmer_solver(sim_dir)
    err, warn, stats = scan_logfile(sim_dir)
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)
    
#ready