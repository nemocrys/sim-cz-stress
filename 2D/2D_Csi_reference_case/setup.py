import os
import numpy as np
from pyelmer import elmerkw as elmer
from pyelmer.execute import run_elmer_solver, run_elmer_grid
from pyelmer.post import scan_logfile
import yaml
import post

from geometry import geometry


def simulation_pyelmer(
    model, config, sim_dir="./simdata", config_mat={}, elmer_config_file="config_elmer.yml"
):


    omega = 2 * np.pi * config["heating_induction"]["frequency"]
    ### SOLVERS SETUP

    sim = elmer.load_simulation("axisymmetric_steady", elmer_config_file)
    sim.settings.update({"Angular Frequency": omega})
    
    solver_mgdyn = elmer.load_solver("MagnetoDynamics2DHarmonic",sim,  elmer_config_file)
    solver_mgdyn.data.update({"Angular Frequency": omega})
    solver_calcfields = elmer.load_solver("MagnetoDynamicsCalcFields",sim,  elmer_config_file)
    solver_calcfields.data.update({"Angular Frequency": omega})

    solver_heat = elmer.load_solver("HeatSolver", sim, elmer_config_file)

    solver_phase_change = elmer.load_solver("SteadyPhaseChange", sim, elmer_config_file)
    solver_phase_change.data["Triple Point Fixed"] = "Logical True"   
    solver_mesh = elmer.load_solver("MeshUpdate", sim, elmer_config_file)
    stress_solver = elmer.load_solver("StressSolver", sim, elmer_config_file)
    flux_solver = elmer.load_solver("FluxSolver", sim, elmer_config_file)
    
    elmer.load_solver("ResultOutputSolver", sim, elmer_config_file)
    elmer.load_solver("SaveLine", sim, elmer_config_file)
    elmer.load_solver("save_scalars", sim, elmer_config_file)

    ### EQUATIONS SETUP

    #solver_calcfields are calculated from the last interation of this run

    eqn_main = elmer.Equation(sim, "eqn_main", [solver_mgdyn ,solver_calcfields, solver_heat , solver_mesh ])
    eqn_phase_change = elmer.Equation(sim, "eqn_phase_change", [solver_phase_change])    
    eqn_crystal = elmer.Equation(sim, "eqn_crystal", [solver_mgdyn ,solver_calcfields, solver_heat,solver_mesh,stress_solver, flux_solver ]) 
    eqn_melt = elmer.Equation(sim, "eqn_melt", [solver_mgdyn ,solver_calcfields, solver_heat , solver_mesh, flux_solver ]) 
    if config["general"]["heat_convection"] : 
        eqn_melt.data.update({"Convection": "Constant"}) #
    eqn_filling = elmer.Equation(sim, "eqn_filling", [solver_mgdyn  , solver_mesh ])



    # crystal boundary export 
    save_line_bnd_crystal = elmer.load_solver("save_line_solver",sim,elmer_config_file)
    save_line_bnd_crystal.data.update(
        {
            "Equation": '"SaveLineCrystal"',
            "Save Mask": 'String "Save Line Crystal"',
            "Filename": '"crystal_T_boundary_side.dat"',
        }
    )

    save_line_bnd_crystal_top = elmer.load_solver("save_line_solver_1",sim,elmer_config_file)
    save_line_bnd_crystal_top.data.update(
        {
            "Equation": '"SaveLineCrystal_top"',
            "Save Mask": 'String "Save Line Crystal Top"',
            "Filename": '"crystal_T_boundary_top.dat"',
        }
    )
    
    # material function
    for function in config_mat["functions"]:
        sim.intro_text += function
    

    # set  growth velocity
    v_pull = config["general"]["v_pull"]
    v_pull /= 6e4 ## to m/s 
    print(f"Growth velocity: {v_pull} m/s")


    # forces
    current_source = elmer.BodyForce(sim, "Current Density")
    

    joule_heat = elmer.BodyForce(sim, "joule_heat")
    joule_heat.joule_heat = True
    if config["general"]["heat_control"] :
        joule_heat.smart_heat_control = True
        if config["smart-heater"]["control-point"]:
            joule_heat.smart_heater_control_point = [
                config["smart_heater"]["x"],
                config["smart_heater"]["y"],
                config["smart_heater"]["z"],
            ]
            joule_heat.smart_heater_T = config["smart_heater"]["T"]


    # add induction heating
    
    for shape in ["inductor_bottom"]:
        
        current_source.current_density = config["heating_induction"]["current"] / model[shape].params.area

        ind= elmer.Body(sim, shape, [model[shape].ph_id])
        material_name = model[shape].params.material
        mat = elmer.Material(sim, material_name, config_mat[material_name])
        ic = elmer.InitialCondition(sim, "T_" + shape, {"Temperature": model[shape].params.T_init})

        ind.equation = eqn_main
        ind.material = mat
        ind.body_force = current_source
        ind.initial_condition = ic




    # add crystal
    crystal = elmer.Body(sim,'crystal', [model["crystal"].ph_id],{ "crystal" : "Logical True"})
    material_name = model["crystal"].params.material
    mat = elmer.Material(sim, material_name, config_mat[material_name])
    melting_point = mat.data["Melting Point"] # set to 1500 °C
    ic = elmer.InitialCondition(sim, "T_crystal", {"Temperature": model["crystal"].params.T_init})

    crystal.equation = eqn_crystal
    crystal.material = mat
    crystal.initial_condition = ic
    crystal.body_force = joule_heat



    # add melt
    melt = elmer.Body(sim, "melt", [model["melt"].ph_id])
    material_name = model["melt"].params.material
    mat = elmer.Material(sim, material_name, config_mat[material_name])
    ic = elmer.InitialCondition(sim, "T_melt", {"Temperature": model["melt"].params.T_init})

    melt.equation = eqn_melt
    melt.material = mat
    melt.initial_condition = ic
    melt.body_force = joule_heat





    # add other bodies

    for shape in [
        "crucible",
        "seedholder",
        "shaft_al2o3",
        "top_axis",
        #"shield",
        #"top_susceptor",
        "outer_insulation",
        "axbot_adapter",
        "bottom_axis",
        "container",
        "inner_insulation_btm",
        "inner_insulation_side",
    ]:
        if shape == "crucible" or shape == "top_susceptor" or shape == "shield" or shape =="shaft_al2o3" or shape == "axbot_adapter" or shape=="bottom_axis" or shape == "container" :
            bdy = elmer.Body(sim, shape, [model[shape].ph_id],{ "joule int " + shape: "Logical True"})
        else:
            bdy = elmer.Body(sim, shape, [model[shape].ph_id])

        material_name = model[shape].params.material
        mat = elmer.Material(sim, material_name, config_mat[material_name])
        ic = elmer.InitialCondition(sim, "T_" + shape, {"Temperature": model[shape].params.T_init})

        bdy.equation = eqn_main
        bdy.material = mat    
        bdy.body_force = joule_heat
        bdy.initial_condition = ic



    atmosphere = elmer.Body(sim, "atmosphere", [model["atmosphere"].ph_id])
    material_name = model["atmosphere"].params.material
    mat = elmer.Material(sim, material_name, config_mat[material_name])
    atmosphere.equation = eqn_filling
    atmosphere.material = mat


    # setup phase change
    melt_crystal_if = elmer.Body(sim, "melt_crystal_if", [model["if_melt_crystal"].ph_id])
    melt_crystal_if.equation = eqn_phase_change
    melt_crystal_if.material = crystal.material

    t0_phase_change = elmer.InitialCondition( sim, "t0_phase_change", {"Temperature": melting_point} )
    t0_phase_change.data = {"PhaseSurface": "Real 0.0"}

    melt_crystal_if.initial_condition = t0_phase_change
    if_melt_crystal = elmer.Boundary(sim,"if_melt_crystal",[model["if_melt_crystal"].ph_id])
    if_melt_crystal.save_line = True
    if_melt_crystal.normal_target_body = crystal
    if_melt_crystal.smart_heater = True
    if_melt_crystal.smart_heater_T = config["smart-heater"]["T"]
    if_melt_crystal.phase_change_steady = True
    if_melt_crystal.phase_change_body = melt_crystal_if
    if_melt_crystal.phase_change_vel = v_pull

    if_melt_crystal.material = crystal.material
    if_melt_crystal.save_scalars = True
    if_melt_crystal.mesh_update = [0, "Equals PhaseSurface"]

    boundries_list =[] #save all boundries for post-processing
    boundries_list.append("melt_crystal_if")

    # boundaries with convection 

        #Crystal bnd

    for bnds in [
        "bnd_crystal_top",
        "bnd_crystal_side",
    ]:
        boundries_list.append(bnds) # for post-pro

        bnd = elmer.Boundary(sim, bnds, [model[bnds].ph_id])
        bnd.radiation = True
        bnd.mesh_update = [0, 0]
        if bnds == "bnd_crystal_top":
            bnd.data.update({"DisplacementCrystal 2": "Real 0.0"})  # Necessary for Stress solver
        bnd.T_ext = config["boundaries"]["crystal"]["T_ext"]
        bnd.heat_transfer_coefficient = config["boundaries"]["crystal"]["htc"]
        bnd.save_scalars = True     


    for bnds in [
        "bnd_melt",
        "bnd_crucible",
        #"bnd_shield" , 
        #"bnd_top_susceptor", 
    ]:
       boundries_list.append(bnds) # for post-pro

       bnd = elmer.Boundary(sim, bnds, [model[bnds].ph_id])
       bnd.radiation = True
       bnd.mesh_update = [0, 0]
       if bnds == "bnd_melt":
            bnd.T_ext = config["boundaries"]["melt"]["T_ext"]
            bnd.heat_transfer_coefficient = config["boundaries"]["melt"]["htc"]
       elif bnds == "bnd_crucible":
            bnd.T_ext = config["boundaries"]["crucible_outside"]["T_ext"]
            bnd.heat_transfer_coefficient = config["boundaries"]["crucible_outside"]["htc"]
       elif bnds == "bnd_shield":
            bnd.T_ext = config["boundaries"]["shield"]["T_ext"]
            bnd.heat_transfer_coefficient = config["boundaries"]["shield"]["htc"]
       elif bnds == "bnd_top_susceptor":
            bnd.T_ext = config["boundaries"]["top_susceptor"]["T_ext"]
            bnd.heat_transfer_coefficient = config["boundaries"]["top_susceptor"]["htc"]
       bnd.save_scalars = True


    # add boundaries with surface-to-surface radiation

    for bnd in [
        "bnd_top_axis",
        "bnd_shaft_al2o3",
        "bnd_seedholder",
        "bnd_outer_insulation",
        "bnd_axbot_adapter",
        "bnd_bottom_axis",
        "bnd_container",
        "bnd_inductor_bottom",
        #"bnd_inductor_top",
        "bnd_inner_insulation_side", 
        #"bnd_inner_insulation_btm", # 
    ]:
        boundries_list.append(bnd) # for post-pro

        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.radiation = True
        bnd.mesh_update = [0, 0]        
        bnd.save_scalars = True
        

    # stationary interfaces
    for bnd in [
        "if_top_axis_container",
        "if_crystal_seedholder",
        "if_seedholder_shaft_al2o3",
        "if_shaft_al2o3_top_axis",
        "if_crucible_melt",
        #"if_crucible_shield",
        #"if_crucible_top_susceptor",
        "if_outer_insulation_axbot_adapter",
        "if_axbot_adapter_bottom_axis",
        "if_bottom_axis_container",
        "if_crucible_inner_insulation_btm", # 
        "if_inner_insulation_btm_outer_insulation", #
        "if_crucible_inner_insulation_side", #
        "if_inner_insulation_side_outer_insulation", #
    ]:
        boundries_list.append(bnd) # for post-pro

        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.mesh_update = [0, 0]
        bnd.save_scalars = True
        

    with open( os.path.join(sim_dir, "boundaries.txt") , 'w') as file:
        for line in boundries_list:
            file.write(line + "\n")


    # add outside boundaries

    for bnd in [
        "bnd_inductor_bottom_inside",
        #"bnd_inductor_top_inside",
    ]:
        bnd = elmer.Boundary(sim, bnd, [model[bnd].ph_id])
        bnd.fixed_temperature = config["boundaries"]["inductor_inside"]["T"]
        bnd.mesh_update = [0, 0]




    bnd = elmer.Boundary(sim, "bnd_outer_container", [model["bnd_outer_container"].ph_id])
    bnd.fixed_temperature = config["boundaries"]["container_outside"]["T"]
    bnd.zero_potential = True
    bnd.mesh_update = [0, 0]
    bnd.save_scalars = True



        # symmetry axis
    bnd = elmer.Boundary(sim, "symmetry_axis", [model["symmetry_axis"].ph_id])
    bnd.mesh_update = [0, None]



    for bnd in ["bnd_crystal_side","if_crystal_seedholder"]:
        sim.boundaries[bnd].data.update({"Save Line Crystal": "Logical True"})

    sim.boundaries["bnd_crystal_top" ].data.update({"Save Line Crystal Top": "Logical True"})

    sim.write_sif(sim_dir)


if __name__ == "__main__":
    sim_dir = "./simdata/10_reference_case_isotropic"
    if os.path.exists(sim_dir):
        raise ValueError("Please remove the old simulation directory.")

    with open("config_geometry.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir)

    with open("config_sim.yml") as f:
        config_sim = yaml.safe_load(f)
    with open("config_mat.yml") as f:
        config_mat = yaml.safe_load(f)


    simulation_pyelmer(model, config_sim, sim_dir, config_mat)

    run_elmer_grid(sim_dir, "mesh.msh")
    run_elmer_solver(sim_dir)
    err, warn, stats = scan_logfile(sim_dir)
    print("Errors:", err)
    print("Warnings:", warn)
    print("Statistics:", stats)

    post.run_post(sim_dir)
    #post.evaluate_heat_flux(sim_dir,config_mat["csi_solid"]["Heat Conductivity"] )
