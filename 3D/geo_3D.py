import os
import numpy as np
from objectgmsh import Model, Shape, MeshControlLinear, MeshControlExponential, cut
import gmsh
import yaml
import matplotlib.pyplot as plt
from my_tools import *
from T_boundary_from2Dto3D import get_crystal_melt_interface

occ = gmsh.model.occ




def geometry(config, sim_dir="./", name="vgf",dir_2D="./", visualize=False, include_atmosphere=True, crystal_height_cut=None):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    model = Model(name)

    l_s = np.array(config["crystal"]["lengths"] )
    d_s = np.array(config["crystal"]["diameters"])
    #-----------------------------------------------------------------------------

    if crystal_height_cut is not None:
        # Use user-supplied value
        pass
    else:
        # Fallback to old default
        crystal_height_cut = np.sum(l_s)  # *(1/2) #  crystal lenght FROM TOP TO BOTTOM ( np.sum(l_s)--> no cut )

    print(f"Crystal height cut used: {crystal_height_cut}")
    #-----------------------------------------------------------------------------

    melt_y_top = (config["melt"]["h"]) # initial h

    # rho_melt = 3.2e+3, rho_crystal = 4.51e+3    CsI
    # rho_melt = 6000, rho_crystal = 5950       Ga2O3
    # rho_melt = 3500, rho_crystal = 3980       Al2O3


    def new_melt_h(h , rho_melt = 3.2e+3, rho_crystal = 4.51e+3  ):  # TODO: make this a function of material
        V_melt_f = V_cylinder(config["melt"]["h"], config["melt"]["r"]) # final
        M_melt_f = rho_melt * V_melt_f
        
        M_crystal_f = rho_crystal * crystal_volume(np.sum(l_s) ,l_s,d_s)
        M_crystal = rho_crystal * crystal_volume(h ,l_s,d_s)
        #### Melt change 
        M_melt = M_melt_f + (M_crystal_f - M_crystal )
        V_melt = M_melt / rho_melt 
        h_new = V_melt / (np.pi* config["melt"]["r"]**2)
        return h_new 

    melt_y_top_new = new_melt_h(crystal_height_cut)

    #melt_y_top_new = melt_y_top  # same level in this simulation
    #-----------------------------------------------------------------------------
    coords_btm = get_crystal_melt_interface(dir_2D + '/results/phase-if.dat')
    #-----------------------------------------------------------------------------

    #crystal = crystal_shape_stress_calc(model,3, l_s , d_s ,coords_btm[:,0],coords_btm[:,1],meniscus_cut=True ,h =crystal_height_cut,starting_point=[0.0, melt_y_top_new, 0.0],name="crystal")

    crystal = crystal_shape_from_2D(model,3, l_s , d_s ,coords_btm[:,0],coords_btm[:,1] ,h =crystal_height_cut,starting_point=[0.0, melt_y_top_new, 0.0],name="crystal")

    crystal.mesh_size = config["crystal"]["mesh_size"]
    crystal.params.material = config["crystal"]["material"]
    crystal.params.T_init = config["crystal"]["T_init"]

    crystal_y_top =crystal.params.h #- crystal_height_cut
    print("-----------------------------------------------------------------------")
    print("Maxmimum height for crystal: " , np.sum(np.array(l_s)))
    print("The crystal height: " , crystal_y_top)
    print("Melt height",melt_y_top_new)
    print("-----------------------------------------------------------------------")


    occ.fragment(
        crystal.dimtags,
        [],
    )
    model.synchronize()

    bnd_crystal_side = Shape(
        model,
        2,
        "bnd_crystal_side",
        crystal.get_boundaries_in_box(
            [-d_s[-1],d_s[-1]],
            [-d_s[-1],d_s[-1]],
            [melt_y_top_new - coords_btm[:,1].max(),melt_y_top_new + np.sum(np.array(l_s)) ],
        ),
    )

    bnd_crystal_top = Shape(
        model,
        2,
        "bnd_crystal_top",
        crystal.get_boundaries_in_box(
            [-d_s[0],d_s[0]],
            [-d_s[0],d_s[0]],
            [melt_y_top_new + crystal_y_top,melt_y_top_new + crystal_y_top],
        ),
    )

    bnd_crystal_btm = Shape(
        model,
        2,
        "bnd_crystal_btm",
        crystal.get_boundaries_in_box(
            [-d_s[-1],d_s[-1]],
            [-d_s[-1],d_s[-1]], 
            [melt_y_top_new - coords_btm[:,1].max(), melt_y_top_new ],
        ),
    )

    bnd_crystal_side -= bnd_crystal_top
    bnd_crystal_side -= bnd_crystal_btm



    

    model.make_physical()
    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()

    model.generate_mesh(**config["mesh"])

    if visualize:
        model.show()
    print(model)
    model.write_msh(sim_dir + "/mesh.msh")

    model.close_gmsh()

    return model


if __name__ == "__main__":


    sim_dir = "./"
    dir_2D = "../2D/Csi_reference_case/simdata/01_case"


    with open("config_geometry.yml") as f:
        config_geo = yaml.safe_load(f)

    model = geometry(config_geo, sim_dir,dir_2D=dir_2D, visualize=True , include_atmosphere=True, crystal_height_cut=None)