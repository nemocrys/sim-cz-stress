import os
import numpy as np
from objectgmsh import Model, Shape, MeshControlLinear, MeshControlExponential, cut
import gmsh
import yaml
import matplotlib.pyplot as plt
from my_tools import *
from my_tools import inductor
from my_tools import crystal_shape
from my_tools import crystal_volume
from my_tools import inductor_filling

occ = gmsh.model.occ

def geometry(config, sim_dir="./", name="vgf", visualize=False, include_atmosphere=True):
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    model = Model(name)

    l_s = np.array(config["crystal"]["lengths"] )
    d_s = np.array(config["crystal"]["diameters"])

    #-----------------------------------------------------------------------------
    crystal_height_cut = np.sum(l_s) #* (1/2)  #  crystal lenght FROM TOP TO BOTTOM ( np.sum(l_s)--> no cut )
    # crystal_height_cut = config["crystal"]["height_cut"]
    #-----------------------------------------------------------------------------


    melt_y_top = (config["melt"]["h"]) # initial h
    melt_y_bottom = (0.0)

    def new_melt_h(h , rho_melt = 3.2e+3, rho_crystal = 4.51e+3 ):
        V_melt_f = V_cylinder(config["melt"]["h"], config["melt"]["r"]) # final
        M_melt_f = rho_melt * V_melt_f
        
        M_crystal_f = rho_crystal * crystal_volume(np.sum(l_s) ,l_s,d_s)
        M_crystal = rho_crystal * crystal_volume(h ,l_s,d_s)
        #### Melt change 
        M_melt = M_melt_f + (M_crystal_f - M_crystal )
        V_melt = M_melt / rho_melt 
        h_new = V_melt / (np.pi* config["melt"]["r"]**2)
        return h_new 
    


    melt_y_top = new_melt_h(crystal_height_cut)
    #--------------------------------------------- MELT --------------------------------------------- #



    melt = Shape(
        model,
        2,
        "melt",
        [
            occ.add_rectangle(
                0,
                0,
                0,
                config["melt"]["r"],  
                melt_y_top,  
            )
        ],
    )
    # we use "params" to save various values we need later
    melt.mesh_size = config["melt"]["mesh_size"]
    melt.params.material = config["melt"]["material"]
    melt.params.T_init = config["melt"]["T_init"]

    #--------------------------------------------- CRYSTAL ---------------------------------------------#
    crystal = crystal_shape(model,2, l_s , d_s ,h =crystal_height_cut,starting_point=[0.0, melt_y_top, 0.0],name="crystal")

    crystal.mesh_size = config["crystal"]["mesh_size"]
    crystal.params.material = config["crystal"]["material"]
    crystal.params.T_init = config["crystal"]["T_init"]

    crystal_y_top =melt_y_top + crystal.params.h #- crystal_height_cut
    crystal.params.r = d_s[-1]

    print("-----------------------------------------------------------------------")
    print("Maxmimum height for crystal: " , np.sum(np.array(l_s)))
    print("-----------------------------------------------------------------------")



    #--------------------------------------------- CRUCIBLE ---------------------------------------------#

    crucible = Shape(
        model,
        2,
        "crucible",
        [
            occ.add_rectangle(
                0,
                config["crucible"]["y0"] ,
                0,
                config["crucible"]["r_out"],
                config["crucible"]["h"],
            )
        ],
    )

    crucible_hole = occ.add_rectangle(
        0,
        config["crucible"]["y0"] + config["crucible"]["t_bt"] ,
        0,
        config["crucible"]["r_in"],
        config["crucible"]["h"],
    )
    occ.cut(crucible.dimtags, [(2, crucible_hole)])

    crucible.mesh_size = config["crucible"]["mesh_size"]
    crucible.params.material = config["crucible"]["material"]
    crucible.params.T_init = config["crucible"]["T_init"]
    crucible.params.h = config["crucible"]["h"]
    crucible.params.r_in = config["crucible"]["r_in"]



    # --------------------------------------------- SEEDHOLDER --------------------------------------------- #


    seedholder_volume = occ.add_rectangle(   
        d_s[0],
        crystal_y_top - config["seedholder"]["overlap"],
        0,
        config["seedholder"]["r_out"] - d_s[0],
        config["seedholder"]["h"] ,
    )
    seedholder = Shape(model, 2, "seedholder", [seedholder_volume])



    seedholder.mesh_size = config["seedholder"]["mesh_size"]
    seedholder.params.material = config["seedholder"]["material"]
    seedholder.params.T_init = config["seedholder"]["T_init"]
    seedholder.params.y_top = crystal_y_top - config["seedholder"]["overlap"]

    # --------------------------------------------- shaft_al2o3 --------------------------------------------- #

    shaft_al2o3 = Shape(
        model,
        2,
        "shaft_al2o3",
        [
            occ.add_rectangle(
                0 ,
                seedholder.params.y_top + config["seedholder"]["h"] - config["shaft_al2o3"]["overlap"] ,
                0,
                config["shaft_al2o3"]["r"],
                0.318 - config["top_axis"]["h"] - (seedholder.params.y_top + config["seedholder"]["h"] - config["shaft_al2o3"]["overlap"] )
                #config["shaft_al2o3"]["h"] ,
            )
        ],
    )   

    shaft_al2o3.mesh_size = config["shaft_al2o3"]["mesh_size"]
    shaft_al2o3.params.material = config["shaft_al2o3"]["material"]
    shaft_al2o3.params.T_init = config["shaft_al2o3"]["T_init"]
    shaft_al2o3.params.y_top= seedholder.params.y_top + config["seedholder"]["h"] - config["shaft_al2o3"]["overlap"] + 0.318 - config["top_axis"]["h"] - (seedholder.params.y_top + config["seedholder"]["h"] - config["shaft_al2o3"]["overlap"] )



    # --------------------------------------------- TOP AXIS --------------------------------------------- #

    top_axis = Shape(
        model,
        2,
        "top_axis",
        [
            occ.add_rectangle(
                0 ,
                shaft_al2o3.params.y_top,
                0,
                config["top_axis"]["r"],
                config["top_axis"]["h"] ,
            )
        ],
    )   

    top_axis.mesh_size = config["top_axis"]["mesh_size"]
    top_axis.params.material = config["top_axis"]["material"]
    top_axis.params.T_init = config["top_axis"]["T_init"]



    #---------------------------------------------  SHIELD --------------------------------------------- #
    dx = - config["shield"]["r"] 
    dy = -(dx) *  np.tan(np.deg2rad(config["shield"]["theta"]))   # adapt the thickness as the rest
    shield_y0 = config["shield"]["y0"]

    points = [
        occ.add_point(config["crucible"]["r_in"],shield_y0 ,0 ),

        occ.add_point(config["crucible"]["r_in"] + dx ,shield_y0 + dy ,0 ),
        occ.add_point(config["crucible"]["r_in"] + dx, shield_y0 + dy + config["shield"]["thickness"]  ,0 ),
        occ.add_point(config["crucible"]["r_in"],shield_y0  +config["shield"]["thickness"] ,0 ),
    ]
    lines = [occ.add_line(points[i - 1], points[i]) for i in range(len(points))]
    loop = occ.add_curve_loop(lines)
    shield_volume = occ.add_surface_filling(loop)

    shield = Shape(model,2,"shield", [shield_volume])

    shield.mesh_size = config["shield"]["mesh_size"] /4
    shield.params.material = config["shield"]["material"]
    shield.params.T_init = config["shield"]["T_init"]

    #---------------------------------------------  TOP SUSCEPTOR --------------------------------------------- #
    top_susceptor_r_out = config["crucible"]["r_in"] + (config["crucible"]["r_out"] - config["crucible"]["r_in"] )
    top_susceptor_r_in = config["crucible"]["r_in"] 
    top_susceptor_y0 = config["crucible"]["h"] - config["crucible"]["t_bt"] 


    top_susceptor = Shape(
        model,
        2,
        "top_susceptor",
        [
            occ.add_rectangle(
                config["top_susceptor"]["r_top"],
                top_susceptor_y0,
                0,
                top_susceptor_r_out - config["top_susceptor"]["r_top"],
                config["top_susceptor"]["h"],
            )
        ],
    )

    top_susceptor_hole = occ.add_rectangle(
        config["top_susceptor"]["r_top"],
        top_susceptor_y0 ,
        0,
        top_susceptor_r_in - config["top_susceptor"]["r_top"],
        config["top_susceptor"]["h"] - config["top_susceptor"]["t_tp"],
    )
    occ.cut(top_susceptor.dimtags, [(2, top_susceptor_hole)])


    top_susceptor.mesh_size = config["top_susceptor"]["mesh_size"] /4
    top_susceptor.params.material = config["top_susceptor"]["material"]
    top_susceptor.params.T_init = config["top_susceptor"]["T_init"]


    #--------------------------------------------- HEATERS ---------------------------------------------#

    inductor_bottom = inductor(model,2, **config["inductor_bottom"],name="inductor_bottom")

    if config["ind_top"]["include"]:
        inductor_top = inductor(model,2, **config["inductor_top"], name="inductor_top")




    #---------------------------------------------INNER INSULATION BOTTOM---------------------------------------------#
    inner_insulation_btm = Shape(
        model,
        2,
        "inner_insulation_btm",
        [
            occ.add_rectangle(
                0,
                config["crucible"]["y0"] - config["inner_insulation_btm"]["h"],
                0,
                config["inner_insulation_btm"]["r"],
                config["inner_insulation_btm"]["h"],
            )
        ],
    )


    inner_insulation_btm.mesh_size = config["inner_insulation_btm"]["mesh_size"]
    inner_insulation_btm.params.material = config["inner_insulation_btm"]["material"]
    inner_insulation_btm.params.T_init = config["inner_insulation_btm"]["T_init"]
    inner_insulation_btm.params.y0 = config["crucible"]["y0"] - config["inner_insulation_btm"]["h"]


    #---------------------------------------------INNER INSULATION SIDE---------------------------------------------#
    inner_insulation_side = Shape(
        model,
        2,
        "inner_insulation_side",
        [
            occ.add_rectangle(
                config["crucible"]["r_out"],
                config["crucible"]["y0"],
                0,
                config["inner_insulation_side"]["r_out"] - config["crucible"]["r_out"],
                config["inner_insulation_side"]["h"],
            )
        ],
    )


    inner_insulation_side.mesh_size = config["inner_insulation_side"]["mesh_size"]
    inner_insulation_side.params.material = config["inner_insulation_side"]["material"]
    inner_insulation_side.params.T_init = config["inner_insulation_side"]["T_init"]
    inner_insulation_side.params.y0 = config["crucible"]["y0"] - config["inner_insulation_side"]["h"]

    #---------------------------------------------- OUTER_INSULATION ----------------------------------------------#
    outer_insulation = Shape(
        model,
        2,
        "outer_insulation",
        [
            occ.add_rectangle(
                0,
                config["outer_insulation"]["y0"] ,
                0,
                config["outer_insulation"]["r_out"] ,
                config["outer_insulation"]["h"] ,
            )
        ]
    )
    outer_insulation_bottom_hole = occ.add_rectangle(
        0,
        inner_insulation_btm.params.y0 ,
        0,
        config["inner_insulation_side"]["r_out"] ,
        config["outer_insulation"]["h"] - config["outer_insulation"]["t_bt"] - config["outer_insulation"]["t_tp"]  ,
    )    

    outer_insulation_top_hole = occ.add_rectangle(
        0,
        config["outer_insulation"]["top_hole_y0"],
        0,
        config["outer_insulation"]["r_in_tp"],
        config["outer_insulation"]["t_tp"],
    )


    occ.cut(outer_insulation.dimtags, [(2, outer_insulation_bottom_hole)])
    occ.cut(outer_insulation.dimtags, [(2, outer_insulation_top_hole)])

    if config["outer_insulation"]["top_half_cut"]:
        top_insulation = occ.add_rectangle(
            0,
            config["outer_insulation"]["y0"] + config["outer_insulation"]["h"]  - config["outer_insulation"]["h_top"]  ,
            0,
            config["outer_insulation"]["r_out"] ,
            config["outer_insulation"]["h"] ,
        )
        occ.cut(outer_insulation.dimtags, [(2, top_insulation)])


    outer_insulation.mesh_size = config["outer_insulation"]["mesh_size"]
    outer_insulation.params.material = config["outer_insulation"]["material"]
    outer_insulation.params.T_init = config["outer_insulation"]["T_init"]

    outer_insulation.params.y0 = inner_insulation_btm.params.y0 -  config["outer_insulation"]["t_bt"] 


    #---------------------------------- AXIS-BT-PLATE---------------------------------- #

    axbot_adapter = Shape(
        model,
        2,
        "axbot_adapter",
        [
            occ.add_rectangle(
                0,
                outer_insulation.params.y0  ,
                0,
                config["axbot_adapter"]["r_out"] , 
                -config["axbot_adapter"]["h"]  ,
            )
        ],
    )

    axbot_adapter.mesh_size = config["axbot_adapter"]["mesh_size"] /2 
    axbot_adapter.params.material = config["axbot_adapter"]["material"]
    axbot_adapter.params.T_init = config["axbot_adapter"]["T_init"]
    axbot_adapter.params.y0 = outer_insulation.params.y0 - config["axbot_adapter"]["h"] 

    #---------------------------------- BOTTOM AXIS ---------------------------------- #

    bottom_axis = Shape(
        model,
        2,
        "bottom_axis",
        [
            occ.add_rectangle(
                config["bottom_axis"]["r_in"] ,
                axbot_adapter.params.y0  ,
                0,
                config["bottom_axis"]["r_out"] - config["bottom_axis"]["r_in"]  , 
                -config["bottom_axis"]["h"]  ,
            )
        ],
    )

    bottom_axis.mesh_size = config["bottom_axis"]["mesh_size"] /2 
    bottom_axis.params.material = config["bottom_axis"]["material"]
    bottom_axis.params.T_init = config["bottom_axis"]["T_init"]

    bottom_axis.params.y0 = axbot_adapter.params.y0 - config["bottom_axis"]["h"] 


    #---------------------------------------CONTAINER ---------------------------------------#   

    container = Shape(
        model,
        2,
        "container",
        [
            occ.add_rectangle(
            0,
            bottom_axis.params.y0 - config["container"]["t"] ,
            0,
            config["container"]["r"]  + config["container"]["t"] ,
            2 *config["container"]["t"] + config["container"]["h"]  ,
            )
        ]
    )

    container_hole = occ.add_rectangle(
        0,
        bottom_axis.params.y0 ,
        0,
        config["container"]["r"] ,
        config["container"]["h"],
    )
    occ.cut(container.dimtags, [(2, container_hole)])

    container.mesh_size = config["container"]["mesh_size"]
    container.params.material = config["container"]["material"]
    container.params.T_init = config["container"]["T_init"]



    #--------------------------------------------- ATMOSPHERE ---------------------------------------------#
    atmosphere = occ.add_rectangle(
        0,
        bottom_axis.params.y0 ,
        0,
        config["container"]["r"] ,
        config["container"]["h"],
    )

    #------------------------------- REMOVE THE INDUCTOR INSIDE ATMOSPHERE ----------------------------------#

    ind_cut_bottom = inductor_filling(model,2, **config["inductor_bottom"],name="ind_cut_bottom")

    if config["ind_top"]["include"]:
        ind_cut_top = inductor_filling(model,2, **config["inductor_top"],name="ind_cut_top")


     #--------------------------------------------------------------------------------------------------------#

    if not config["inner_insulation_btm"]["include"]:
        model.remove_shape(inner_insulation_btm)
    if not config["inner_insulation_side"]["include"]:
        model.remove_shape(inner_insulation_side)

    if not config["shield"]["include"]:
        model.remove_shape(shield)

    if not config["top_susceptor"]["include"]:
        model.remove_shape(top_susceptor)

    

    model.synchronize()
    

    # atmosphere, will be used to determine surfaces
    # it is just a helper shape and will be removed later
    shapes = model.get_shapes(2)
    atmosphere = Shape(model, 2, "atmosphere", [atmosphere])
    atmosphere.mesh_size = config["atmosphere"]["mesh_size"]
    atmosphere.params.material = config["atmosphere"]["material"]
    atmosphere.params.T_init = config["atmosphere"]["T_init"]

    for shape in shapes:
        atmosphere.geo_ids = cut(atmosphere.dimtags, shape.dimtags, remove_tool=False)
    
    # set interfaces between shapes, this removes duplicate lines and ensures a consistent mesh
    #print(atmosphere.geo_ids)
    occ.fragment(
        melt.dimtags
        + crystal.dimtags
        + crucible.dimtags
        + seedholder.dimtags
        + shaft_al2o3.dimtags
        + top_axis.dimtags
        + shield.dimtags
        + top_susceptor.dimtags
        + inductor_top.dimtags
        + inductor_bottom.dimtags

        + inner_insulation_btm.dimtags
        + inner_insulation_side.dimtags

        + outer_insulation.dimtags
        + axbot_adapter.dimtags
        + bottom_axis.dimtags
        + container.dimtags
        + ind_cut_bottom.dimtags
        + ind_cut_top.dimtags
        + atmosphere.dimtags,
        [],
    )

    model.synchronize()

    #----------------------------------- extract phase interface ------------------------------------- # 
    if_melt_crystal = Shape(model, 1, "if_melt_crystal", melt.get_interface(crystal))

    #----------------------------------- extract interfaces ------------------------------------------ # 

    if_crystal_seedholder = Shape(model, 1, "if_crystal_seedholder", seedholder.get_interface(crystal))
    if_seedholder_shaft_al2o3 = Shape(model, 1, "if_seedholder_shaft_al2o3", seedholder.get_interface(shaft_al2o3))
    if_shaft_al2o3_top_axis = Shape(model, 1, "if_shaft_al2o3_top_axis", shaft_al2o3.get_interface(top_axis))
    if_top_axis_container = Shape(model, 1, "if_top_axis_container", top_axis.get_interface(container))

    if_crucible_melt = Shape(model, 1, "if_crucible_melt", crucible.get_interface(melt))
    if_outer_insulation_axbot_adapter = Shape(model, 1, "if_outer_insulation_axbot_adapter", outer_insulation.get_interface(axbot_adapter))
    if_axbot_adapter_bottom_axis = Shape(model, 1, "if_axbot_adapter_bottom_axis", axbot_adapter.get_interface(bottom_axis))
    if_bottom_axis_container = Shape(model, 1, "if_bottom_axis_container", bottom_axis.get_interface(container))



    if  config["shield"]["include"]:
        if_shield_top_susceptor = Shape(model, 1, "if_shield_top_susceptor", shield.get_interface(top_susceptor))
    
    if  config["top_susceptor"]["include"]:
        if_crucible_top_susceptor = Shape(model, 1, "if_crucible_top_susceptor", crucible.get_interface(top_susceptor))

    if config["inner_insulation_btm"]["include"]:
        if_crucible_inner_insulation_btm = Shape(model, 1, "if_crucible_inner_insulation_btm", crucible.get_interface(inner_insulation_btm))
        if_inner_insulation_btm_outer_insulation = Shape(model, 1, "if_inner_insulation_btm_outer_insulation", inner_insulation_btm.get_interface(outer_insulation))

    if config["inner_insulation_side"]["include"]:
        if_crucible_inner_insulation_side= Shape(model, 1, "if_crucible_inner_insulation_side", crucible.get_interface(inner_insulation_side))
        if_inner_insulation_side_outer_insulation = Shape(model, 1, "if_inner_insulation_side_outer_insulation", inner_insulation_side.get_interface(outer_insulation))


     #----------------------------------- extract radiation boundaries ------------------------------------------ # 

    bnd_top_axis = Shape(model, 1, "bnd_top_axis", top_axis.get_interface(atmosphere))
    bnd_shaft_al2o3 = Shape(model, 1, "bnd_shaft_al2o3", shaft_al2o3.get_interface(atmosphere))
    bnd_seedholder = Shape(model, 1, "bnd_seedholder", seedholder.get_interface(atmosphere))
    bnd_crystal_side = Shape(model, 1, "bnd_crystal_side", crystal.get_interface(atmosphere))
    bnd_melt = Shape(model, 1, "bnd_melt", melt.get_interface(atmosphere))
    bnd_crucible = Shape(model, 1, "bnd_crucible", crucible.get_interface(atmosphere))
    bnd_outer_insulation = Shape(model, 1, "bnd_outer_insulation", outer_insulation.get_interface(atmosphere))
    bnd_inductor_bottom = Shape(model, 1, "bnd_inductor_bottom", inductor_bottom.get_interface(atmosphere))
    bnd_axbot_adapter= Shape(model, 1, "bnd_axbot_adapter", axbot_adapter.get_interface(atmosphere))
    bnd_bottom_axis = Shape(model, 1, "bnd_bottom_axis", bottom_axis.get_interface(atmosphere))
    bnd_container = Shape(model, 1, "bnd_container", container.get_interface(atmosphere))
    bnd_outer_container = Shape(model, 1, "bnd_outer_container", [container.bottom_boundary,container.top_boundary, container.right_boundary])

    bnd_crystal_top = Shape(model, 1, "bnd_crystal_top", crystal.get_boundaries_in_box(
            [0, config["seedholder"]["r_in"]],
            [
                crystal_y_top,
                crystal_y_top,
            ],
        ),
    )
    bnd_crystal_side -= bnd_crystal_top


    if  config["shield"]["include"]:
        bnd_shield = Shape(model, 1, "bnd_shield", shield.get_interface(atmosphere))

    if  config["top_susceptor"]["include"]:
        bnd_top_susceptor = Shape(model, 1, "bnd_top_susceptor", top_susceptor.get_interface(atmosphere))

    if config["ind_top"]["include"]:
        bnd_inductor_top = Shape(model, 1, "bnd_inductor_top", inductor_top.get_interface(atmosphere))


    if config["inner_insulation_side"]["include"]:
        bnd_inner_insulation_side = Shape(model, 1, "bnd_inner_insulation_side", inner_insulation_side.get_interface(atmosphere))

    if config["inner_insulation_btm"]["include"] and not(config["inner_insulation_side"]["include"]):
        bnd_inner_insulation_btm = Shape(model, 1, "bnd_inner_insulation_btm", inner_insulation_btm.get_interface(atmosphere))



    bnd_inductor_bottom_inside = Shape(model, 1, "bnd_inductor_bottom_inside", inductor_bottom.get_interface(ind_cut_bottom))
    if config["ind_top"]["include"]:
        bnd_inductor_top_inside = Shape(model, 1, "bnd_inductor_top_inside", inductor_top.get_interface(ind_cut_top))


    # symmetry axis
    symmetry_axis = Shape(model, 1, "symmetry_axis", model.symmetry_axis)


    if not include_atmosphere:
        model.remove_shape(atmosphere)

    model.remove_shape(ind_cut_bottom)
    if config["ind_top"]["include"]:
        model.remove_shape(ind_cut_top)
    model.make_physical()

    #------------------------------------------------------------------------ #

    #------------------------------------ mesh settings ------------------------------------ #



    #---refine mesh for Iridium and steel parts to 3 elements per sking depth δ 

    def delta(el_cond,f):
        m0 = 4 * np.pi * 1e-7 
        return np.sqrt(1/(np.pi * m0 * el_cond* f))

    frequency = 14.2e+3 

    steel_res = delta(1.37e+6 , frequency) /3 
    graphite_CZ3R6300_res = delta(5.88e+4 , frequency) /3 
    insulator_material_res =  delta(1.95e-4  , frequency) /3 
    graphite_FU8957_res = delta(7.14e+4 , frequency) /3 
    iridium_res = delta(2380952 , frequency) /3 

    print("steel_skin_depth :",delta(1.37e+6 , frequency) )
    print("Insulator_skin_depth :",delta(1.95e-4  , frequency) )
    print("graphite_skin_depth :",delta(7.14e+4 , frequency) )
    print("iridium_skin_depth :",delta(2380952 , frequency))

    shapes = model.get_shapes(2)
    for shape in shapes:
        if shape.params.material == "air" :
            gmsh.model.setColor(shape.dimtags, 173, 216, 230)  # light blue atmosphere
        elif shape.params.material == "graphite-FU8957":
            gmsh.model.setColor(shape.dimtags, 50,50,50)  # dark grey
            shape.mesh_size = graphite_FU8957_res 
        elif shape.params.material == "csi_liquid":
            gmsh.model.setColor(shape.dimtags, 101, 67, 33) # deep brown
        elif shape.params.material == "csi_solid":
            gmsh.model.setColor(shape.dimtags, 184, 66, 33) # deep red 
        elif shape.params.material == "al2o3":
            gmsh.model.setColor(shape.dimtags, 0,128,0) # green
        elif shape.params.material == "steel-1.4541":
            gmsh.model.setColor(shape.dimtags, 0,0,255) # blue
            shape.mesh_size = steel_res
        elif shape.params.material == "copper-inductor":
            gmsh.model.setColor(shape.dimtags, 255,165,0) # blue
        elif shape.params.material == "insulation":
            gmsh.model.setColor(shape.dimtags, 255, 219, 88) #deep yellow
            # shape.mesh_size = insulator_material_res
        elif shape.params.material == "Iridium":
            gmsh.model.setColor(shape.dimtags, 50,50,50)  # dark grey
            shape.mesh_size = iridium_res

    model.deactivate_characteristic_length()
    model.set_const_mesh_sizes()

    # add linear mesh control to ensure smooth transition in mesh sizes
    max_meshsize = config["atmosphere"]["mesh_size"]
    for shape in model.get_shapes(2):
        MeshControlLinear(model, shape, shape.mesh_size, max_meshsize)


    # mesh_size = char_length + fact * distance ^ exp

        #---at melt ---#
    MeshControlExponential(model, if_melt_crystal, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(model, bnd_melt, melt.mesh_size / 5, exp=1.6, fact=3)
    MeshControlExponential(model, if_crucible_melt, melt.mesh_size / 5, exp=1.6, fact=3)


    #MeshControlExponential(model, bnd_outer_insulation_bottom, outer_insulation_bottom.mesh_size / 7, exp=1.2, fact=3)
    #---at crucible ---#
    MeshControlExponential(model, bnd_crucible, melt.mesh_size / 5, exp=1.6, fact=3)
     #---at inductors ---#   
    MeshControlExponential(model, inductor_bottom, inductor_bottom.mesh_size*3, exp=1.4 )

    if config["ind_top"]["include"]:
        MeshControlExponential(model, inductor_top, inductor_top.mesh_size*3, exp=1.4)

    if  config["top_susceptor"]["include"]:
        MeshControlExponential(model, bnd_top_susceptor, top_susceptor.mesh_size / 2, exp=1.2, fact=3)

    if  config["shield"]["include"]:
        MeshControlExponential(model, bnd_shield, shield.mesh_size / 2, exp=1.2, fact=3)






    model.generate_mesh(**config["mesh"])



    
    if visualize:
        model.show()
    print(model)
    model.write_msh(sim_dir + "/mesh.msh")

    model.close_gmsh()
    return model


if __name__ == "__main__":
    sim_dir = "./"

    with open("config_geometry.yml") as f:
        config_geo = yaml.safe_load(f)
    model = geometry(config_geo, sim_dir, visualize=True , include_atmosphere=True)