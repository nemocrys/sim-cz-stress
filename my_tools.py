import yaml
from objectgmsh import *
import numpy as np

occ = gmsh.model.occ


def inductor( 
    model,
    dim,
    d,
    d_in,
    X0,
    g=0,
    n=1,
    char_l=0,
    T_init=273.15,
    material="",
    name="inductor",
):
    """2D inductor, defined as a couple of circles
 ### https://github.com/nemocrys/test-cz-induction#overview ### 
    Args:
        model (Model): objectgmsh model
        dim (int): dimension
        d (float): diameter of windings
        d_in (flaot): inner diameter of windings (internal cooling)
        X0 (float): origin, center of bottom winding
        g (float, optional): gap between windings. Defaults to 0.
        n (int, optional): number of windings. Defaults to 1.
        char_l (int, optional): mesh size characteristic length. Defaults to 0.
        T_init (float, optional): initial temperature. Defaults to 273.15.
        material (str, optional): material name. Defaults to "".
        name (str, optional): shape name. Defaults to "inductor".

    Returns:
        Shape: objectgmsh shape
    """
    # X0: center of bottom winding
    ind = Shape(model, dim, name)
    ind.params.d = d
    ind.params.d_in = d_in
    ind.params.g = g
    ind.params.n = n
    ind.params.X0 = X0
    ind.params.T_init = T_init
    ind.params.material = material
    ind.params.area = np.pi * (d ** 2 - d_in ** 2) / 4
    if char_l == 0:
        ind.mesh_size = d / 10
    else:
        ind.mesh_size = char_l

    x = X0[0]
    y = X0[1]
    for _ in range(n):
        circle_1d = factory.addCircle(x, y, 0, d / 2)
        circle = factory.addSurfaceFilling(factory.addCurveLoop([circle_1d]))
        hole_1d = factory.addCircle(x, y, 0, d_in / 2)
        hole = factory.addSurfaceFilling(factory.addCurveLoop([hole_1d]))
        factory.synchronize()
        factory.cut([(2, circle)], [(2, hole)])
        if dim == 3:
            circle = rotate(circle)
        ind.geo_ids.append(circle)
        y += g + d

    return ind

def inductor_filling( 
    model,
    dim,
    d,
    d_in,
    X0,
    g=0,
    n=1,
    char_l=0,
    T_init=273.15,
    material="",
    name="inductor_filling",
):
    # X0: center of bottom winding
    filling = Shape(model, dim, name)
    filling.params.d = d
    filling.params.d_in = d_in
    filling.params.g = g
    filling.params.n = n
    filling.params.X0 = X0
    filling.params.T_init = T_init
    filling.params.material = material
    
    x = X0[0]
    y = X0[1]
    for _ in range(n):
        #circle_1d = factory.addCircle(x, y, 0, d / 2)
        #circle = factory.addSurfaceFilling(factory.addCurveLoop([circle_1d]))
        hole_1d = factory.addCircle(x, y, 0, d_in / 2)
        hole = factory.addSurfaceFilling(factory.addCurveLoop([hole_1d]))
        factory.synchronize()
        if dim == 3:
            circle = rotate(hole)
        filling.geo_ids.append(hole)
        y += g + d


    return filling



def crystal_shape(
    model,
    dim,  
    l_s,
    d_s,
    h = 0.0,
    starting_point = [0.0,0.0,0.0],
    name= "crystal", material="", mesh_size=0.05, T_init=1511):


    """starts the crystal shape from top left point. 
    model (Model): objectgmsh model
    dim (int): dimension
    d_s (float): diameters from top to bottom
    l_s (float): lenghts from top to bottom
    starting_point (float): position of the crystal
    h (float) : the cut height position FROM TOP TO BOTTOM
    name (str, optional): shape name. Defaults to "crystal"
    material (str, optional): material name. Defaults to ""
    T_init (float, optional): initial temperature. Defaults to 1511"""


    crystal = Shape(model, dim, name)
    points = []
    
    starting_h = sum(l_s)
    d_s = np.array(d_s) 
    l_s = np.array(l_s) 

    
    crystal.params.h = starting_h
    dx = starting_point[0]
    dy = starting_point[1]


    if h != 0.0:  # adjust to melt surface
        crystal.params.h =  h
        dy -= starting_h - h
        
     
    points.append(occ.add_point(dx, dy,0)) #starting point
    points.append(occ.add_point(dx, starting_h +dy, 0)) # backbone
    y= 0
    # remaining points
    points.append(occ.add_point(d_s[0]+dx, starting_h+dy , 0))

    for i in range(len(l_s)):
        x = d_s[i] 
        y = starting_h - np.sum(l_s[:i+1])
        z = 0
        points.append(occ.add_point(x +dx, y + dy, z))

    lines = [occ.add_line(points[i - 1], points[i]) for i in range(len(points))]
    loop = occ.add_curve_loop(lines)
    volume = occ.add_surface_filling(loop)


    if h != 0.0 and h < starting_h:
        if h> starting_h:
            raise ValueError("The given height excess the crystal size")
        
        crystal_hole = occ.add_rectangle(
            dx,
            dy,
            0.0,
            np.max(d_s),
            starting_h - h,
        )
        factory.cut([(2, volume)] , [(2, crystal_hole)])
    

    crystal.params.material = material
    crystal.mesh_size = mesh_size
    crystal.params.T_init = T_init
    crystal.params.seed_r = d_s[0]/2
    



    if dim==3:
        volume = rotate(volume)

    crystal.geo_ids = [volume]

    return crystal







def crystal_shape_stress_calc(
    model,
    dim,  
    l_s,
    d_s,
    btm_x ,
    btm_y ,
    meniscus_cut = False,
    h = 0.0,
    starting_point = [0.0,0.0,0.0],
    name= "crystal", material="", mesh_size=0.05, T_init=1511):


    """starts the crystal shape from top left point. 
    model (Model): objectgmsh model
    dim (int): dimension
    d_s (float): diameters from top to bottom
    l_s (float): lenghts from top to bottom
    starting_point (float): position of the crystal
    h (float) : the cut height position FROM TOP TO BOTTOM
    name (str, optional): shape name. Defaults to "crystal"
    material (str, optional): material name. Defaults to ""
    T_init (float, optional): initial temperature. Defaults to 1511"""


    crystal = Shape(model, dim, name)
    points = []
    starting_h = sum(l_s)
    d_s = np.array(d_s) 
    l_s = np.array(l_s) 


    
    crystal.params.h = starting_h
    dx = starting_point[0]
    dy = starting_point[1]


    if h != 0.0:  # adjust to melt surface
        crystal.params.h =  h
        dy -= starting_h - h
        
   
    
    points = []
    #btm_left = occ.add_point(0, dy, 0)

    top_left = occ.add_point(dx, 0, starting_h +dy)

    points.append(occ.add_point(dx, 0, starting_h +dy)) # backbone
    # remaining points
    points.append(occ.add_point(d_s[0]+dx,0 ,  starting_h+dy))

    for i in range(len(l_s)):
        x = d_s[i] 
        z = starting_h - np.sum(l_s[:i+1])
        y = 0
        points.append(occ.add_point(x +dx, y, z+ dy))


    lines = [occ.add_line(points[i ], points[i+1]) for i in range(len(points) - 1 )]



    meniscus_start = occ.add_point(btm_x[0], 0, btm_y[0] )


    spline_points =  [occ.add_point(btm_x[i], 0, btm_y[i]) for i in range(len(btm_x))]
    meniscus_spline = occ.add_spline(spline_points)

    left = occ.add_line(meniscus_start, top_left)


    loop = factory.addCurveLoop([meniscus_spline,left] + lines)
    
    volume = occ.add_surface_filling(loop)   

    

    if h != 0.0 and h < starting_h:
        if h> starting_h:
            raise ValueError("The given height excess the crystal size")
        
        crystal_hole = occ.add_rectangle(
            dx,
            dy,
            0.0,
            np.max(d_s),
            starting_h - h,
        )
        factory.cut([(2, volume)] , [(2, crystal_hole)])
    

    crystal.params.material = material
    crystal.mesh_size = mesh_size
    crystal.params.T_init = T_init
    crystal.params.seed_r = d_s[0]/2
    


    if dim==3:
        #volume = rotate(volume)
        factory.revolve([(2, volume)], 0, 0, 0, 0, 0, 1, 2 * np.pi)
    crystal.geo_ids = [volume]

    return crystal



#------------------------ updated ---------------------------- #


def crystal_shape_from_2D(
    model,
    dim,  
    l_s,
    d_s,
    btm_x ,
    btm_y ,
    h = 0.0,
    starting_point = [0.0,0.0,0.0],
    name= "crystal", material="", mesh_size=0.05, T_init=1511):


    """starts the crystal shape from top left point. 
    model (Model): objectgmsh model
    dim (int): dimension
    d_s (float): diameters from top to bottom
    l_s (float): lenghts from top to bottom
    starting_point (float): position of the crystal
    h (float) : the cut height position FROM TOP TO BOTTOM
    name (str, optional): shape name. Defaults to "crystal"
    material (str, optional): material name. Defaults to ""
    T_init (float, optional): initial temperature. Defaults to 1511"""


    crystal = Shape(model, dim, name)
    points = []
    total_length = sum(l_s)
    d_s = np.array(d_s) 
    l_s = np.array(l_s) 
    crystal.params.h = total_length
    dx = starting_point[0]
    dz = starting_point[1]

    if h != 0.0:  # adjust to melt surface
        crystal.params.h =  h
        dz -= total_length - h



    top_left = occ.add_point(dx, 0, total_length +dz)
    points.append(occ.add_point(dx, 0, total_length +dz)) # top left
    # remaining points
    points.append(occ.add_point(d_s[0]+dx,0 ,  total_length+dz)) #top right

    for i in range(len(l_s)): # vale mia if 
        x = d_s[i] 
        z = total_length - np.sum(l_s[:i+1])
        if (z+ dz) < starting_point[1] :
            points.append(occ.add_point(btm_x[-1], 0, btm_y[-1] ))
            break
        y = 0
        points.append(occ.add_point(x +dx, y, z+ dz))



    lines = [occ.add_line(points[i ], points[i+1]) for i in range(len(points) - 1 )] # connects the outline points
    meniscus_start = occ.add_point(btm_x[0], 0, btm_y[0] )
    spline_points =  [occ.add_point(btm_x[i], 0, btm_y[i]) for i in range(len(btm_x))]
    meniscus_spline = occ.add_spline(spline_points)
    left = occ.add_line(meniscus_start, top_left)
    loop = factory.addCurveLoop([meniscus_spline,left] + lines)

    volume = occ.add_surface_filling(loop)   



    if dim==3:
        #volume = rotate(volume)
        factory.revolve([(2, volume)], 0, 0, 0, 0, 0, 1, 2 * np.pi)
    crystal.geo_ids = [volume]

    crystal.params.material = material
    crystal.mesh_size = mesh_size
    crystal.params.T_init = T_init

    return crystal




def V_trunc_cone(h,r,R): # 2D
    return  np.pi *  h / 3 * (R**2 + R*r + r**2) 
def V_cylinder(h,r): # 2D
    return np.pi *h*r**2


 
def crystal_volume(h, l_s, d_s):
    total_length = np.sum(l_s)
    total_volume = 0.0
    cum_height = 0.0

    for i in range(len(l_s)):
        segment_top = cum_height + l_s[i]

        if h >= segment_top:  # full segment is included
            if i == 0 or d_s[i] == d_s[i-1]:
                total_volume += V_cylinder(l_s[i], d_s[i])
            else:
                total_volume += V_trunc_cone(l_s[i], d_s[i-1], d_s[i])
        elif h > cum_height:  # partial segment
            remaining = h - cum_height
            if i == 0 or d_s[i] == d_s[i-1]:
                total_volume += V_cylinder(remaining, d_s[i])
            else:
                # Interpolate diameter at the cutoff height
                d_start = d_s[i-1]
                d_end = d_s[i]
                slope = (d_end - d_start) / l_s[i]
                d_cut = d_start + slope * remaining
                total_volume += V_trunc_cone(remaining, d_start, d_cut)
            break
        else:
            break  # no more volume to add

        cum_height = segment_top

    return total_volume
