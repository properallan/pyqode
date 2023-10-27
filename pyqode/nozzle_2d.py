from pyqode.utils import inch_to_meter, inch2_to_meter2, area_to_radius
import numpy as np
import os
import subprocess
from pyqode.meshing import Gmsh
from pathlib import Path

GMSH = Path(__file__).parent / 'src' / 'gmsh'

def nasa_cdv_nozzle_verification_geometry(npts):
    X = np.linspace(0,10,npts)
    area = np.zeros_like(X)

    for i, x in enumerate(X):
        if ( x < 5.0 ):
            area[i] = 1.75 - 0.75 * np.cos( ( 0.2 * x - 1.0 ) * np.pi )
        else:
            area[i] = 1.25 - 0.25 * np.cos( ( 0.2 * x - 1.0 ) * np.pi )
        
    X = inch_to_meter(X)
    area = inch2_to_meter2(area)
    
    return X, area

def gen_wall_mesh(xn, Sn, thickness, Nx, Ny, inflation_rate, output_file, gmshsolver, symmetry=False):
    x = xn
    r = area_to_radius(Sn)

    geo = Gmsh()

    geo.writeLine('cl__1 = 1;\n')

    wp = []
    wp_2 = []
    for xi,ri in zip(x,r):
        wp.append(  geo.addPoint(xi, ri, 0) )
        if not symmetry:
            wp_2.append(  geo.addPoint(xi, -ri, 0) )

    l1 = geo.addSpline(wp)
    geo.transfiniteLine(l1, Nx, 1)

    # add solid wall
    dyWall = thickness
    l2, l3, l4 = geo.extrudeMesh([0, dyWall, 0], l1, layers=Ny)

    #geo.gfile += 'Physical Surface(6) = {5};\n'

    if not symmetry:
        l1_2 = geo.addSpline(wp_2)
        geo.transfiniteLine(l1_2, Nx, 1)
        l2_2, l3_2, l4_2 = geo.extrudeMesh([0, -dyWall, 0], l1_2, layers=Ny)


    geo.addPhysicalLine("INNERWALL",l1)
    geo.addPhysicalLine("OUTERWALL",l2)
    geo.addPhysicalLine("INWALL",l3)
    geo.addPhysicalLine("OUTWALL",l4)

    if not symmetry:
        geo.addPhysicalLine("INNERWALL_2",l1_2)
        geo.addPhysicalLine("OUTERWALL_2",l2_2)
        geo.addPhysicalLine("INWALL_2",l3_2)
        geo.addPhysicalLine("OUTWALL_2",l4_2)#

    #geo.gfile += 'Recombine Surface {5};\nPhysical Surface(6) = {5};\n'
    if not symmetry:
        geo.gfile += 'Physical Surface(10) = {5, 9};\n'
    else:
        geo.gfile += 'Physical Surface(10) = {5};\n'
    
    output_file = Path(output_file)

    os.makedirs(output_file.parent, exist_ok=True)

    geofilename = output_file.with_suffix('.geo')

    geo.writeFile(geofilename)
    
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','su2','-o', output_file.with_suffix('.su2')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.msh')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.vtk')])
    p.wait()
            
    return (output_file).resolve().with_suffix('.msh').__str__()


def gen_nozzle_mesh_from_config(config, gmshsolver=GMSH):
    xn = config['domain_x']
    rn = area_to_radius(config['domain_area'])
    symmetry = config['domain_axisymmetric']
    Nx = config['domain_Nx']
    Ny = config['domain_Ny']
    inflation_rate = config['domain_inflation_rate']
    using='Progression'
    hwall_n=config['domain_hwall_n']
    beta=config['domain_bl_beta']
    nblayers=config['domain_bl_nblayers']
    output_file = Path(config['working_dir']).resolve() / config['domain_mesh']

    x0 = xn[0]
    L = xn[-1]-xn[0]
    x = xn
    r = rn

    geo = Gmsh()

    geo.writeLine('cl__1 = 1;\n')

    wp = []
    for xi,ri in zip(x,r):
        wp.append(  geo.addPoint(xi, ri, 0) )

    upper_wall_points = wp
    upper_wall = geo.addLine(wp)
    geo.transfiniteLine(upper_wall, Nx, 1)

    if symmetry == 'YES':
        bottom_wall_points = []
        for xi in x:
            bottom_wall_points.append(  geo.addPoint(xi, 0, 0) )
        bottom_wall = geo.addLine(bottom_wall_points)
    elif symmetry == 'NO':
        bottom_wall_points = [i + len(wp) for i in wp]
        bwp = bottom_wall_points
        bottom_wall = geo.addSymmetryCurve([0,1,0], upper_wall)
        
    geo.transfiniteLine(bottom_wall, Nx, 1)

    inflow = geo.addLine([upper_wall_points[0], bottom_wall_points[0]])
    geo.transfiniteLine(inflow, Ny, inflation_rate, using=using)

    outflow = geo.addLine([upper_wall_points[-1], bottom_wall_points[-1]])
    geo.transfiniteLine(outflow, Ny, inflation_rate, using=using)

    loop1 = geo.addLoop([upper_wall,outflow,-bottom_wall,-inflow])
    surface1 = geo.addPlaneSurface([loop1])

    geo.transfiniteSurface(surface1, [upper_wall_points[0],upper_wall_points[-1],bottom_wall_points[-1],bottom_wall_points[0]])
    geo.recombineSurface(surface1)

    # add boundary
    geo.addPhysicalLine("UPPER_WALL",upper_wall)
    if symmetry == 'YES':
        geo.addPhysicalLine("SYMMETRY",bottom_wall)
    elif symmetry == 'NO':
        geo.addPhysicalLine("BOTTOM_WALL",bottom_wall)
    geo.addPhysicalLine("INFLOW",inflow)
    geo.addPhysicalLine("OUTFLOW",outflow)

    geo.addPhysicalSurface(surface1)

    output_file = Path(output_file)

    os.makedirs(output_file.parent, exist_ok=True)

    geofilename = output_file.with_suffix('.geo')
    
    if symmetry == 'NO' and nblayers>0:    
        geo.addBoundaryLayer(edges=[upper_wall], nodes=[wp[0],wp[-1]], hwall_n=hwall_n, ratio=0, thickness=0, beta=beta, nblayers=nblayers)
        geo.addBoundaryLayer(edges=[bottom_wall], nodes=[bwp[0],bwp[-1]], hwall_n=hwall_n, ratio=0, thickness=0, beta=beta, nblayers=nblayers)
    elif symmetry == 'YES' and nblayers>0:
        geo.addBoundaryLayer(edges=[upper_wall], nodes=[wp[0],wp[-1]], hwall_n=hwall_n, ratio=0, thickness=0, beta=beta, nblayers=nblayers)

    geo.writeFile(geofilename)

    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','su2','-o', output_file.with_suffix('.su2')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.msh')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.vtk')])
    p.wait()

    return output_file.resolve().with_suffix('.su2').__str__()


def gen_nozzle_mesh(xn, rn, Nx, Ny, inflation_rate, output_file, gmshsolver, symmetry=False, using='Progression', hwall_n=0.00001, thickness=0.01, ratio=1.2, beta=1.001, nblayers=10):
    x0 = xn[0]
    L = xn[-1]-xn[0]
    x = xn
    r = rn

    geo = Gmsh()

    geo.writeLine('cl__1 = 1;\n')

    wp = []
    for xi,ri in zip(x,r):
        wp.append(  geo.addPoint(xi, ri, 0) )

    upper_wall_points = wp
    upper_wall = geo.addLine(wp)
    geo.transfiniteLine(upper_wall, Nx, 1)

    if symmetry:
        bottom_wall_points = []
        for xi in x:
            bottom_wall_points.append(  geo.addPoint(xi, 0, 0) )
        bottom_wall = geo.addLine(bottom_wall_points)
    else:
        bottom_wall_points = [i + len(wp) for i in wp]
        bottom_wall = geo.addSymmetryCurve([0,1,0], upper_wall)
        
    geo.transfiniteLine(bottom_wall, Nx, 1)

    inflow = geo.addLine([upper_wall_points[0], bottom_wall_points[0]])
    geo.transfiniteLine(inflow, Ny, inflation_rate, using=using)

    outflow = geo.addLine([upper_wall_points[-1], bottom_wall_points[-1]])
    geo.transfiniteLine(outflow, Ny, inflation_rate, using=using)

    loop1 = geo.addLoop([upper_wall,outflow,-bottom_wall,-inflow])
    surface1 = geo.addPlaneSurface([loop1])

    geo.transfiniteSurface(surface1, [upper_wall_points[0],upper_wall_points[-1],bottom_wall_points[-1],bottom_wall_points[0]])
    geo.recombineSurface(surface1)

    # add boundary
    geo.addPhysicalLine("UPPER_WALL",upper_wall)
    if symmetry:
        geo.addPhysicalLine("SYMMETRY",bottom_wall)
    else:
        geo.addPhysicalLine("BOTTOM_WALL",bottom_wall)
    geo.addPhysicalLine("INFLOW",inflow)
    geo.addPhysicalLine("OUTFLOW",outflow)

    geo.addPhysicalSurface(surface1)

    output_file = Path(output_file)

    os.makedirs(output_file.parent, exist_ok=True)

    geofilename = output_file.with_suffix('.geo')

    geo.addBoundaryLayer(edges=[upper_wall], nodes=[wp[0],wp[-1]], hwall_n=hwall_n, ratio=ratio, thickness=thickness, beta=beta, nblayers=nblayers)

    geo.writeFile(geofilename)

    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','su2','-o', output_file.with_suffix('.su2')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.msh')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.vtk')])
    p.wait()

    return output_file.resolve().with_suffix('.msh').__str__()

def gen_nozzle_mesh_3d(xn, rn, Nx, Ny, inflation_rate, output_file, gmshsolver, symmetry=False):
    x0 = xn[0]
    L = xn[-1]-xn[0]
    x = xn
    r = rn

    geo = Gmsh()

    geo.writeLine('cl__1 = 1;\n')

    wp = []
    for xi,ri in zip(x,r):
        wp.append(  geo.addPoint(xi, ri, 0) )

    upper_wall_points = wp
    upper_wall = geo.addLine(wp)
    geo.transfiniteLine(upper_wall, Nx, 1)

    if symmetry:
        bottom_wall_points = []
        for xi in x:
            bottom_wall_points.append(  geo.addPoint(xi, 0, 0) )
        bottom_wall = geo.addLine(bottom_wall_points)
    else:
        bottom_wall_points = [i + len(wp) for i in wp]
        bottom_wall = geo.addSymmetryCurve([0,1,0], upper_wall)
        
    geo.transfiniteLine(bottom_wall, Nx, 1)

    inflow = geo.addLine([upper_wall_points[0], bottom_wall_points[0]])
    geo.transfiniteLine(inflow, Ny, inflation_rate-1)

    outflow = geo.addLine([upper_wall_points[-1], bottom_wall_points[-1]])
    geo.transfiniteLine(outflow, Ny, inflation_rate-1)

    loop1 = geo.addLoop([upper_wall,outflow,-bottom_wall,-inflow])
    surface1 = geo.addPlaneSurface([loop1])

    geo.transfiniteSurface(surface1, [upper_wall_points[0],upper_wall_points[-1],bottom_wall_points[-1],bottom_wall_points[0]])
    geo.recombineSurface(surface1)

    # add boundary
    geo.addPhysicalLine("UPPER_WALL",upper_wall)
    if symmetry:
        geo.addPhysicalLine("SYMMETRY",bottom_wall)
    else:
        geo.addPhysicalLine("BOTTOM_WALL",bottom_wall)
    geo.addPhysicalLine("INFLOW",inflow)
    geo.addPhysicalLine("OUTFLOW",outflow)

    geo.addPhysicalSurface(surface1)

    output_file = Path(output_file)

    os.makedirs(output_file.parent, exist_ok=True)

    geofilename = output_file.with_suffix('.geo')

    geo.writeFile(geofilename)

    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','su2','-o', output_file.with_suffix('.su2')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.msh')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.vtk')])
    p.wait()

    return output_file.resolve().with_suffix('.msh').__str__()

def gen_nozzle_mesh_with_solid(xn, Sn, Nx, Ny, inflation_rate, output_file, gmshsolver):
    x0 = xn[0]
    L = xn[-1]-xn[0]
    x = xn
    r = area_to_radius(Sn)

    geo = Gmsh()

    geo.writeLine('cl__1 = 1;\n')

    wp = []
    for xi,ri in zip(x,r):
        wp.append(  geo.addPoint(xi, ri, 0) )

    upper_wall_points = wp
    upper_wall = geo.addLine(wp)
    geo.transfiniteLine(upper_wall, Nx, 1)

    bottom_wall_points = [i + len(wp) for i in wp]
    bottom_wall = geo.addSymmetryCurve([0,1,0], upper_wall)
    geo.transfiniteLine(bottom_wall, Nx, 1)

    inflow = geo.addLine([upper_wall_points[0], bottom_wall_points[0]])
    geo.transfiniteLine(inflow, Ny, inflation_rate-1)

    outflow = geo.addLine([upper_wall_points[-1], bottom_wall_points[-1]])
    geo.transfiniteLine(outflow, Ny, inflation_rate-1)

    loop1 = geo.addLoop([upper_wall,outflow,-bottom_wall,-inflow])
    surface1 = geo.addPlaneSurface([loop1])

    geo.transfiniteSurface(surface1, [upper_wall_points[0],upper_wall_points[-1],bottom_wall_points[-1],bottom_wall_points[0]])
    geo.recombineSurface(surface1)

    # add boundary
    geo.addPhysicalLine("UPPER_WALL",upper_wall)
    geo.addPhysicalLine("BOTTOM_WALL",bottom_wall)
    geo.addPhysicalLine("INFLOW",inflow)
    geo.addPhysicalLine("OUTFLOW",outflow)

    geo.addPhysicalSurface(surface1)

    output_file = Path(output_file)

    os.makedirs(output_file.parent, exist_ok=True)

    geofilename = output_file.with_suffix('.geo')

    geo.writeFile(geofilename)

    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','su2','-o', output_file.with_suffix('.su2')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.msh')])
    p.wait()
    p = subprocess.Popen([f'{gmshsolver}', geofilename,'-0','-2','-format','msh','-o', output_file.with_suffix('.vtk')])
    p.wait()

def gen_multizone_mesh(input_meshes, outputfile, gmshsolver):
        import subprocess
        #-------------------------------------------------------------------#
        # enter the names of .geo files here, without .geo -ending
        #geo_files = ['channel_only','top_block_only','bottom_block_only']

        geo_files = input_meshes

        # output filename
        #out_mesh_name = 'multizone.su2'

        # mesh dimension
        ndim = 2
      
        #-------------------------------------------------------------------#
        # merge .su2 meshes of all zones to single mesh
        nzones = len(geo_files)

        with open(Path(outputfile).with_suffix('.su2'),'w') as out_mesh:
            out_mesh.write('NDIME=' + str(ndim) +  '\n')
            out_mesh.write('NZONE=' + str(nzones) + '\n')

            for idx,file in enumerate(geo_files):
                out_mesh.write('\nIZONE=' + str(idx+1) + '\n\n')
                #input.pipe
                with open(Path(geo_files[idx]).with_suffix('.su2').resolve() ,'r') as in_mesh:
                    out_mesh.write(in_mesh.read())

        
        return Path(outputfile).resolve().__str__()