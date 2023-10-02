import numpy as np
import os
import subprocess
import shutil
import timeit
import pandas as pd
import meshio
import pyvista as pv

def saveVtkMesh(vtkfile):
    #vtkfile = '/home/ppiper/Dropbox/local/thermal_nozzle/src_air/yeom/outputs/fluid.vtk'
    vtk = pv.read(vtkfile)
    vtk.clear_data()
    dot_idx = vtkfile.rfind('.')
    vtkfile_mesh = vtkfile[:dot_idx]+'_mesh.vtk'
    vtk.save(vtkfile_mesh)
    return 0

class gmsh:
    pcount = 0
    lcount = 0
    loopcount = 0
    scount = 0
    phycount = 0
    gfile = ''

    def addPoint(self, x, y, z):
        self.pcount += 1
        self.gfile += 'Point('+str(self.pcount)+') = {'+str(x)+','+str(y)+','+str(z)+', cl__1};\n'
        return self.pcount

    def extrudeMesh(self, direction, line, layers):
        d = direction
        self.lcount += 3
        self.gfile += 'Extrude { ' + str(d[0]) + ', ' + str(d[1]) + ', ' + str(d[2]) + '} { Curve{'+ str(line) +'}; Layers {' + str(layers) + '}; Recombine;}\n'

        return self.lcount-2, self.lcount-1, self.lcount

    def addLine(self, points):
        self.lcount += 1
        self.gfile += 'Line('+str(self.lcount)+') = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.lcount

    def addSpline(self, points):
        self.lcount += 1
        self.gfile += 'Spline('+str(self.lcount)+') = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.lcount

    def addLoop(self, lines):
        self.loopcount += 1
        self.gfile += 'Line Loop('+str(self.loopcount)+') = {'
        for i,l in enumerate(lines):
            self.gfile += str(l)
            if i < np.size(lines)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.loopcount

    def addPlaneSurface(self, loops):
        self.scount += 1
        self.gfile += 'Plane Surface('+str(self.scount)+') = {'
        for i,l in enumerate(loops):
            self.gfile += str(l)
            if i < np.size(loops)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
        return self.scount
            
    def transfiniteLine(self, line, points, progression=1):
        self.gfile += 'Transfinite Line {' + str(line) + '} = ' + str(points) + ' Using Bump ' + str(progression) + ';\n'

    def transfiniteSurface(self, surface, points):
        self.gfile += 'Transfinite Surface {' + str(surface) + '} = {'
        for i,p in enumerate(points):
            self.gfile += str(p)
            if i < np.size(points)-1:
                self.gfile += ', '
        self.gfile += '};\n'
        
    def recombineSurface(self, surface):
        self.gfile += 'Recombine Surface {' + str(surface) + '};\n'

    def addPhysicalLine(self, boundName, line):
        self.phycount +=1
        self.gfile += 'Physical Line(\"' + boundName +'\") = {' + str(line) + '};\n'

    def addPhysicalSurface(self, surface):
        self.phycount +=1
        self.gfile += 'Physical Surface(' + str(self.phycount) +') = {' + str(surface) + '};\n'
        
    def writeLine(self, s):
        self.gfile += s

    def print(self):
        print(self.gfile)

    def writeFile(self, filename):
        f = open(filename,'w')
        f.write(self.gfile)
        f.close()

    def addSymmetryCurve(self, direction, line):
        self.lcount += 1
        d = direction
        self.gfile += f"Symmetry {{{d[0]}, {d[1]}, {d[2]}, 0}} {{ Duplicata {{ Curve{{{line}}}; }} }}\n"

        return self.lcount

def getHLine(field,prop, idx, N, Ny):
        centers = field.cell_centers()
        field = field.probe(centers)
        # the field at the idx-th hvline
        return field[prop][idx:N+idx:Ny], field.points[idx:N+idx:Ny]

def getVLine(field, prop, idx, N, Ny):
    # get field at the idx-th vline
    # N is the 
    centers = field.cell_centers()
    field = field.probe(centers)
    return field[prop][idx*Ny:idx*Ny + Ny], field.points[idx*Ny:idx*Ny + Ny]


class nozzle:
    def __init__(self):
        self.mesh = {}
        self.vtk = {}
        self.hprop = {}
        self.vprop = {}

    def setX(self, x):
        self.xn = x

    def setS(self, s):
        self.Sn = s
        
    def setBC(self, p0in, T0in, Min, pb):
        self.p0in = p0in
        self.T0in = T0in
        self.pb = pb
        self.Min = Min

    def setFluid(self, R, gamma):
        self.R = R
        self.gamma = gamma

    def setupQ1DSolver(self, itamx, itprint, CFL, tol,
                  tscheme, fscheme, dttype, dim):
        self.itmax = itamx
        self.itprint = itprint
        self.CFL = CFL
        self.tol = tol
        self.tscheme = tscheme
        self.fscheme = fscheme
        self.dttype = dttype
        self.dim = dim

    def setupQ1D(self, filepath, filename):
        self.filepath = filepath
        infilepath = filepath + 'inputs/'
        self.q1dinfilepath = infilepath

        self.setupfile = infilepath + filename
        os.makedirs(infilepath, exist_ok=True)

        outfilepath = filepath + 'outputs/'
        self.q1doutfilepath = outfilepath
        
        os.makedirs(outfilepath, exist_ok=True)

        # interpolate domain to NxQ1D accuracy
        self.xnQ1D = np.linspace(self.xn[0], self.xn[-1], self.NxQ1D)
        self.SnQ1D = np.interp(self.xnQ1D, self.xn, self.Sn)

        np.savetxt(infilepath+'xn.txt', self.xnQ1D)
        np.savetxt(infilepath+'Sn.txt', self.SnQ1D)

        f = open(infilepath + filename,'w')
        f.write('# Domain x-coordinates at cell faces (no need for ghost cells)\n')
        f.write(infilepath + 'xn.txt' + "\n")
        f.write('# Area distributuin at cell faces (no need for ghost cells)\n')
        f.write(infilepath + 'Sn.txt' + "\n")
        f.write('# Inlet Total Pressure [Pa]\n')
        f.write(str(self.p0in) + "\n")
        f.write('# Inlet Total Temperature [K]\n')
        f.write(str(self.T0in) + "\n")
        f.write('# Inlet Mach Number\n')
        f.write(str(self.Min) + "\n")
        f.write('# Outlet Static Pressure [Pa]\n')
        f.write(str(self.pb) + "\n")
        f.write('# Gas constant [J/kg/K]\n')
        f.write(str(self.R) + "\n")
        f.write('# Specific heat ratio\n')
        f.write(str(self.gamma) + "\n")
        f.write('# Maximum number of iterations \n')
        f.write(str(self.itmax) + "\n")
        f.write('# Interval to print iterations \n')
        f.write(str(self.itprint) + "\n")
        f.write('# CFL number \n')
        f.write(str(self.CFL) + "\n")
        f.write('# Convergence criteria \n')
        f.write(str(self.tol) + "\n")
        f.write('# Time integration scheme (Euler, RK4) \n')
        f.write(str(self.tscheme) + "\n")
        f.write('# Flux scheme (Basic, Roe, Split, SplitTVD, VanLeer, LaxFriedrichs, StegerWarming, AUSM) \n')
        f.write(str(self.fscheme) + "\n")
        f.write('# Timestep type (Global, Local) \n')
        f.write(str(self.dttype) + "\n")
        f.write('# Form to be solved (Dimensional, Dimensionless) \n')
        f.write(str(self.dim) + "\n")
        f.write('# Set output directory \n')
        f.write( outfilepath + "\n")
        f.close()
        
    def solveQ1D(self, solver):
        self.runtimeQ1D = timeit.default_timer()
        print([solver, self.setupfile])
        p = subprocess.Popen([solver, self.setupfile])
        p.wait()
        self.runtimeQ1D = timeit.default_timer() - self.runtimeQ1D
        print('runtime: ' + str(self.runtimeQ1D) + 's')
        #input()
        

    def genMultizone(self, outputpath, fnames):
        import subprocess
        #-------------------------------------------------------------------#
        # enter the names of .geo files here, without .geo -ending
        #geo_files = ['channel_only','top_block_only','bottom_block_only']

        geo_files = fnames

        # output filename
        #out_mesh_name = 'multizone.su2'

        # mesh dimension
        ndim = 2
        create_zone_meshes = 1

        #-------------------------------------------------------------------#
        # Create .su2 meshes out of .geo files with gmsh
        if create_zone_meshes:
            for file in geo_files:
                #create_mesh = ['./gmsh', '-3', '-smooth', '5' ,'-optimize', outputpath +file+'.geo', '-2', '-f', 'su2', '-o', outputpath+file+'.su2']
                create_mesh = ['./gmsh' , outputpath +file+'.geo', '-2', '-f', 'su2', '-o', outputpath+file+'.su2']
                
                #create_mesh = 'gmsh ' + file+'.geo' + ' -2' + ' -f' + ' su2' + ' -saveall' + ' -o ' + file+'.su2' 
                process = subprocess.check_output(create_mesh, stderr=subprocess.STDOUT)
    
            print("Meshes for each zone created.")

        #-------------------------------------------------------------------#
        # merge .su2 meshes of all zones to single mesh
        nzones = len(geo_files)

        with open(self.multizonefile,'w') as out_mesh:
            out_mesh.write('NDIME=' + str(ndim) +  '\n')
            out_mesh.write('NZONE=' + str(nzones) + '\n')

            for idx,file in enumerate(geo_files):
                out_mesh.write('\nIZONE=' + str(idx+1) + '\n\n')
                #input.pipe
                with open(outputpath + geo_files[idx] + '.su2','r') as in_mesh:
                    out_mesh.write(in_mesh.read())

        #-------------------------------------------------------------------#
        print("Input files:")
        for idx,file in enumerate(geo_files): print(str(idx+1) + ". " + file)
        print("Multi-zone mesh written to file: " + str(self.multizonefile))
 
        

    def genGmshWall(self, outputpath, fname, thickness):
        x = self.xn
        r = (self.Sn/np.pi)**0.5
        Nx = self.NxWall
        Ny = self.NyWall

        geo = gmsh()

        geo.writeLine('cl__1 = 1;\n')

        wp = []
        for xi,ri in zip(x,r):
            wp.append(  geo.addPoint(xi, ri, 0) )

        l1 = geo.addSpline(wp)
        geo.transfiniteLine(l1, Nx, 1)

        # add solid wall
        dyWall = thickness
        l2, l3, l4 = geo.extrudeMesh([0, dyWall, 0], l1, layers=Ny)


        #loop1 = geo.addLoop([l4,l3,-l1,-l2])
        #surface1 = geo.addPlaneSurface([loop1])

        #geo.transfiniteSurface(surface1, [wp[0],wp[-1],p2,p1])
        #geo.recombineSurface(surface1)

        geo.addPhysicalLine("INNERWALL",l1)
        geo.addPhysicalLine("OUTERWALL",l2)
        geo.addPhysicalLine("INWALL",l3)
        geo.addPhysicalLine("OUTWALL",l4)

        #geo.gfile += 'Recombine Surface {5};\nPhysical Surface(6) = {5};\n'
        geo.gfile += 'Physical Surface(6) = {5};\n'

        os.makedirs(outputpath, exist_ok=True)

        fname = fname.split('.cfg')[0]
        geofilename = outputpath + fname + '.geo'

        geo.writeFile(geofilename)
        
        p = subprocess.Popen(['./gmsh', '-3', '-smooth', '5' ,'-optimize', geofilename,'-0','-2','-format','su2','-o', outputpath +fname+'.su2'])
        p.wait()
        p = subprocess.Popen(['./gmsh', '-3', '-smooth', '5' ,'-optimize', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.msh'])
        p.wait()
        p = subprocess.Popen(['./gmsh', '-3', '-smooth', '5' ,'-optimize', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.vtk'])
        p.wait()
                
        #self.meshfilename = outputpath +fname+'.su2'

    def genGmshWall_2(self, outputpath, fname, thickness):
        x = self.xn
        r = (self.Sn/np.pi)**0.5
        Nx = self.NxWall
        Ny = self.NyWall

        geo = gmsh()

        geo.writeLine('cl__1 = 1;\n')

        wp_2 = []
        for xi,ri in zip(x,r):
            wp_2.append(  geo.addPoint(xi, -ri, 0) )

        # add solid wall
        dyWall = thickness

        l1_2 = geo.addSpline(wp_2)
        geo.transfiniteLine(l1_2, Nx, 1)
        l2_2, l3_2, l4_2 = geo.extrudeMesh([0, -dyWall, 0], l1_2, layers=Ny)
            
        #loop1 = geo.addLoop([l4_2,l3_2,-l1_2,-l2_2])
        #surface1 = geo.addPlaneSurface([loop1])

        #geo.transfiniteSurface(surface1, [wp_2[0],wp_2[-1],p2,p1])
        #geo.recombineSurface(surface1)

        geo.addPhysicalLine("INNERWALL_2",l1_2)
        geo.addPhysicalLine("OUTERWALL_2",l2_2)
        geo.addPhysicalLine("INWALL_2",l3_2)
        geo.addPhysicalLine("OUTWALL_2",l4_2)

        #geo.gfile += 'Recombine Surface {5};\nPhysical Surface(6) = {5};\n'
        geo.gfile += 'Physical Surface(6) = {5};\n'

        os.makedirs(outputpath, exist_ok=True)

        fname = fname.split('.cfg')[0]
        geofilename = outputpath + fname + '.geo'

        geo.writeFile(geofilename)
        
        p = subprocess.Popen(['./gmsh', '-3', '-smooth', '5' ,'-optimize', geofilename,'-0','-2','-format','su2','-o', outputpath +fname+'.su2'])
        p.wait()
        p = subprocess.Popen(['./gmsh', '-3', '-smooth', '5' ,'-optimize', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.msh'])
        p.wait()
        p = subprocess.Popen(['./gmsh', '-3', '-smooth', '5' ,'-optimize', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.vtk'])
        p.wait()
                
        #self.meshfilename = outputpath +fname+'.su2'


    def genGmshNozzle(self, outputpath, fname):
        x0 = self.xn[0]
        L = self.xn[-1]-self.xn[0]
        x = self.xn
        r = (self.Sn/np.pi)**0.5
        Nx = self.NxSU2
        Ny = self.NySU2
        inflationRate = self.inflationRate

        geo = gmsh()

        geo.writeLine('cl__1 = 1;\n')

        #p1 = geo.addPoint(x0, 0, 0)
        #p2 = geo.addPoint(x0+L, 0, 0)
        # write wall countour
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
        geo.transfiniteLine(inflow, Ny, inflationRate-1)

        outflow = geo.addLine([upper_wall_points[-1], bottom_wall_points[-1]])
        geo.transfiniteLine(outflow, Ny, inflationRate-1)
        #l1 = geo.addLine([p1,p2])
        #geo.transfiniteLine(l1, Nx, 1)

        # inlet
        #l2 = geo.addLine([wp[0],p1])
        #geo.transfiniteLine(l2, Ny, inflationRate)

        # outlet
        #l3 = geo.addLine([wp[-1],p2])
        #geo.transfiniteLine(l3, Ny, inflationRate)

        # upper wall
        #l4 = geo.addSpline(wp)
        #geo.transfiniteLine(l4, Nx, 1)
                        
        loop1 = geo.addLoop([upper_wall,outflow,-bottom_wall,-inflow])
        surface1 = geo.addPlaneSurface([loop1])

        geo.transfiniteSurface(surface1, [upper_wall_points[0],upper_wall_points[-1],bottom_wall_points[-1],bottom_wall_points[0]])
        geo.recombineSurface(surface1)


        # add boundary
        geo.addPhysicalLine("UPPER_WALL",upper_wall)
        geo.addPhysicalLine("BOTTOM_WALL",bottom_wall)
        geo.addPhysicalLine("INFLOW",inflow)
        geo.addPhysicalLine("OUTFLOW",outflow)
        #geo.addPhysicalLine("SYMMETRY",l1)

        geo.addPhysicalSurface(surface1)
        
        #geo.gfile += 'Recombine Surface {5};\nPhysical Surface(6) = {5};\n'

        self.geo = geo

        os.makedirs(outputpath, exist_ok=True)

        fname = fname.split('.cfg')[0]
        geofilename = outputpath + fname + '.geo'

        geo.writeFile(geofilename)
        
        p = subprocess.Popen(['./gmsh', geofilename,'-0','-2','-format','su2','-o', outputpath +fname+'.su2'])
        p.wait()
        p = subprocess.Popen(['./gmsh', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.msh'])
        p.wait()
        p = subprocess.Popen(['./gmsh', geofilename,'-0','-2','-format','msh','-o', outputpath +fname+'.vtk'])
        p.wait()

        self.meshfilename = outputpath +fname+'.su2'

    def setupSU2Solver(self, itmax, itprint, CFL, tol, tscheme, fscheme, dim):
        self.itmaxSU2 = itmax
        self.itprintSU2 = itprint
        self.CFLSU2 = CFL
        self.tolSU2 = tol
        self.tschemeSU2 = tscheme
        self.fschemeSU2 = fscheme
        self.dimSU2 = dim

    def solveSU2CHT(self, solver, cores=None):
        oldcfg = self.su2cfgfile
        self.su2cfgfile = self.su2infilepath + 'cht_setupSU2.cfg'
        self.solveSU2(solver, cores)
        self.su2cfgfile = oldcfg

    def solveSU2(self, solver, cores=None):
        self.runtimeSU2 = timeit.default_timer()
        
    
        if cores is not None: 
            p = subprocess.Popen(['mpirun', '-n', f'{cores}', solver, self.su2cfgfile]) 
        else: 
            p = subprocess.Popen([solver, self.su2cfgfile])
        p.wait()
        # move results to output file
        os.makedirs(self.su2outfilepath, exist_ok=True)
        try:
            shutil.move("./restart_0.csv", self.su2outfilepath +  "restart_fluid.csv")
            shutil.move("./restart_1.csv", self.su2outfilepath +  "restart_solid.csv")
            shutil.move("./fluid_0.vtk", self.su2outfilepath +  "fluid.vtk")
            shutil.move("./solid_1.vtk", self.su2outfilepath +  "solid.vtk")
            
            saveVtkMesh(self.su2outfilepath +  "fluid.vtk")
            saveVtkMesh(self.su2outfilepath +  "solid.vtk")
        except:
            print("No solution files exists")

        self.runtimeSU2 = timeit.default_timer() - self.runtimeSU2
        print('runtime: ' + str(self.runtimeSU2) + 's')
        
    def getMesh(self, meshfile):
        mesh = meshio.read(self.su2infilepath + meshfile)
        Ny = mesh.cell_sets_dict['INWALL']['line'].size
        Nx = mesh.cell_sets_dict['OUTERWALL']['line'].size
        N = Ny*Nx
        try:
            self.mesh[meshfile]
        except:
            self.mesh[meshfile] = {}
        self.mesh[meshfile]['mesh'] = mesh
        self.mesh[meshfile]['Ny'] = Ny
        self.mesh[meshfile]['Nx'] = Nx
        self.mesh[meshfile]['N'] = N
        return self.mesh[meshfile]['mesh']

    
    def getVtk(self, vtkfile):
        solid = pv.read(self.su2outfilepath + vtkfile)
        centers = solid.cell_centers()
        sampled = solid.probe(centers)
        self.vtk[vtkfile] = solid
        return self.vtk

    def getHprop(self, meshfile, vtkfile, prop_str):
        Ny = self.mesh[meshfile]['Ny']
        Nx = self.mesh[meshfile]['Nx']
        N = self.mesh[meshfile]['N']
        vtk = self.vtk[vtkfile]
        
        prop = np.zeros((Ny, Nx))
        coord = np.zeros((Ny, Nx, 3))

        for i in range(Ny):
            p, c = getHLine(vtk, prop_str, i, N, Ny)
            prop[i,:] = np.array(p)
            coord[i,:,:] = np.array(c)

        try:
            self.hprop[vtkfile]
        except:
            self.hprop[vtkfile]={}

        self.hprop[vtkfile][prop_str] = prop
        self.hprop[vtkfile]['coordinates'] = coord
        
        return prop, coord

    def getVprop(self, meshfile, vtkfile, prop_str):
        Ny = self.mesh[meshfile]['Ny']
        Nx = self.mesh[meshfile]['Nx']
        N = self.mesh[meshfile]['N']
        vtk = self.vtk[vtkfile]
        
        prop = np.zeros((Nx, Ny))
        coord = np.zeros((Nx, Ny, 3))

        for i in range(Nx):
            p, c = getVLine(vtk, prop_str, i, N, Ny)
            prop[i,:] = np.array(p)
            coord[i,:,:] = np.array(c)

        try:
            self.vprop[vtkfile]
        except:
            self.vprop[vtkfile]={}

        self.vprop[vtkfile][prop_str] = prop
        self.vprop[vtkfile]['coordinates'] = coord
        
        return prop, coord

    def setupSU2(self, su2filepath, su2filename):
        self.su2infilepath = su2filepath
        self.su2filename = su2filename
        self.su2infilepath = su2filepath + 'inputs/'
        self.su2outfilepath = su2filepath + 'outputs/'
        self.su2cfgfile = os.path.abspath(self.su2infilepath + su2filename)        
        #self.su2cfgfile_solid = os.path.abspath(self.su2infilepath + f"solid_{su2filename}")
        #self.su2cfgfile_cht = os.path.abspath(self.su2infilepath + f"cht_{su2filename}")
        #self.multizonefile = self.su2infilepath + 'multizone.su2'
        

        # Generate the mesh files
        self.genGmshNozzle(self.su2infilepath, su2filename)

        # Setup the SU2 .cfg file
        meshfile = self.meshfilename
        solfile = 'solution.dat'

        with open(self.su2cfgfile, 'w') as f:
            newcfg = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% SU2 configuration file                                                       %
% Case description: Non-ideal compressible fluid flow in a converging-         %
%                   diverging supersonic nozzle for siloxane fluid MDM         %
% Author: Alberto Guardone                                                     %
% Institution: Politecnico di Milano                                           %
% Date: 2019.05.03                                                             %
% File Version 6.2.0 "Falcon"                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               FEM_EULER, FEM_NAVIER_STOKES, FEM_RANS, FEM_LES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)
SOLVER= RANS
%AXISYMMETRIC= YES
%
% Specify turbulence model (NONE, SA, SA_NEG, SST, SA_E, SA_COMP, SA_E_COMP)
KIND_TURB_MODEL= {self.turbulenceModel}
%
% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)
MATH_PROBLEM= DIRECT
%
% Restart solution (NO, YES)
RESTART_SOL= NO
%
% System of measurements (SI, US)
% International system of units (SI): ( meters, kilograms, Kelvins,
%                                       Newtons = kg m/s^2, Pascals = N/m^2,
%                                       Density = kg/m^3, Speed = m/s,
%                                       Equiv. Area = m^2 )
% United States customary units (US): ( inches, slug, Rankines, lbf = slug ft/s^2,
%                                       psf = lbf/ft^2, Density = slug/ft^3,
%                                       Speed = ft/s, Equiv. Area = ft^2 )
SYSTEM_MEASUREMENTS= SI

%
% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%
%
% Mach number (non-dimensional, based on the free-stream values)
MACH_NUMBER= {self.Min}
%
% Angle of attack (degrees, only for compressible flows)
AOA= 0.0
%
% Side-slip angle (degrees, only for compressible flows)
SIDESLIP_ANGLE= 0.0
%
% Init option to choose between Reynolds (default) or thermodynamics quantities
% for initializing the solution (REYNOLDS, TD_CONDITIONS)
INIT_OPTION= TD_CONDITIONS
%
% Free-stream option to choose between density and temperature (default) for
% initializing the solution (TEMPERATURE_FS, DENSITY_FS)
FREESTREAM_OPTION= TEMPERATURE_FS
%
% Free-stream pressure (101325.0 N/m^2, 2116.216 psf by default)
FREESTREAM_PRESSURE = {self.p0in}

% Free-stream temperature (288.15 K, 518.67 R by default)
FREESTREAM_TEMPERATURE= {self.T0in}
%
% Compressible flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,
%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)
REF_DIMENSIONALIZATION= {self.dimSU2}

% ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------%
%
% Fluid model (STANDARD_AIR, IDEAL_GAS, VW_GAS, PR_GAS,
%              CONSTANT_DENSITY, INC_IDEAL_GAS, INC_IDEAL_GAS_POLY)
FLUID_MODEL= {self.fluidModel}
%
% Ratio of specific heats (1.4 default and the value is hardcoded
%                          for the model STANDARD_AIR, compressible only)
GAMMA_VALUE= {self.gamma}
%
% Specific gas constant (287.058 J/kg*K default and this value is hardcoded
%                        for the model STANDARD_AIR, compressible only)
GAS_CONSTANT= {self.R}
%
% Critical Temperature (131.00 K by default)
CRITICAL_TEMPERATURE= {self.criticalTemperature}
%
% Critical Pressure (3588550.0 N/m^2 by default)
CRITICAL_PRESSURE= {self.criticalPressure}
%
% Acentric factor (0.035 (air))
ACENTRIC_FACTOR= {self.acentricFactor}

% --------------------------- VISCOSITY MODEL ---------------------------------%
%
% Viscosity model (SUTHERLAND, CONSTANT_VISCOSITY).
VISCOSITY_MODEL= {self.viscosityModel}
%
% Sutherland Viscosity Ref (1.716E-5 default value for AIR SI)
MU_REF= {self.sutherlandViscosity}
%
% Sutherland Temperature Ref (273.15 K default value for AIR SI)
MU_T_REF= {self.sutherlandTemperature}
%
% Sutherland constant (110.4 default value for AIR SI)
SUTHERLAND_CONSTANT= {self.sutherlandConstant}

% --------------------------- THERMAL CONDUCTIVITY MODEL ----------------------%
%
% Conductivity model (CONSTANT_CONDUCTIVITY, CONSTANT_PRANDTL).
CONDUCTIVITY_MODEL= {self.conductivityModel}
%
% Laminar Prandtl number (0.72 (air), only for CONSTANT_PRANDTL)
PRANDTL_LAM= {self.laminarPrandtl}
%
% Turbulent Prandtl number (0.9 (air), only for CONSTANT_PRANDTL)
PRANDTL_TURB= {self.turbulentPrandtl}

% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%
%
% Navier-Stokes (no-slip), constant heat flux wall  marker(s) (NONE = no marker)
% Format: ( marker name, constant heat flux (J/m^2), ... )
MARKER_HEATFLUX= ( WALL, 0.0 )
%
% Symmetry boundary marker(s) (NONE = no marker)
MARKER_SYM= ( SYMMETRY )
%
% Riemann boundary marker(s) (NONE = no marker)
% Format: (marker, data kind flag, list of data)
MARKER_RIEMANN= ( INFLOW, TOTAL_CONDITIONS_PT, {self.p0in}, {self.T0in}, 1.0, 0.0, 0.0 , OUTFLOW, STATIC_PRESSURE, {self.pb}, 0.0, 0.0, 0.0, 0.0 )

% MARKER_SUPERSONIC_OUTLET = ( OUTFLOW )

% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%
%
% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)
NUM_METHOD_GRAD= GREEN_GAUSS
%
% CFL number (initial value for the adaptive CFL number)
CFL_NUMBER= {self.CFLSU2}
%
% Adaptive CFL number (NO, YES)
CFL_ADAPT= YES
%
% Parameters of the adaptive CFL number (factor down, factor up, CFL min value,
%                                        CFL max value )
CFL_ADAPT_PARAM= ( 0.5, 2.0, 0.1 , 100000.0 )

%
% Maximum Delta Time in local time stepping simulations
MAX_DELTA_TIME= 1E6

% ----------- SLOPE LIMITER AND DISSIPATION SENSOR DEFINITION -----------------%
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the flow equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_FLOW= YES
%
% Slope limiter (NONE, VENKATAKRISHNAN, VENKATAKRISHNAN_WANG,
%                BARTH_JESPERSEN, VAN_ALBADA_EDGE)
SLOPE_LIMITER_FLOW= NONE
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the turbulence equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_TURB= NO

% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%
%
% Linear solver or smoother for implicit formulations (BCGSTAB, FGMRES, SMOOTHER_JACOBI,
%                                                      SMOOTHER_ILU, SMOOTHER_LUSGS,
%                                                      SMOOTHER_LINELET)
LINEAR_SOLVER= FGMRES
%
% Preconditioner of the Krylov linear solver (ILU, LU_SGS, LINELET, JACOBI)
LINEAR_SOLVER_PREC= ILU
%
% Linael solver ILU preconditioner fill-in level (0 by default)
LINEAR_SOLVER_ILU_FILL_IN= 0
%
% Minimum error of the linear solver for implicit formulations
LINEAR_SOLVER_ERROR= 1E-6
%
% Max number of iterations of the linear solver for the implicit formulation
LINEAR_SOLVER_ITER= 10

% -------------------------- MULTIGRID PARAMETERS -----------------------------%
%
% Multi-grid levels (0 = no multi-grid)
MGLEVEL= 3

% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%
%
% Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, AUSMPLUSUP, AUSMPLUSUP2, HLLC,
%                              TURKEL_PREC, MSW, FDS)
CONV_NUM_METHOD_FLOW= {self.fschemeSU2}
%
% Entropy fix coefficient (0.0 implies no entropy fixing, 1.0 implies scalar
%                          artificial dissipation)
ENTROPY_FIX_COEFF= 0.1
%
% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)
TIME_DISCRE_FLOW= {self.tschemeSU2}

% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%
%
% Convective numerical method (SCALAR_UPWIND)
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
%
% Time discretization (EULER_IMPLICIT)
TIME_DISCRE_TURB= EULER_IMPLICIT
%
% Reduction factor of the CFL coefficient in the turbulence problem
CFL_REDUCTION_TURB= 1.0

% --------------------------- CONVERGENCE PARAMETERS --------------------------%
%
% Number of total iterations
TIME_ITER= {self.itmaxSU2}
%
% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= {self.tolSU2}
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10

% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%
%
% Mesh input file
MESH_FILENAME= {self.meshfilename}
%
% Mesh input file format (SU2, CGNS)
MESH_FORMAT= SU2
%
% Mesh output file
%MESH_OUT_FILENAME= mesh_out.su2
%
% Restart flow input file
SOLUTION_FILENAME= solution_fluid.dat
%
% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW, PARAVIEW_BINARY,
%                     FIELDVIEW, FIELDVIEW_BINARY)
OUTPUT_FILES = (CSV, PARAVIEW_ASCII)
OUTPUT_WRT_FREQ= {self.itmaxSU2}
TABULAR_FORMAT= CSV
%
% Output file convergence history (w/o extension)
CONV_FILENAME= history_fluid
%
% Output file restart flow
%RESTART_FILENAME= restart_fluid.csv
%
% Output file flow (w/o extension) variables
VOLUME_FILENAME= fluid
%
% Output file surface flow coefficient (w/o extension)
SURFACE_FILENAME= fluid_surface
%
%HISTORY_OUTPUT=OUTER_ITER, RMS_DENSITY, RMS_ENERGY, RMS_MOMENTUM-X
"""
            f.write(newcfg)
        
        
    def setupSU2CHT(self, su2filepath, su2filename):
        self.su2infilepath = su2filepath
        self.su2filename = su2filename
        self.su2infilepath = su2filepath + 'inputs/'
        self.su2outfilepath = su2filepath + 'outputs/'
        self.su2cfgfile = os.path.abspath(self.su2infilepath + su2filename)        
        self.su2cfgfile_solid = os.path.abspath(self.su2infilepath + f"solid_{su2filename}")
        self.su2cfgfile_cht = os.path.abspath(self.su2infilepath + f"cht_{su2filename}")
        self.multizonefile = self.su2infilepath + 'multizone.su2'
        

        # Generate the mesh files
        self.genGmshNozzle(self.su2infilepath, su2filename)

        # Setup the SU2 .cfg file
        meshfile = self.meshfilename
        solfile = 'solution.dat'

        with open(self.su2cfgfile, 'w') as f:
            newcfg = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% SU2 configuration file                                                       %
% Case description: Non-ideal compressible fluid flow in a converging-         %
%                   diverging supersonic nozzle for siloxane fluid MDM         %
% Author: Alberto Guardone                                                     %
% Institution: Politecnico di Milano                                           %
% Date: 2019.05.03                                                             %
% File Version 6.2.0 "Falcon"                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               FEM_EULER, FEM_NAVIER_STOKES, FEM_RANS, FEM_LES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)
SOLVER= RANS
%AXISYMMETRIC= YES
%
% Specify turbulence model (NONE, SA, SA_NEG, SST, SA_E, SA_COMP, SA_E_COMP)
KIND_TURB_MODEL= {self.turbulenceModel}
%
% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)
MATH_PROBLEM= DIRECT
%
% Restart solution (NO, YES)
RESTART_SOL= NO
%
% System of measurements (SI, US)
% International system of units (SI): ( meters, kilograms, Kelvins,
%                                       Newtons = kg m/s^2, Pascals = N/m^2,
%                                       Density = kg/m^3, Speed = m/s,
%                                       Equiv. Area = m^2 )
% United States customary units (US): ( inches, slug, Rankines, lbf = slug ft/s^2,
%                                       psf = lbf/ft^2, Density = slug/ft^3,
%                                       Speed = ft/s, Equiv. Area = ft^2 )
SYSTEM_MEASUREMENTS= SI

%
% -------------------- COMPRESSIBLE FREE-STREAM DEFINITION --------------------%
%
% Mach number (non-dimensional, based on the free-stream values)
MACH_NUMBER= {self.Min}
%
% Angle of attack (degrees, only for compressible flows)
AOA= 0.0
%
% Side-slip angle (degrees, only for compressible flows)
SIDESLIP_ANGLE= 0.0
%
% Init option to choose between Reynolds (default) or thermodynamics quantities
% for initializing the solution (REYNOLDS, TD_CONDITIONS)
INIT_OPTION= TD_CONDITIONS
%
% Free-stream option to choose between density and temperature (default) for
% initializing the solution (TEMPERATURE_FS, DENSITY_FS)
FREESTREAM_OPTION= TEMPERATURE_FS
%
% Free-stream pressure (101325.0 N/m^2, 2116.216 psf by default)
FREESTREAM_PRESSURE = {self.p0in}

% Free-stream temperature (288.15 K, 518.67 R by default)
FREESTREAM_TEMPERATURE= {self.T0in}
%
% Compressible flow non-dimensionalization (DIMENSIONAL, FREESTREAM_PRESS_EQ_ONE,
%                              FREESTREAM_VEL_EQ_MACH, FREESTREAM_VEL_EQ_ONE)
REF_DIMENSIONALIZATION= {self.dimSU2}

% ---- IDEAL GAS, POLYTROPIC, VAN DER WAALS AND PENG ROBINSON CONSTANTS -------%
%
% Fluid model (STANDARD_AIR, IDEAL_GAS, VW_GAS, PR_GAS,
%              CONSTANT_DENSITY, INC_IDEAL_GAS, INC_IDEAL_GAS_POLY)
FLUID_MODEL= {self.fluidModel}
%
% Ratio of specific heats (1.4 default and the value is hardcoded
%                          for the model STANDARD_AIR, compressible only)
GAMMA_VALUE= {self.gamma}
%
% Specific gas constant (287.058 J/kg*K default and this value is hardcoded
%                        for the model STANDARD_AIR, compressible only)
GAS_CONSTANT= {self.R}
%
% Critical Temperature (131.00 K by default)
CRITICAL_TEMPERATURE= {self.criticalTemperature}
%
% Critical Pressure (3588550.0 N/m^2 by default)
CRITICAL_PRESSURE= {self.criticalPressure}
%
% Acentric factor (0.035 (air))
ACENTRIC_FACTOR= {self.acentricFactor}

% --------------------------- VISCOSITY MODEL ---------------------------------%
%
% Viscosity model (SUTHERLAND, CONSTANT_VISCOSITY).
VISCOSITY_MODEL= {self.viscosityModel}
%
% Sutherland Viscosity Ref (1.716E-5 default value for AIR SI)
MU_REF= {self.sutherlandViscosity}
%
% Sutherland Temperature Ref (273.15 K default value for AIR SI)
MU_T_REF= {self.sutherlandTemperature}
%
% Sutherland constant (110.4 default value for AIR SI)
SUTHERLAND_CONSTANT= {self.sutherlandConstant}

% --------------------------- THERMAL CONDUCTIVITY MODEL ----------------------%
%
% Conductivity model (CONSTANT_CONDUCTIVITY, CONSTANT_PRANDTL).
CONDUCTIVITY_MODEL= {self.conductivityModel}
%
% Laminar Prandtl number (0.72 (air), only for CONSTANT_PRANDTL)
PRANDTL_LAM= {self.laminarPrandtl}
%
% Turbulent Prandtl number (0.9 (air), only for CONSTANT_PRANDTL)
PRANDTL_TURB= {self.turbulentPrandtl}

% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%
%
% Navier-Stokes (no-slip), constant heat flux wall  marker(s) (NONE = no marker)
% Format: ( marker name, constant heat flux (J/m^2), ... )
% MARKER_HEATFLUX= ( WALL, 0.0 )
%
% Symmetry boundary marker(s) (NONE = no marker)
% MARKER_SYM= ( SYMMETRY )
%
% Riemann boundary marker(s) (NONE = no marker)
% Format: (marker, data kind flag, list of data)
MARKER_RIEMANN= ( INFLOW, TOTAL_CONDITIONS_PT, {self.p0in}, {self.T0in}, 1.0, 0.0, 0.0 , OUTFLOW, STATIC_PRESSURE, {self.pb}, 0.0, 0.0, 0.0, 0.0 )

%MARKER_SUPERSONIC_OUTLET = ( OUTFLOW )

% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%
%
% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)
NUM_METHOD_GRAD= GREEN_GAUSS
%
% CFL number (initial value for the adaptive CFL number)
CFL_NUMBER= {self.CFLSU2}
%
% Adaptive CFL number (NO, YES)
CFL_ADAPT= YES
%
% Parameters of the adaptive CFL number (factor down, factor up, CFL min value,
%                                        CFL max value )
CFL_ADAPT_PARAM= ( 0.5, 2.0, 0.1 , 100000.0 )

%
% Maximum Delta Time in local time stepping simulations
MAX_DELTA_TIME= 1E6

% ----------- SLOPE LIMITER AND DISSIPATION SENSOR DEFINITION -----------------%
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the flow equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_FLOW= YES
%
% Slope limiter (NONE, VENKATAKRISHNAN, VENKATAKRISHNAN_WANG,
%                BARTH_JESPERSEN, VAN_ALBADA_EDGE)
SLOPE_LIMITER_FLOW= NONE
%
% Monotonic Upwind Scheme for Conservation Laws (TVD) in the turbulence equations.
%           Required for 2nd order upwind schemes (NO, YES)
MUSCL_TURB= NO

% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%
%
% Linear solver or smoother for implicit formulations (BCGSTAB, FGMRES, SMOOTHER_JACOBI,
%                                                      SMOOTHER_ILU, SMOOTHER_LUSGS,
%                                                      SMOOTHER_LINELET)
LINEAR_SOLVER= FGMRES
%
% Preconditioner of the Krylov linear solver (ILU, LU_SGS, LINELET, JACOBI)
LINEAR_SOLVER_PREC= ILU
%
% Linael solver ILU preconditioner fill-in level (0 by default)
LINEAR_SOLVER_ILU_FILL_IN= 0
%
% Minimum error of the linear solver for implicit formulations
LINEAR_SOLVER_ERROR= 1E-6
%
% Max number of iterations of the linear solver for the implicit formulation
LINEAR_SOLVER_ITER= 10

% -------------------------- MULTIGRID PARAMETERS -----------------------------%
%
% Multi-grid levels (0 = no multi-grid)
MGLEVEL= 3

% -------------------- FLOW NUMERICAL METHOD DEFINITION -----------------------%
%
% Convective numerical method (JST, LAX-FRIEDRICH, CUSP, ROE, AUSM, AUSMPLUSUP, AUSMPLUSUP2, HLLC,
%                              TURKEL_PREC, MSW, FDS)
CONV_NUM_METHOD_FLOW= {self.fschemeSU2}
%
% Entropy fix coefficient (0.0 implies no entropy fixing, 1.0 implies scalar
%                          artificial dissipation)
ENTROPY_FIX_COEFF= 0.1
%
% Time discretization (RUNGE-KUTTA_EXPLICIT, EULER_IMPLICIT, EULER_EXPLICIT)
TIME_DISCRE_FLOW= {self.tschemeSU2}

% -------------------- TURBULENT NUMERICAL METHOD DEFINITION ------------------%
%
% Convective numerical method (SCALAR_UPWIND)
CONV_NUM_METHOD_TURB= SCALAR_UPWIND
%
% Time discretization (EULER_IMPLICIT)
TIME_DISCRE_TURB= EULER_IMPLICIT
%
% Reduction factor of the CFL coefficient in the turbulence problem
CFL_REDUCTION_TURB= 1.0

% --------------------------- CONVERGENCE PARAMETERS --------------------------%
%
% Number of total iterations
TIME_ITER= {self.itmaxSU2}
%
% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= {self.tolSU2}
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10

% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%
%
% Mesh input file
%MESH_FILENAME= {self.meshfilename}
%
% Mesh input file format (SU2, CGNS)
%MESH_FORMAT= SU2
%
% Mesh output file
%MESH_OUT_FILENAME= mesh_out.su2
%
% Restart flow input file
SOLUTION_FILENAME= solution_fluid.dat
%
% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW, PARAVIEW_BINARY,
%                     FIELDVIEW, FIELDVIEW_BINARY)
OUTPUT_FILES = (CSV, PARAVIEW_ASCII)

TABULAR_FORMAT= CSV
%
% Output file convergence history (w/o extension)
CONV_FILENAME= history_fluid
%
% Output file restart flow
%RESTART_FILENAME= restart_fluid.csv
%
% Output file flow (w/o extension) variables
VOLUME_FILENAME= fluid
%
% Output file surface flow coefficient (w/o extension)
%SURFACE_FILENAME= surface_flow
%
%HISTORY_OUTPUT=OUTER_ITER, RMS_DENSITY, RMS_ENERGY, RMS_MOMENTUM-X
"""
            f.write(newcfg)
        
        with open(self.su2cfgfile_solid, 'w') as f:
            solidcfg = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% SU2 configuration file                                                       %
% Case description: Steady incompressible laminar flow around heated cylinders %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)                           
SOLVER= HEAT_EQUATION
%AXISYMMETRIC= YES
%
% Restart solution (NO, YES)
RESTART_SOL= NO

% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%
%
MARKER_ISOTHERMAL= ( OUTERWALL, {self.outerTemperature})
%, INWALL, {self.outerTemperature}, OUTWALL, {self.outerTemperature})
MARKER_HEATFLUX= ( INWALL, 0.0, OUTWALL, 0.0)
%
% Marker(s) of the surface to be plotted or designed
MARKER_PLOTTING= ( INNERWALL )
%
% Marker(s) of the surface where the functional (Cd, Cl, etc.) will be evaluated
MARKER_MONITORING= ( INNERWALL )

% ---------------- (SOLIDS) CONDUCTION CONDITION DEFINITION -------------------%
%
% We should keep the dimensionalization of the coupled flow solver
INC_NONDIM= DIMENSIONAL
%
% Temperature initialization value
FREESTREAM_TEMPERATURE= {self.outerTemperature}
%
% Nettis case: hollow cylinder (air w/ 4x the conductivity)
%
% Solid density (kg/m^3)
MATERIAL_DENSITY= {self.solidDensity}

% Solid specific heat (J/kg*K)
SPECIFIC_HEAT_CP= {self.solidHeatCP}
%
% Solid thermal conductivity (W/m*K)
THERMAL_CONDUCTIVITY_CONSTANT= {self.solidConductivity}

% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%
%
% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)
NUM_METHOD_GRAD= GREEN_GAUSS
%
% Courant-Friedrichs-Lewy condition of the finest grid
%CFL_NUMBER= 10.0
CFL_NUMBER= {self.CFLSU2}
%
% Adaptive CFL number (NO, YES)
CFL_ADAPT= NO
%
% Runge-Kutta alpha coefficients
RK_ALPHA_COEFF= ( 0.66667, 0.66667, 1.000000 )

% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%
%
% Linear solver or smoother for implicit formulations (BCGSTAB, FGMRES, SMOOTHER_JACOBI,
%                                                      SMOOTHER_ILU, SMOOTHER_LUSGS,
%                                                      SMOOTHER_LINELET)
LINEAR_SOLVER= FGMRES
%
% Preconditioner of the Krylov linear solver (ILU, LU_SGS, LINELET, JACOBI)
LINEAR_SOLVER_PREC= ILU
%
% Linear solver ILU preconditioner fill-in level (0 by default)
LINEAR_SOLVER_ILU_FILL_IN= 0
%
% Minimum error of the linear solver for implicit formulations
LINEAR_SOLVER_ERROR= 1E-15
%
% Max number of iterations of the linear solver for the implicit formulation
LINEAR_SOLVER_ITER= 5

% -------------------- HEAT NUMERICAL METHOD DEFINITION -----------------------%
%
TIME_DISCRE_HEAT= EULER_IMPLICIT

% --------------------------- CONVERGENCE PARAMETERS --------------------------%
%
% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= -19
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10
%
% Number of elements to apply the criteria
CONV_CAUCHY_ELEMS= 100
%
% Epsilon to control the series convergence
CONV_CAUCHY_EPS= 1E-6

% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%
%
% Restart flow input file
SOLUTION_FILENAME= solution_solid.dat
%
% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW, PARAVIEW_BINARY,
%                     FIELDVIEW, FIELDVIEW_BINARY)
OUTPUT_FILES = (CSV, PARAVIEW_ASCII)
%
TABULAR_FORMAT= CSV
%
% Output file convergence history (w/o extension)
CONV_FILENAME= history_solid
%
% Output file with the forces breakdown
%BREAKDOWN_FILENAME= forces_breakdown.dat
%
% Output file restart flow
%RESTART_FILENAME= restart_solid.dat
%
%
% Output file flow (w/o extension) variables
VOLUME_FILENAME= solid
%
% Output file surface flow coefficient (w/o extension)
% SURFACE_FILENAME= surface_flow
%
% History output
%HISTORY_OUTPUT=OUTER_ITER, RMS_TEMPERATURE

% Screen output
%SCREEN_OUTPUT=OUTER_ITER, RMS_TEMPERATURE
"""
            f.write(solidcfg)


        with open(self.su2cfgfile_solid.replace('.cfg','_2.cfg'), 'w') as f:
            solidcfg = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% SU2 configuration file                                                       %
% Case description: Steady incompressible laminar flow around heated cylinders %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- DIRECT, ADJOINT, AND LINEARIZED PROBLEM DEFINITION ------------%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)                           
SOLVER= HEAT_EQUATION
%AXISYMMETRIC= YES
%
% Restart solution (NO, YES)
RESTART_SOL= NO

% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%
%
MARKER_ISOTHERMAL= ( OUTERWALL_2, {self.outerTemperature})
%, INWALL, {self.outerTemperature}, OUTWALL, {self.outerTemperature})
MARKER_HEATFLUX= ( INWALL_2, 0.0, OUTWALL_2, 0.0)
%
% Marker(s) of the surface to be plotted or designed
MARKER_PLOTTING= ( INNERWALL_2 )
%
% Marker(s) of the surface where the functional (Cd, Cl, etc.) will be evaluated
MARKER_MONITORING= ( INNERWALL_2 )

% ---------------- (SOLIDS) CONDUCTION CONDITION DEFINITION -------------------%
%
% We should keep the dimensionalization of the coupled flow solver
INC_NONDIM= DIMENSIONAL
%
% Temperature initialization value
FREESTREAM_TEMPERATURE= {self.outerTemperature}
%
% Nettis case: hollow cylinder (air w/ 4x the conductivity)
%
% Solid density (kg/m^3)
MATERIAL_DENSITY= {self.solidDensity}

% Solid specific heat (J/kg*K)
SPECIFIC_HEAT_CP= {self.solidHeatCP}
%
% Solid thermal conductivity (W/m*K)
THERMAL_CONDUCTIVITY_CONSTANT= {self.solidConductivity}

% ------------- COMMON PARAMETERS DEFINING THE NUMERICAL METHOD ---------------%
%
% Numerical method for spatial gradients (GREEN_GAUSS, WEIGHTED_LEAST_SQUARES)
NUM_METHOD_GRAD= GREEN_GAUSS
%
% Courant-Friedrichs-Lewy condition of the finest grid
%CFL_NUMBER= 10.0
CFL_NUMBER= {self.CFLSU2}
%
% Adaptive CFL number (NO, YES)
CFL_ADAPT= NO
%
% Runge-Kutta alpha coefficients
RK_ALPHA_COEFF= ( 0.66667, 0.66667, 1.000000 )

% ------------------------ LINEAR SOLVER DEFINITION ---------------------------%
%
% Linear solver or smoother for implicit formulations (BCGSTAB, FGMRES, SMOOTHER_JACOBI,
%                                                      SMOOTHER_ILU, SMOOTHER_LUSGS,
%                                                      SMOOTHER_LINELET)
LINEAR_SOLVER= FGMRES
%
% Preconditioner of the Krylov linear solver (ILU, LU_SGS, LINELET, JACOBI)
LINEAR_SOLVER_PREC= ILU
%
% Linear solver ILU preconditioner fill-in level (0 by default)
LINEAR_SOLVER_ILU_FILL_IN= 0
%
% Minimum error of the linear solver for implicit formulations
LINEAR_SOLVER_ERROR= 1E-15
%
% Max number of iterations of the linear solver for the implicit formulation
LINEAR_SOLVER_ITER= 5

% -------------------- HEAT NUMERICAL METHOD DEFINITION -----------------------%
%
TIME_DISCRE_HEAT= EULER_IMPLICIT

% --------------------------- CONVERGENCE PARAMETERS --------------------------%
%
% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= -19
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10
%
% Number of elements to apply the criteria
CONV_CAUCHY_ELEMS= 100
%
% Epsilon to control the series convergence
CONV_CAUCHY_EPS= 1E-6

% ------------------------- INPUT/OUTPUT INFORMATION --------------------------%
%
% Restart flow input file
SOLUTION_FILENAME= solution_solid.dat
%
% Output file format (TECPLOT, TECPLOT_BINARY, PARAVIEW, PARAVIEW_BINARY,
%                     FIELDVIEW, FIELDVIEW_BINARY)
OUTPUT_FILES = (CSV, PARAVIEW_ASCII)
%
TABULAR_FORMAT= CSV
%
% Output file convergence history (w/o extension)
CONV_FILENAME= history_solid
%
% Output file with the forces breakdown
%BREAKDOWN_FILENAME= forces_breakdown.dat
%
% Output file restart flow
%RESTART_FILENAME= restart_solid.dat
%
%
% Output file flow (w/o extension) variables
VOLUME_FILENAME= solid
%
% Output file surface flow coefficient (w/o extension)
% SURFACE_FILENAME= surface_flow
%
% History output
%HISTORY_OUTPUT=OUTER_ITER, RMS_TEMPERATURE

% Screen output
%SCREEN_OUTPUT=OUTER_ITER, RMS_TEMPERATURE
"""
            f.write(solidcfg)

        with open(self.su2cfgfile_cht, 'w') as f:

            chtcfg = f"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                               %
% SU2 configuration file                                                        %
% Case description: 2D cylinder array with CHT couplings                        %
% Author: O. Burghardt, T. Economon                                             %
% Institution: Chair for Scientific Computing, TU Kaiserslautern                %
% Date: August 8, 2019                                                          %
% File Version 7.1.1 "Blackbird"                                                %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Physical governing equations (EULER, NAVIER_STOKES,
%                               WAVE_EQUATION, HEAT_EQUATION, FEM_ELASTICITY,
%                               POISSON_EQUATION)             
SOLVER= MULTIPHYSICS
%AXISYMMETRIC= YES
%
% Mathematical problem (DIRECT, CONTINUOUS_ADJOINT, DISCRETE_ADJOINT)
MATH_PROBLEM= DIRECT
%
%

CONFIG_LIST = ({self.su2cfgfile}, {self.su2cfgfile_solid},  {self.su2cfgfile_solid.replace('.cfg','_2.cfg')} )
%
%
MARKER_ZONE_INTERFACE= (UPPER_WALL, INNERWALL, BOTTOM_WALL, INNERwALL_2)
%
%
MARKER_CHT_INTERFACE= (UPPER_WALL, INNERWALL, BOTTOM_WALL, INNERWALL_2)
%
%
%MARKER_SYM= ( SYMMETRY )
%
%
CHT_COUPLING_METHOD= DIRECT_TEMPERATURE_ROBIN_HEATFLUX
%
%
TIME_DOMAIN = NO
%
% Number of total iterations (15000 for suitable results)
%OUTER_ITER = 10000
OUTER_ITER = {self.itmaxSU2}


% Min value of the residual (log10 of the residual)
CONV_RESIDUAL_MINVAL= {self.tolSU2}
%
% Start convergence criteria at iteration number
CONV_STARTITER= 10
%
CONV_FILENAME= history_multizone
%
% Mesh input file
MESH_FILENAME= {self.multizonefile}

% Mesh input file format (SU2, CGNS, NETCDF_ASCII)
MESH_FORMAT= SU2
%
% Output file format
OUTPUT_FILES= (CSV, PARAVIEW_ASCII)

TABULAR_FORMAT=CSV

%CONV_FILENAME= multizone

% Data written to history file
HISTORY_OUTPUT=OUTER_ITER, RMS_DENSITY[0], RMS_ENERGY[0], RMS_MOMENTUM-X[0], RMS_TEMPERATURE[1]
% Writing frequency for history output
HISTORY_WRT_FREQ_INNER= 1
%
HISTORY_WRT_FREQ_OUTER= 1
% 
HISTORY_WRT_FREQ_TIME= 1


% Screen output
SCREEN_OUTPUT=OUTER_ITER, RMS_DENSITY[0], RMS_ENERGY[0], RMS_MOMENTUM-X[0], RMS_TEMPERATURE[1]

% Writing solution file frequency
OUTPUT_WRT_FREQ= {self.itmaxSU2}

% Writing frequency for screen output
SCREEN_WRT_FREQ_OUTER= {self.itprintSU2}

%CONV_FIELD = RMS_DENSITY[0], RMS_ENERGY[0], RMS_MOMENTUM-X[0], RMS_TEMPERATURE[1]
CONV_FIELD = RMS_ENERGY[0], RMS_ENERGY[0], RMS_TEMPERATURE[1]

"""
            f.write(chtcfg)

        self.genGmshWall(self.su2infilepath, 'solid_setupSU2.cfg', self.thickness)
        self.genGmshWall_2(self.su2infilepath, 'solid_setupSU2_2.cfg', self.thickness)
        self.genMultizone(self.su2infilepath, ['setupSU2','solid_setupSU2', 'solid_setupSU2_2'])
        

        