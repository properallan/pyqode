import numpy as np
import os
import timeit
import subprocess

class eulerQ1D:
    def __init__(self):
        self.property_names = { 'SOUND_SPEED' : 'c.txt',
                                'ENERGY' : 'e.txt',
                                'MACH' : 'M.txt',
                                'PRESSURE' : 'p.txt',
                                'RESIDUE_ENERGY' : 'rese.txt',
                                'RESIDUE_MASS' : 'resrho.txt',
                                'RESIDUE_MOMENTUM' : 'resrhou.txt',
                                'DENSITY' : 'rho.txt',
                                'AREA' : 'S.txt',
                                'TEMPERATURE' : 'T.txt',
                                'VELOCITY' : 'u.txt',
                                'X_COORDINATE' : 'x.txt'}

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

        self._get_results()
    
    def _get_results(self):
        for property, file_name in self.property_names.items():
            setattr(self, property, np.loadtxt(self.q1doutfilepath + file_name))

def sod_shock():
    Nx = 100
    x = np.linspace(0, 1, Nx)
    A = np.ones_like(x)*1.0

    sod = eulerQ1D()

    sod.setX(x)
    sod.setS(A)

    
    R = 287.0
    gamma = 1.4
    rho_in = 1.0
    p0in = 1.0
    T0in = p0in/rho_in/R
    
    pe = 0.1
    Min = 0.0
    sod.setBC(p0in, T0in, Min, pe)

    sod.R = R
    sod.gamma = gamma

    ## run quasi-1D model
    itmaxQ1D = 50000
    NxQ1D = 100
    itprintQ1D = 1000
    CFLQ1D = 0.15
    tolQ1D = 1e-6
    tschemeQ1D = 'RK4'
    fschemeQ1D = 'AUSM'
    dttypeQ1D = 'Global'
    dimQ1D = 'Dimensionless'

    sod.NxQ1D = 100
    sod.setupQ1DSolver(itmaxQ1D, itprintQ1D, CFLQ1D, tolQ1D, tschemeQ1D, fschemeQ1D, dttypeQ1D, dimQ1D)
    rootfile = './sod_shock/'
    sod.setupQ1D(rootfile,'setupQ1D.txt')
    sod.solveQ1D('./eulerQ1D')

    with open(sod.q1doutfilepath+'solver.log','w') as f:
        f.write(f"runtime={sod.runtimeQ1D}\n")
        f.write(f"p0in={sod.p0in}\n")
        f.write(f"T0in={sod.T0in}\n")

    sod.runtimeQ1D

    return sod

def sample_run():
    xr = np.loadtxt('/home/ppiper/Dropbox/local/thermal_nozzle/data/nozzle_air/nozzle.txt', delimiter=',')
    x = xr[:,0]
    r = xr[:,1]
    A = np.pi*r**2

    n = eulerQ1D()
    n.setX(x)
    n.setS(A)
    
    # Boundary conditions
    T0in = 842.78
    p0in = 11.91e5
    pr = 101e3/p0in
    Min = 0.01
    n.setBC(p0in, T0in, Min, pr*p0in)

    # Fluid properties
    fluid = 'Air'
    gamma = 1.4
    R = 287.0
    n.R = R
    n.gamma = gamma

    ## run quasi-1D model
    itmaxQ1D = 50000
    NxQ1D = 100
    itprintQ1D = 1000
    CFLQ1D = 0.15
    tolQ1D = 1e-6
    tschemeQ1D = 'RK4'
    fschemeQ1D = 'AUSM'
    dttypeQ1D = 'Global'
    dimQ1D = 'Dimensionless'
    
    n.NxQ1D = 100
    n.setupQ1DSolver(itmaxQ1D, itprintQ1D, CFLQ1D, tolQ1D, tschemeQ1D, fschemeQ1D, dttypeQ1D, dimQ1D)
    rootfile = './sample_run/'
    n.setupQ1D(rootfile,'setupQ1D.txt')
    n.solveQ1D('./eulerQ1D')

    with open(n.q1doutfilepath+'solver.log','w') as f:
        f.write(f"runtime={n.runtimeQ1D}\n")
        f.write(f"p0in={n.p0in}\n")
        f.write(f"T0in={n.T0in}\n")

    n.runtimeQ1D

    return n

if __name__== "__main__":
    sample_run()