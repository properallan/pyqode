import timeit
import subprocess
import os
from pathlib import Path
import numpy as np
from scipy import interpolate
from typing import Union
from pyqode.utils import parse_su2_config
from numpy import interp

PYQODE_SOLVER = Path(__file__).parent / 'src' / 'eulerQ1D'
SU2_SOLVER = Path(__file__).parent / 'src' / 'SU2'

class PostProcessor:
    def __init__(self, path: str|Path):
        for file_ in Path(path).glob('*.txt'):
            print(file_.stem)
            self.__setattr__(file_.stem, np.loadtxt(file_))

class Solver:
    def __init__(self, 
            config: dict,
            executable : Union[str,Path] = PYQODE_SOLVER):
        self.config = config
        self.executable = executable

        self.setup(self.config)

    def interpolate_domain(self, config):
        x = config['domain_x']
        area = config['domain_area']
        size = config['domain_size']

        new_x = np.linspace(x[0], x[-1], size)
        try:
            spl = interpolate.InterpolatedUnivariateSpline(x, area)
            new_s = spl(new_x)
        except:
            new_s = interp(new_x, x, area)

        return new_x, new_s


    def setup(self, config):
        working_dir = Path(config['working_dir']).resolve()
        setup_file = working_dir / 'setup.txt'
        input_path = setup_file.parent / 'inputs'
        os.makedirs(input_path, exist_ok=True)
        output_path = setup_file.parent / 'outputs'
        os.makedirs(output_path, exist_ok=True)
        domain_x_file = input_path / 'xn.txt'
        domain_s_file = input_path / 'sn.txt'

        x, area = self.interpolate_domain(config)

        np.savetxt(domain_x_file, x)
        np.savetxt(domain_s_file, area)

        with open(setup_file, 'w') as f:
            f.write('# Domain x-coordinates at cell faces (no need for ghost cells)\n')
            f.write(f"{domain_x_file}\n")
            f.write('# Area distributuin at cell faces (no need for ghost cells)\n')
            f.write(f"{domain_s_file}\n")
            f.write('# Inlet Total Pressure [Pa]\n')
            f.write(f"{config['bc_p0']}\n")
            f.write('# Inlet Total Temperature [K]\n')
            f.write(f"{config['bc_T0']}\n")
            f.write('# Inlet Mach Number\n')
            f.write(f"{config['bc_M']}\n")
            f.write('# Outlet Static Pressure [Pa]\n')
            f.write(f"{config['bc_pb']}\n")
            f.write('# Gas constant [J/kg/K]\n')
            f.write(f"{config['fluid_R']}\n")
            f.write('# Specific heat ratio\n')
            f.write(f"{config['fluid_gamma']}\n")
            f.write('# Maximum number of iterations \n')
            f.write(f"{config['solver_itmax']}\n")
            f.write('# Interval to print iterations \n')
            f.write(f"{config['solver_itprint']}\n")
            f.write('# CFL number \n')
            f.write(f"{config['solver_CFL']}\n")
            f.write('# Convergence criteria \n')
            f.write(f"{config['solver_tol']}\n")
            f.write('# Time integration scheme (Euler, RK4) \n')
            f.write(f"{config['solver_tscheme']}\n")
            f.write('# Flux scheme (Basic, Roe, Split, SplitTVD, VanLeer, LaxFriedrichs, StegerWarming, AUSM) \n')
            f.write(f"{config['solver_fscheme']}\n")
            f.write('# Timestep type (Global, Local) \n')
            f.write(f"{config['solver_dttype']}\n")
            f.write('# Form to be solved (Dimensional, Dimensionless) \n')
            f.write(f"{config['solver_dim']}\n")
            f.write('# Set output directory \n')
            f.write(f"{output_path}/\n")
        self.setup_file = setup_file

    def run(self, setupfile: str|Path=None):
        if setupfile is None:
            setupfile = self.setup_file

        self.runtime = timeit.default_timer()
        p = subprocess.Popen([self.executable, setupfile])
        p.wait()
        self.runtime = timeit.default_timer() - self.runtime
        print('runtime: ' + str(self.runtime) + 's')


def gen_su2_setup(template, config, output_file): 
    template = Path(template)
    os.makedirs(Path(output_file).parent, exist_ok=True)

    with open(output_file, 'w') as f:
        f.write(parse_su2_config(template,config))

    return Path(output_file).resolve().__str__()

def gen_multizone_su2_setup(template, config, output_file, config_list):
    config['config_list'] = config_list
    output_file = gen_su2_setup(template, config, output_file)
    return Path(output_file).resolve().__str__()

class SU2Solver:
    def __init__(self, 
        config_file: Union[str, Path],
        executable : Union[str,Path] = SU2_SOLVER):
        self.config_file = Path(config_file).resolve()
        self.executable = Path(executable).resolve()

    def run(self, config_file: Union[str,Path]=None, cwd: Union[str,Path]=None):
        if config_file is None:
            config_file = self.config_file
        else:
            config_file = Path(config_file).resolve()
        if cwd is None:
            cwd = Path(self.config_file).parent

        self.runtime = timeit.default_timer()
        p = subprocess.Popen([self.executable, config_file.name], cwd=cwd)
        p.wait()
        self.runtime = timeit.default_timer() - self.runtime
        print('runtime: ' + str(self.runtime) + 's')