import numpy as np

def ft_to_meter(ft):
    """ Converts ft to meter

    :param ft: ft
    :type ft: float
    :return: meter
    :rtype: float
    """
    return ft * 0.3048

def inch_to_meter(inch):
    """ Converts inch to meter

    :param inch: inch
    :type inch: float
    :return: meter
    :rtype: float
    """
    return inch * 0.0254

def inch2_to_meter2(inch2):
    """ Converts inch^2 to meter^2

    :param inch2: inch^2
    :type inch2: float
    :return: meter^2
    :rtype: float
    """
    return inch2 * 0.0254**2

def meter_to_inch(meter):
    """ Converts meter to inch

    :param meter: meter
    :type meter: float
    :return: inch
    :rtype: float
    """
    return meter / 0.0254

def psi_to_pa(psi):
    """ Converts psi to Pa

    :param psi: psi
    :type psi: float
    :return: Pa
    :rtype: float
    """
    return psi * 6894.76

def r_to_k(r):
    """ Converts R to K

    :param r: R
    :type r: float
    :return: K
    :rtype: float
    """
    return r / 1.8

def area_to_radius(area):
    """ Converts area to radius
    
    :param area: area
    :type area: float
    :return: radius
    :rtype: float
    """
    return np.sqrt(area / np.pi)

def radius_to_area(radius):
    """ Converts radius to area

    :param radius: radius
    :type radius: float
    :return: area
    :rtype: float
    """

    return np.pi * radius**2

def parse_su2_config(template, config):
   with open(template, 'r') as f:
      fstring = f.read()
   return eval(fstring, config)

def pa_to_psi(pa):
    """ Converts Pa to psi

    :param pa: Pa
    :type pa: float
    :return: psi
    :rtype: float
    """
    return pa / 6894.76
