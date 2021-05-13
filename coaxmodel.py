#!/usr/bin/python3

import numpy as np
from math import atanh
from rsclib.capacitance import c

ft_per_m = 0.3048

""" Transmission Line Properties of Cables
    Either from manufacturer data or from measured parameters.
    Inspired by two articles by Frank Witt [3,4] in an old ARRL Antenna
    Compendium [5]. Unfortunately Witt uses a now obsolete software and
    doesn't publish all formulas in the paper. Most, however, are in
    Chipman [1] and Walter Maxwell [2] has a Basic implementation from
    which some of the formulas can be reverse engineered (the Basic
    dialect used doesn't support complex numbers so everything there
    looks very bare-footed from todays perspective).
    [1] Robert A. Chipman. Theory and Problems of Transmission Lines.
        Schaums Outline. McGraw-Hill, 1968. Section 7.5 "Determination
        of transmission line characteristics from impedance
        measurements", p.  134.
    [2] M. Walter Maxwell. Reflections III, Transmission Lines
        and Antennas. CQ Communications, Inc., Hicksville, NY, third
        edition, 2010. In particular chapter 15 introducing computer
        programs (in Basic), but this becomes much easier and is
        closer to Chipman with direct use of complex numbers. I'm
        using Maxwells results to test this implementation.
    [3] Frank Witt. Transmission line properties from manufacturer’s
        data. In Straw [5], pages 179–183.
    [4] Frank Witt. Transmission line properties from measured data.
        In Straw [5], pages 184–188.
    [5] R. Dean Straw, editor. The ARRL Antenna Compendium, volume 6.
        American Radio Relay League (ARRL), 1999.
"""

class Manufacturer_Data_Cable :
    """ From manufacturer data we can interpolate cable parameters. [3]
        Z0: Characteristic impedance of the cable
        vf: velocity factor of the cable
        Cpl: Capacity per length (in F/m, Farad per meter)
        The product of Z0, vf and Cpl is a constant according to Witt [3].
        Consulting Chipman [1] what he calls the "high frequency
        approximation" we have (p. 49) (Note that Chipman uses C for what
        we called Cpl above and L is the inductivity per length):
        (1) v_phf = 1 / sqrt (LC)
        (2) Z0 = sqrt (L/C)
        solving (2) for L and replacing in (1) we get:
        (3) Z0 = 1 / (v_phf C)
        v_phf is the velocity of propagation on the line under the
        high-frequency approximation. So this is
        v_phf = c * vf
        So replacing this in (3) we get
        (4) Z0 = 1 / (vf c C)
        And Witts constant is obtained as 1/c
        expressing 1/c in µs / ft correctly resolves to Witts formula.
        Like Witt [3] we correct vf if all three parameters are given.

    >>> print ("%.3f" % (1.0 / c * 1e12 * ft_per_m))
    1016.703
    >>> Cpl = 30.8e-12 / ft_per_m
    >>> cable = Manufacturer_Data_Cable (50, .66)
    >>> print ("%.1f pF/ft" % (cable.Cpl * 1e12 * ft_per_m))
    30.8 pF/ft
    >>> cable = Manufacturer_Data_Cable (50, .66, Cpl)
    >>> print ("vf = %.6f" % cable.vf)
    vf = 0.660197

    """

    def __init__ (self, Z0, vf, Cpl = None) :
        self.Z0 = Z0
        self.vf = vf
        if Cpl is None :
            self.Cpl = 1.0 / (c * Z0 * vf)
        else :
            self.Cpl = Cpl
            self.vf  = 1.0 / (c * Z0 * Cpl)
    # end def __init__

    def fit (self, loss_data) :
        """ Gets a list of frequency/loss pairs
            Note that the loss is in dB per 100m (not ft)
        """
        pass
    # end def fit

# end class Manufacturer_Data_Cable

class Measured_Cable :
    """ By measuring cable impedance (both R and X) for open and closed
        circuit we can compute other cable paramters. The measurements
        should be made on a line with length close to lambda/8 or an odd
        multiple of it. [4,2]
        Doctests are from Maxwell [2] p.15-3

    >>> mc = Measured_Cable (1.9e6, 1, 2.0, 0.6, 50.0, -48.5)
    >>> print ("%.3f" % abs (mc.zc))
    49.266
    >>> print ("%.3f" % (np.angle (mc.zc) / np.pi * 180))
    -0.791

    >>> mc = Measured_Cable (1.9e6, 1, 2.0 + 50.0j, 0.6 - 48.5j)
    >>> print ("%.3f" % abs (mc.zc))
    49.266
    >>> print ("%.3f" % (np.angle (mc.zc) / np.pi * 180))
    -0.791

#    >>> print ("%.5f" % (mc.phi / np.pi * 180))
#    45.44678

    #>>> "%.5f" % (mc.attenuation_db * ft_per_m)
    #>>> mc.attenuation_neper
    #>>> mc.len_per_mhz

    >>> mc = Measured_Cable (5.65e6, 1, 2.1, 2.58, -50.177, 47.778)
    >>> print ("%.3f" % abs (mc.zc))
    49.020
    >>> print ("%.3f" % (np.angle (mc.zc) / np.pi * 180))
    -0.347

    >>> mc = Measured_Cable (5.7e6, 1, 1.9, 2.58, -48.42, 49.47)
    >>> print ("%.3f" % abs (mc.zc))
    48.994
    >>> print ("%.3f" % (np.angle (mc.zc) / np.pi * 180))
    -0.369
    """

    def __init__ (self, f0, length, r_sc, r_oc, x_sc = 0, x_oc = 0) :
        self.f0     = f0
        self.length = length
        self.z_sc   = r_sc + 1j * x_sc
        self.z_oc   = r_oc + 1j * x_oc
        # See Chipman formula (7.27)
        self.zc     = np.sqrt (self.z_sc * self.z_oc)
        # Seems atanh can't deal with complex numbers
        #self.gamma  = atanh (np.sqrt (self.z_sc / self.z_oc) / length)
        #self.beta   = self.gamma.real
        #self.lamda  = 2 * np.pi / self.beta
        #self.phi    = length / self.lamda * 2 * np.pi
        #alpha       = gamma.
    # end def __init__

# end class Measured_Cable
