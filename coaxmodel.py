#!/usr/bin/python3

import numpy as np
from math import atanh
from scipy.optimize import curve_fit
from rsclib.capacitance import c

m_per_ft = 0.3048

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
        Note that with our fit method we obtain a better solution than
        Witt who requires two frequencies matching exactly.
    >>> print ("%.3f" % (1.0 / c * 1e12 * m_per_ft))
    1016.703
    >>> Cpl = 30.8e-12 / m_per_ft
    >>> cable = Manufacturer_Data_Cable (50, .66)
    >>> print ("%.1f pF/ft" % (cable.Cpl * 1e12 * m_per_ft))
    30.8 pF/ft
    >>> cable = Manufacturer_Data_Cable (50, .66, Cpl)
    >>> print ("vf = %.6f" % cable.vf)
    vf = 0.660197
    >>> l = []
    >>> l.append ((1e6,    0.44 / m_per_ft))
    >>> l.append ((10e6,   1.4  / m_per_ft))
    >>> l.append ((50e6,   3.3  / m_per_ft))
    >>> l.append ((100e6,  4.9  / m_per_ft))
    >>> l.append ((200e6,  7.3  / m_per_ft))
    >>> l.append ((400e6, 11.5  / m_per_ft))
    >>> l.append ((700e6, 17.0  / m_per_ft))
    >>> l.append ((900e6, 20.0  / m_per_ft))
    >>> l.append ((1e9,   21.5  / m_per_ft))
    >>> cable.fit (l)
    >>> for x, y in l :
    ...    y = y * m_per_ft
    ...    a = cable.loss (x) * m_per_ft
    ...    print ("%4.0f %5.2f %5.2f %6.2f" % (x / 1e6, y, a, (a - y)))
       1  0.44  0.41  -0.03
      10  1.40  1.35  -0.05
      50  3.30  3.26  -0.04
     100  4.90  4.88  -0.02
     200  7.30  7.42   0.12
     400 11.50 11.55   0.05
     700 17.00 16.82  -0.18
     900 20.00 20.02   0.02
    1000 21.50 21.57   0.07
    >>> print ("f0:  %6.2f MHz" % (cable.f0 / 1e6))
    f0:  500.50 MHz
    >>> print ("a0r (f0): %.3f dB/100ft"  % (cable.a0r * m_per_ft))
    a0r (f0): 8.887 dB/100ft
    >>> print ("a0g (f0): %.3f dB/100ft" % (cable.a0g * m_per_ft))
    a0g (f0): 4.511 dB/100ft
    >>> print ("g:   %.4f" % cable.g)
    g: 0.9990
    >>> for x in (1, 5, 10, 50, 100, 500, 500.5, 1000, 10000) :
    ...    x *= 1e6
    ...    a = cable.loss (x) * m_per_ft
    ...    lr = cable.loss_r (x) * m_per_ft
    ...    lg = cable.loss_g (x) * m_per_ft
    ...    print ("%7.1f %6.2f %5.2f %5.2f" % (x / 1e6, a, lr, lg))
        1.0   0.41  0.40  0.01
        5.0   0.93  0.89  0.05
       10.0   1.35  1.26  0.09
       50.0   3.26  2.81  0.45
      100.0   4.88  3.97  0.90
      500.0  13.39  8.88  4.51
      500.5  13.40  8.89  4.51
     1000.0  21.57 12.56  9.01
    10000.0 129.57 39.72 89.84
    >>> print ("%7.3f" % (cable.fx / 1e6))
    1948.522
    >>> print ("%5.2f" % cable.loss_g (cable.fx))
    57.53
    >>> print ("%5.2f" % cable.loss_r (cable.fx))
    57.53

    Try fitting to Witts computed values should get us his value for g
    >>> l = []
    >>> l.append ((1e6,    0.43 / m_per_ft))
    >>> l.append ((10e6,   1.39 / m_per_ft))
    >>> l.append ((50e6,   3.30 / m_per_ft))
    >>> l.append ((100e6,  4.89 / m_per_ft))
    >>> l.append ((200e6,  7.39 / m_per_ft))
    >>> l.append ((400e6, 11.46 / m_per_ft))
    >>> l.append ((700e6, 16.74 / m_per_ft))
    >>> l.append ((900e6, 20.00 / m_per_ft))
    >>> l.append ((1e9,   21.58 / m_per_ft))
    >>> cable.fit (l)
    >>> print ("%.2f" % cable.g)
    1.10
    >>> for x, y in l :
    ...    a = cable.loss (x) * m_per_ft
    ...    print ("%4.0f %5.2f" % (x / 1e6, a))
       1  0.43
      10  1.39
      50  3.30
     100  4.89
     200  7.39
     400 11.46
     700 16.74
     900 20.00
    1000 21.58
    >>> print ("%7.3f" % (cable.fx / 1e6))
    2280.569
    >>> print ("%5.2f" % cable.loss_g (cable.fx))
    66.41
    >>> print ("%5.2f" % cable.loss_r (cable.fx))
    66.41
    """

    def __init__ (self, Z0, vf, Cpl = None) :
        self.Z0 = Z0
        self.vf = vf
        if Cpl is None :
            self.Cpl = 1.0 / (c * Z0 * vf)
        else :
            self.Cpl = Cpl
            self.vf  = 1.0 / (c * Z0 * Cpl)
        self.f0 = self.a0r = self.a0g = self.g = None
    # end def __init__

    def loss_r (self, f, a0r = None) :
        """ Loss due to skin effect (R loss)
        """
        a0r = a0r or self.a0r
        return a0r * np.sqrt (f / self.f0)
    # end def loss_r

    def loss_g (self, f, a0g = None, g = None) :
        """ Dieelectric loss
        """
        a0g = a0g or self.a0g
        g   = g   or self.g
        return a0g * (f / self.f0) ** g
    # end def loss_g

    def loss (self, f, a0r = None, a0g = None, g = None) :
        return self.loss_r (f, a0r) + self.loss_g (f, a0g, g)
    # end def loss 

    def fit (self, loss_data) :
        """ Gets a list of frequency/loss pairs
            Note that the loss is in dB per 100m (not ft)
        """
        self.loss_data = np.array (loss_data)
        self.f0 = (self.loss_data [0][0] + self.loss_data [-1][0]) / 2.0
        x = self.loss_data [:,0]
        y = self.loss_data [:,1]
        popt, pcov = curve_fit (self.loss, x, y)
        self.a0r, self.a0g, self.g = popt
    # end def fit

    @property
    def fx (self) :
        """ Frequency where the impedance is real.
            This is the point where loss_g and loss_r are equal.
        """
        return self.f0 * (self.a0r / self.a0g) ** (1.0 / (self.g - .5))
    # end def fx

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

    #>>> "%.5f" % (mc.attenuation_db * m_per_ft)
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
