#!/usr/bin/python3

import numpy as np
from math import atanh
from scipy.optimize import curve_fit
from rsclib.capacitance import c
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

m_per_ft = 0.3048

eps = 1e-4

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
    [6] William E. Sabin. Computer modeling of coax cable circuits.
        QEX, pages 3–10, August 1996.
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
    >>> z0f = cable.z0f (3e9)
    >>> print ("%.3f %+.3fj" % (z0f.real, z0f.imag))
    50.000 +0.008j
    >>> z0f = cable.z0f (cable.fx)
    >>> print ("%.3f %+.3fj" % (z0f.real, z0f.imag))
    50.000 -0.000j

    # Witt Table 3
    >>> f = 28.8e6
    >>> l = 50 * m_per_ft
    >>> z0f = cable.z0f (f)
    >>> print ("%.3f %+.3fj" % (z0f.real, z0f.imag))
    50.002 -0.436j
    >>> print ("%.3fe-3" % (cable.alpha (f) * m_per_ft * 1000))
    2.809e-3
    >>> print ("%.2f" % (cable.phi (f, l) / np.pi * 180))
    798.34
    >>> print (cable.summary (f, l, metric = False))
            Electrical length 798.34°
                  Wavelengths 2.218λ
    Interesting Lengths at 28.8 MHz:
              Half Wavelength 11.27 feet
           Quarter Wavelength  5.64 feet
            Eighth Wavelength  2.82 feet
    Interesting Frequencies for 50.00 feet:
              Half Wavelength 6.494 MHz
           Quarter Wavelength 3.247 MHz
            Eighth Wavelength 1.623 MHz
     Characteristic Impedance 50.002 -0.436j Ohm
         Attenuation Constant 2.809e-3 nepers/foot
               Phase Constant 0.279 rad/foot
              Resistance/foot 0.2619 Ohm
              Inductance/foot 0.077 μH
             Conductance/foot 7.616 μS
             Capacitance/foot 30.800 pF
                 Matched Loss 1.220 dB
    >>> f, l = 14e6, 100 * m_per_ft
    >>> sm = cable.summary_match
    >>> print (sm (f, l, 1500, z_l = 50 -500j, metric = False))
    100.00 feet at 14.00 MHz with 1500 W applied
               Load impedance 50.000 -500.000j Ohm
              Input impedance 12.374 -25.607j Ohm
                 Matched Loss 1.661 dB
                   Total Loss 13.148 dB
             abs(rho) at load 0.981
                 VSWR at load 101.990
            abs(rho) at input 0.675
                VSWR at input 5.154

    >>> print (sm (f, l, 1500, z_i = 12.374351 -25.607388j, metric = False))
    100.00 feet at 14.00 MHz with 1500 W applied
               Load impedance 50.000 -500.000j Ohm
              Input impedance 12.374 -25.607j Ohm
                 Matched Loss 1.661 dB
                   Total Loss 13.148 dB
             abs(rho) at load 0.981
                 VSWR at load 101.990
            abs(rho) at input 0.675
                VSWR at input 5.154


    # Sabin [6] example worksheet
    >>> mpf = 0.3047 # Wrong value conversion from foot by Sabin (sic)
    >>> l = 50 * mpf
    >>> f = 7.15e6
    >>> z_l = 43 + 30j
    >>> p = 100
    >>> cable = Manufacturer_Data_Cable (50, .66)
    >>> cable.set_loss_constants (f, 0.57 / mpf)
    >>> cable.set_freq_params (f, l, p, z_l)
    >>> print ("alpha: %.3f * 10^-3" % (cable.alpha (f) * 1000))
    alpha: 2.154 * 10^-3
    >>> z0f = cable.z0f (f)
    >>> print ("%.3f %+.3fj" % (z0f.real, z0f.imag))
    50.002 -0.474j
    >>> zin = cable.z_i

    # Looks like here values from Sabin differ slightly
    >>> print ("%.3f %+.3fj" % (zin.real, zin.imag))
    65.695 +31.901j
    >>> print ("%.3f" % abs (zin))
    73.031
    >>> u = cable.U_i
    >>> print ("%.2f" % u)
    90.10
    >>> u_l = cable.U_l

    # Here it differs from Sabin in that the real part is negative
    >>> print ("%.3f %+.3fj" % (u_l.real, u_l.imag))
    -75.341 +15.475j
    >>> print ("%.3f" % abs (u_l))
    76.914
    >>> p_l = cable.P_l
    >>> print ("%.3f" % p_l)
    92.535
    >>> print ("%.3f" % cable.combined_loss)
    0.337
    """

    def __init__ (self, Z0, vf, Cpl = None, use_sabin = False) :
        self.Z0 = Z0
        self.vf = vf
        if Cpl is None :
            self.Cpl = 1.0 / (c * Z0 * vf)
        else :
            self.Cpl = Cpl
            self.vf  = 1.0 / (c * Z0 * Cpl)
        self.f0 = self.a0r = self.a0g = self.g = None
        self.use_sabin = use_sabin
    # end def __init__

    def alpha (self, f) :
        """ Attenuation in nepers / m
        """
        return self.loss (f) * np.log (10) / 20 / 100
    # end def alpha

    def beta (self, f) :
        """ Phase constant in rad / m
        """
        return 2 * np.pi * f / (c * self.vf)
    # end def beta

    def gamma (self, f = None) :
        """ Propagation constant = alpha + j beta
        """
        if f is None :
            f = self.f
        return self.alpha (f) + 1j * self.beta (f)
    # end def gamma

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

    def set_loss_constants (self, f0, a0r, a0g = 0) :
        """ Alternative to fit if we have less data, allows to specify
            a0r (and optionally a0g) for a single frequency.
        """
        self.f0  = f0
        self.a0r = a0r
        self.a0g = a0g
        self.g   = 1.0
    # end def set_loss_constants

    def set_freq_params (self, f, l, p, z_l = None, z_i = None) :
        self.f = f
        self.l = l
        self.p = p
        assert z_l or z_i
        assert not (z_l and z_i)
        if z_l is not None :
            self.z_l = z_l
            self.z_i = self.z_d (f, l, z_l)
        elif z_i is not None :
            self.z_i = z_i
            self.z_l = self.z_load (l, f, z_i)
    # end def set_freq_params

    def freq (self, l) :
        """ Frequency at given length lambda
        """
        return c * self.vf / l
    # end def freq

    @property
    def fx (self) :
        """ Frequency where the impedance is real.
            This is the point where loss_g and loss_r are equal.
        """
        return self.f0 * (self.a0r / self.a0g) ** (1.0 / (self.g - .5))
    # end def fx

    @property
    def L (self) :
        return self.Cpl * self.Z0 ** 2
    # end def L

    @property
    def P_l (self) :
        return abs (self.U_l) ** 2 * (1.0 / self.z_l).real
    # end def P_l

    @property
    def combined_loss (self) :
        return 10 * np.log (self.p / self.P_l) / np.log (10)
    # end def combined_loss

    @property
    def U_i (self) :
        return np.sqrt (self.p / (1.0 / self.z_i).real)
    # end def U_i

    @property
    def U_l (self) :
        u = self.U_i
        z0 = self.z0f   ()
        gm = self.gamma ()
        zi = self.z_i
        return u * (np.cosh (gm * self.l) - z0 / zi * np.sinh (gm * self.l))
    # end def U_l

    def lamda (self, f = None) :
        if f is None :
            f = self.f
        return 2 * np.pi / self.beta (f)
    # end def lamda

    def phi (self, f = None, l = None) :
        """ Electrical length at given frequency (in rad)
        """
        if f is None :
            f = self.f
        if l is None :
            l = self.l
        return self.beta (f) * l
    # end def phi

    def conductance (self, f, z0 = None) :
        ln10 = np.log (10)
        if z0 is None :
            z0 = self.Z0
        g = self.loss_g (f) * ln10 * z0.real / (1000.0 * abs (z0) ** 2)
        return g
    # end def conductance

    def resistance (self, f, z0 = None) :
        ln10 = np.log (10)
        if z0 is None :
            z0 = self.Z0
        r = self.loss_r (f) * ln10 * z0.real / 1000.0
        return r
    # end def resistance

    def summary (self, f, l, metric = True) :
        unit = units = 'm'
        cv   = 1.0
        if not metric :
            unit  = 'foot'
            units = 'feet'
            cv    = m_per_ft
        ohm = '\u3a09'
        ohm = 'Ohm'
        lm  = '\u03bb'
        mu  = '\u03bc'
        r = []
        z0f = self.z0f (f)
        phi = self.phi (f, l)
        r.append ('%25s %.2f°' % ('Electrical length', phi / np.pi * 180))
        ll = l / self.lamda (f)
        r.append ('%25s %.3f%s' % ('Wavelengths', ll, lm))
        r.append ('Interesting Lengths at %.1f MHz:' % (f * 1e-6))
        r.append \
            ( '%25s %5.2f %s'
            % ('Half Wavelength',    (self.lamda (f) / 2 / cv), units)
            )
        r.append \
            ( '%25s %5.2f %s'
            % ('Quarter Wavelength', (self.lamda (f) / 4 / cv), units)
            )
        r.append \
            ( '%25s %5.2f %s'
            % ('Eighth Wavelength',  (self.lamda (f) / 8 / cv), units)
            )
        
        freq = self.freq (l)
        r.append ('Interesting Frequencies for %.2f %s:' % (l / cv, units))
        r.append ('%25s %.3f MHz' % ('Half Wavelength',    (freq / 2e6)))
        r.append ('%25s %.3f MHz' % ('Quarter Wavelength', (freq / 4e6)))
        r.append ('%25s %.3f MHz' % ('Eighth Wavelength',  (freq / 8e6)))
        r.append \
            ( '%25s %.3f %+.3fj %s'
            % ('Characteristic Impedance', z0f.real, z0f.imag, ohm)
            )
        r.append \
            ( '%25s %.3fe-3 nepers/%s'
            % ('Attenuation Constant', self.alpha (f) * cv * 1000, unit)
            )
        r.append \
            ( '%25s %.3f rad/%s'
            % ('Phase Constant', self.beta (f) * cv, unit)
            )
        r.append \
            ( '%25s %.4f %s'
            % ('Resistance/%s' % unit, self.resistance (f) * cv, ohm)
            )
        r.append \
            ( '%25s %.3f %sH'
            % ('Inductance/%s' % unit, self.L * cv * 1e6, mu)
            )
        r.append \
            ( '%25s %.3f %sS'
            % ('Conductance/%s' % unit, self.conductance (f) * cv * 1e6, mu)
            )
        r.append \
            ( '%25s %.3f pF'
            % ('Capacitance/%s' % unit, self.Cpl * cv * 1e12)
            )
        r.append \
            ( '%25s %.3f dB'
            % ('Matched Loss', self.loss (f) / 100 * l)
            )

        return '\n'.join (r)
    # end def summary

    def summary_match (self, f, l, p, z_l = None, z_i = None, metric = True) :
        """ Return summary of matching parameters given f, l, power p
            and load impedance z_l or input impedance z_i
        """
        self.set_freq_params (f, l, p, z_l, z_i)
        z_l = self.z_l
        z_i = self.z_i
        unit = units = 'm'
        cv   = 1.0
        if not metric :
            unit  = 'foot'
            units = 'feet'
            cv    = m_per_ft
        ohm = 'Ohm'
        r = []
        r.append \
            ( '%.2f %s at %.2f MHz with %.0f W applied'
            % (l / cv, units, f * 1e-6, p)
            )
        r.append \
            ( '%25s %.3f %+.3fj %s'
            % ('Load impedance', z_l.real, z_l.imag, ohm)
            )
        r.append \
            ( '%25s %.3f %+.3fj %s'
            % ('Input impedance', z_i.real, z_i.imag, ohm)
            )
        r.append \
            ( '%25s %.3f dB'
            % ('Matched Loss', self.loss (f) / 100 * l)
            )
        u = self.U_i
        r.append ('%25s %.3f dB' % ('Total Loss', self.combined_loss))
        rho_l  = (z_l - self.Z0) / (z_l + self.Z0)
        vswr_l = (1 + abs (rho_l)) / (1 - abs (rho_l))
        rho_i  = (z_i - self.Z0) / (z_i + self.Z0)
        vswr_i = (1 + abs (rho_i)) / (1 - abs (rho_i))
        r.append ('%25s %.3f' % ('abs(rho) at load', abs (rho_l)))
        r.append ('%25s %.3f' % ('VSWR at load', vswr_l))
        r.append ('%25s %.3f' % ('abs(rho) at input', abs (rho_i)))
        r.append ('%25s %.3f' % ('VSWR at input', vswr_i))
        return '\n'.join (r)
    # end def summary_match

    def z0f_witt (self, f, z0, r, g) :
        """ From Witt [3]
        """
        cpart = 2j * np.pi * f * self.Cpl
        return np.sqrt ((r + cpart * z0 ** 2) / (g + cpart))
    # end def z0f_witt

    def z0f_sabin (self, f, z0, r, g) :
        """ From Sabin [6] p.4:
        """
        l = self.Cpl * z0 ** 2
        o2 = 4 * np.pi * f
        return z0.real * (1 - 1j * (r / (o2 * l) - g / (o2 * self.Cpl)))
    # end def z0f_sabin

    def z0f (self, f = None, z0 = None) :
        """ Characteristic impedance X_0f for a given frequency
            *without* the high-frequency approximation. We have from
            Chipman [1] Formula 4.12 p.32:
            Z_0 = sqrt ((R + j omega L) / (G + j omega C))
            And from the high-frequency approximation we have
            Z0 = sqrt (L/C)
            Solving for L:
            L = C * Z0 ** 2
            We get
            (1) Z_0 = sqrt ((R + j omega C * Z0 ** 2) / (G + j omega C))
            From Eq 5.31 p. 59 [1] we have:
            alpha = R / (2R_0) + (G * abs (Z_0) ** 2) / (2R_0)
            which we can separate into the losses due to R and due to G:
            (2) alpha_R = R / (2R_0)
            (3) alpha_G = (G * abs (Z_0) ** 2) / (2R_0)
            These are in nepers/m, so when converting to dB/100m and
            solving for R and G we get:
            R = 2 * alpha_R R_0 = A_R ln (10) R_0 / 1000
            G = 2 * alpha_G R_0 / (abs (Z_0) ** 2)
              = R_0 A_G ln (10) / (1000 * (abs (Z_0) ** 2))
            Z_0 is a function of f, Z_0 (f) R_0 (f) is the real part.
            A_R is self.loss_r (f)
            A_G is self.loss_g (f)
            Witt tells us [3] can use this to compute Z_0 (f)
            iteratively, this apparently is not true, it diverges on
            repeated iteration.
            Another formula from Sabin [6] (which is non-iterative)
            yields almost the same result as the non-iterated formula of
            Witt [3]
        """
        if f is None :
            f = self.f
        if z0 is None :
            z0 = self.Z0
        r = self.resistance  (f, z0)
        g = self.conductance (f, z0)
        if self.use_sabin :
            return self.z0f_sabin (f, z0, r, g)
        return self.z0f_witt  (f, z0, r, g)
    # end def z0f

    def plot_z0f (self, f) :
        """ Since Witt [3] tells us that the z0_f formula above can be
            iterated, this plots the absolute value of the difference of
            an input z0 with an output z0. It clearly shows this diverges.
        """
        x = np.arange (-1.5, 1.5, 0.02)
        y = np.arange (49-0.5, 51+0.5, 0.02)
        z = []
        for vx in x :
            zz = []
            z.append (zz)
            for vy in y :
                z0  = 50 + 0j
                z0f = self.z0f (f, vy + vx*1j)
                zz.append (abs (z0 - z0f))
                #zz.append ((z0 - z0f).real)
                #zz.append ((z0 - z0f).imag)
        z    = np.array (z)
        x, y = np.meshgrid (x, y)
        fig  = plt.figure ()
        ax   = fig.add_subplot (111, projection = '3d')
        ax.plot_wireframe (x, y, z)
        plt.show ()
    # end def plot_z0f

    def z_d (self, f, d, z_l = None) :
        """ Compute impedance at distance d from load, given load
            impedance z_l, Eq 24 from [6] p.7
            Special case z_l = None is open circuit.
        """
        z0 = self.z0f   (f)
        gm = self.gamma (f)
        ep = np.e ** (gm * d)
        em = np.e ** (-gm * d)
        zz = z_l / z0
        if z_l is None :
            return z0 * (ep + em) / (ep - em)
        return z0 * ( (ep * (zz + 1) + (zz - 1) / ep)
                    / (ep * (zz + 1) - (zz - 1) / ep)
                    )
    # end def z_d

    def z_d_open (self, f, d) :
        """ Same as z_d but for open circuit (infinite z_l)
        """
        return self.z_d (f, d, z_l = None)
    # end def z_d_open

    def z_d_short (self, f, d) :
        """ Same as z_d but for short circuit (z_l = 0)
        """
        return self.z_d (f, d, z_l = 0.0)
    # end def z_d_short

    def z_load (self, d, f = None, z_i = None) :
        """ Compute impedance at load from input impedance given length
        """
        if f is None :
            f = self.f
        if z_i is None :
            z_i = self.z_i
        z0 = self.z0f (f)
        gm = self.gamma (f)
        ep = np.e ** (gm * d)
        return -((z0 * z_i - z0**2) * ep**2 + z0 * z_i + z0**2) \
                / ((z_i - z0) * ep**2 - z_i - z0)
    # end def z_load

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
