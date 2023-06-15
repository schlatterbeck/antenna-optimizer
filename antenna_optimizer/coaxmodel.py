#!/usr/bin/python3

import numpy as np
from math import atanh
from scipy.optimize import curve_fit
from argparse import ArgumentParser
from rsclib.capacitance import c
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

m_per_ft = 0.3048

eps = 1e-6

""" Transmission Line Properties of Cables
    Either from manufacturer data or from measured parameters.
    Inspired by two articles by Frank Witt [3,4] in an old ARRL Antenna
    Compendium [5]. Unfortunately Witt uses a now obsolete software and
    doesn't publish all formulas in the paper. Most, however, are in
    Chipman [1] and Walter Maxwell [2] has a Basic implementation from
    which some of the formulas can be reverse engineered (the Basic
    dialect used doesn't support complex numbers so everything there
    looks very bare-footed from todays perspective).
    [1]  Robert A. Chipman. Theory and Problems of Transmission Lines.
         Schaums Outline. McGraw-Hill, 1968. Section 7.5 "Determination
         of transmission line characteristics from impedance
         measurements", p.  134.
    [2]  M. Walter Maxwell. Reflections III, Transmission Lines
         and Antennas. CQ Communications, Inc., Hicksville, NY, third
         edition, 2010. In particular chapter 15 introducing computer
         programs (in Basic), but this becomes much easier and is
         closer to Chipman with direct use of complex numbers. I'm
         using Maxwells results to test this implementation.
    [3]  Frank Witt. Transmission line properties from manufacturer’s
         data. In Straw [5], pages 179–183.
    [4]  Frank Witt. Transmission line properties from measured data.
         In Straw [5], pages 184–188.
    [5]  R. Dean Straw, editor. The ARRL Antenna Compendium, volume 6.
         American Radio Relay League (ARRL), 1999.
    [6]  William E. Sabin. Computer modeling of coax cable circuits.
         QEX, pages 3–10, August 1996.
    [7]  Walter C. Johnson. Transmission Lines and Networks.
         McGraw-Hill, international student edition, 1950.
    [8]  Frederick Emmons Terman. Radio Engineers’ Handbook.
         McGraw-Hill, first edition, 1943.
    [9]  Frank Witt. The coaxial resonator match. In Hall [10], pages 110–118.
    [10] Gerald L. (Jerry) Hall, editor. The ARRL Antenna Compendium,
         volume 2. American Radio Relay League (ARRL), 1989.


"""

class Manufacturer_Data_Cable:
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
    >>> print (cable.summary_loss_data (metric = False))
    f (MHz)   Manu    Fit   Diff   (in dB/100 feet)
          1   0.44   0.41  -0.03
         10   1.40   1.35  -0.05
         50   3.30   3.26  -0.04
        100   4.90   4.88  -0.02
        200   7.30   7.42   0.12
        400  11.50  11.55   0.05
        700  17.00  16.82  -0.18
        900  20.00  20.02   0.02
       1000  21.50  21.57   0.07
    >>> print ("f0:  %6.2f MHz" % (cable.f0 / 1e6))
    f0:  500.50 MHz
    >>> print ("a0r (f0): %.3f dB/100ft"  % (cable.a0r * m_per_ft))
    a0r (f0): 8.887 dB/100ft
    >>> print ("a0g (f0): %.3f dB/100ft" % (cable.a0g * m_per_ft))
    a0g (f0): 4.511 dB/100ft
    >>> print ("g:   %.4f" % cable.g)
    g:   0.9990
    >>> for x in (1, 5, 10, 50, 100, 500, 500.5, 1000, 10000):
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
    >>> print (cable.summary_loss_data (metric = False))
    f (MHz)   Manu    Fit   Diff   (in dB/100 feet)
          1   0.43   0.43  -0.00
         10   1.39   1.39   0.00
         50   3.30   3.30   0.00
        100   4.89   4.89  -0.00
        200   7.39   7.39  -0.00
        400  11.46  11.46   0.00
        700  16.74  16.74  -0.00
        900  20.00  20.00  -0.00
       1000  21.58  21.58   0.00
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
     Characteristic Impedance 50.002 -0.436j Ω
         Attenuation Constant 2.809e-3 nepers/foot
               Phase Constant 0.279 rad/foot
              Resistance/foot 0.2619 Ω
              Inductance/foot 0.077 μH
             Conductance/foot 7.616 μS
             Capacitance/foot 30.800 pF
                 Matched Loss 1.220 dB
    >>> print (cable.summary (f, l))
            Electrical length 798.34°
                  Wavelengths 2.218λ
    Interesting Lengths at 28.8 MHz:
              Half Wavelength  3.44 m
           Quarter Wavelength  1.72 m
            Eighth Wavelength  0.86 m
    Interesting Frequencies for 15.24 m:
              Half Wavelength 6.494 MHz
           Quarter Wavelength 3.247 MHz
            Eighth Wavelength 1.623 MHz
     Characteristic Impedance 50.002 -0.436j Ω
         Attenuation Constant 9.217e-3 nepers/m
               Phase Constant 0.914 rad/m
                 Resistance/m 0.8593 Ω
                 Inductance/m 0.253 μH
                Conductance/m 24.987 μS
                Capacitance/m 101.050 pF
                 Matched Loss 1.220 dB
    >>> f, l = 14e6, 100 * m_per_ft
    >>> print (cable.summary (f, l))
            Electrical length 776.16°
                  Wavelengths 2.156λ
    Interesting Lengths at 14.0 MHz:
              Half Wavelength  7.07 m
           Quarter Wavelength  3.53 m
            Eighth Wavelength  1.77 m
    Interesting Frequencies for 30.48 m:
              Half Wavelength 3.247 MHz
           Quarter Wavelength 1.623 MHz
            Eighth Wavelength 0.812 MHz
     Characteristic Impedance 50.005 -0.642j Ω
         Attenuation Constant 6.273e-3 nepers/m
               Phase Constant 0.444 rad/m
                 Resistance/m 0.5991 Ω
                 Inductance/m 0.253 μH
                Conductance/m 11.304 μS
                Capacitance/m 101.050 pF
                 Matched Loss 1.661 dB

    >>> cbl = belden_8295
    >>> sm = cbl.summary_match
    >>> cbl.set_freq_params (28e6, 15, 100, z_l = 77+15j)
    >>> print (sm ())
    15.00 m at 28.00 MHz with 100 W applied
               Load impedance 77.000 +15.000j Ω
              Input impedance 54.963 -19.658j Ω
                 Matched Loss 1.159 dB
                   Total Loss 1.252 dB
             abs(rho) at load 0.242
                 VSWR at load 1.637
            abs(rho) at input 0.190
                VSWR at input 1.469
              Maximum Voltage 85.69 V RMS
              Maximum Current 1.71 A RMS
    Inductive stub with short circuit at end:
                Stub attached 1.22967 m from load
                  Stub length 1.25091 m
          Resulting impedance 50.00 +0.00j
    Inductive stub with open circuit at end:
                Stub attached 1.20655 m from load
                  Stub length 3.02279 m
          Resulting impedance 50.00 -0.00j
    Capacitive stub with short circuit at end:
                Stub attached 2.77576 m from load
                  Stub length 2.25368 m
          Resulting impedance 50.00 -0.00j
    Capacitive stub with open circuit at end:
                Stub attached 2.75264 m from load
                  Stub length 0.49129 m
          Resulting impedance 50.00 -0.00j

    >>> sm = cable.summary_match
    >>> cable.set_freq_params (f, l, 1500, z_l = 50 -500j)
    >>> print (sm (metric = False))
    100.00 feet at 14.00 MHz with 1500 W applied
               Load impedance 50.000 -500.000j Ω
              Input impedance 12.374 -25.607j Ω
                 Matched Loss 1.661 dB
                   Total Loss 13.148 dB
             abs(rho) at load 0.981
                 VSWR at load 101.990
            abs(rho) at input 0.675
                VSWR at input 5.154
              Maximum Voltage 621.73 V RMS
              Maximum Current 12.43 A RMS
    Inductive stub with short circuit at end:
                Stub attached 9.56111 feet from load
                  Stub length 1.33886 feet
          Resulting impedance 50.00 +0.00j
    Inductive stub with open circuit at end:
                Stub attached 9.14802 feet from load
                  Stub length 13.31534 feet
          Resulting impedance 50.00 -0.00j
    Capacitive stub with short circuit at end:
                Stub attached 12.91373 feet from load
                  Stub length 21.14709 feet
          Resulting impedance 50.00 +0.00j
    Capacitive stub with open circuit at end:
                Stub attached 12.57213 feet from load
                  Stub length 9.85829 feet
          Resulting impedance 50.00 -0.00j


    >>> d = cable.d_voltage_min ()
    >>> print ("%.5f" % d)
    3.31225
    >>> print ("%.5f" % cable.d_voltage_min (cable.vswr_l * cable.Z0))
    6.84657
    >>> print ("%.5f" % cable.d_voltage_min (1 / cable.vswr_l * cable.Z0))
    3.31225
    >>> zd = 0.6309966224159409-26.78963638394159j
    >>> print ("%.5f" % cable.d_voltage_min (zd))
    2.20542
    >>> print ("%.5f" % cable.d_voltage_min (50 -500j))
    -0.00000
    >>> print ("%.5f" % cable.d_voltage_min (12.374351 -25.607388j))
    2.20171
    >>> z_l = 50 -500j
    >>> zd = cable.z_d (cable.f, 30.48 - 2*cable.lamda (), z_l)
    >>> print ("%.5f" % cable.d_voltage_min (zd))
    2.20464
    >>> d, y = cable.stub_match (capacitive = True)
    >>> z = (1 / (y.imag * 1j)).imag
    >>> print ("%.6f" % d)
    2.943836
    >>> print ("%.5f" % (d / cable.lamda ()))
    0.20823
    >>> print ("%.6f" % z)
    -8.519958
    >>> print ("%.5f" % cable.d_voltage_min ())
    3.31225
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.6f" % l_sc)
    0.379754
    >>> l_sc = cable.stub_short_iter (-z)
    >>> print ("%.6f" % l_sc)
    0.379512
    >>> z_m = cable.stub_impedance (cable.f, d, l_sc, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    43.115 -0.000j

    # Now try iterative stub matching
    >>> d, l = cable.stub_match_iterative ()
    >>> print ("%.6f" % d)
    2.914225
    >>> print ("%.6f" % l)
    0.408086
    >>> z_m = cable.stub_impedance (cable.f, d, l, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    # This is the result of nec simulation
    >>> z_l = 50 -500j
    >>> d, l = 2.9175946668, 0.4033092359
    >>> z_m = cable.stub_impedance (cable.f, d, l, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    49.142 +1.021j

    # Same as above but with open circuit match
    >>> d, l = 2.7938663047, 4.0508900506
    >>> z_m = cable.stub_impedance (cable.f, d, l, z_l, shortcircuit = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    48.784 +0.837j

    # Same as above but at capacitive impedance at length l, short circuit
    >>> d, l = 3.9268108472, 6.4538558619
    >>> z_m = cable.stub_impedance (cable.f, d, l, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    48.669 +0.323j

    # Same as above but at capacitive impedance at length l, open circuit
    >>> d, l = 3.8244235289, 3.0105401976
    >>> z_m = cable.stub_impedance (cable.f, d, l, z_l, shortcircuit = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    48.721 +0.680j

    # Another result of a nec simulation
    >>> z_l  = 500.
    >>> d, l = 2.7777599127, 0.8415883400
    >>> z_m = cable.stub_impedance (cable.f, d, l, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    49.941 +0.314j

    >>> f, l = 14e6, 100 * m_per_ft
    >>> z_i = 12.374351 -25.607388j
    >>> cable.set_freq_params (f, l, 1500, z_l = None, z_i = z_i)
    >>> print (sm (metric = False))
    100.00 feet at 14.00 MHz with 1500 W applied
               Load impedance 50.000 -500.000j Ω
              Input impedance 12.374 -25.607j Ω
                 Matched Loss 1.661 dB
                   Total Loss 13.148 dB
             abs(rho) at load 0.981
                 VSWR at load 101.990
            abs(rho) at input 0.675
                VSWR at input 5.154
              Maximum Voltage 621.73 V RMS
              Maximum Current 12.43 A RMS
    Inductive stub with short circuit at end:
                Stub attached 9.56111 feet from load
                  Stub length 1.33886 feet
          Resulting impedance 50.00 +0.00j
    Inductive stub with open circuit at end:
                Stub attached 9.14802 feet from load
                  Stub length 13.31534 feet
          Resulting impedance 50.00 -0.00j
    Capacitive stub with short circuit at end:
                Stub attached 12.91373 feet from load
                  Stub length 21.14709 feet
          Resulting impedance 50.00 -0.00j
    Capacitive stub with open circuit at end:
                Stub attached 12.57213 feet from load
                  Stub length 9.85829 feet
          Resulting impedance 50.00 -0.00j

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

    # Example with lossless 50 Ohm line terminated in 150 Ohm
    >>> cable = Manufacturer_Data_Cable (50, 1.0)
    >>> f = 3.5e6
    >>> cable.set_loss_constants (f, 0.0)
    >>> z_l = 150.0
    >>> cable.set_freq_params (f, 1, 1, z_l)
    >>> print ("%.2f" % (cable.lamda () / 2.0))
    42.83
    >>> print ("%.2f" % (cable.lamda () / 4.0))
    21.41

    >>> d, y = cable.stub_match (capacitive = True)
    >>> z = (1 / (y.imag * 1j)).imag
    >>> print ("%.2f" % (cable.lamda () / 6))
    14.28
    >>> print ("%.2f" % d)
    14.28
    >>> print ("%.3f" % z)
    -43.301
    >>> print ("%.3f" % (1.0 / (z / cable.Z0)))
    -1.155
    >>> l_oc = cable.stub_open (-z)
    >>> print ("%.2f" % l_oc)
    31.14
    >>> l_oc = cable.stub_open_iter (-z)
    >>> print ("%.2f" % l_oc)
    31.14
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    9.73
    >>> l_sc = cable.stub_short_iter (-z)
    >>> print ("%.2f" % l_sc)
    9.73
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, shortcircuit = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    >>> d, y = cable.stub_match (capacitive = False)
    >>> z = (1 / (y.imag * 1j)).imag
    >>> print ("%.2f" % (cable.lamda () / 3))
    28.55
    >>> print ("%.2f" % d)
    28.55
    >>> print ("%.3f" % z)
    43.301
    >>> print ("%.3f" % (1.0 / (z / cable.Z0)))
    1.155
    >>> l_oc = cable.stub_open (-z)
    >>> print ("%.2f" % l_oc)
    11.68
    >>> l_oc = cable.stub_open_iter (-z)
    >>> print ("%.2f" % l_oc)
    11.68
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    33.10
    >>> l_sc = cable.stub_short_iter (-z)
    >>> print ("%.2f" % l_sc)
    33.10
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, shortcircuit = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    # Example with 50 Ohm line terminated in succeptance 2 -2j Chart Siemens
    # that's 0.25+0.25j chart Ohm or 12.5+12.5j Ohm
    >>> z_l = 12.5 +12.5j
    >>> cable.set_freq_params (f, 1, 1, 12.5 +12.5j)

    >>> d, y = cable.stub_match (capacitive = False)
    >>> z = (1 / (y.imag * 1j)).imag
    >>> print ("%.4f" % d)
    2.6085
    >>> print ("%.4f" % (d / cable.lamda ()))
    0.0305
    >>> print ("%.3f" % (1.0 / (z / cable.Z0)))
    1.581
    >>> l_oc = cable.stub_open (-z)
    >>> print ("%.2f" % l_oc)
    13.73
    >>> l_oc = cable.stub_open_iter (-z)
    >>> print ("%.2f" % l_oc)
    13.73
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    35.14
    >>> l_sc = cable.stub_short_iter (-z)
    >>> print ("%.2f" % l_sc)
    35.14
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, shortcircuit = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    >>> d, y = cable.stub_match (capacitive = True)
    >>> z = (1 / (y.imag * 1j)).imag
    >>> print ("%.4f" % d)
    33.1418
    >>> print ("%.4f" % (d / cable.lamda ()))
    0.3869
    >>> print ("%.3f" % (1.0 / (z / cable.Z0)))
    -1.581
    >>> l_oc = cable.stub_open (-z)
    >>> print ("%.2f" % l_oc)
    29.10
    >>> l_oc = cable.stub_open_iter (-z)
    >>> print ("%.2f" % l_oc)
    29.10
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    7.69
    >>> l_sc = cable.stub_short_iter (-z)
    >>> print ("%.2f" % l_sc)
    7.69
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, shortcircuit = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, shortcircuit = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    # Example of series / parallel resonators at 21 MHz Witt [3]
    # Tables 6, 7
    >>> cable = belden_8295
    >>> cable.set_freq_params (21e6, 1, 1, 50)
    >>> print (cable.summary_resonator (metric = False))
    Series Resonator at 21.00 MHz
    Open circuited quarter wave:
      Length: 7.730 feet
     abs (z): 0.895 Ω
           Q: 47.114
    Short circuited half wave:
      Length: 15.461 feet
     abs (z): 1.789 Ω
           Q: 47.114
    Parallel Resonator at 21.00 MHz
    Short circuited quarter wave:
      Length: 7.730 feet
     abs (z): 2795.210 Ω
           Q: 47.114
    Open circuited half wave:
      Length: 15.461 feet
     abs (z): 1398.052 Ω
           Q: 47.114
    >>> cable.set_freq_params (21e6, 1, 1, 50)
    >>> l = cable.lamda () / 4.0
    >>> print ("%.3f" % (l / m_per_ft))
    7.730
    >>> z = cable.z_d_open (cable.f, l)
    >>> y11 = cable.y11 (cable.f, l)
    >>> print ("%.3f" % abs (1 / y11))
    2795.210
    >>> y22 = cable.y22 (cable.f, l, 50)
    >>> y12 = cable.y12 (cable.f, l)
    >>> print ("%.3f %.3f" % (abs (z), abs (1/z)))
    0.895 1.118
    >>> print ("%.3f %+.3fj" % (z.real, z.imag))
    0.895 -0.008j
    >>> print ("%.3f %+.3fj" % ((1/z).real, (1/z).imag))
    1.118 +0.010j
    >>> print ("%.3f" % abs (1 / y11))
    2795.210
    >>> print ("%.3f" % abs (1 / y22))
    49.121
    >>> y11 = cable.y11 (cable.f, l)
    >>> y22 = cable.y22 (cable.f, l, None)
    >>> print ("%.3f" % abs (1 / y22))
    2795.210
    >>> print ("%.3f" % abs (1 / y12))
    50.013

    # Short-circuit should give 1/y11
    >>> z_in = 1 / admittance (y11, y12, y22, 0)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    2795.091 -25.819j
    >>> z_in = 1 / admittance (y11, y12, y22, 50)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    50.002 -0.908j
    >>> z_in = 1 / admittance (y11, y12, y22, 25)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    97.413 -1.761j
    >>> z_in = 1 / admittance (y11, y12, y22, 100)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    25.666 -0.464j
    >>> y22 = cable.y22 (cable.f, l, 50)
    >>> print ("%.3f" % abs (1 / y22))
    49.121
    >>> print ("%.3f" % abs (1 / y12))
    50.013
    >>> z_in = 1 / admittance (y11, y12, y22, 1e40)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    50.002 -0.908j

    # Short-circuit should give 1/y11
    >>> z_in = 1 / admittance (y11, y12, y22, 0)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    2795.091 -25.819j

    # Output impedance is 50, another 50 in parallel should give 25
    # and lambda/4 transformation should give 100
    >>> z_in = 1 / admittance (y11, y12, y22, 50)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    97.413 -1.761j
    >>> z_in = 1 / admittance (y11, y12, y22, 25)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    143.215 -2.571j

    # Now check lambda / 2
    >>> l = cable.lamda () / 2.0
    >>> print ("%.3f" % (l / m_per_ft))
    15.461

    # open circuit case
    >>> y11 = cable.y11 (cable.f, l)
    >>> print ("%.3f" % abs (1 / y11))
    1.789
    >>> y22 = cable.y22 (cable.f, l)
    >>> print ("%.3f" % abs (1 / y22))
    1.789
    >>> y12 = cable.y12 (cable.f, l)
    >>> print ("%.3f" % abs (1 / y12))
    1.790
    >>> z_in = 1 / admittance (y11, y12, y22, 0)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    1.789 -0.017j
    >>> z_in = 1 / admittance (y11, y12, y22, 25)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    26.318 -0.021j
    >>> z_in = 1 / admittance (y11, y12, y22, 50)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    50.000 -0.032j
    >>> z_in = 1 / admittance (y11, y12, y22, 100)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    94.994 -0.074j
    >>> z_in = 1 / admittance (y11, y12, y22, None)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    1397.993 -12.914j
    >>> z_in = 1 / admittance (y11, y12, y22, 1e30)
    >>> print ("%.3f %+.3fj" % (z_in.real, z_in.imag))
    1397.993 -12.914j

    >>> cable.set_freq_params (21e6, 1, 1, 50)
    >>> print (cable.summary_stub (21e6,  100, metric = False))
    Inductive impedance 100.00 Ω
    Closed-Circuited
         Length: 5.44861 feet
             Zi: 4.07 +99.90j Ω
    Effective L: 0.757 µH
              Q: 34.04
    >>> print (cable.summary_stub (21e6, -100, metric = False))
    Capacitive impedance -100.00 Ω
    Open-Circuited
         Length: 2.28175 feet
             Zi: 0.40 -100.00j Ω
    Effective C: 75.785 pF
              Q: 81.29

    """

    def __init__ \
        ( self
        , Z0
        , vf
        , Cpl       = None
        , use_sabin = False
        , name      = None
        , metric    = True
        ):
        self.Z0     = Z0
        self.vf     = vf
        self.metric = metric
        self.name   = name or 'custom cable'
        if Cpl is None:
            self.Cpl = 1.0 / (c * Z0 * vf)
        else:
            self.Cpl = Cpl
            self.vf  = 1.0 / (c * Z0 * Cpl)
        self.f0 = self.a0r = self.a0g = self.g = None
        self.use_sabin = use_sabin
    # end def __init__

    def _units (self, metric = True):
        unit = units = 'm'
        cv   = 1.0
        if metric is None:
            metric = self.metric
        if not metric:
            unit  = 'foot'
            units = 'feet'
            cv    = m_per_ft
        return unit, units, cv
    # end def _units

    def alpha (self, f):
        """ Attenuation in nepers / m
        """
        return self.loss (f) * np.log (10) / 20 / 100
    # end def alpha

    def beta (self, f):
        """ Phase constant in rad / m
        """
        return 2 * np.pi * f / (c * self.vf)
    # end def beta

    def gamma (self, f = None):
        """ Propagation constant = alpha + j beta
        """
        if f is None:
            f = self.f
        return self.alpha (f) + 1j * self.beta (f)
    # end def gamma

    def loss_r (self, f, a0r = None):
        """ Loss due to skin effect (R loss)
        """
        a0r = a0r or self.a0r
        return a0r * np.sqrt (f / self.f0)
    # end def loss_r

    def loss_g (self, f, a0g = None, g = None):
        """ Dieelectric loss
        """
        a0g = a0g or self.a0g
        g   = g   or self.g
        return a0g * (f / self.f0) ** g
    # end def loss_g

    def loss (self, f, a0r = None, a0g = None, g = None):
        return self.loss_r (f, a0r) + self.loss_g (f, a0g, g)
    # end def loss

    def fit (self, loss_data):
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

    def resonator_q (self, f = None):
        """ Resonator Q
            loss_dB = a * l / 100 = a * (lambda / 4) / 100
            loss_factor = R * I**2 / R0 * I**2
            loss_factor = 10 ** (loss_dB / 10)
            R = R0 * loss_factor
            lambda = c * vf / F0
            -> R = R0 * (1 - (10 ** (-a * (lambda) / 100 / 10)))
                 = R0 * (1 - (10 ** (-a * c * vf / (F0 * 100) / 10))
            Q = R0 / R = 1 / (1 - (10 ** (-a * c * vf / (F0 * 100) / 10)))
            This produces Q=44.7 for A=0.345 dB/100ft in the original formula
            The original formula in Witt [9] p.113 is
            (2.774 * F0) / (A * vf)
        >>> f = 3670830.9685955304
        >>> cable = Manufacturer_Data_Cable (50, 0.66)
        >>> cable.set_loss_constants (f, 0.345 / m_per_ft)
        >>> cable.set_freq_params (f, 1, 100, 50)
        >>> print ("%.4f" % cable.loss (f))
        1.1319
        >>> print ("%.4f" % cable.resonator_q ())
        47.9411
        """
        if f is None:
            f = self.f
        return 2 * np.pi \
            / (1 - 10 ** (-self.loss (f) * c * self.vf / (1000 * f)))
    # end def resonator_q

    def resonator_q_approx (self, f = None):
        """ Approximate formula from Chipman [1] p. 222 Eq 10.20
        >>> f = 3670830.9685955304
        >>> cable = Manufacturer_Data_Cable (50, 0.66)
        >>> cable.set_loss_constants (f, 0.345 / m_per_ft)
        >>> cable.set_freq_params (f, 1, 100, 50)
        >>> print ("%.1f" % cable.resonator_q_approx (f))
        44.7
        """
        if f is None:
            f = self.f
        return self.beta (f) / (2 * self.alpha (f))
    # end def resonator_q_approx

    def reactance_q (self, f = None, l = None):
        """ This is from Terman [8] p. 193 (73)
        """
        if f is None:
            f = self.f
        if l is None:
            l = self.l
        phi = 2 * np.pi * (l / self.lamda (f))
        return self.resonator_q (f) * np.sin (2 * phi) / phi
    # end def reactance_q

    def set_loss_constants (self, f0, a0r, a0g = 0):
        """ Alternative to fit if we have less data, allows to specify
            a0r (and optionally a0g) for a single frequency.
        """
        self.f0  = f0
        self.a0r = a0r
        self.a0g = a0g
        self.g   = 1.0
    # end def set_loss_constants

    def set_freq_params (self, f, l, p, z_l = None, z_i = None):
        self.f = f
        self.l = l
        self.p = p
        assert z_l or z_i
        assert not (z_l and z_i)
        if z_l is not None:
            self.z_l = complex (z_l)
            self.z_i = self.z_d (f, l, z_l)
        elif z_i is not None:
            self.z_i = complex (z_i)
            self.z_l = self.z_load (l, f, z_i)
    # end def set_freq_params

    def freq (self, l):
        """ Frequency at given length lambda
        """
        return c * self.vf / l
    # end def freq

    @property
    def fx (self):
        """ Frequency where the impedance is real.
            This is the point where loss_g and loss_r are equal.
        """
        return self.f0 * (self.a0r / self.a0g) ** (1.0 / (self.g - .5))
    # end def fx

    @property
    def L (self):
        return self.Cpl * self.Z0 ** 2
    # end def L

    @property
    def P_l (self):
        return abs (self.U_l) ** 2 * (1.0 / self.z_l).real
    # end def P_l

    @property
    def combined_loss (self):
        return 10 * np.log (self.p / self.P_l) / np.log (10)
    # end def combined_loss

    @property
    def U_i (self):
        return np.sqrt (self.p / (1.0 / self.z_i).real)
    # end def U_i

    @property
    def U_l (self):
        u = self.U_i
        z0 = self.z0f   ()
        gm = self.gamma ()
        zi = self.z_i
        return u * (np.cosh (gm * self.l) - z0 / zi * np.sinh (gm * self.l))
    # end def U_l

    @property
    def U_max (self):
        """ Maximum voltage given p
        """
        return np.sqrt (self.p * self.Z0) * np.sqrt (self.vswr_i)
    # end def U_max

    @property
    def I_max (self):
        """ Maximum current given p
        """
        return np.sqrt (self.p * self.vswr_i / self.Z0)
    # end def I_max

    @property
    def rho_l (self):
        return (self.z_l - self.Z0) / (self.z_l + self.Z0)
    # end def rho_l

    @property
    def vswr_l (self):
        return (1 + abs (self.rho_l)) / (1 - abs (self.rho_l))
    # end def vswr_l

    @property
    def rho_i (self):
        return (self.z_i - self.Z0) / (self.z_i + self.Z0)
    # end def rho_i

    @property
    def vswr_i (self):
        return (1 + abs (self.rho_i)) / (1 - abs (self.rho_i))
    # end def vswr_i

    def lamda (self, f = None):
        if f is None:
            f = self.f
        return 2 * np.pi / self.beta (f)
    # end def lamda

    def phi (self, f = None, l = None):
        """ Electrical length at given frequency (in rad)
        """
        if f is None:
            f = self.f
        if l is None:
            l = self.l
        return self.beta (f) * l
    # end def phi

    def conductance (self, f = None, z0 = None):
        if f is None:
            f = self.f
        ln10 = np.log (10)
        if z0 is None:
            z0 = self.Z0
        g = self.loss_g (f) * ln10 * z0.real / (1000.0 * abs (z0) ** 2)
        return g
    # end def conductance

    def resistance (self, f = None, z0 = None):
        if f is None:
            f = self.f
        ln10 = np.log (10)
        if z0 is None:
            z0 = self.Z0
        r = self.loss_r (f) * ln10 * z0.real / 1000.0
        return r
    # end def resistance

    def summary (self, f, l, metric = True):
        unit, units, cv = self._units (metric)
        ohm = '\u3a09'
        ohm = '\u2126'
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

    def summary_loss_data (self, metric = True):
        """ Return summary of loss parameters (from match and from
            manufacturer data)
        """
        r = []
        unit, units, cv = self._units (metric)
        r.append ("f (MHz)   Manu    Fit   Diff   (in dB/100 %s)" % units)
        for x, y in self.loss_data:
            y = y * cv
            a = self.loss (x) * cv
            r.append ("   %4.0f %6.2f %6.2f %6.2f" % (x / 1e6, y, a, (a - y)))
        return '\n'.join (r)
    # end def summary_loss_data

    def summary_loss (self, metric = True):
        r = []
        unit, units, cv = self._units (metric)
        r.append ("f (MHz)  Loss  (in dB/100 %s)" % units)
        frq = [ 10e6, 20e6, 50e6, 100e6, 200e6, 300e6, 500e6, 800e6, 1e9
              , 2.4e9, 5e9, 10e9, 20e9
              ]
        for f in frq:
            a = self.loss (f) * cv
            r.append ("  %5.0f %8.2f" % (f / 1e6, a))
        return '\n'.join (r)
    # end def summary_loss

    def summary_match \
        ( self
        , f      = None
        , l      = None
        , p      = None
        , z_l    = None
        , z_i    = None
        , metric = True
        ):
        """ Return summary of matching parameters, this asumes that
            self.f, self.l, self.p and either self.z_l or self.z_i have
            been set.
        """
        if f is None:
            f = self.f
        if l is None:
            l = self.l
        if p is None:
            p = self.p
        z_l = self.z_l
        z_i = self.z_i
        unit, units, cv = self._units (metric)
        ohm = '\u2126'
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
        r.append ('%25s %.3f' % ('abs(rho) at load', abs (self.rho_l)))
        r.append ('%25s %.3f' % ('VSWR at load', self.vswr_l))
        r.append ('%25s %.3f' % ('abs(rho) at input', abs (self.rho_i)))
        r.append ('%25s %.3f' % ('VSWR at input', self.vswr_i))
        r.append ('%25s %.2f V RMS' % ('Maximum Voltage', self.U_max))
        r.append ('%25s %.2f A RMS' % ('Maximum Current', self.I_max))
        # Too low
        #r.append \
        #    ( '%25s %.2f W per %s'
        #    % ('Maximum Power Stress', (self.p - self.P_l) / (l * cv), unit)
        #    )
        # Thats not per foot
        #r.append \
        #    ( '%25s %.2f W per %s'
        #    % ('Maximum Power Stress'
        #      , self.I_max ** 2 * self.resistance ()
        #      , unit
        #      )
        #    )
        for is_cap in True, False:
            for is_short in True, False:
                sd, sl = self.stub_match_iterative (is_cap, is_short)
                zi     = self.stub_impedance (f, sd, sl, self.z_l, is_short)
                type   = 'Inductive'
                circ   = 'open'
                if not is_cap:
                    type = 'Capacitive'
                if is_short:
                    circ = 'short'
                r.append ('%s stub with %s circuit at end:' % (type, circ))
                r.append \
                    ( '%25s %.5f %s from load'
                    % ('Stub attached', sd / cv, units)
                    )
                r.append ('%25s %.5f %s' % ('Stub length', sl / cv, units))
                r.append \
                    ( '%25s %.2f %+.2fj'
                    % ('Resulting impedance', zi.real, zi.imag)
                    )

        return '\n'.join (r)
    # end def summary_match

    def summary_resonator (self, f = None, metric = True):
        unit, units, cv = self._units (metric)
        if f is None:
            f = self.f
        r = []
        l = self.lamda (f)
        r.append ("Series Resonator at %.2f MHz" % (self.f / 1e6))
        r.append ("Open circuited quarter wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 4 / cv, units))
        r.append \
            ("%8s: %.3f \u2126" % ('abs (z)', abs (self.z_d_open (f, l/4))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        r.append ("Short circuited half wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 2 / cv, units))
        r.append \
            ("%8s: %.3f \u2126" % ('abs (z)', abs (self.z_d_short (f, l/2))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        r.append ("Parallel Resonator at %.2f MHz" % (self.f / 1e6))
        r.append ("Short circuited quarter wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 4 / cv, units))
        r.append \
            ("%8s: %.3f \u2126" % ('abs (z)', abs (self.z_d_short (f, l/4))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        r.append ("Open circuited half wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 2 / cv, units))
        r.append \
            ("%8s: %.3f \u2126" % ('abs (z)', abs (self.z_d_open (f, l/2))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        return '\n'.join (r)
    # end def summary_resonator

    def summary_stub (self, f = None, reactance = None, metric = True):
        if f is None:
            f = self.f
        if reactance is None:
            reactance = self.reactance
        z = reactance
        unit, units, cv = self._units (metric)
        self.f = f
        r = []
        if z.imag:
            z = z.imag
        if z < 0:
            r.append ('Capacitive impedance %.2f \u2126' % z)
            r.append ('Open-Circuited')
            l = self.stub_open (z)
            z = self.z_d_open (f, l)
            a = 'C'
            u = 'pF'
            v = 1 / (2 * np.pi * f * (-z.imag)) * 1e12
        else:
            r.append ('Inductive impedance %.2f \u2126' % z)
            r.append ('Closed-Circuited')
            l = self.stub_short (z)
            z = self.z_d_short (f, l)
            a = 'L'
            u = 'µH'
            v = z.imag / (2 * np.pi * f) * 1e6
        r.append ("%11s: %.5f %s" % ('Length', l / cv, units))
        r.append ("%11s: %.2f %+.2fj \u2126" % ('Zi', z.real, z.imag))
        r.append ("%11s: %.3f %s" % ('Effective %s' % a, v, u))
        r.append ("%11s: %.2f" % ('Q', self.reactance_q (f, l)))
        return '\n'.join (r)
    # end def summary_stub

    # The Y-Parameters (Admittance Parameters)
    # These are needed for simulating a lossy cable in NEC

    def y11 (self, f, l):
        """ Admittance Parameter Y11: This is simply the reciprocal
            value of z_d_short for the given f and l. Note that for an
            unterminated cable Y11 = Y22. See Wikipedia
            https://en.wikipedia.org/wiki/Admittance_parameters
        """
        return 1.0 / self.z_d_short (f, l)
    # end def y11

    def y22 (self, f, l, z_l = None):
        """ Admittance Parameter Y22: This simply adds 1/z_l in parallel
            to y11. Default value for z_l is open circuit.
        """
        if z_l is None:
            return self.y11 (f, l)
        elif z_l == 0:
            return 1e50
        return 1.0 / z_l + self.y11 (f, l)
    # end def y22

    def y12 (self, f, l):
        """ Admittance Parameter Y12: Since this is a reciprocal network
            Y12 = Y21.
            Y_in = Y11 - Y12 ** 2 / (Y22 + Y_L)
            Where Y_L is the load admittance and Y_in is the admittance to
            be measured at the input. We can compute Y_in from Y_L using
            self.z_d, from this we can compute Y12.
            Y12 = +/- sqrt (Y11 * (Y_L + Y22) - Y_in * (Y_L + Y22))
            We use the result with a positive real part.
            See Wikipedia https://en.wikipedia.org/wiki/Admittance_parameters
            Note that we set Y_L to 0 (open circuit) to compute the
            cable parameter y12 and we set y11 = y22.
        """
        y_in = 1 / self.z_d_open (f, l)
        y11  = self.y11 (f, l)
        r    = np.sqrt (y11 ** 2 - y_in * (y11))
        if r.real < 0:
            return -r
        return r
    # end def y12

    def z0f_witt (self, f, z0, r, g):
        """ From Witt [3]
        """
        cpart = 2j * np.pi * f * self.Cpl
        return np.sqrt ((r + cpart * z0 ** 2) / (g + cpart))
    # end def z0f_witt

    def z0f_sabin (self, f, z0, r, g):
        """ From Sabin [6] p.4:
        """
        l = self.Cpl * z0 ** 2
        o2 = 4 * np.pi * f
        return z0.real * (1 - 1j * (r / (o2 * l) - g / (o2 * self.Cpl)))
    # end def z0f_sabin

    def z0f (self, f = None, z0 = None):
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
        if f is None:
            f = self.f
        if z0 is None:
            z0 = self.Z0
        r = self.resistance  (f, z0)
        g = self.conductance (f, z0)
        if self.use_sabin:
            return self.z0f_sabin (f, z0, r, g)
        return self.z0f_witt  (f, z0, r, g)
    # end def z0f

    def plot_z0f (self, f):
        """ Since Witt [3] tells us that the z0_f formula above can be
            iterated, this plots the absolute value of the difference of
            an input z0 with an output z0. It clearly shows this diverges.
        """
        x = np.arange (-1.5, 1.5, 0.02)
        y = np.arange (49-0.5, 51+0.5, 0.02)
        z = []
        for vx in x:
            zz = []
            z.append (zz)
            for vy in y:
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

    def _stub_match_iter (self, d, y, goal = None):
        """ Find better approximation of d with a binary search
        """
        if goal is None:
            goal = 1.0 / self.Z0
        # Preconditions for binary search, current value is valid
        assert y.imag != 0
        v = abs ((1.0 / self.z_d (self.f, d, self.z_l)).real - y.real)
        assert v < eps ** 2
        dir = np.sign (y.real - goal)
        du  = d + self.lamda () / 50.
        dl  = d - self.lamda () / 50.
        zu  = (1.0 / self.z_d (self.f, du, self.z_l)).real
        zl  = (1.0 / self.z_d (self.f, dl, self.z_l)).real
        yr  = y.real
        if np.sign (zu - y.real) == dir:
            du = d
            zu = y.real
        else:
            assert np.sign (zl - y.real) == dir
            dl = d
            zl = y.real
        for k in range (100):
            if abs (yr - goal) / goal < 1e-3:
                break
            d  = (dl + du) / 2.
            y  = 1.0 / self.z_d (self.f, d, self.z_l)
            yr = y.real
            if yr > goal:
                du = d
                zu = yr
            else:
                dl = d
                zl = yr
        #print ("  corrected: d: %.3f y:%.6f %+.6f" % (d, y.real, y.imag))
        return d, y
    # end def _stub_match_iter

    def stub_match (self, capacitive = True):
        """ Distance *from load* of stub match point and impedance to be
            matched.
            We first compute the point from the voltage minimum and then
            add the distance from the load to the voltage minimum.
            In the example by Chipman [1] p. 178 (problem 8.3) the matching
            point nearer to the load is the shorter one with
            short-circuit termination. We don't count on this, we
            compute *both* succeptances and use the one with a positive
            imaginary part, thats the capacitive one (succeptance!).
            In fact if the load reactance is inside the unit reactance
            circle then the inductive matching point is nearer to the
            load.
            We then further optimize the matching point taking loss into
            account. There doesn't seem to be a closed-form method to
            compute this with loss directly.
            Note that the capacitive flag means the to-be-compensated
            impedance. So the default is to match with a short
            short-circuit stub.
            We return the admittance of which the imaginary part (the
            susceptance) has to be compensated with a stub.
        """
        d_v = self.lamda () * abs (np.arccos (abs (self.rho_l))) / (4 * np.pi)
        # Compute both points and determine which one is capacitive
        # unless capacitive is False in which case we return the other
        # one.
        dmin = self.d_voltage_min ()
        #print ("d_v:", d_v)
        #print ("wl d_v:", d_v / self.lamda ())
        #print ("dmin:", dmin)
        D = []
        for x in (dmin + d_v, dmin - d_v):
            if x / self.lamda () > 0.5:
                D.append (x - self.lamda () / 2.0)
            elif x < 0:
                D.append (x + self.lamda () / 2.0)
            else:
                D.append (x)
        assert len (D) == 2
        #print (D)
        Y = [1.0 / self.z_d (self.f, x, self.z_l) for x in D]
        #print (Y)
        d = D [0]
        for n, y in enumerate (Y):
            #print (y)
            if capacitive:
                if y.imag > 0:
                    d = D [n]
                    break
            else:
                if y.imag < 0:
                    d = D [n]
                    break
        #print ("uncorrected: d: %.3f y:%.6f %+.6f" % (d, y.real, y.imag))
        return self._stub_match_iter (d, y)
    # end def stub_match

    def _stub_match_step (self, d, shortcircuit):
        z = self.z_d (self.f, d, self.z_l)
        l = self.stub_short_open_iter (-z, shortcircuit = shortcircuit)
        y = 1.0 / self.stub_impedance \
            (self.f, d, l, self.z_l, shortcircuit = shortcircuit)
        return z, l, y
    # end def _stub_match_step

    def stub_match_iterative (self, capacitive = True, shortcircuit = True):
        """ Perform iterative matching: Compute the best match with
            y.real = 1/Z0, then compute the stub length and the
            resulting impedance. Iterate if the goal is not reached.
        """
        goal = (1.0 / self.Z0).real
        dmin = self.d_voltage_min ()
        d, ym = self.stub_match (capacitive = capacitive)
        l  = self.stub_short_open_iter (-1/ym, shortcircuit = shortcircuit)
        y  = 1.0 / self.stub_impedance \
            (self.f, d, l, self.z_l, shortcircuit = shortcircuit)
        dl  = d - self.lamda () / 30
        if d > dmin and dl < dmin:
            dl = dmin + eps
        zl, ll, yl = self._stub_match_step (dl, shortcircuit)
        du  = d + self.lamda () / 30
        if d < dmin and du > dmin:
            du = dmin - eps
        zu, lu, yu = self._stub_match_step (du, shortcircuit)
        # Both values on the same side of y?
        # Code currently not triggered by any test.
        # But it triggers when dmin-limit code above is disabled
        if np.sign (y.real - yl.real) == np.sign (y.real - yu.real):
            # Sort-of lim (yl.real -> y.real) and lim (yu.real -> y.real)
            dleps = d - eps
            zleps, lleps, yleps = self._stub_match_step (dleps, shortcircuit)
            ylr = yleps.real
            dueps = d + eps
            zueps, lueps, yueps = self._stub_match_step (dueps, shortcircuit)
            yur = yueps.real
            assert np.sign (y.real - ylr) != np.sign (y.real - yur)
            if np.sign (y.real - ylr) != np.sign (y.real - yl.real):
                s  = np.sign (y.real - yl.real)
                dn = dl
                vn = 'l'
            else:
                assert np.sign (y.real - yur) != np.sign (y.real - yu.real)
                s  = np.sign (y.real - yu.real)
                dn = du
                vn = 'u'
            for i in range (20):
                dn = (d + dn) / 2
                zn, ln, yn = self._stub_match_step (dn, shortcircuit)
                if np.sign (y.real - yn) != s:
                    break
            else:
                assert 0
            if vn == 'l':
                dl, zl, ll, yl = dn, zn, ln, yn
            else:
                du, zu, lu, yu = dn, zn, ln, yn

        if yl.real > y.real:
            yl, yu = yu, yl
            ll, lu = lu, ll
            dl, du = du, dl

        for i in range (200):
            if abs (y.real - goal) < eps:
                break
            if y.real > goal:
                yu = y
                lu = l
                du = d
            else:
                yl = y
                ll = l
                dl = d
            d = (dl + du) / 2.0
            z = self.z_d (self.f, d, self.z_l)
            l = self.stub_short_open_iter (-z, shortcircuit = shortcircuit)
            y = 1.0 / self.stub_impedance \
                (self.f, d, l, self.z_l, shortcircuit = shortcircuit)
        return d, l
    # end def stub_match_iterative

    def stub_short (self, z, f = None):
        """ Compute length of a short-circuited (closed) stub that has
            the given reactance as the imaginary part of a complex number.
            phi = 2*pi*d / lamda
            Formula from Chipman [1] p.131 (7.21)
            from Johnson [7] p.155 (6.21) and Terman [8] p.192 (72)
            z = j*z0 * tan (phi)
            j * tan (phi) = z / z0

        >>> cable = Manufacturer_Data_Cable (73, .66)
        >>> f = 200e6
        >>> cable.set_loss_constants (f, 10.6)
        >>> cable.set_freq_params (f, 1, 100, 73)
        >>> print ("%.5f" % (cable.lamda () / 2.0))
        0.49466

        # This should be 0
        >>> print ("%.5f" % cable.stub_short (0.0))
        0.00000

        # This should be lamda/8
        >>> print ("%.5f" % (cable.lamda () / 8))
        0.12366
        >>> print ("%.5f" % cable.stub_short (z = 73.0))
        0.12366

        # This should be 3*lamda/8
        >>> print ("%.5f" % (3 * cable.lamda () / 8))
        0.37099
        >>> print ("%.5f" % cable.stub_short (z = -73.0))
        0.37099

        >>> print ("%.5f" % (cable.lamda () * 0.4))
        0.39573
        >>> print ("%.5f" % (cable.lamda () * 0.5))
        0.49466
        >>> print ("%.5f" % cable.stub_short (z = -53.4))
        0.39522
        >>> print ("%.5f" % cable.stub_short (z = -53.4j))
        0.39522
        >>> print ("%.5f" % cable.stub_short (z = 53.4))
        0.09944
        >>> print ("%.5f" % cable.stub_short (z = 53.4j))
        0.09944
        """
        if f is None:
            f = self.f
        if z.imag:
            z = z.imag
        z /= self.Z0
        phi = np.arctan (z)
        if phi < 0:
            phi += np.pi
        return self.lamda (f) * phi / (2 * np.pi)
    # end def stub_short

    def stub_short_open_iter (self, z, f = None, shortcircuit = True):
        """ Iterative method of computing the matching impedance
            We search for the stub length where the imaginary part most
            closely matches the given admittance 1/(z*j).
        """
        stubmethod = self.stub_open
        z_d_method = self.z_d_open
        if shortcircuit:
            stubmethod = self.stub_short
            z_d_method = self.z_d_short
        if z.imag:
            goal = (1/z).imag
        else:
            goal = (1.0 / (z * 1j)).imag

        l = stubmethod (-1/goal, f)
        f = f or self.f
        assert 0 <= l <= self.lamda (f) / 2
        y   = 1.0 / z_d_method (f, l)
        ll  = l - self.lamda (f) / 50.0
        lu  = l + self.lamda (f) / 50.0
        if 0 <= l < self.lamda (f) / 4:
            if ll < 0:
                ll = eps
            if lu >= self.lamda (f) / 4.0:
                lu = self.lamda (f) / 4.0 - eps
        else:
            assert self.lamda (f) / 4.0 <= l <= self.lamda (f) / 2.0
            if ll <= self.lamda (f) / 4.0:
                ll = self.lamda (f) / 4.0 + eps
            if lu >= self.lamda (f) / 2.0:
                lu = self.lamda (f) / 2.0 - eps
        yl = 1.0 / z_d_method (f, ll)
        yu = 1.0 / z_d_method (f, lu)
        if yl.imag > goal:
            yl, yu = yu, yl
            ll, lu = lu, ll
        for i in range (200):
            if abs (y.imag - goal) < eps**2:
                break
            if y.imag > goal:
                lu = l
                yu = y
            else:
                ll = l
                yl = y
            l = (ll + lu) / 2.0
            y = 1.0 / z_d_method (f, l)
        return l
    # end def stub_short_open_iter

    def stub_short_iter (self, z, f = None):
        """ Note: The iterative variant tries to optimize the susceptance
            part of y = 1/z (y.imag) to be as exactly as possible.
            This will sometimes yield worse results if we're aiming for
            a close *reactance* (z.imag) match. The iterative version is
            mostly used for (parallel) stub matching where we want to
            eliminate the susceptance part as closely as possible.
        """
        return self.stub_short_open_iter (z, f, shortcircuit = True)
    # end def stub_short_iter

    def stub_open (self, z, f = None):
        """ Compute length of an open-circuit stub that has the given
            reactance as the imaginary part of a complex number.
            phi = 2*pi*d / lamda
            Formula from Chipman [1] p.131 (7.22)
            z = -j*z0 * cot phi
            from Johnson [7] p.160 (6.24), same from Terman [8] p.192 (72)
            z = -j*z0 / tan phi
            tan phi / -j = z0 / z
            j * tan phi = z0 / z

        >>> cable = Manufacturer_Data_Cable (73, .66)
        >>> f = 200e6
        >>> cable.set_loss_constants (f, 10.6)
        >>> cable.set_freq_params (f, 1, 100, 73)
        >>> print ("%.5f" % (cable.lamda () / 2.0))
        0.49466

        # This should be lamda/8
        >>> print ("%.5f" % (cable.lamda () / 8))
        0.12366
        >>> print ("%.5f" % cable.stub_open (z = -73.0))
        0.12366

        # This should be lamda/4
        >>> print ("%.5f" % (cable.lamda () / 4))
        0.24733
        >>> print ("%.5f" % cable.stub_open (z = 0.0))
        0.24733

        # This should be 3*lamda/8
        >>> print ("%.5f" % (3 * cable.lamda () / 8))
        0.37099
        >>> print ("%.5f" % cable.stub_open (z = 73.0))
        0.37099

        >>> print ("%.5f" % cable.stub_open (z = 73**2.0/53.4))
        0.39522
        >>> print ("%.5f" % cable.stub_open (z = -73**2.0/53.4j))
        0.39522
        >>> print ("%.5f" % cable.stub_open (z = -73**2.0/53.4))
        0.09944
        >>> print ("%.5f" % cable.stub_open (z = 73**2.0/53.4j))
        0.09944

        # [1] Example 7.5 p.132
        >>> cable = Manufacturer_Data_Cable (68.0336051416609, .66, 52.5e-12)
        >>> print ("%.5f" % cable.vf)
        0.93389
        >>> f = 500e6
        >>> cable.set_loss_constants (f, 10.6)
        >>> cable.set_freq_params (f, 1, 100, 73)
        >>> print ("%.5f" % cable.stub_open (z = -1 / 0.025))
        0.09262
        >>> print ("%.5f" % cable.stub_open (z = 1 / 0.025j))
        0.09262
        """
        if f is None:
            f = self.f
        if z.imag:
            z = z.imag
        z /= -self.Z0
        if z == 0:
            return self.lamda () / 4
        phi = np.arctan (1.0 / z)
        if phi < 0:
            phi += np.pi
        return self.lamda (f) * phi / (2 * np.pi)
    # end def stub_open

    def stub_open_iter (self, z, f = None):
        """ Note: The iterative variant tries to optimize the susceptance
            part of y = 1/z (y.imag) to be as exactly as possible.
            This will sometimes yield worse results if we're aiming for
            a close *reactance* (z.imag) match. The iterative version is
            mostly used for (parallel) stub matching where we want to
            eliminate the susceptance part as closely as possible.
        """
        return self.stub_short_open_iter (z, f, shortcircuit = False)
    # end def stub_open_iter

    def stub_impedance (self, f, stub_d, stub_l, z_l, shortcircuit = True):
        """ Compute stub impedance with distance from load stub_d and
            stub length stub_l and load impedance z_l. By default a
            closed (short circuit) stub is assumed.
        """
        #print (z_l)
        z_i = self.z_d (f, stub_d, z_l)
        #print (1/z_i)
        method = self.z_d_short
        if not shortcircuit:
            method = self.z_d_open
        z_s = method (f, stub_l)
        #print (1/z_s)
        return 1.0 / ((1.0 / z_i) + (1.0 / z_s))
    # end def stub_impedance

    def d_voltage_min (self, zd = None):
        """ Compute the (approximate) distance d from load where the
            impedance is real. We use vswr_l for this and postulate that
            the magnitude is the same as at the load (so we're assuming a
            lossless line here). At this point everything is real (no
            imaginary part). There are two solutions, one at voltage
            maximum, one at current maximum.  We return the one at
            current maximum (voltage minimum). Optionally we allow to
            compute the point of a given impedance.
        """
        z0  = self.Z0
        zl  = self.z_l
        zd  = zd or 1 / self.vswr_l * z0 # voltage minimum
        nom = (zl - z0) * zd + z0 * zl - z0 ** 2
        den = (zl + z0) * zd - z0 * zl - z0 ** 2
        ep  = np.sqrt (nom / den)
        a1  = np.angle (ep)
        a2  = np.angle (-ep)
        if a1 < 0:
            a1 += 2 * np.pi
        if a2 < 0:
            a2 += 2 * np.pi
        a = a1 if a1 < a2 else a2
        return a / 2 / np.pi * self.lamda ()
    # end def d_voltage_min

    def z_d (self, f, d, z_l):
        """ Compute impedance at distance d from load, given load
            impedance z_l, Eq 24 from [6] p.7 or better (exponential
            form) from 7.15 [1] p.130
            Special case z_l = None is open circuit.
            Note that we may return None for an open circuit if too near
            the load.
        """
        z0 = self.z0f   (f)
        gm = self.gamma (f)
        ep = np.e ** (gm * d)
        if abs (d) < 1e-20:
            return z_l
        if z_l is None:
            return z0 * (ep + 1/ep) / (ep - 1/ep)
        zz = z_l / z0
        return z0 * ( (ep * (zz + 1) + (zz - 1) / ep)
                    / (ep * (zz + 1) - (zz - 1) / ep)
                    )
    # end def z_d

    def z_d_open (self, f, d):
        """ Same as z_d but for open circuit (infinite z_l)
        """
        return self.z_d (f, d, z_l = None)
    # end def z_d_open

    def z_d_short (self, f, d):
        """ Same as z_d but for short circuit (z_l = 0)
        """
        return self.z_d (f, d, z_l = 0.0)
    # end def z_d_short

    def z_load (self, d, f = None, z_i = None):
        """ Compute impedance at load from input impedance given length
        """
        if f is None:
            f = self.f
        if z_i is None:
            z_i = self.z_i
        z0 = self.z0f (f)
        gm = self.gamma (f)
        ep = np.e ** (gm * d)
        return -((z0 * z_i - z0**2) * ep**2 + z0 * z_i + z0**2) \
                / ((z_i - z0) * ep**2 - z_i - z0)
    # end def z_load

# end class Manufacturer_Data_Cable

class Measured_Cable:
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

    def __init__ (self, f0, length, r_sc, r_oc, x_sc = 0, x_oc = 0):
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

def admittance (y11, y12, y22, z_l):
    """ Given the admittance matrix for a cable (note that y12 = y21)
        compute z_i from z_l.
    """
    if z_l == 0:
        return y11
    if z_l is None:
        y_l = 0
    else:
        y_l = 1.0 / z_l
    return y11 - (y12 ** 2) / (y22 + y_l)
# end def admittance

# Definining new cable parameters:
# We make a list of attenuation per 100m
# Then we call the constructor of Measured_Cable with
# - Characteristic impedance Z0 (typically 50.0 Ohm)
# - Velocity factor vf
# - Capacitance per meter Cpl
# All other constructor parameters are typically left at their defaults

belden_8295_data = \
    [ (1e6,    0.44 / m_per_ft)
    , (10e6,   1.4  / m_per_ft)
    , (50e6,   3.3  / m_per_ft)
    , (100e6,  4.9  / m_per_ft)
    , (200e6,  7.3  / m_per_ft)
    , (400e6, 11.5  / m_per_ft)
    , (700e6, 17.0  / m_per_ft)
    , (900e6, 20.0  / m_per_ft)
    , (1e9,   21.5  / m_per_ft)
    ]

belden_8295 = Manufacturer_Data_Cable \
    (50, .66, 30.8e-12 / m_per_ft, name = 'belden_8295')
belden_8295.fit (belden_8295_data)

sytronic_RG_213_UBX_data = \
    [ ( 10e6,  2.0)
    , ( 20e6,  3.0)
    , ( 50e6,  4.8)
    , (100e6,  7.8)
    , (200e6, 10.6)
    , (300e6, 13.4)
    , (500e6, 17.2)
    , (800e6, 24.0)
    , (  1e9, 27.5)
    ]
sytronic_RG_213_UBX = Manufacturer_Data_Cable \
    (50, .66, 103e-12, name = 'SYTRONIC RG 213 UBX (also RG 8/U)')
sytronic_RG_213_UBX.fit (sytronic_RG_213_UBX_data)

sytronic_RG_213_U_data = \
    [ ( 10e6,  1.8)
    , ( 20e6,  2.8)
    , ( 50e6,  4.4)
    , (100e6,  6.8)
    , (200e6,  9.7)
    , (300e6, 12.3)
    , (500e6, 16.2)
    , (800e6, 21.5)
    , (  1e9, 24.5)
    ]
sytronic_RG_213_U = Manufacturer_Data_Cable \
    (50, .66, 103e-12, name = 'SYTRONIC RG 213/U')
sytronic_RG_213_U.fit (sytronic_RG_213_U_data)

sytronic_RG_58_C_U_data = \
    [ ( 10e6,  4.7)
    , ( 20e6,  7.2)
    , ( 50e6, 10.7)
    , (100e6, 15.3)
    , (200e6, 22.8)
    , (300e6, 28.3)
    , (500e6, 37.0)
    , (800e6, 48.8)
    , (  1e9, 55.5)
    ]
sytronic_RG_58_C_U = Manufacturer_Data_Cable \
    (50, .66, 103e-12, name = 'SYTRONIC RG 58 C/U')
sytronic_RG_58_C_U.fit (sytronic_RG_58_C_U_data)

sytronic_RG_174_A_U_data = \
    [ (10e6,    9.6)
    , (20e6,   13.7)
    , (50e6,   21.8)
    , (100e6,  31.1)
    , (200e6,  44.5)
    , (300e6,  50.3)
    , (500e6,  72.7)
    , (800e6,  91.3)
    , (1e9,   106.1)
    ]
sytronic_RG_174_A_U = Manufacturer_Data_Cable \
    (50, 0.66, 103e-12, name = 'SYTRONIC RG 174 A/U (Reichelt RG 174-50)')
sytronic_RG_174_A_U.fit (sytronic_RG_174_A_U_data)

rs_222_8610_RG174A_U_data = \
    [ (200e6, 42.0)
    , (400e6, 60.0)
    , (3e9,  220.0)
    ]
rs_222_8610_RG174A_U = Manufacturer_Data_Cable \
    (50, 0.659, 106e-12, name = 'RS 222-8610 RG174A/U')
rs_222_8610_RG174A_U.fit (rs_222_8610_RG174A_U_data)

sytronic_RG_316_B_U_data = \
    [ (50e6,   19.2)
    , (100e6,  28.7)
    , (1e9,   104.8)
    , (3e9,   209.2)
    ]
sytronic_RG_316_B_U = Manufacturer_Data_Cable \
    (50, 0.7, 91e-12, name = 'SYTRONIC RG 316 B/U')
sytronic_RG_316_B_U.fit (sytronic_RG_316_B_U_data)

sytronic_RG_178_B_U_data = \
    [ (50e6,   38.0)
    , (100e6,  52.5)
    , (200e6,  65.3)
    , (300e6,  81.0)
    , (500e6, 120.7)
    , (1e9,   170.0)
    , (3e9,   308.0)
    ]
sytronic_RG_178_B_U = Manufacturer_Data_Cable \
    (50, 0.7, 94e-12, name = 'SYTRONIC RG 178 B/U')
sytronic_RG_178_B_U.fit (sytronic_RG_178_B_U_data)

# This does not work, curve-fit is not possible with this data.
# Note that the low-frequency components on the data sheet that have
# presumably been inserted by someone else (Seller?) don't make it
# better.
#bkl_RG_178_B_U_data = \
#    [ (1e9,   171.0)
#    , (1.8e9, 269.0)
#    , (2.4e9, 312.0)
#    , (3e9,   358.0)
#    , (5.2e9, 472.0)
#    , (6e9,   518.0)
#    ]
#bkl_RG_178_B_U = Manufacturer_Data_Cable \
#    (50, 0.7, 94e-12, name = 'BKL Electronic 1511006 RG 178 B/U')
#bkl_RG_178_B_U.fit (bkl_RG_178_B_U_data)

sytronic_RG_179_B_U_data = \
    [ (1e6,     3.0)
    , (5e6,    10.0)
    , (10e6,   12.0)
    , (50e6,   15.0)
    , (100e6,  21.0)
    , (300e6,  41.0)
    , (500e6,  58.0)
    , (800e6,  78.0)
    , (1e9,    90.0)
    ]
sytronic_RG_179_B_U = Manufacturer_Data_Cable \
    (75, 0.7, 102e-12, name = 'SYTRONIC RG 179 B/U')
sytronic_RG_179_B_U.fit (sytronic_RG_179_B_U_data)

# These are from Funkamateur Taschenkalender 2023
Airborne5 = Manufacturer_Data_Cable \
    (50, 0.85, name = 'Airborne5')
Airborne5_data = \
    [ (  10e6,  2.9)
    , (  14e6,  3.8)
    , (  28e6,  5.4)
    , (  50e6,  7.0)
    , ( 100e6,  9.4)
    , ( 144e6, 11)
    , ( 200e6, 12.9)
    , ( 432e6, 19)
    , (1296e6, 34.5)
    ]
Airborne5.fit (Airborne5_data)
Airborne10 = Manufacturer_Data_Cable \
    (50, 0.87, name = 'Airborne10')
Airborne10_data = \
    [ (  10e6,  1.2)
    , (  14e6,  1.4)
    , (  28e6,  1.9)
    , (  50e6,  2.4)
    , ( 100e6,  3.5)
    , ( 144e6,  4.2)
    , ( 200e6,  5.0)
    , ( 432e6,  7.6)
    , (1296e6, 13.6)
    ]
Airborne10.fit (Airborne10_data)
Aircell5 = Manufacturer_Data_Cable \
    (50, 0.82, name = 'Aircell5')
Aircell5_data = \
    [ (  10e6,  2.93)
    , (  28e6,  4.9)
    , (  50e6,  6.61)
    , ( 100e6,  9.4)
    , ( 144e6, 11.3)
    , ( 200e6, 13.4)
    , ( 432e6, 19.9)
    , (1296e6, 35.7)
    ]
Aircell5.fit (Aircell5_data)
Aircell7 = Manufacturer_Data_Cable \
    (50, 0.83, name = 'Aircell7')
Aircell7_data = \
    [ (  10e6,  2.2)
    , (  28e6,  3.7)
    , (  50e6,  4.5)
    , ( 100e6,  6.3)
    , ( 144e6,  7.6)
    , ( 200e6,  9.04)
    , ( 432e6, 13.6)
    , (1296e6, 24.8)
    ]
Aircell7.fit (Aircell7_data)
Aircom = Manufacturer_Data_Cable \
    (50, 0.85, name = 'Aircom Prem.')
Aircom_data = \
    [ (  10e6,  1.1)
    , (  50e6,  2.7)
    , ( 100e6,  3.3)
    , ( 144e6,  4.5)
    , ( 432e6,  8.5)
    , (1296e6, 12.5)
    ]
Aircom.fit (Aircom_data)
Ecoflex10 = Manufacturer_Data_Cable \
    (50, 0.85, name = 'Ecoflex10')
Ecoflex10_data = \
    [ (  10e6,  1.2)
    , (  14e6,  1.6)
    , (  28e6,  2.1)
    , (  50e6,  2.8)
    , ( 100e6,  4.0)
    , ( 144e6,  4.9)
    , ( 200e6,  5.8)
    , ( 432e6,  8.9)
    , (1296e6, 16.5)
    ]
Ecoflex10.fit (Ecoflex10_data)
Ecoflex10plus = Manufacturer_Data_Cable \
    (50, 0.85, name = 'Ecoflex10+')
Ecoflex10plus_data = \
    [ (  10e6,  1.3)
    , (  30e6,  2.3)
    , (  50e6,  2.9)
    , ( 100e6,  4.1)
    , ( 144e6,  5.0)
    , ( 432e6,  8.9)
    , (1296e6, 16.2)
    ]
Ecoflex10plus.fit (Ecoflex10plus_data)
Ecoflex15 = Manufacturer_Data_Cable \
    (50, 0.86, name = 'Ecoflex15')
Ecoflex15_data = \
    [ (  10e6,  0.86)
    , (  14e6,  1.0)
    , (  28e6,  1.4)
    , (  50e6,  1.96)
    , ( 100e6,  2.8)
    , ( 144e6,  3.4)
    , ( 200e6,  4.05)
    , ( 432e6,  6.1)
    , (1296e6, 11.4)
    ]
Ecoflex15.fit (Ecoflex15_data)
Ecoflex15plus = Manufacturer_Data_Cable \
    (50, 0.86, name = 'Ecoflex15+')
Ecoflex15plus_data = \
    [ (  10e6,  0.83)
    , (  50e6,  1.87)
    , ( 100e6,  2.67)
    , ( 144e6,  3.23)
    , ( 432e6,  5.8)
    , (1296e6, 10.5)
    ]
Ecoflex15plus.fit (Ecoflex15plus_data)
EcoflexMulti = Manufacturer_Data_Cable \
    (50, 0.85, name = 'Ecoflex Multi.')
EcoflexMulti_data = \
    [ (  10e6,  2.9)
    , ( 100e6,  9.4)
    , ( 500e6, 21.6)
    , (1000e6, 31.1)
    ]
EcoflexMulti.fit (EcoflexMulti_data)
H155 = Manufacturer_Data_Cable \
    (50, 0.81, name = 'H155')
H155_data = \
    [ (  10e6,  3.0)
    , (  50e6,  6.5)
    , ( 100e6,  9.3)
    , ( 144e6, 11.2)
    , ( 200e6, 13.2)
    , ( 432e6, 19.8)
    , (1296e6, 34.9)
    ]
H155.fit (H155_data)
H2000_Flex = Manufacturer_Data_Cable \
    (50, 0.85, name = 'H2000-Flex')
H2000_Flex_data = \
    [ (  10e6,  1.2)
    , (  14e6,  1.4)
    , (  28e6,  2.0)
    , (  50e6,  2.7)
    , ( 100e6,  3.9)
    , ( 144e6,  4.8)
    , ( 432e6,  8.5)
    , (1296e6, 15.7)
    ]
H2000_Flex.fit (H2000_Flex_data)
H2005 = Manufacturer_Data_Cable \
    (50, 0.85, name = 'H2005')
H2005_data = \
    [ (  10e6,  2.98)
    , (  14e6,  3.83)
    , (  28e6,  5.37)
    , (  50e6,  6.98)
    , ( 144e6, 11.0)
    , ( 432e6, 19.5)
    , (1296e6, 33.8)
    ]
H2005.fit (H2005_data)

coax_models  = dict \
    ( belden_8295          = belden_8295
    , sytronic_RG_174_A_U  = sytronic_RG_174_A_U
    , sytronic_RG_58_C_U   = sytronic_RG_58_C_U
    , sytronic_RG_213_U    = sytronic_RG_213_U
    , sytronic_RG_213_UBX  = sytronic_RG_213_UBX
    , sytronic_RG_178_B_U  = sytronic_RG_178_B_U
    , sytronic_RG_179_B_U  = sytronic_RG_179_B_U
    , sytronic_RG_316_B_U  = sytronic_RG_316_B_U
    , rs_222_8610_RG174A_U = rs_222_8610_RG174A_U
    )

def main ():
    cmd = ArgumentParser ()
    actions = ['loss', 'loss-data', 'match', 'resonator', 'stub']
    cmd.add_argument \
        ( 'action'
        , help = "Action to perform, one of %s" % ', '.join (actions)
        )
    cmd.add_argument \
        ( '-c', '--coaxmodel'
        , help    = "Coax to model, one of %s, default=%%(default)s"
                  % ', '.join (coax_models)
        , default = 'belden_8295'
        )
    cmd.add_argument \
        ( '-f', '--frequency'
        , type    = float
        , help    = "Frequency to match transmission line (Hz), "
                    "default=%(default)s"
        , default = 3.5e6
        )
    cmd.add_argument \
        ( '-I', '--imperial'
        , help    = "Use imperial length (in feet) instead of metric length"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '-l', '--length'
        , type    = float
        , help    = "Length of feed line (m), default=%(default)s"
        , default = 100 * m_per_ft
        )
    eo = " (either specify load *or* input impedance)"
    cmd.add_argument \
        ( '-z', '--z-load'
        , help    = "Impedance at load to be matched" + eo
                    + " Default=50 Ohm if no impedance given"
        , type    = complex
        )
    cmd.add_argument \
        ( '-i', '--z-input'
        , help    = "Impedance at input to be matched" + eo
        , type    = complex
        )
    cmd.add_argument \
        ( '-p', '--power'
        , help    = "Power applied to the cable, default=%(default)s"
        , type    = float
        , default = 100.0
        )
    cmd.add_argument \
        ( '-x', '--reactance'
        , help    = "Reactance for stub report, default=%(default)s"
        , type    = float
        , default = -100.0
        )
    args = cmd.parse_args ()
    cable = coax_models [args.coaxmodel]
    cable.reactance = args.reactance
    if args.z_load is None and args.z_input is None:
        args.z_load = 50.0
    cable.set_freq_params \
        (args.frequency, args.length, args.power, args.z_load, args.z_input)
    method = getattr (cable, 'summary_' + args.action.replace ('-', '_'))
    print (method (metric = not args.imperial))
# end def main

if __name__ == '__main__':
    main ()
