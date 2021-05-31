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
    >>> print (cable.summary_loss (metric = False))
    f (MHz)  Manu   Fit   Diff   (in dB/100 feet)
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
    >>> print (cable.summary_loss (metric = False))
    f (MHz)  Manu   Fit   Diff   (in dB/100 feet)
          1  0.43  0.43  -0.00
         10  1.39  1.39   0.00
         50  3.30  3.30   0.00
        100  4.89  4.89  -0.00
        200  7.39  7.39  -0.00
        400 11.46 11.46   0.00
        700 16.74 16.74  -0.00
        900 20.00 20.00  -0.00
       1000 21.58 21.58   0.00
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
     Characteristic Impedance 50.002 -0.436j Ohm
         Attenuation Constant 9.217e-3 nepers/m
               Phase Constant 0.914 rad/m
                 Resistance/m 0.8593 Ohm
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
     Characteristic Impedance 50.005 -0.642j Ohm
         Attenuation Constant 6.273e-3 nepers/m
               Phase Constant 0.444 rad/m
                 Resistance/m 0.5991 Ohm
                 Inductance/m 0.253 μH
                Conductance/m 11.304 μS
                Capacitance/m 101.050 pF
                 Matched Loss 1.661 dB

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
              Maximum Voltage 621.73 V RMS
              Maximum Current 12.43 A RMS

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
    >>> zd = cable.z_d (cable.f, 30.48 - 2*cable.lamda (), 50 -500j)
    >>> print ("%.5f" % cable.d_voltage_min (zd))
    2.20464
    >>> d, z = cable.stub_match (capacitive = True)
    >>> print ("%.4f" % d)
    2.9438
    >>> print ("%.5f" % (d / cable.lamda ()))
    0.20823
    >>> print ("%.6f" % z)
    -8.519958
    >>> print ("%.5f" % cable.d_voltage_min ())
    3.31225

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
              Maximum Voltage 621.73 V RMS
              Maximum Current 12.43 A RMS

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

    >>> d, z = cable.stub_match (capacitive = True)
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
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    9.73
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, closed = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, closed = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    >>> d, z = cable.stub_match (capacitive = False)
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
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    33.10
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, closed = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, closed = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    # Example with 50 Ohm line terminated in succeptance 2 -2j Chart Siemens
    # that's 0.25+0.25j chart Ohm or 12.5+12.5j Ohm
    >>> z_l = 12.5 +12.5j
    >>> cable.set_freq_params (f, 1, 1, 12.5 +12.5j)

    >>> d, z = cable.stub_match (capacitive = False)
    >>> print ("%.4f" % d)
    2.6085
    >>> print ("%.4f" % (d / cable.lamda ()))
    0.0305
    >>> print ("%.3f" % (1.0 / (z / cable.Z0)))
    1.581
    >>> l_oc = cable.stub_open (-z)
    >>> print ("%.2f" % l_oc)
    13.73
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    35.14
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, closed = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, closed = True)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j

    >>> d, z = cable.stub_match (capacitive = True)
    >>> print ("%.4f" % d)
    33.1418
    >>> print ("%.4f" % (d / cable.lamda ()))
    0.3869
    >>> print ("%.3f" % (1.0 / (z / cable.Z0)))
    -1.581
    >>> l_oc = cable.stub_open (-z)
    >>> print ("%.2f" % l_oc)
    29.10
    >>> l_sc = cable.stub_short (-z)
    >>> print ("%.2f" % l_sc)
    7.69
    >>> z_m = cable.stub_impedance (f, d, l_oc, z_l, closed = False)
    >>> print ("%.3f %+.3fj" % (z_m.real, z_m.imag))
    50.000 +0.000j
    >>> z_m = cable.stub_impedance (f, d, l_sc, z_l, closed = True)
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
     abs (z): 0.895 Ohm
           Q: 47.114
    Short circuited half wave:
      Length: 15.461 feet
     abs (z): 1.789 Ohm
           Q: 47.114
    Parallel Resonator at 21.00 MHz
    Short circuited quarter wave:
      Length: 7.730 feet
     abs (z): 2795.210 Ohm
           Q: 47.114
    Open circuited half wave:
      Length: 15.461 feet
     abs (z): 1398.052 Ohm
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
    Inductive
    Closed-Circuited
         Length: 5.45 feet
             Zi: 4.07 +99.90j Ohm
    Effective L: 0.757 µH
              Q: 34.04
    >>> print (cable.summary_stub (21e6, -100, metric = False))
    Capacitive
    Open-Circuited
         Length: 2.28 feet
             Zi: 0.40 -100.00j Ohm
    Effective C: 75.785 pF
              Q: 81.29

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

    def _units (self, metric) :
        unit = units = 'm'
        cv   = 1.0
        if not metric :
            unit  = 'foot'
            units = 'feet'
            cv    = m_per_ft
        return unit, units, cv
    # end def _units

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

    def resonator_q (self, f = None) :
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
        if f is None :
            f = self.f
        return 2 * np.pi \
            / (1 - 10 ** (-self.loss (f) * c * self.vf / (1000 * f)))
    # end def resonator_q

    def resonator_q_approx (self, f = None) :
        """ Approximate formula from Chipman [1] p. 222 Eq 10.20
        >>> f = 3670830.9685955304
        >>> cable = Manufacturer_Data_Cable (50, 0.66)
        >>> cable.set_loss_constants (f, 0.345 / m_per_ft)
        >>> cable.set_freq_params (f, 1, 100, 50)
        >>> print ("%.1f" % cable.resonator_q_approx (f))
        44.7
        """
        if f is None :
            f = self.f
        return self.beta (f) / (2 * self.alpha (f))
    # end def resonator_q_approx

    def reactance_q (self, f = None, l = None) :
        """ This is from Terman [8] p. 193 (73)
        """
        if f is None :
            f = self.f
        if l is None :
            l = self.l
        phi = 2 * np.pi * (l / self.lamda (f))
        return self.resonator_q (f) * np.sin (2 * phi) / phi
    # end def reactance_q

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
            self.z_l = complex (z_l)
            self.z_i = self.z_d (f, l, z_l)
        elif z_i is not None :
            self.z_i = complex (z_i)
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

    @property
    def U_max (self) :
        """ Maximum voltage given p
        """
        return np.sqrt (self.p * self.Z0) * np.sqrt (self.vswr_i)
    # end def U_max

    @property
    def I_max (self) :
        """ Maximum current given p
        """
        return np.sqrt (self.p * self.vswr_i / self.Z0)
    # end def I_max

    @property
    def rho_l (self) :
        return (self.z_l - self.Z0) / (self.z_l + self.Z0)
    # end def rho_l

    @property
    def vswr_l (self) :
        return (1 + abs (self.rho_l)) / (1 - abs (self.rho_l))
    # end def vswr_l

    @property
    def rho_i (self) :
        return (self.z_i - self.Z0) / (self.z_i + self.Z0)
    # end def rho_i

    @property
    def vswr_i (self) :
        return (1 + abs (self.rho_i)) / (1 - abs (self.rho_i))
    # end def vswr_i

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

    def conductance (self, f = None, z0 = None) :
        if f is None :
            f = self.f
        ln10 = np.log (10)
        if z0 is None :
            z0 = self.Z0
        g = self.loss_g (f) * ln10 * z0.real / (1000.0 * abs (z0) ** 2)
        return g
    # end def conductance

    def resistance (self, f = None, z0 = None) :
        if f is None :
            f = self.f
        ln10 = np.log (10)
        if z0 is None :
            z0 = self.Z0
        r = self.loss_r (f) * ln10 * z0.real / 1000.0
        return r
    # end def resistance

    def summary (self, f, l, metric = True) :
        unit, units, cv = self._units (metric)
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

    def summary_loss (self, metric = True) :
        """ Return summary of loss parameters (from match and from
            manufacturer data)
        """
        r = []
        unit, units, cv = self._units (metric)
        r.append ("f (MHz)  Manu   Fit   Diff   (in dB/100 %s)" % units)
        for x, y in self.loss_data :
           y = y * cv
           a = self.loss (x) * cv
           r.append ("   %4.0f %5.2f %5.2f %6.2f" % (x / 1e6, y, a, (a - y)))
        return '\n'.join (r)
    # end def summary_loss

    def summary_match (self, f, l, p, z_l = None, z_i = None, metric = True) :
        """ Return summary of matching parameters given f, l, power p
            and load impedance z_l or input impedance z_i
        """
        self.set_freq_params (f, l, p, z_l, z_i)
        z_l = self.z_l
        z_i = self.z_i
        unit, units, cv = self._units (metric)
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
        return '\n'.join (r)
    # end def summary_match

    def summary_resonator (self, f = None, metric = True) :
        unit, units, cv = self._units (metric)
        if f is None :
            f = self.f
        r = []
        l = self.lamda (f)
        r.append ("Series Resonator at %.2f MHz" % (self.f / 1e6))
        r.append ("Open circuited quarter wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 4 / cv, units))
        r.append ("%8s: %.3f Ohm" % ('abs (z)', abs (self.z_d_open (f, l/4))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        r.append ("Short circuited half wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 2 / cv, units))
        r.append ("%8s: %.3f Ohm" % ('abs (z)', abs (self.z_d_short (f, l/2))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        r.append ("Parallel Resonator at %.2f MHz" % (self.f / 1e6))
        r.append ("Short circuited quarter wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 4 / cv, units))
        r.append ("%8s: %.3f Ohm" % ('abs (z)', abs (self.z_d_short (f, l/4))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        r.append ("Open circuited half wave:")
        r.append ("%8s: %.3f %s" % ('Length', l / 2 / cv, units))
        r.append ("%8s: %.3f Ohm" % ('abs (z)', abs (self.z_d_open (f, l/2))))
        r.append ("%8s: %.3f" % ('Q', self.resonator_q (f)))
        return '\n'.join (r)
    # end def summary_resonator

    def summary_stub (self, f, z, metric = True) :
        unit, units, cv = self._units (metric)
        self.f = f
        r = []
        if z.imag :
            z = z.imag
        if z < 0 :
            r.append ('Capacitive')
            r.append ('Open-Circuited')
            l = self.stub_open (z)
            z = self.z_d_open (f, l)
            a = 'C'
            u = 'pF'
            v = 1 / (2 * np.pi * f * (-z.imag)) * 1e12
        else :
            r.append ('Inductive')
            r.append ('Closed-Circuited')
            l = self.stub_short (z)
            z = self.z_d_short (f, l)
            a = 'L'
            u = 'µH'
            v = z.imag / (2 * np.pi * f) * 1e6
        r.append ("%11s: %.2f %s" % ('Length', l / cv, units))
        r.append ("%11s: %.2f %+.2fj Ohm" % ('Zi', z.real, z.imag))
        r.append ("%11s: %.3f %s" % ('Effective %s' % a, v, u))
        r.append ("%11s: %.2f" % ('Q', self.reactance_q (f, l)))
        return '\n'.join (r)
    # end def summary_stub

    # The Y-Parameters (Admittance Parameters)
    # These are needed for simulating a lossy cable in NEC

    def y11 (self, f, l) :
        """ Admittance Parameter Y11: This is simply the reciprocal
            value of z_d_short for the given f and l. Note that for an
            unterminated cable Y11 = Y22. See Wikipedia
            https://en.wikipedia.org/wiki/Admittance_parameters
        """
        return 1.0 / self.z_d_short (f, l)
    # end def y11

    def y22 (self, f, l, z_l = None) :
        """ Admittance Parameter Y22: This simply adds 1/z_l in parallel
            to y11. Default value for z_l is open circuit.
        """
        if z_l is None :
            return self.y11 (f, l)
        elif z_l == 0 :
            return 1e50
        return 1.0 / z_l + self.y11 (f, l)
    # end def y22

    def y12 (self, f, l) :
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
        if r.real < 0 :
            return -r
        return r
    # end def y12

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

    def stub_match (self, capacitive = True) :
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
            Note that the capacitive flag means to to-be-compensated
            impedance. So the default is to match with a short closed
            stub.
        """
        d_v = self.lamda () * abs (np.arccos (abs (self.rho_l))) / (4 * np.pi)
        # Compute both points and determine which one is capacitive
        # unless capacitive is False in which case we return the other
        # one.
        dmin = self.d_voltage_min ()
        #print ("d_v:", d_v)
        #print ("wl d_v:", d_v / self.lamda ())
        #print ("dmin:", dmin)
        d = []
        for x in (dmin + d_v, dmin - d_v) :
            if x / self.lamda () > 0.5 :
                d.append (x - self.lamda () / 2.0)
            elif x < 0 :
                d.append (x + self.lamda () / 2.0)
            else :
                d.append (x)
        assert len (d) == 2
        #print (d)
        Y = [1.0 / self.z_d (self.f, x, self.z_l) for x in d]
        #print (Y)
        r = d [0]
        for n, y in enumerate (Y) :
            #print (y)
            if capacitive :
                if y.imag > 0 :
                    r = d [n]
                    break
            else :
                if y.imag < 0 :
                    r = d [n]
                    break
        #print ("uncorrected: r: %.3f y:%.6f %+.6f" % (r, y.real, y.imag))
        goal = 1.0 / self.Z0
        # Find better approximation with a binary search
        dir = np.sign (y.real - goal)
        u   = r + self.lamda () / 50.
        l   = r - self.lamda () / 50.
        gu  = (1.0 / self.z_d (self.f, u, self.z_l)).real
        gl  = (1.0 / self.z_d (self.f, l, self.z_l)).real
        g   = y
        gr  = y.real
        if np.sign (gu - y.real) == dir :
            u  = r
            gu = y.real
        else :
            assert np.sign (gl - y.real) == dir
            l  = r
            gl = y.real
        for k in range (10) :
            if abs (gr.real - goal) / goal < 1e-3 :
                break
            r  = (l + u) / 2.
            g  = 1.0 / self.z_d (self.f, r, self.z_l)
            gr = g.real
            if gr > goal :
                u  = r
                gu = gr
            else :
                l  = r
                gl = gr
        #print ("  corrected: r: %.3f y:%.6f %+.6f" % (r, g.real, g.imag))
        return r, (1 / (g.imag * 1j)).imag
    # end def stub_match

    def stub_short (self, z, f = None) :
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
        if f is None :
            f = self.f
        if z.imag :
            z = z.imag
        z /= self.Z0
        phi = np.arctan (z)
        if phi < 0 :
            phi += np.pi
        return self.lamda (f) * phi / (2 * np.pi)
    # end def stub_short

    def stub_open (self, z, f = None) :
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
        if f is None :
            f = self.f
        if z.imag :
            z = z.imag
        z /= -self.Z0
        if z == 0 :
            return self.lamda () / 4
        phi = np.arctan (1.0 / z)
        if phi < 0 :
            phi += np.pi
        return self.lamda (f) * phi / (2 * np.pi)
    # end def stub_open

    def stub_impedance (self, f, stub_d, stub_l, z_l, closed = True) :
        """ Compute stub impedance with distance from load stub_d and
            stub length stub_l and load impedance z_l. By default a
            closed (short circuit) stub is assumed.
        """
        #print (z_l)
        z_i = self.z_d (f, stub_d, z_l)
        #print (1/z_i)
        method = self.z_d_short
        if not closed :
            method = self.z_d_open
        z_s = method (f, stub_l)
        #print (1/z_s)
        return 1.0 / ((1.0 / z_i) + (1.0 / z_s))
    # end def stub_impedance

    def d_voltage_min (self, zd = None) :
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
        if a1 < 0 :
            a1 += 2 * np.pi
        if a2 < 0 :
            a2 += 2 * np.pi
        a = a1 if a1 < a2 else a2
        return a / 2 / np.pi * self.lamda ()
    # end def d_voltage_min

    def z_d (self, f, d, z_l) :
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
        if abs (d) < 1e-20 :
            return z_l
        if z_l is None :
            return z0 * (ep + 1/ep) / (ep - 1/ep)
        zz = z_l / z0
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

def admittance (y11, y12, y22, z_l) :
    """ Given the admittance matrix for a cable (note that y12 = y21)
        compute z_i from z_l.
    """
    if z_l == 0 :
        return y11
    if z_l is None :
        y_l = 0
    else :
        y_l = 1.0 / z_l
    return y11 - (y12 ** 2) / (y22 + y_l)
# end def admittance

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

belden_8295 = Manufacturer_Data_Cable (50, .66, 30.8e-12 / m_per_ft)
belden_8295.fit (belden_8295_data)

sytronic_RG213UBX_data = \
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
sytronic_RG213UBX = Manufacturer_Data_Cable (50, .66, 103e-12)
sytronic_RG213UBX.fit (sytronic_RG213UBX_data)

sytronic_RG_58_CU_data = \
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
sytronic_RG_58_CU = Manufacturer_Data_Cable (50, .66, 103e-12)
sytronic_RG_58_CU.fit (sytronic_RG_58_CU_data)
