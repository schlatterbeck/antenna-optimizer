#!/usr/bin/python3
from __future__ import print_function
from math       import pi, sqrt, acosh, sin, cos, log

c0        = 2.99792458e8
mu0       = 4e-7 * pi
eps0      = 1 / (mu0 * c0 ** 2)
# alternatively: z0 = mu0 * c0
z0        = sqrt (mu0 / eps0)
eps_r_air = 1.00054

def transmission_line_z (wire_dia, wire_dist, eps_r = eps_r_air):
    """ Impedance of transmission line
        https://hamwaves.com/zc.circular/en/
        Default eps_r is for air.
        The wire_dist is the distance (center to center) of the two
        wires. Both need to have the same dimension (e.g. mm or cm or in).
    """
    zc   = z0 / (pi * sqrt (eps_r)) * acosh (wire_dist / wire_dia)
    return zc
# end def transmission_line_z

def transmission_line_z_square (wire_dia, wire_dist):
    """ Impedance of transmission line with two square conductors
        https://www.owenduffy.net/calc/tstl.htm
        https://hamwaves.com/zc.square/en/
        For now no eps_r can be specified.
    """
    a=0.539774145266
    b=0.404050444546
    c=0.009504588299
    zc = log ((wire_dist / wire_dia - a) / b) / c
    return zc
# end transmission_line_z_square

class Z_Interpolation (object):
    """ Z-Interpolation by Hartwig Harm, DH2MIC
        http://dh2mic.darc.de/tlc/tlc.pdf
        Note that instead of the inner diameter Harm uses the distance
        from the center of the inner conductor to the outer conductor.
        These are more like radii and he calls them a and b. So
        wire_dist in the other formulae above is 2 * a.  This means for
        the two-wire cases in free space, the distance is 2 * a.

        So we have the following cases, the names are the names of the
        respective Z_Interpolation objects:
        - z_round: round wire in round conductor: a = radius of inner
          wall of outer conductor
        - z_rectangular: round wire in square conductor: a = distance
          from middle of inner conductor to nearest wall
        - z_square: round wire in rectangular conductor: a = distance to nearer
          wall, b = distance to far wall
        - z_u_shaped: round wire in U-shaped conductor: b = distance to
          bottom, a = distance from middle of inner conductor to a side
          wall, both side walls equal distance
        - z_l_shaped: round wire in L-shaped conductor: a = distance to
          nearer wall, b = distance to far wall
        - z_two_planes_equal: two parallel planes equally spaced from
          wire: a is the distance to wall
        - z_two_planes_unequal: two parallel planes unequally spaced
          from wire: a is the distance to the nearer wall, b is the
          distance to the far wall
        - z_round_single_wire_plane: single round wire above plane: a is
          the distance to the plane
        - z_square_single_wire_plane: single square wire above plane: a
          is the distance to the plane
        - z_two_wire_line_round: two round conductors: 2 * a is the
          distance between the centers of the conductors
        - z_two_wire_line_square: two square conductors: 2 * a is the
          distance between the centers of the conductors

        We instatiate with two parameters, kmin and kmax. If kmax = 0
        we interprete this as having only a single k factor without
        interpolation (e.g. for the two parallel wire case). If only a
        single k is in use, no second distance (b) is used.

        Note that for the 2-wire case (for the round version of which we
        have the exact version above) we need to take the resulting Z
        with factor 2, this is the class parameter mult. For other
        structures the resulting Z is used as-is. The z_interpolation is
        off by 0.07% (exactly the same) in all cases, not yet clear why.

    >>> zi  = Z_Interpolation ("two round conductors", mult = 2.0)
    >>> dia = 10.0
    >>> for i in (15, 30, 50, 65, 80):
    ...     t = transmission_line_z (10.0, 1.0 * i)
    ...     z = zi.z_interpolation  (10.0, 0.5 * i)
    ...     err = abs (t - z) / t * 100
    ...     print ("%3.3f %7.3f %7.3f %7.5f%%" % (1.0 * i, t, z, err))
    15.000 115.380 115.460 0.06923%
    30.000 211.326 211.473 0.06923%
    50.000 274.827 275.018 0.06923%
    65.000 306.782 306.995 0.06923%
    80.000 331.920 332.149 0.06923%
    >>> for i in (15, 30, 50, 65, 80):
    ...     t = transmission_line_z   (10.0, 1.0 * i)
    ...     z = z_two_wire_line_round (10.0, 0.5 * i)
    ...     err = abs (t - z) / t * 100
    ...     print ("%3.3f %7.3f %7.3f %7.5f%%" % (1.0 * i, t, z, err))
    15.000 115.380 115.460 0.06923%
    30.000 211.326 211.473 0.06923%
    50.000 274.827 275.018 0.06923%
    65.000 306.782 306.995 0.06923%
    80.000 331.920 332.149 0.06923%
    >>> for x in range (3, 17):
    ...     v1 = transmission_line_z_square (10.0, 5.0 * x)
    ...     v2 = z_two_wire_line_square     (10.0, 2.5 * x)
    ...     err = abs (v1 - v2) / v1 * 100.0
    ...     print ("%3.3f %7.3f %7.3f %7.5f%%" % (5.0 * x, v1, v2, err))
    15.000  91.075  96.916 6.41324%
    20.000 135.178 137.222 1.51278%
    25.000 166.159 166.306 0.08858%
    30.000 190.063 189.360 0.36970%
    35.000 209.528 208.542 0.47066%
    40.000 225.949 225.000 0.41983%
    45.000 240.149 239.427 0.30062%
    50.000 252.658 252.277 0.15091%
    55.000 263.838 263.866 0.01081%
    60.000 273.942 274.422 0.17510%
    65.000 283.160 284.115 0.33709%
    70.000 291.636 293.077 0.49425%
    75.000 299.479 301.412 0.64539%
    80.000 306.778 309.202 0.79002%
    >>> def f (fun, a, b, c=0):
    ...     print ("%6.2f" % fun (a, b, c))
    >>> f (z_round,                    3.0, 10.0)
    113.80
    >>> f (z_rectangular,              3.0, 10.0, 14.0)
    125.36
    >>> f (z_rectangular,              3.0, 10.0, 23.0)
    127.84
    >>> f (z_rectangular,              3.0, 10.0, 10.0)
    118.38
    >>> f (z_rectangular,              3.0,  5.0, 14.0)
     86.09
    >>> f (z_square,                   3.0, 10.0)
    118.38
    >>> f (z_u_shaped,                 3.0, 10.0, 12.0)
    125.47
    >>> f (z_u_shaped,                 3.0, 10.0, 23.0)
    127.95
    >>> f (z_u_shaped,                 3.0,  5.0, 12.0)
     86.04
    >>> f (z_u_shaped,                 3.0,  5.0, 7.0)
     84.74
    >>> f (z_l_shaped,                 3.0, 10.0, 12.0)
    138.98
    >>> f (z_l_shaped,                 3.0,  5.0,  7.0)
    100.14
    >>> f (z_l_shaped,                 3.0, 10.0, 40.0)
    152.91
    >>> f (z_l_shaped,                 3.0,  6.0, 22.0)
    121.34
    >>> f (z_two_planes_equal,         3.0, 10.0)
    128.17
    >>> f (z_two_planes_unequal,       3.0, 10.0, 23.0)
    146.49
    >>> f (z_two_planes_unequal,       3.0, 10.0, 10.0)
    128.17
    >>> f (z_round_single_wire_plane,  3.0, 10.0)
    155.03
    >>> f (z_square_single_wire_plane, 3.0, 10.0)
    143.59
    >>> f (z_two_wire_line_round,      3.0, 10.0)
    310.07
    >>> f (z_two_wire_line_square,     3.0, 10.0)
    287.18
    """

    def __init__ (self, name, kmin = 2.0, kmax = 0, exp = 0, mult = 1.0):
        self.name  = name
        self.kmin  = kmin
        self.kmax  = kmax
        if kmax == 0:
            self.k     = kmin
            self.log2k = log (self.k) / log (2.0)
        self.mult  = mult
        self.exp   = exp
    # end def __init__

    def k_factor (self, a, b):
        """ Getting the object from a function is different for python3
            and python2. So we add a function for retrieving the object
            first.
        >>> def get_obj (f):
        ...     try:
        ...         return f.__self__
        ...     except AttributeError:
        ...         pass
        ...     return f.im_self
        >>> obj = get_obj (z_rectangular)
        >>> a = z_rectangular (3.0, 10.0, 14.0)
        >>> print ("%6.4f" % obj.k)
        1.2145
        >>> a = z_rectangular (3.0, 10.0, 23.0)
        >>> print ("%6.4f" % obj.k)
        1.2663
        >>> a = z_rectangular (3.0, 10.0, 10.0)
        >>> print ("%6.4f" % obj.k)
        1.0800
        >>> a = z_rectangular (3.0,  5.0, 14.0)
        >>> print ("%6.4f" % obj.k)
        1.2703
        >>> obj = get_obj (z_u_shaped)
	>>> a = z_u_shaped (3.0, 10.0, 12.0)
        >>> print ("%6.4f" % obj.k)
	1.2167
	>>> a = z_u_shaped (3.0, 10.0, 23.0)
        >>> print ("%6.4f" % obj.k)
        1.2686
	>>> a = z_u_shaped (3.0,  5.0, 12.0)
        >>> print ("%6.4f" % obj.k)
        1.2693
	>>> a = z_u_shaped (3.0,  5.0, 7.0)
        >>> print ("%6.4f" % obj.k)
        1.2412
        >>> obj = get_obj (z_l_shaped)
        >>> a = z_l_shaped (3.0, 10.0, 12.0)
        >>> print ("%6.4f" % obj.k)
        1.5270
        >>> a = z_l_shaped (3.0,  5.0,  7.0)
        >>> print ("%6.4f" % obj.k)
        1.6188
        >>> a = z_l_shaped (3.0, 10.0, 40.0)
        >>> print ("%6.4f" % obj.k)
        1.9299
        >>> a = z_l_shaped (3.0,  6.0, 22.0)
        >>> print ("%6.4f" % obj.k)
        1.9186
        >>> obj = get_obj (z_two_planes_unequal)
        >>> a = z_two_planes_unequal (3.0, 10.0, 23.0)
        >>> print ("%6.4f" % obj.k)
        1.7325
        >>> a = z_two_planes_unequal (3.0, 10.0, 10.0)
        >>> print ("%6.4f" % obj.k)
        1.2732
        """
        m = self.kmax - 1.0
        i = self.kmin - 1.0
        p = (m - i) / (m + i)
        self.p = p
        q = ( (1.0 - p * (1.0 * a / b) ** self.exp)
            / (1.0 + p * (1.0 * a / b) ** self.exp)
            )
        self.q = q
        self.k = 1.0 + m * q
        self.log2k = log (self.k) / log (2.0)
    # end def k_factor

    def z_interpolation (self, wire_dia, a, b = 0, eps_r = eps_r_air):
        if self.kmax:
            if a > b:
                a, b = b, a
            self.k_factor (a, b)
        factor = 2.0 * a / wire_dia
        zc = \
            ( 60.0 / sqrt (eps_r)
            * ( log (factor)
              + self.log2k * log (1.0 + sqrt (1.0 - (factor ** -2)))
              )
            )
        return zc * self.mult
    # end def z_interpolation
    
# end class Z_Interpolation

z_round                    = Z_Interpolation \
    ( "round wire in round outer conductor"
    , kmin = 1.0
    ).z_interpolation

z_rectangular              = Z_Interpolation \
    ( "round wire in rectangular outer conductor"
    , kmin = 1.08, kmax = 4.0 / pi, exp = 4.5
    ).z_interpolation

z_square                   = Z_Interpolation \
    ( "round wire in square outer conductor"
    , kmin = 1.08
    ).z_interpolation

z_u_shaped                 = Z_Interpolation \
    ( "round wire in U-shaped conductor"
    , kmin = 1.167754622, kmax = 4 / pi, exp = 4.0
    ).z_interpolation

z_l_shaped                 = Z_Interpolation \
    ( "round wire in L-shaped conductor"
    , kmin = 1.4, kmax = 2.0, exp = 1.78
    ).z_interpolation

z_two_planes_equal         = Z_Interpolation \
    ( "two parallel planes equally spaced from wire"
    , kmin = 4.0 / pi
    ).z_interpolation

z_two_planes_unequal       = Z_Interpolation \
    ( "two parallel planes unequally spaced from wire"
    , kmin = 4.0 / pi, kmax = 2.0, exp = 1.57
    ).z_interpolation

z_round_single_wire_plane  = Z_Interpolation \
    ( "single round wire above plane"
    ).z_interpolation

z_square_single_wire_plane = Z_Interpolation \
    ( "single square wire above plane"
    , kmin = 1.65
    ).z_interpolation

z_two_wire_line_round      = Z_Interpolation \
    ( "two round conductors"
    , mult = 2.0
    ).z_interpolation

z_two_wire_line_square     = Z_Interpolation \
    ( "two square conductors"
    , kmin = 1.65, mult = 2.0
    ).z_interpolation


def wire_L_from_Z (z, eps_r = eps_r_air):
    """ Compute L (per length) from Zo
    """
    return z * sqrt (eps_r) / c0
# end def wire_L_from_Z

def wire_C_from_Z (z, eps_r = eps_r_air):
    """ Compute C (per length) from Zo
    """
    return sqrt (eps_r) / (z * c0)
# end def wire_C_from_Z

def phase_shift (f, length, vf = 1):
    """ Compute phase shift in radians for a phasing stub with the given
        length and the given frequency f (Hz).
    >>> print ("%1.3f" % phase_shift (435000000, .35))
    3.191
    """
    lamb  = c0 / f
    lvf   = lamb * vf
    phase = (length % lvf) / lvf * 2 * pi
    return phase
# end def phase_shift

def complex_voltage (phi, u = 1):
    """ Given a phase phi (in rad) and a voltage (with only real part)
        return the complex voltage (real, imag)
    >>> print ("%2.3f %2.3f" % complex_voltage (pi))
    -1.000 0.000
    >>> print ("%2.3f %2.3f" % complex_voltage (2 * pi))
    1.000 -0.000
    >>> print ("%2.3f %2.3f" % complex_voltage (pi / 2))
    0.000 1.000
    >>> print ("%2.3f %2.3f" % complex_voltage (3 * pi / 2))
    -0.000 -1.000
    >>> print ("%2.3f %2.3f" % complex_voltage (pi / 6))
    0.866 0.500
    >>> print ("%2.3f %2.3f" % complex_voltage (pi + pi / 6))
    -0.866 -0.500
    """
    return (u * cos (phi), u * sin (phi))
# end def complex_voltage
