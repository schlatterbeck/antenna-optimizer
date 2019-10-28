#!/usr/bin/python
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from pga import PGA, PGA_STOP_TOOSIMILAR, PGA_STOP_MAXITER, \
                PGA_STOP_NOCHANGE, PGA_REPORT_STRING, PGA_POPREPL_BEST
from rsclib.autosuper import autosuper
from argparse import ArgumentParser

import math
import PyNEC

class Nec_File (object) :
    """ Implements the methods of geometry and nec context so we can
        output an nec file given the called methods
    """
    def __init__ (self) :
        self.repr = []
        self.repr.append ("CE")
    # end def __init__

    def arc (self, tag, segs, rad, a1, a2, r) :
        self.repr.append ("GA %d %d %g %g %g %g" % (tag, segs, rad, a1, a2, r))
    # end def arc

    def move (self, rox, roy, roz, xs, ys, zs, its, nrpt, itgi) :
        self.repr.append \
            ( "GM %d %d %g %g %g %g %g %g %d"
            % (itgi, nrpt, rox, roy, roz, xs, ys, zs, its)
            )
    # end def move

    def wire (self, tag, segs, x1, y1, z1, x2, y2, z2, r, rdel, rrad) :
        self.repr.append \
            ( "GW %d %d %g %g %g %g %g %g %g"
            % (tag, segs, x1, y1, z1, x2, y2, z2, r)
            )
        if rdel != 1 or rrad != 1 :
            self.repr.append ("GC 0 0 %g %g %g" % (rdel, rad, rrad * rad))
    # end def wire

    def geometry_complete (self, gnd) :
        self.repr.append ("GE %d" % gnd)
    # end def geometry_complete

    def set_extended_thin_wire_kernel (self, flag) :
        self.repr.append ("EK %d" % int (bool (flag)))
    # end def set_extended_thin_wire_kernel

    def ld_card (self, ldtyp, ldtag, ldtagf, ldtagt, zlr, zli, zlc) :
        self.repr.append \
            ( "LD %d %d %d %d %g %g %d"
            % (ldtyp, ldtag, ldtagf, ldtagt, zlr, zli, zlc)
            )
    # end def ld_card

    def ex_card (self, *val) :
        # Only the 10-value variant is implemented here
        assert 10 <= len (val) <= 10
        self.repr.append ( "EX %d %d %d %d %g %g %g %g %g %g" % val)
    # end def ex_card


    def fr_card (self, ifrq, nfrq, freq_hz, del_freq) :
        self.repr.append \
            ("FR %d %d 0 0 %g %g" % (ifrq, nfrq, freq_hz, del_freq))
    # end def fr_card

    def rp_card \
        ( self
        , calc_mode
        , n_theta
        , n_phi
        , output_format
        , normalization
        , D
        , A
        , theta0
        , phi0
        , delta_theta
        , delta_phi
        , radial_distance
        , gain_norm
        ) :
        self.repr.append \
            ( "RP %d %d %d %d %g %g %g %g %g %g"
            % ( calc_mode, n_theta, n_phi
              , output_format % 1000 + normalization * 100 + D * 10 + A
              , theta0, phi0, delta_theta, delta_phi, radial_distance, gain_norm
              )
            )
    # end def rp_card

    def __repr__ (self) :
        return '\n'.join (self.repr + ['EN'])
    # end def __repr__

# end class Nec_File

class Folded_Dipole (object) :

    wire_radius   = 1.5e-3 / 2.0
    dipole_radius = 0.010
    lambda_4      = .146
    segs_dipole   = 19
    segs_arc      = 17
    segs_boom     =  5
    reflector     = 0.2
    impedance     = 50.0
    frqidxmax     = 201
    frqinc        = 0.05
    frqstart      = 430 # MHz
    frqend        = 440
    phi_inc       = 5
    theta_inc     = 5
    theta_max     = 180 / theta_inc + 1
    phi_max       = 360 / phi_inc   + 1

    def __init__ \
        ( self
        , refl_dist     = .01
        , dipole_radius = dipole_radius
        , lambda_4      = lambda_4
        , reflector     = reflector
        , frqidxmax     = 201
        ) :
        self.refl_dist     = refl_dist
        self.dipole_radius = dipole_radius
        self.lambda_4      = lambda_4
        self.reflector     = reflector
        self.frqidxmax     = frqidxmax
        self.frqinc        = (self.frqend - self.frqstart) / (frqidxmax - 1.0)
        self.nec           = PyNEC.nec_context ()
        self.geometry   ()
        self.nec_params ()
    # end def __init__

    @property
    def frqidxrange (self) :
        return range (self.frqidxmax)
    # end def frqidxrange

    def as_nec (self) :
        n = Nec_File ()
        self.geometry (n)
        self.nec_params (n)
        self._compute (n)
        return repr (n)
    # end def as_nec

    def _compute (self, nec = None) :
        if nec is None :
            nec = self.nec
        nec.rp_card \
            (0, self.theta_max, self.phi_max, 0, 0, 0, 0, 0, 0
            , self.theta_inc, self.phi_inc, 0, 0
            )
    # end def _compute

    def compute (self, frqidx = None) :
        self._compute ()
        if frqidx is None :
            frqidx = self.frqidxmax // 2
        self.rp = self.nec.get_radiation_pattern (frqidx)
    # end def compute

    def geometry (self, geo = None) :
        self.tag = 1
        self.ex  = None
        if geo is None :
            geo = self.nec.get_geometry ()
        a = ((-90, 90), (90, 270))
        for n, z in enumerate ((1, -1)) :
            a1, a2 = a [n]
            geo.arc \
                ( self.tag
                , self.segs_arc
                , self.dipole_radius
                , a1
                , a2
                , self.wire_radius
                )
            geo.move (0, 0, 0, z * self.lambda_4, 0, 0, self.tag, 0, 0)
            self.tag += 1
        for x in (-self.lambda_4, self.lambda_4) :
            for z in (self.dipole_radius, -self.dipole_radius) :
                geo.wire \
                    ( self.tag
                    , self.segs_dipole
                    , 0, 0, z
                    , x, 0, z
                    , self.wire_radius
                    , 1, 1
                    )
                self.tag += 1
        self.ex = self.tag - 1
        # first part of boom across folded part
        geo.wire \
            ( self.tag
            , self.segs_boom
            , 0, 0,  self.dipole_radius
            , 0, 0, -self.dipole_radius
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # second part of boom from dipole to reflector
        geo.wire \
            ( self.tag
            , self.segs_boom
            , 0, 0, -self.dipole_radius
            , 0, 0, -(self.refl_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Reflector
        geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , self.reflector, 0, -(self.refl_dist + self.dipole_radius)
            , 0,              0, -(self.refl_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , -self.reflector, 0, -(self.refl_dist + self.dipole_radius)
            , 0,               0, -(self.refl_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Turn around Y by 270 deg,
        # move everything up by max (reflector length, lambda_4 + r)
        up = max (self.reflector, self.lambda_4 + self.dipole_radius)
        geo.move (0, 270, 0, 0, 0, up, 0, 0, 0)
    # end def geometry

    def nec_params (self, nec = None) :
        if nec is None :
            nec = self.nec
        nec.geometry_complete (0)
        nec.set_extended_thin_wire_kernel (True)
        nec.ld_card (5, 0, 0, 0, 37735849, 0, 0)
        nec.ex_card (0, self.ex, 1, 0, 1, 0, 0, 0, 0, 0)
        nec.fr_card (0, self.frqidxmax, self.frqstart, self.frqinc)
    # end def nec_params

    def max_f_r_gain (self) :
        """ Maximum forward and backward gain
        """
        gains = self.rp.get_gain ()
        n1max = n2max = -1
        gmax  = None
        for n1, ga in enumerate (gains) :
            for n2, g in enumerate (ga) :
                if gmax is None or g > gmax :
                    gmax = g
                    n1max = n1
                    n2max = n2
        # The other side of gmax: Find the maximum 30 deg around the
        # opposite side but with the same theta angle (also +- 30 deg)
        rmax = None
        phi2 = self.phi_max   // 2
        the2 = self.theta_max // 2
        p30  = 30 // self.phi_inc
        t30  = 30 // self.theta_inc
        tm   = None
        pm   = None
        for t in range (n1max - t30, n1max + t30 + 1) :
            theta = t
            if theta < 0 :
                theta += self.theta_max - 1
            if theta >= self.theta_max :
                theta -= self.theta_max - 1
            assert 0 <= theta < self.theta_max
            for p in range (n2max - p30 - phi2, n2max + p30 + 1 - phi2) :
                phi = p
                if phi < 0 :
                    phi += self.phi_max - 1
                if phi >= self.phi_max :
                    phi -= self.phi_max - 1
                assert 0 <= phi < self.phi_max
                if rmax is None or gains [theta][phi] > rmax :
                    rmax = gains [theta][phi]
                    pm   = phi
                    tm   = theta
        return gmax, rmax
    # end def max_f_r_gain

    def plot (self, frqidx = 100) :
        if not self.rp :
            self.compute (frqidx)

        # 0: linear, 1: right, 2: left
        print (self.rp.get_pol_sense_index ())
        print (self.rp.get_pol_tilt ())
        print (self.rp.get_pol_axial_ratio ())

        gains  = self.rp.get_gain ()
        gains  = 10.0 ** (gains / 10.0)
        # Thetas are upside down (count from top)
        thetas = self.rp.get_theta_angles () * np.pi / 180.0
        thetas = -thetas + np.pi
        phis   = self.rp.get_phi_angles ()   * np.pi / 180.0

        P, T = np.meshgrid (phis, thetas)

        X = np.cos (P) * np.sin (T) * gains
        Y = np.sin (P) * np.sin (T) * gains
        Z = np.cos (T) * gains

        fig = plt.figure ()
        ax  = fig.gca (projection='3d')

        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array \
            ( [ X.max () - X.min ()
              , Y.max () - Y.min ()
              , Z.max () - Z.min ()
              ]
            ).max() / 2.0

        mid_x = (X.max () + X.min ()) * 0.5
        mid_y = (Y.max () + Y.min ()) * 0.5
        mid_z = (Z.max () + Z.min ()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        ax.plot_wireframe (X, Y, Z, color = 'r')
        plt.show ()
    # end def plot

    def swr_plot (self) :
        fun   = self.nec.get_radiation_pattern
        frqs  = list (fun (i).get_frequency () for i in self.frqidxrange)
        vswrs = list (self.vswr (i) for i in self.frqidxrange)
        fig   = plt.figure ()
        ax    = fig.add_subplot (111)
        ax.plot (frqs, vswrs)
        plt.show ()
    # end def swr_plot

    def vswr (self, frqidx):
        ipt = self.nec.get_input_parameters (frqidx)
        z   = ipt.get_impedance ()
        gamma = np.abs ((z - self.impedance) / (z + self.impedance))
        return ((1. + gamma) / (1. - gamma)) [0]
    # end def vswr

# end class Folded_Dipole

class Dipole_Optimizer (PGA, autosuper) :
    """ Optimize given folded dipole
        Length are encoded as integers with a resolution of .5mm
        We use:                             min max
        * 8mm  <= dipole_radius <= 10cm  ->  16 200
        * 8mm  <= refl_dist     <= 10cm  ->  16 200
        * 10cm <= reflector     <= 40cm  -> 200 800
        * 10cm <= lambda_4      <= 20cm  -> 200 400
    """

    resolution = 0.5e-3 # 0.5 mm in meter

    def __init__ (self, srand = 42, verbose = False, random_seed = 42) :
        self.verbose     = verbose
        self.random_seed = random_seed
        stop_on = [PGA_STOP_NOCHANGE, PGA_STOP_MAXITER, PGA_STOP_TOOSIMILAR]
        PGA.__init__ \
            ( self
            , int
            , 4
            , init = [(16, 200), (16, 200), (200, 800), (200, 400)]
            , maximize            = False
            , pop_size            = 100
            , num_replace         = 50
            , random_seed         = self.random_seed
            , print_options       = [PGA_REPORT_STRING]
            , stopping_rule_types = stop_on
            )
    # end def __init__

    def to_meter (self, l) :
        return (l) * self.resolution
    # end def to_meter

    def evaluate (self, p, pop) :
        dipole_radius = self.to_meter (self.get_allele (p, pop, 0))
        refl_dist     = self.to_meter (self.get_allele (p, pop, 1))
        reflector     = self.to_meter (self.get_allele (p, pop, 2))
        lambda_4      = self.to_meter (self.get_allele (p, pop, 3))
        f = Folded_Dipole \
            ( dipole_radius = dipole_radius
            , refl_dist     = refl_dist
            , reflector     = reflector
            , lambda_4      = lambda_4
            , frqidxmax     = 3
            )
        f.compute ()
        if self.verbose :
            print \
                ( "r = %1.4f refd = %1.4f refl = %1.4f l/4 = %1.4f"
                % (dipole_radius, refl_dist, reflector, lambda_4)
                )
        vswrs = list (f.vswr (i) for i in f.frqidxrange)
        eval  = sum ([v, v ** 2][v > 1.8] for v in vswrs)
        eval *= 1 + sum (2 * bool (v > 1.8) for v in vswrs)
        if abs (vswrs [0] - vswrs [-1]) > 0.2 :
            eval *= 10
        gmax, rmax = f.max_f_r_gain ()

        if rmax < -10 :
            rmax = -10.0
        eval *= (10 - gmax) + (rmax / 2 + 5)
        assert eval > 0
        if self.verbose :
            print \
                ( "VSWR: %s\nGMAX: %s, RMAX: %s"
                % (vswrs, gmax, rmax)
                )
            print ("Eval: %3.2f" % eval)
        return eval
    # end def evaluate

    def print_string (self, file, p, pop) :
        verbose      = self.verbose
        self.verbose = True
        self.evaluate (p, pop)
        self.verbose = verbose
        return self.__super.print_string (file, p, pop)
    # end def print_string

# end class Dipole_Optimizer


if __name__ == '__main__' :
    actions = ['optimize', 'print', 'swr', 'gain']
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( 'action'
        , help = "Action to perform, one of %s" % actions
        )
    cmd.add_argument \
        ( '-4', '--lambda-len'
        , type = float
        , help = "(Half) Length of the dipole without rounded part"
        , default = 0.146
        )
    cmd.add_argument \
        ( '-d', '--reflector-distance'
        , type    = float
        , help    = "Distance of the reflector from nearest dipole part"
        , default = 0.01
        )
    cmd.add_argument \
        ( '-l', '--reflector-length'
        , type = float
        , help = "(Half) Length of the reflector"
        , default = 0.2
        )
    cmd.add_argument \
        ( '-R', '--random-seed'
        , type    = int
        , help    = "Random number seed for optimizer, default=%(default)s"
        , default = 42
        )
    cmd.add_argument \
        ( '-r', '--dipole-radius'
        , type    = float
        , help    = "Radius of the rounded corner of the folded dipole"
        , default = 0.01
        )
    cmd.add_argument \
        ( '-v', '--verbose'
        , help    = "Verbose reporting in every generation"
        , action  = 'store_true'
        )
    args = cmd.parse_args ()
    if args.action == 'optimize' :
        do = Dipole_Optimizer \
            (verbose = args.verbose, random_seed = args.random_seed)
        do.run ()
    else :
        f = Folded_Dipole \
            ( dipole_radius = args.dipole_radius
            , refl_dist     = args.reflector_distance
            , reflector     = args.reflector_length
            , lambda_4      = args.lambda_len
            )
        if args.action == 'print' :
            print (f.as_nec ())
        elif args.action == 'swr' :
            f.compute ()
            f.swr_plot ()
        elif args.action == 'gain' :
            f.compute ()
            f.plot ()
        elif args.action == 'frgain' :
            f.compute ()
            f, b = f.max_f_r_gain ()
            print ("forward: %2.2g backward: %2.2g" % (f, b))

