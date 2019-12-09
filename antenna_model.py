#!/usr/bin/python
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from pga import PGA, PGA_STOP_TOOSIMILAR, PGA_STOP_MAXITER, \
                PGA_STOP_NOCHANGE, PGA_REPORT_STRING, PGA_POPREPL_BEST, \
                PGA_NEWPOP
from rsclib.autosuper import autosuper
from argparse import ArgumentParser

import sys
import math
import PyNEC

class Nec_File (object) :
    """ Implements the methods of geometry and nec context so we can
        output an nec file given the called methods
    """
    def __init__ (self, comment) :
        self.repr = []
        if comment :
            if not isinstance (comment, type ([])) :
                comment = [comment]
            for c in comment :
                self.add_parameter_comment (c)
        self.end_comments ()
    # end def __init__

    def arc (self, tag, segs, rad, a1, a2, r) :
        self.repr.append ("GA %d %d %g %g %g %g" % (tag, segs, rad, a1, a2, r))
    # end def arc

    def add_parameter_comment (self, comment) :
        self.repr.append ("CM %s" % comment)
    # end def add_parameter_comment

    def end_comments (self) :
        self.repr.append ("CE")
    # end def end_comments

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

class Antenna_Model (autosuper) :

    wire_radius   = 1.5e-3 / 2.0
    boom_radius   = wire_radius
    impedance     = 50.0
    frqidxmax     = 201
    frqinc        = 0.05
    frqstart      = 430 # MHz
    frqend        = 440
    phi_inc       = 5
    theta_inc     = 5
    theta_max     = int (180 / theta_inc + 1)
    phi_max       = int (360 / phi_inc   + 1)

    def __init__ \
        ( self
        , frqidxmax     = 201
        , frqidxnec     = 201 # only for necout
        , wire_radius   = wire_radius
        , boom_radius   = boom_radius
        ) :
        self.wire_radius   = wire_radius
        self.boom_radius   = boom_radius
        self.frqidxmax     = frqidxmax
        self.frqidxnec     = frqidxnec
        self.frqinc        = (self.frqend - self.frqstart) / (frqidxmax - 1.0)
        self.frqincnec     = (self.frqend - self.frqstart) / (frqidxnec - 1.0)
        self.nec           = PyNEC.nec_context ()
        self.rp            = {}
        self.geometry   ()
        self.nec_params ()
    # end def __init__

    def as_nec (self) :
        c = [self.cmdline ()]
        if not self.rp :
            self._compute ()
        c.extend (self.show_gains ())
        n = Nec_File (c)
        self.geometry (n)
        self.nec_params (n)
        self._compute (n)
        return repr (n)
    # end def as_nec

    def cmdline (self) :
        """ This should be overridden in derived class to print out
            command-line parameters to regenerate this model.
        """
        return ''
    # end def cmdline

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
        if frqidx not in self.rp :
            self.rp [frqidx] = self.nec.get_radiation_pattern (frqidx)
    # end def compute

    def frqidxrange (self, step = 1) :
        return range (0, self.frqidxmax, step)
    # end def frqidxrange

    def geometry (self, geo = None) :
        raise NotImplemented ("Derived class must implement 'geometry'")
    # end def geometry

    def nec_params (self, nec = None) :
        frqidxmax = self.frqidxnec
        frqinc    = self.frqincnec
        if nec is None :
            nec       = self.nec
            frqidxmax = self.frqidxmax
            frqinc    = self.frqinc
        nec.geometry_complete (0)
        nec.set_extended_thin_wire_kernel (True)
        nec.ld_card (5, 0, 0, 0, 37735849, 0, 0)
        nec.ex_card (0, self.ex, self.exseg, 0, 1, 0, 0, 0, 0, 0)
        nec.fr_card (0, frqidxmax, self.frqstart, frqinc)
    # end def nec_params

    def max_f_r_gain (self, frqidx = None) :
        """ Maximum forward and backward gain
        """
        if frqidx is None :
            frqidx = self.frqidxmax // 2
        if frqidx not in self.rp :
            self.rp [frqidx] = self.nec.get_radiation_pattern (frqidx)
        gains = self.rp [frqidx].get_gain ()
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

    def show_gains (self, prefix = '') :
        r = []
        step = self.frqidxmax // 2
        for idx in self.frqidxrange (step) :
            f, b = self.max_f_r_gain (idx)
            frq = self.rp [idx].get_frequency ()
            frq = frq / 1e6
            r.append ("%sFRQ: %3.2f fw: %2.2f bw: %2.2f" % (prefix, frq, f, b))
        vswrs = list (self.vswr (i) for i in self.frqidxrange (step))
        r.append ("SWR: %1.2f %1.2f %1.2f" % tuple (vswrs))
        return r
    # end def show_gains

    def plot (self, frqidx = 100) :
        if not self.rp :
            self.compute (frqidx)
        if frqidx not in self.rp :
            self.rp [frqidx] = self.nec.get_radiation_pattern (frqidx)

        # 0: linear, 1: right, 2: left
        #print (self.rp [frqidx].get_pol_sense_index ())
        #print (self.rp [frqidx].get_pol_tilt ())
        #print (self.rp [frqidx].get_pol_axial_ratio ())

        gains  = self.rp [frqidx].get_gain ()
        gains  = 10.0 ** (gains / 10.0)
        # Thetas are upside down (count from top)
        thetas = self.rp [frqidx].get_theta_angles () * np.pi / 180.0
        thetas = -thetas + np.pi
        phis   = self.rp [frqidx].get_phi_angles ()   * np.pi / 180.0

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
        frqs  = list (fun (i).get_frequency () for i in self.frqidxrange ())
        vswrs = list (self.vswr (i) for i in self.frqidxrange ())
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

# end class Antenna_Model

class Antenna_Optimizer (PGA, autosuper) :
    """ Optimize given antenna, needs to be subclassed.
    """

    resolution = 0.5e-3 # 0.5 mm in meter

    def __init__ \
        (self, random_seed = 42, verbose = False, wire_radius = 0.00075) :
        self.verbose     = verbose
        self.random_seed = random_seed
        self.wire_radius = wire_radius
        stop_on = [PGA_STOP_NOCHANGE, PGA_STOP_MAXITER, PGA_STOP_TOOSIMILAR]
        # Determine number of bits needed from minmax,
        # we need at least self.resolution precision.
        self.nbits = []
        for l, u in self.minmax :
            n = (u - l) / self.resolution
            self.nbits.append (int (math.ceil (math.log (n) / math.log (2))))
        self.bitidx = []
        l   = 0
        for b in self.nbits :
            u = l + b - 1
            self.bitidx.append ((l, u))
            l = u + 1
        PGA.__init__ \
            ( self
            , bool
            , sum (self.nbits)
            , maximize            = True
            , pop_size            = 100
            , num_replace         = 50
            , random_seed         = self.random_seed
            , print_options       = [PGA_REPORT_STRING]
            , stopping_rule_types = stop_on
            , pop_replace_type    = PGA_POPREPL_BEST
            , print_frequency     = 10
            )
        self.cache = {}
        self.cache_hits = 0
        self.nohits     = 0
        self.file       = sys.stdout
    # end def __init__

    def get_parameter (self, p, pop, i) :
        return self.get_real_from_gray_code \
            (p, pop, *(self.bitidx [i] + self.minmax [i]))
    # end def get_parameter

    def set_parameter (self, p, pop, i, val) :
        self.encode_real_as_gray_code \
            (p, pop, *(self.bitidx [i] + self.minmax [i] + (val,)))
    # end def set_parameter

    def evaluate (self, p, pop) :
        cache_key = 0
        for k in range (len (self)) :
            cache_key <<= 1
            cache_key |= int (self.get_allele (p, pop, k))
        if cache_key in self.cache and not self.verbose :
            self.cache_hits += 1
            return self.cache [cache_key]
        else :
            self.nohits += 1

        antenna = self.compute_antenna (p, pop)
        antenna.compute ()
        if self.verbose :
            print (antenna.cmdline (), file = self.file)
        vswrs = list (antenna.vswr (i) for i in antenna.frqidxrange ())
        swr_eval  = sum (v for v in vswrs) / 3.0
        swr_eval *= 1 + sum (6 * bool (v > 1.8) for v in vswrs)
        diff = abs (vswrs [0] - vswrs [-1])
        if diff > 0.2 :
            swr_eval *= 1.0 + 15 * diff
        swr_eval **= (1./2)
        gmax, rmax = antenna.max_f_r_gain ()

        gmax = max (gmax, -20.0)
        rmax = max (rmax, -20.0)
        eval = ( 50.0
               + gmax ** 3.0
               - rmax * 4
               ) / swr_eval
        assert eval > 0
        if self.verbose :
            print \
                ( "VSWR: %s\nGMAX: %s, RMAX: %s"
                % (vswrs, gmax, rmax)
                , file = self.file
                )
            ch = self.cache_hits
            cn = self.nohits + self.cache_hits
            print \
                ( "Cache hits: %s/%s %2.2f%%"
                % (ch, cn, 100.0 * ch / cn)
                , file = self.file
                )
            print ("Eval: %3.2f" % eval, file = self.file)
        self.cache [cache_key] = eval
        return eval
    # end def evaluate

    def endofgen (self) :
        """ Simple hill-climb: Loop over all individuums, select one of
            the four alleles by random and try to inc/dec (randomly).
            If the inc/decremented gene is better, update allele and
            evaluation.
        """
        pop  = PGA_NEWPOP
        l    = len (self.nbits) # number of bits per float
        bidx = self.get_best_index (pop)
        best = self.get_evaluation (bidx, pop)
        calc = False
        for p in range (self.pop_size) :
            ev  = self.get_evaluation (p, pop)
            assert self.get_evaluation_up_to_date (p, pop)
            idx = self.random_interval (0, l - 1)
            val = self.get_parameter (p, pop, idx)
            if self.random_flip (0.5) :
                if val + self.resolution <= self.minmax [idx][1] :
                    self.set_parameter (p, pop, idx, val + self.resolution)
            else :
                if val - self.resolution >= self.minmax [idx][0] :
                    self.set_parameter (p, pop, idx, val - self.resolution)
            evnew = self.evaluate (p, pop)
            if evnew <= ev :
                # undo
                self.set_parameter (p, pop, idx, val)
            else :
                self.set_evaluation (p, pop, evnew)
                if evnew > best :
                    calc = True
        # Re-calculate fitness values if the best index changed
        if calc :
            self.fitness (pop)
    # end def endofgen

    def print_string (self, file, p, pop) :
        f            = self.file
        self.file    = file
        verbose      = self.verbose
        self.verbose = True
        self.evaluate (p, pop)
        self.verbose = verbose
        self.file.flush ()
        self.file    = f
        x = self.__super.print_string (file, p, pop)
        return x
    # end def print_string

# end class Antenna_Optimizer
