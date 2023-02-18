#!/usr/bin/python
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import pga
from mpl_toolkits.mplot3d import Axes3D
from rsclib.autosuper import autosuper
from argparse import ArgumentParser
from math import ceil, log, isnan

import sys
import numbers
import PyNEC

class Excitation (object):
    """ An excitation of the antenna, stores the element tag and segment
        and the voltage (both the real and the imag part
    """
    def __init__ (self, tag, segment, u_real = 1.0, u_imag = 0.0):
        self.tag     = tag
        self.segment = segment
        self.u_real  = u_real
        self.u_imag  = u_imag
    # end def __init__

# end class Excitation

class Nec_File (object):
    """ Implements the methods of geometry and nec context so we can
        output an nec file given the called methods
    """
    def __init__ (self, comment):
        self.repr = []
        if comment:
            if not isinstance (comment, type ([])):
                comment = [comment]
            for c in comment:
                self.add_parameter_comment (c)
        self.end_comments ()
    # end def __init__

    def arc (self, tag, segs, rad, a1, a2, r):
        self.repr.append ("GA %d %d %g %g %g %g" % (tag, segs, rad, a1, a2, r))
    # end def arc

    def add_parameter_comment (self, comment):
        self.repr.append ("CM %s" % comment)
    # end def add_parameter_comment

    def end_comments (self):
        self.repr.append ("CE")
    # end def end_comments

    def move (self, rox, roy, roz, xs, ys, zs, its, nrpt, itgi):
        self.repr.append \
            ( "GM %d %d %g %g %g %g %g %g %d"
            % (itgi, nrpt, rox, roy, roz, xs, ys, zs, its)
            )
    # end def move

    # Patches

    def arbitrary_shaped_patch (self, x1, y1, z1, elev, azim, area):
        self.repr.append \
            ( "SP 0 0 %g %g %g %g %g %g" % (x1, y1, z1, elev, azim, area))
    # end def arbitrary_shaped_patch

    def rectangular_patch (self, x1, y1, z1, x2, y2, z2, x3, y3, z3):
        self.repr.append \
            ( "SP 0 1 %g %g %g %g %g %g" % (x1, y1, z1, x2, y2, z2))
        self.repr.append ( "SC 0 0 %g %g %g" % (x3, y3, z3))
    # end def rectangular_patch

    def triangular_patch (self, x1, y1, z1, x2, y2, z2, x3, y3, z3):
        self.repr.append \
            ( "SP 0 2 %g %g %g %g %g %g" % (x1, y1, z1, x2, y2, z2))
        self.repr.append ( "SC 0 0 %g %g %g" % (x3, y3, z3))
    # end def triangular_patch

    def quadrilateral_patch \
        (self, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
        self.repr.append \
            ( "SP 0 3 %g %g %g %g %g %g" % (x1, y1, z1, x2, y2, z2))
        self.repr.append \
            ( "SC 0 0 %g %g %g %g %g %g" % (x3, y3, z3, x4, y4, z4))
    # end def quadrilateral_patch

    def multiple_patch (self, nx, ny, x1, y1, z1, x2, y2, z2, x3, y3, z3, *x):
        self.repr.append \
            ( "SM %d %d %g %g %g %g %g %g" % (nx, ny, x1, y1, z1, x2, y2, z2))
        self.repr.append ( "SC 0 0 %g %g %g" % (x3, y3, z3))
    # end def multiple_patch

    def helix (self, tag, segs, pitch, h, a1, b1, a2, b2, r):
        self.repr.append \
            ( "GH %d %d %g %g %g %g %g %g %g"
            % (tag, segs, pitch, h, a1, b1, a2, b2, r)
            )
    # end def helix

    def wire (self, tag, segs, x1, y1, z1, x2, y2, z2, r, rdel, rrad):
        self.repr.append \
            ( "GW %d %d %g %g %g %g %g %g %g"
            % (tag, segs, x1, y1, z1, x2, y2, z2, r)
            )
        if rdel != 1 or rrad != 1:
            self.repr.append ("GC 0 0 %g %g %g" % (rdel, rad, rrad * rad))
    # end def wire

    def geometry_complete (self, gnd):
        self.repr.append ("GE %d" % gnd)
    # end def geometry_complete

    def get_geometry (self):
        return self
    # end def get_geometry

    def set_extended_thin_wire_kernel (self, flag):
        self.repr.append ("EK %d" % int (bool (flag)))
    # end def set_extended_thin_wire_kernel

    def ld_card (self, ldtyp, ldtag, ldtagf, ldtagt, zlr, zli, zlc):
        self.repr.append \
            ( "LD %d %d %d %d %g %g %g"
            % (ldtyp, ldtag, ldtagf, ldtagt, zlr, zli, zlc)
            )
    # end def ld_card

    def ex_card (self, *val):
        # Only the 10-value variant is implemented here
        assert 10 <= len (val) <= 10
        self.repr.append ( "EX %d %d %d %d %g %g %g %g %g %g" % val)
    # end def ex_card

    def fr_card (self, ifrq, nfrq, freq_hz, del_freq):
        self.repr.append \
            ("FR %d %d 0 0 %g %g" % (ifrq, nfrq, freq_hz, del_freq))
    # end def fr_card

    def gn_card (self, type, nwires, dieel, cond, t1=0, t2=0, t3=0, t4=0):
        self.repr.append \
            ( "GN %d %d 0 0 %g %g %g %g %g %g"
            % (type, nwires, dieel, cond, t1, t2, t3, t4)
            )
    # end def gn_card

    def nt_card \
        ( self, tag1, seg1, tag2, seg2
        , y11r, y11i, y12r, y12i, y22r, y22i
        ):
        self.repr.append \
            ( "NT %d %d %d %d %.10g %.10g %.10g %.10g %.10g %.10g"
            % ( tag1, seg1, tag2, seg2
              , y11r, y11i, y12r, y12i, y22r, y22i
              )
            )
    # end def nt_card

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
        ):
        self.repr.append \
            ( "RP %d %d %d %d %g %g %g %g %g %g"
            % ( calc_mode, n_theta, n_phi
              , output_format * 1000 + normalization * 100 + D * 10 + A
              , theta0, phi0, delta_theta, delta_phi, radial_distance, gain_norm
              )
            )
    # end def rp_card

    def tl_card \
        ( self, tag1, seg1, tag2, seg2
        , impedance, length
        , shunt_r1, shunt_i1, shunt_r2, shunt_i2
        ):
        self.repr.append \
            ( "TL %d %d %d %d %.10g %.10g %.10g %.10g %.10g %.10g"
            % ( tag1, seg1, tag2, seg2
              , impedance, length
              , shunt_r1, shunt_i1, shunt_r2, shunt_i2
              )
            )
    # end def tl_card

    def __repr__ (self):
        return '\n'.join (self.repr + ['EN'])
    # end def __repr__

# end class Nec_File

class Antenna_Model (autosuper):

    name          = 'Antenna Model'
    wire_radius   = 1.5e-3 / 2.0
    boom_radius   = wire_radius
    impedance     = 50.0
    frq_step_max  = 201
    frq_ranges    = [(430, 440)]
    phi_inc       = 5
    theta_inc     = 5
    phi_range     = 360
    theta_range   = 180

    def __init__ \
        ( self
        , frq_step_max     = 201
        , frq_step_nec     = 201 # only for necout
        , wire_radius      = wire_radius
        , boom_radius      = boom_radius
        , avg_gain         = False
        , force_horizontal = False
        , force_forward    = False
        , force_backward   = False
        , force_same_theta = False
        , copper_loading   = True
        , frq_min          = None
        , frq_max          = None
        ):
        self.theta_max     = int (self.theta_range / self.theta_inc + 1)
        self.phi_max       = int (self.phi_range   / self.phi_inc   + 1)
        self.theta_horz    = 90 / self.theta_inc
        self.phi_mid       = self.phi_max   // 2
        self.wire_radius   = wire_radius
        self.boom_radius   = boom_radius
        self.frq_step_max  = frq_step_max
        self.frq_step_nec  = frq_step_nec
        self.frqinc        = []
        self.frqincnec     = []
        if frq_min:
            self.frq_ranges = []
            for fl, fh in zip (frq_min, frq_max):
                self.frq_ranges.append ((fl, fh))
        for lo, hi in self.frq_ranges:
            self.frqinc.append    ((hi - lo) / (frq_step_max - 1.0))
            self.frqincnec.append ((hi - lo) / (frq_step_nec - 1.0))
        self.force_horizontal = force_horizontal
        # force_forward and force_forward are *not* mutually exclusive,
        # if both are set we use the forward *or* the backward gain
        # whatever is more.
        self.force_forward    = force_forward
        self.force_backward   = force_backward
        self.force_same_theta = force_same_theta
        self.copper_loading   = copper_loading
        # This can be set by register_frequency_callback and is called
        # for each frequency. We can implement frequency dependent
        # network cards where the admittance is different for each
        # frequency.
        self.tl_by_frq        = None
        self.nec              = PyNEC.nec_context ()
        self.avg_gain         = avg_gain
        self.rp               = {}
        self.rp_avg_gain      = {}
        self.geometry          ()
        self.geometry_complete ()
        self.nec_params        ()
        self.transmission_line ()
        self.handle_frequency  ()
    # end def __init__

    def as_nec (self, compute = True):
        c = self.cmdline ().split ('\n')
        if compute:
            if not self.rp:
                if self.avg_gain:
                    self._compute (avgain = True)
                self._compute ()
            for frq_idx in range (len (self.frq_ranges)):
                c.extend (self.show_gains (frq_idx))
        n = Nec_File (c)
        self.geometry          (n)
        self.geometry_complete (n)
        self.nec_params        (n)
        self.transmission_line (n)
        self.handle_frequency  (n)
        if self.avg_gain:
            self._compute (n, avgain = True)
        self._compute          (n)
        return repr            (n)
    # end def as_nec

    def cmdline (self):
        """ This should be overridden in derived class to print out
            command-line parameters to regenerate this model.
        """
        return ''
    # end def cmdline

    def _compute (self, nec = None, avgain = False):
        if nec is None:
            nec = self.nec
        if avgain:
            self.nec_params_avg_gain (nec)
        else:
            self.nec_params_compute (nec)
        for n, (lo, hi) in enumerate (self.frq_ranges):
            if callable (self.tl_by_frq):
                for i in range (self.frq_max_idx):
                    f = lo + i * self.frq_inc [n]
                    self.tl_by_frq (nec, f)
                    nec.fr_card (0, 1, f, 0)
                    nec.rp_card \
                        ( 0, self.theta_max, self.phi_max
                        , 0, 0, 0, int (self.avg_gain), 0, 0
                        , self.theta_inc, self.phi_inc, 0, 0
                        )
            else:
                nec.fr_card (0, self.frq_max_idx, lo, self.frq_inc [n])
                nec.rp_card \
                    ( 0, self.theta_max, self.phi_max
                    , 0, 0, 0, int (self.avg_gain), 0, 0
                    , self.theta_inc, self.phi_inc, 0, 0
                    )
    # end def _compute

    def compute (self, frq_step = None, avgain = False):
        self._compute (avgain = avgain)
        rp  = self.rp
        off = self.avg_offset
        if avgain:
            rp  = self.rp_avg_gain
            off = 0
        # If index isn't explicitly given, use the middle
        if frq_step is None:
            frq_step = self.frq_step_max // 2
        for n, (lo, hi) in enumerate (self.frq_ranges):
            idx = n * self.frq_step_max + frq_step
            if idx not in rp:
                rp [idx] = self.nec.get_radiation_pattern (idx + off)
    # end def compute

    def frq_step_range (self, step = 1):
        return range (0, self.frq_step_max, step)
    # end def frq_step_range

    def geometry (self, nec):
        """ Derived class *must no longer* call geometry_complete!
        """
        raise NotImplemented ("Derived class must implement 'geometry'")
    # end def geometry

    def geometry_complete (self, nec = None):
        """ Currently no ground model
            Derived classes may implement a 'ground' method, if this
            exists we set the ground flag on nec.geometry_complete.
            Note that derived classes should *not* apply ground in
            geometry_complete anymore but should implement their own
            ground method.
        """
        if nec is None:
            nec = self.nec
        ground = 0
        if hasattr (self, 'ground'):
            ground = 1
        nec.geometry_complete (ground)
    # end def geometry_complete

    def transmission_line (self, nec = None):
        """ This optionally implements transmission-line (TL) or network
            cards. These come *after* the excitation (EX) card.
            Note that this can -- instead of immediately generating the
            TL or NT cards -- call register_frequency_callback to
            produce frequency-specific NT cards instead.
        """
    # end def transmission_line

    def handle_frequency (self, nec = None):
        frq_step_max = self.frq_step_nec
        frqinc       = self.frqincnec
        if nec is None:
            nec          = self.nec
            frq_step_max = self.frq_step_max
            frqinc       = self.frqinc
        self.frq_max_idx = frq_step_max
        self.frq_inc     = frqinc
        self.avg_offset  = 0
        if self.avg_gain:
            self.avg_offset = self.frq_max_idx * len (self.frq_ranges)
    # end def handle_frequency

    def nec_params (self, nec = None):
        if nec is None:
            nec = self.nec
        nec.set_extended_thin_wire_kernel (True)
        if isinstance (self.ex, Excitation):
            self.ex = [self.ex]
        for ex in self.ex:
            nec.ex_card \
                (0, ex.tag, ex.segment, 0, ex.u_real, ex.u_imag, 0, 0, 0, 0)
    # end def nec_params

    def nec_params_compute (self, nec = None):
        """ NEC cards to set when doing the *real* computation, if
            average gain computation has been specified these must not
            be set at first (because we compute average gain with a
            perfect ground (if there *is* ground) and without copper or
            other resistive loading.
        """
        if nec is None:
            nec = self.nec
        if hasattr (self, 'ground'):
            self.ground (nec)
        if self.copper_loading:
            nec.ld_card (5, 0, 0, 0, 37735849, 0, 0)
    # end def nec_params_compute

    def nec_params_avg_gain (self, nec = None):
        """ NEC cards to set when doing the *average gain* computation.
            Only called when average gain has been requested.
        """
        if nec is None:
            nec = self.nec
        if hasattr (self, 'ground'):
            # Perfect ground
            nec.gn_card (1, 0, 0, 0, 0, 0, 0, 0)
    # end def nec_params_avg_gain

    def max_f_r_gain (self, frq = 0, frq_step = None):
        """ Maximum forward and backward gain
            If we have requested average gain computation, this corrects
            the gain by the average gain.
        """
        if frq_step is None:
            frq_step = self.frq_step_max // 2
        idx = frq * self.frq_step_max + frq_step
        if idx not in self.rp:
            self.rp [idx] = self.nec.get_radiation_pattern \
                (idx + self.avg_offset)
        if self.avg_gain and idx not in self.rp_avg_gain:
            self.rp_avg_gain [idx] = self.nec.get_radiation_pattern (idx)
        gains = self.rp [idx].get_gain ()
        n1max = n2max = -1
        gmax  = None
        # First loop is theta, second is phi
        for n1, ga in enumerate (gains):
            for n2, g in enumerate (ga):
                if not self.force_horizontal or n1 == self.theta_horz:
                    if  (  (not self.force_forward and not self.force_backward)
                        or (self.force_forward  and n2 ==  0)
                        or (self.force_backward and n2 == self.phi_mid)
                        ):
                        if gmax is None or g > gmax:
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
        trange = range (n1max - t30, n1max + t30 + 1)
        if self.force_same_theta:
            trange = [n1max]
        for t in trange:
            theta = t
            if theta < 0:
                theta += self.theta_max - 1
            if theta >= self.theta_max:
                theta -= self.theta_max - 1
            assert 0 <= theta < self.theta_max
            for p in range (n2max - p30 - phi2, n2max + p30 + 1 - phi2):
                phi = p
                if phi < 0:
                    phi += self.phi_max - 1
                if phi >= self.phi_max:
                    phi -= self.phi_max - 1
                assert 0 <= phi < self.phi_max
                if rmax is None or gains [theta][phi] > rmax:
                    rmax = gains [theta][phi]
                    pm   = phi
                    tm   = theta
        if self.avg_gain:
            avg  = self.rp_avg_gain [idx].get_average_power_gain ()
            # Seems to happen for ill-conditioned antennas, obviously
            # the average gain should be always positive.
            # We make it very large to subtract a high amount from the
            # gain
            if avg <= 0:
                avg = 1e12
            exp  = 1.0
            if hasattr (self, 'ground'):
                assert self.theta_range == 90
                exp = 2.0
            avdb = 10 * (log (exp / avg) / log (10))
            # We don't want ill-conditioned antennas where we correct
            # the gain *up*. Only down corrections of the gain happen.
            if avdb < 0:
                gmax += avdb
                rmax += avdb
        return gmax, rmax
    # end def max_f_r_gain

    def show_gains (self, frq_idx = 0, prefix = ''):
        r = []
        step = self.frq_step_max // 2
        r.append ('FRQ Range: %.2f-%.2f' % self.frq_ranges [frq_idx])
        for frqstep in self.frq_step_range (step):
            idx = frq_idx * self.frq_step_max + frqstep
            f, b = self.max_f_r_gain (frq_idx, frqstep)
            frq = self.rp [idx].get_frequency ()
            frq = frq / 1e6
            rr = "%sFRQ: %3.2f fw: %2.2f bw: %2.2f" % (prefix, frq, f, b)
            if self.avg_gain:
                rpa = self.rp_avg_gain [idx]
                rr += " average gain: %.5f solid angle: %.4f" \
                    % ( rpa.get_average_power_gain ()
                      , rpa.get_average_power_solid_angle ()
                      )
            r.append (rr)
        vswrs = list \
            (self.vswr (frq_idx, i) for i in self.frq_step_range (step))
        r.append ("SWR: %1.2f %1.2f %1.2f" % tuple (vswrs))
        return r
    # end def show_gains

    def plot (self, frq_idx = 0, frq_step = None):
        if frq_step is None:
            frq_step = self.frq_step_max // 2
        idx = self.frq_step_max * frq_idx + frq_step
        if not self.rp or idx not in self.rp:
            if self.avg_gain:
                self.compute (idx, avgain = True)
            self.compute (idx)

        # 0: linear, 1: right, 2: left
        #print (self.rp [frq_step].get_pol_sense_index ())
        #print (self.rp [frq_step].get_pol_tilt ())
        #print (self.rp [frq_step].get_pol_axial_ratio ())

        f = self.rp [idx].get_frequency ()
        gains  = self.rp [idx].get_gain ()
        gains  = 10.0 ** (gains / 10.0)
        # Display max gain in dBi

        thetas = self.rp [idx].get_theta_angles () * np.pi / 180.0
        phis   = self.rp [idx].get_phi_angles ()   * np.pi / 180.0

        P, T = np.meshgrid (phis, thetas)

        X = np.cos (P) * np.sin (T) * gains
        Y = np.sin (P) * np.sin (T) * gains
        Z = np.cos (T) * gains

        fig = plt.figure ()
        ax  = fig.gca (projection='3d')
        t   = "%s %.2f MHz" % (self.name, f / 1e6)
        ax.set_title (t)
        ax.set_xlabel ('X')
        ax.set_ylabel ('Y')
        ax.set_zlabel ('Z')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])

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

        ax.plot_wireframe (X, Y, Z, color = 'r', linewidth = 0.5)
        plt.show ()
    # end def plot

    def swr_plot (self):
        """ If we have several frequency ranges we do a plot for each
        """
        fun   = self.nec.get_radiation_pattern
        for frq in range (len (self.frq_ranges)):
            offset = frq * self.frq_step_max + self.avg_offset
            frqs  = []
            vswrs = []
            for i in self.frq_step_range ():
                frqs.append  (fun (i + offset).get_frequency ())
                vswrs.append (self.vswr (frq, i))
            fig = plt.figure ()
            ax  = fig.add_subplot (111)
            ax.plot (frqs, vswrs)
            ax.set_title ('Freq range: %.2f - %.2f MHz' % self.frq_ranges [frq])
            plt.show ()
    # end def swr_plot

    def vswr (self, frq_idx, frq_step):
        off = frq_idx * self.frq_step_max + self.avg_offset
        ipt = self.nec.get_input_parameters (off + frq_step)
        z   = ipt.get_impedance ()
        rho = np.abs ((z - self.impedance) / (z + self.impedance))
        return ((1. + rho) / (1. - rho)) [0]
    # end def vswr

    def register_frequency_callback (self, method):
        self.tl_by_frq = method
    # end def register_frequency_callback

# end class Antenna_Model

class Antenna_Phenotype (autosuper):
    def __init__ (self, optimizer, antenna, frq_idx):
        self.optimizer = optimizer
        self.antenna   = antenna
        self.frq_idx   = frq_idx
        self.offset    = frq_idx * antenna.frq_step_max
        self.vswrs     = vswrs = list \
            (antenna.vswr (frq_idx, i) for i in antenna.frq_step_range ())
        # Looks like NEC sometimes computes negative SWR
        # We set the SWR to something very high in that case
        for swr in vswrs:
            if swr < 0:
                swr_eval = swr_med = 1e6
                gmax, rmax = (-20.0, 0.0)
                break
        else:
            swr_eval  = sum (vswrs) / 3.0
            swr_med   = swr_eval
            swr_eval *= 1 + sum (6 * bool (v > optimizer.maxswr) for v in vswrs)
            diff = abs (vswrs [0] - vswrs [-1])
            if diff > 0.2:
                swr_eval *= 1.0 + 15 * diff
            # If relax_swr is given, do not bother as long a swr is
            # below optimizer.maxswr
            if optimizer.relax_swr and max (vswrs) <= optimizer.maxswr:
                swr_eval = 1.0

            # We take the *minimum* gain over all frequencies
            # and the *maximum* rear gain over all frequencies
            gmax = None
            rmax = None
            for idx in antenna.frq_step_range ():
                f, b = antenna.max_f_r_gain (frq_idx, idx)
                if gmax is None or gmax > f:
                    gmax = f
                if rmax is None or rmax < b:
                    rmax = b
        mid = antenna.frq_step_range () [len (antenna.frq_step_range ()) // 2]
        self.gmid, self.rmid = antenna.max_f_r_gain (frq_idx, mid)
        if optimizer.nofb:
            rmax = 0.0
        swr_eval **= (1./2)
        self.gmax = gmax = max (gmax, -20.0)
        self.rmax = rmax = min (rmax,  20.0) # Don't allow too large values
        self.swr_eval = swr_eval
        self.swr_med  = swr_med
    # end def __init__
# end class Antenna_Phenotype

class Antenna_Optimizer (pga.PGA, autosuper):
    """ Optimize given antenna, needs to be subclassed.
    """

    resolution = 0.5e-3 # 0.5 mm in meter
    ant_cls    = Antenna_Model

    def __init__ \
        ( self
        , random_seed      = 42
        , verbose          = False
        , wire_radius      = 0.00075
        , nofb             = False
        , maxswr           = 1.8
        , use_rtr          = True
        , randselect       = False
        , use_de           = True
        , popsize          = 100
        , force_horizontal = False
        , force_forward    = False
        , force_backward   = False
        , force_same_theta = False
        , relax_swr        = False
        , copper_loading   = True
        , multiobjective   = False
        , nsga_iii         = None
        , title            = None
        , stagnation_max   = 100
        , avg_gain         = False
        , min_gain         = 0.0
        , min_fb           = 0.0
        , use_mid          = False
        , frq_min          = None
        , frq_max          = None
        , ** kw
        ):
        self.verbose          = verbose
        self.wire_radius      = wire_radius
        self.nofb             = nofb
        self.maxswr           = maxswr
        self.use_rtr          = use_rtr
        self.use_de           = use_de
        self.popsize          = popsize
        self.force_horizontal = force_horizontal
        self.force_forward    = force_forward
        self.force_backward   = force_backward
        self.force_same_theta = force_same_theta
        self.relax_swr        = relax_swr
        self.stag_count       = 0
        self.copper_loading   = copper_loading
        self.multiobjective   = multiobjective
        self.nsga_iii         = nsga_iii
        self.title            = title
        self.stagnation_max   = stagnation_max
        self.avg_gain         = avg_gain
        self.min_gain         = min_gain
        self.min_fb           = min_fb
        self.use_mid          = use_mid
        if frq_min:
            self.frq_ranges = []
            for fl, fh in zip (frq_min, frq_max):
                self.frq_ranges.append ((fl, fh))
        else:
            self.frq_ranges   = self.ant_cls.frq_ranges
        evc                   = self.get_eval_and_constraints ()
        refpoints             = None
        stop_on               = \
            [ pga.PGA_STOP_NOCHANGE
            , pga.PGA_STOP_MAXITER
            , pga.PGA_STOP_TOOSIMILAR
            ]
        num_replace           = popsize // 2
        pop_replace_type      = pga.PGA_POPREPL_BEST
        typ                   = bool
        if use_de:
            num_replace       = popsize
            pop_replace_type  = pga.PGA_POPREPL_PAIRWISE_BEST
            length            = len (self.minmax)
            typ               = float
        elif use_rtr:
            num_replace       = popsize
            pop_replace_type  = pga.PGA_POPREPL_RTR
        if multiobjective:
            if nsga_iii:
                pop_replace_type = pga.PGA_POPREPL_NSGA_III
                dim = evc ['num_eval'] - evc ['num_constraint']
                refpoints = pga.das_dennis (dim, nsga_iii)
                lref = len (refpoints)
                # Round up to next even number
                if lref > popsize:
                    popsize = int ((lref / 2 + .5) * 2)
            else:
                pop_replace_type = pga.PGA_POPREPL_NSGA_II
            num_replace       = popsize
        # Determine number of bits needed from minmax,
        # we need at least self.resolution precision.
        if not use_de:
            self.nbits = []
            for l, u in self.minmax:
                n = (u - l) / self.resolution
                self.nbits.append (int (ceil (log (n) / log (2))))
            self.bitidx = []
            l   = 0
            for b in self.nbits:
                u = l + b - 1
                self.bitidx.append ((l, u))
                l = u + 1
            length = sum (self.nbits)
        args = dict \
            ( maximize            = True
            , pop_size            = popsize
            , num_replace         = num_replace
            , random_seed         = random_seed
            , print_options       = [pga.PGA_REPORT_STRING]
            , stopping_rule_types = stop_on
            , pop_replace_type    = pop_replace_type
            , print_frequency     = 10
            , randomize_select    = bool (randselect)
            )
        if refpoints is not None:
            args ['reference_points']     = refpoints
        if self.use_de:
            args ['init']                 = self.minmax
            args ['select_type']          = pga.PGA_SELECT_LINEAR
            args ['mutation_bounce_back'] = True
            args ['mutation_only']        = True
            args ['mutation_type']        = pga.PGA_MUTATION_DE
            args ['DE_variant']           = pga.PGA_DE_VARIANT_BEST
            args ['DE_crossover_prob']    = 0.2
            args ['DE_jitter']            = 0.001
            args ['DE_scale_factor']      = 0.85 - (popsize * 0.0005)
            args ['DE_crossover_type']    = pga.PGA_DE_CROSSOVER_BIN
        if 'DE_variant' in kw:
            if kw ['DE_variant'] == 'rand':
                args ['DE_variant'] = pga.PGA_DE_VARIANT_RAND
            if kw ['DE_variant'] == 'either-or':
                args ['DE_variant'] = pga.PGA_DE_VARIANT_EITHER_OR
        if 'DE_crossover_prob' in kw:
            args ['DE_crossover_prob'] = kw ['DE_crossover_prob']
        if 'DE_jitter' in kw:
            args ['DE_jitter'] = kw ['DE_jitter']
        if 'DE_dither' in kw:
            args ['DE_dither'] = kw ['DE_dither']
        if 'epsilon_generation' in kw:
            args ['epsilon_generation'] = kw ['epsilon_generation']
        args.update (evc)
        # Always sum up constraint violations. Alternative seems to be buggy
        if multiobjective:
            args ['sum_constraints'] = True
        pga.PGA.__init__ (self, typ, length, ** args)
        self.last_best = [float ('nan')] * (self.num_eval - self.num_constraint)
        if self.title is None:
            self.title = "%s %s" % (self.__class__.__name__, self.random_seed)
        self.cache = {}
        self.cache_hits = 0
        self.nohits     = 0
        self.file       = sys.stdout
    # end def __init__

    @property
    def antenna_args (self):
        """ Optimizer arguments passed to the constructed antenna during
            optimization.
        """
        d = dict \
            ( copper_loading   = self.copper_loading
            , force_backward   = self.force_backward
            , force_forward    = self.force_forward
            , force_horizontal = self.force_horizontal
            , force_same_theta = self.force_same_theta
            , wire_radius      = self.wire_radius
            , frq_step_max     = 3
            )
        return d
    # end def antenna_args

    @property
    def nfreq (self):
        return len (self.frq_ranges)
    # end def nfreq

    def get_eval_and_constraints (self):
        """ This computes the number of evaluations and constraints.
            By default there are no constraints *and* the number of
            evals is the default of 1 (no multiobjective eval). For
            multiobjective evaluation we have one constraint (the VSWR
            must be below a threshold), in addition we may want a
            minimum gain and a mimimum forward/backward ratio as a
            constraint which allows us to easier find a good solution.

            Of course if we optimize something differently, e.g., for
            multiple frequencies, we want to return something completely
            different here.
        """
        args = {}
        if self.multiobjective:
            to_add = bool (self.min_gain) + bool (self.min_fb)
            args ['num_eval']       = (3 + to_add) * self.nfreq
            args ['num_constraint'] = (1 + to_add) * self.nfreq
        return args
    # end def get_eval_and_constraints

    def get_parameter (self, p, pop, i):
        """ Get floating-point value from encoded allele
            We tried gray code but now use binary (BCD) encoding.
        """
        if self.use_de:
            return self.get_allele (p, pop, i)
        return self.get_real_from_binary \
            (p, pop, *(self.bitidx [i] + self.minmax [i]))
    # end def get_parameter

    def set_parameter (self, p, pop, i, val):
        """ set gene from floating-point value
            We tried gray code but now use binary (BCD) encoding.
            There is obviously a rounding error sometimes when *reading*
            the value (see get_parameter above), so we limit val to
            minmax [0] <= val <= minmax [1]
        """
        if self.use_de:
            self.set_allele (p, pop, i, val)
        else:
            if val > self.minmax [i][1]:
                print ('Oops')
                val = self.minmax [i][1]
            if val < self.minmax [i][0]:
                print ('Oops')
                val = self.minmax [i][0]
            self.encode_real_as_binary \
                (p, pop, *(self.bitidx [i] + self.minmax [i] + (val,)))
    # end def set_parameter

    def cache_key (self, p, pop):
        if self.use_de:
            return tuple \
                (self.get_allele (p, pop, k) for k in range (len (self)))
        ck = 0
        for k in range (len (self)):
            ck <<= 1
            ck |= int (self.get_allele (p, pop, k))
        return ck
    # end def cache_key

    def pre_eval (self, pop):
        # Do not run before very first eval
        if pop != pga.PGA_NEWPOP:
            return
        for p in range (self.pop_size):
            if self.get_evaluation_up_to_date (p, pop):
                continue
            ck = self.cache_key (p, pop)
            if ck in self.cache:
                self.cache_hits += 1
                self.set_evaluation (p, pop, self.cache [ck])
                self.set_evaluation_up_to_date (p, pop, True)
            else:
                self.nohits += 1
    # end def pre_eval

    def phenotype (self, p, pop):
        antenna = self.compute_antenna (p, pop)
        antenna.compute ()
        pheno = []
        for n, frq in enumerate (antenna.frq_ranges):
            pheno.append (Antenna_Phenotype (self, antenna, n))
        return pheno
    # end def phenotype

    def evaluate (self, p, pop):
        phenos = self.phenotype (p, pop)
        if self.multiobjective:
            retval = []
            for pheno in phenos:
                swr_max = max (pheno.vswrs)
                if self.use_mid:
                    retval.append \
                        ([ pheno.gmid
                         , pheno.gmid - pheno.rmid
                         , swr_max - self.maxswr
                        ])
                else:
                    retval.append \
                        ([ pheno.gmax
                         , pheno.gmax - pheno.rmax
                         , swr_max - self.maxswr
                        ])
                if self.min_gain:
                    retval [-1].append (self.min_gain - pheno.gmax)
                if self.min_fb:
                    retval [-1].append (self.min_fb + pheno.rmax - pheno.gmax)
            return tuple (np.array (retval).T.flatten ())

        # If not multiobjective we do currently not support multiple
        # frequencies
        assert len (phenos) == 1
        pheno = phenos [0]
        egm   = pheno.gmax ** 3.0
        # Don't use gmax ** 3 if too much swr:
        if pheno.swr_med > 3.0:
            egm = pheno.gmax

        eval = ( 100.0
               + egm
               - pheno.rmax * 4
               ) / pheno.swr_eval
        if eval < 0:
            eval = 0.0
        assert eval >= 0
        return eval
    # end def evaluate

    def endofgen (self):
        """ Simple hill-climb: Loop over all individuums, select one of
            the four alleles by random and try to inc/dec (randomly).
            If the inc/decremented gene is better, update allele and
            evaluation.
        """
        # No hillclimbing for Differential Evolution
        if self.use_de:
            return
        pop  = pga.PGA_NEWPOP
        l    = len (self.nbits) # number of bits per float
        bidx = self.get_best_index (pop)
        best = self.get_evaluation (bidx, pop)
        calc = False
        for p in range (self.pop_size):
            ck  = self.cache_key (p, pop)
            ev  = self.get_evaluation (p, pop)
            assert self.get_evaluation_up_to_date (p, pop)
            if ck not in self.cache:
                self.cache [ck] = ev
            idx = self.random_interval (0, l - 1)
            val = self.get_parameter (p, pop, idx)
            if self.random_flip (0.5):
                if val + self.resolution <= self.minmax [idx][1]:
                    self.set_parameter (p, pop, idx, val + self.resolution)
            else:
                if val - self.resolution >= self.minmax [idx][0]:
                    self.set_parameter (p, pop, idx, val - self.resolution)
            ck  = self.cache_key (p, pop)
            if ck in self.cache:
                self.cache_hits += 1
                evnew = self.cache [ck]
            else:
                self.nohits += 1
                self.set_evaluation_up_to_date (p, pop, False)
                evnew = self.evaluate (p, pop)
                self.cache [ck] = evnew
                self.set_evaluation_up_to_date (p, pop, True)

            if evnew <= ev:
                # undo
                self.set_parameter (p, pop, idx, val)
            else:
                self.set_evaluation (p, pop, evnew)
                if evnew > best:
                    calc = True
        # Re-calculate fitness values if the best index changed
        if calc:
            self.fitness (pop)
    # end def endofgen

    def get_gain_fw_maxswr (self, sub_eval):
        """ From the relevant parameters for a single frequency-range
            get the gain, the f/b ratio (in dB) and the max vswr.
            Needs to return a tuple.
        """
        sub_eval [-1] += self.maxswr
        return tuple (sub_eval)
    # end def get_gain_fw_maxswr

    def print_string (self, file, p, pop):
        f         = self.file
        self.file = file
        antenna   = self.compute_antenna (p, pop)
        print ("Title: %s" % self.title, file = file)
        print (antenna.cmdline (), file = file)
        if self.verbose:
            phenos = self.phenotype (p, pop)
            for n, pheno in enumerate (phenos):
                fr = pheno.antenna.frq_ranges [n]
                print \
                    ( "F:%g-%g VSWR: %s\nGMAX: %.2f, RMAX: %.2f"
                    % (fr [0], fr [1], pheno.vswrs, pheno.gmax, pheno.rmax)
                    , file = file
                    )
        if self.multiobjective:
            arg  = self.get_eval_and_constraints ()
            l    = arg ['num_eval'] // self.nfreq
            eval = np.array (self.get_evaluation (p, pop))
            eval = np.reshape (eval, (l, self.nfreq)).T
            for n in range (self.nfreq):
                fr = self.frq_ranges [n]
                ev = fr + self.get_gain_fw_maxswr (eval [n][:3])
                assert isinstance (ev, tuple)
                print ( "F:%g-%g Gain: %.2f dBi, "
                        "F/B ratio: %.2f dB, max VSWR: %.2f"
                      % ev, file = file
                      )
        ch = self.cache_hits
        cn = self.nohits + self.cache_hits
        print \
            ( "Cache hits: %s/%s %2.2f%%" % (ch, cn, 100.0 * ch / cn)
            , file = file
            )
        print \
            ( "Iter: %s Evals: %s Stag: %s"
            % (self.GA_iter, self.eval_count, self.stag_count)
            , file = file
            )
        file.flush ()
        self.file = f
        x = self.__super.print_string (file, p, pop)
        return x
    # end def print_string

    def stop_cond (self):
        """ Experimental early stopping when stagnating
            for Differential Evolution
        """
        if self.use_de:
            num_f = self.num_eval - self.num_constraint
            for k in range (num_f):
                best_ev = self.get_best_report (pga.PGA_OLDPOP, k)
                lbest   = self.last_best [k]
                if isnan (lbest) or abs (lbest - best_ev) >= abs (lbest) / 500:
                    self.stag_count = 0
                    break
            else:
                self.stag_count += 1
                if self.stag_count >= self.stagnation_max:
                    return True
            for k in range (num_f):
                self.last_best [k] = self.get_best_report (pga.PGA_OLDPOP, k)
        return self.check_stopping_conditions ()
    # end def stop_cond

# end class Antenna_Optimizer

class Arg_Handler:

    """ Encapsulate options that occur in (almost) every antenna
        or optimizer for an antenna.
    """
    actions = ['optimize', 'necout', 'swr', 'gain', 'frgain']

    def __init__ (self, **default):
        self.default = default
        self.cmd = cmd = ArgumentParser ()
        cmd.add_argument \
            ( 'action'
            , help = "Action to perform, one of %s" % ', '.join (self.actions)
            )
        cmd.add_argument \
            ( '-a', '--average-gain'
            , action  = "store_true"
            , help    = "Output average gain in nec file and use during"
                        " optimization"
                        " (unsupported by xnec2c)"
            )
        cmd.add_argument \
            ( '--epsilon-generation'
            , help    = "Use epsilon constraints until this generation"
            , type    = int
            , default = 0
            )
        cmd.add_argument \
            ( '--frq-max'
            , help    = "Add a maximum frequency, must be matched with "
                        "frq-min option, can be specified multiple times"
            , action  = 'append'
            , type    = float
            , default = []
            )
        cmd.add_argument \
            ( '--frq-min'
            , help    = "Add a minimum frequency, must be matched with "
                        "frq-max option, can be specified multiple times"
            , action  = 'append'
            , type    = float
            , default = []
            )
        cmd.add_argument \
            ( '--force-horizontal'
            , help    = "Consider gain only in horizontal plane"
            , action  = 'store_true'
            )
        cmd.add_argument \
            ( '--force-forward'
            , help    = "Consider only forward gain"
            , action  = 'store_true'
            )
        cmd.add_argument \
            ( '--force-backward'
            , help    = "Consider only backward gain"
            , action  = 'store_true'
            )
        cmd.add_argument \
            ( '--force-same_theta'
            , help    = "Use same theta in backward direction as in forward"
            , action  = 'store_true'
            )
        cmd.add_argument \
            ( '-i', '--frq-step-max'
            , type = int
            , help = "Number of frequency steps, default=%(default)s"
            , default = self.default.get ('frq_step_max', 21)
            )
        cmd.add_argument \
            ( '--max-swr'
            , type    = float
            , help    = "Maximum SWR to considered good when optimizing"
                        ", default=%(default)g"
            , default = self.default.get ('max_swr', 1.8)
            )
        cmd.add_argument \
            ( '--no-copper-loading'
            , help    = "Do not insert an LD card defining material copper"
            , action  = 'store_false'
            , dest    = 'copper_loading'
            , default = self.default.get ('copper_loading', True)
            )
        cmd.add_argument \
            ( '-P', '--popsize'
            , type    = int
            , help    = "Population size, default=%(default)s"
            , default = self.default.get ('popsize', 100)
            )
        cmd.add_argument \
            ( '-R', '--random-seed'
            , type    = int
            , help    = "Random number seed for optimizer, default=%(default)s"
            , default = self.default.get ('random_seed', 42)
            )
        cmd.add_argument \
            ( '--stagnation-max'
            , type    = int
            , help    = "Maximum stagnation generations, default=%(default)s"
            , default = self.default.get ('stagnation_max', 100)
            )
        cmd.add_argument \
            ( '-w', '--wire-radius'
            , type    = float
            , help    = "Radius of the wire, default=%(default)g"
            , default = self.default.get ('wire_radius', 0.00075)
            )
        cmd.add_argument \
            ( '-v', '--verbose'
            , help    = "Verbose reporting in every generation"
            , action  = 'store_true'
            )
        cmd.add_argument \
            ( '--no-rtr'
            , help    = "Do not use restricted tournament replacement"
            , dest    = "use_rtr"
            , default = True
            , action  = "store_false"
            )
        cmd.add_argument \
            ( '--no-de'
            , help    = "Do not use differential evolution"
            , dest    = "use_de"
            , default = True
            , action  = "store_false"
            )
        cmd.add_argument \
            ( '--multiobjective'
            , help    = "Use multi-objective optimization"
            , dest    = "multiobjective"
            , action  = "store_true"
            , default = self.default.get ('multiobjective', False)
            )
        cmd.add_argument \
            ( '-N', '--nsga-iii'
            , help    = "Use NSGA-III for multi-objective optimization,"
                        " specify number of Das/Dennis partitions"
            , type    = int
            )
        cmd.add_argument \
            ( '--randomize-select'
            , help    = "Randomize select again for backward compatibility"
            , action  = "store_true"
            )
        cmd.add_argument \
            ( '--relax-swr'
            , help    = "Don't optimize SWR below max_swr"
            , action  = "store_true"
            )
        cmd.add_argument \
            ( '--DE-variant'
            , help    = "Differential Evolution variant, one of rand, either_or"
                        ", default='%(default)s'"
            , default = 'best'
            )
        cmd.add_argument \
            ( '--DE-crossover-prob'
            , help    = "Differential Evolution crossover probability, "
                        "default=%(default)s"
            , type    = float
            , default = 0.9
            )
        cmd.add_argument \
            ( '--DE-jitter'
            , help    = "Differential Evolution jitter, default=%(default)s"
            , type    = float
            , default = 0.001
            )
        cmd.add_argument \
            ( '--DE-dither'
            , help    = "Differential Evolution dither, default=%(default)s"
            , type    = float
            , default = 0.2
            )
        cmd.add_argument \
            ( '--min-gain'
            , help    = "Minimum gain as a constraint"
            , type    = float
            , default = 0.0
            )
        cmd.add_argument \
            ( '--min-fb'
            , help    = "Minimum forward/backward ratio as a constraint"
            , type    = float
            , default = 0.0
            )
        cmd.add_argument \
            ( '--use-mid'
            , help    = "Use middle frequency for gain and f/b ratio"
                        " (Default is to use maximum gain and minimum"
                        " f/b ratio over all frequencies)"
            , action  = 'store_true'
            )
    # end def __init__

    def __getattr__ (self, name):
        """ Delegate to self.cmd, cache result """
        # This will raise AttributeError if name doesn't exist in self.cmd
        setattr (self, name, getattr (self.cmd, name))
        return getattr (self, name)
    # end def __getattr__

    @property
    def default_optimization_args (self):
        d = dict \
            ( avg_gain           = self.args.average_gain
            , random_seed        = self.args.random_seed
            , randselect         = self.args.randomize_select
            , use_rtr            = self.args.use_rtr
            , verbose            = self.args.verbose
            , wire_radius        = self.args.wire_radius
            , use_de             = self.args.use_de
            , popsize            = self.args.popsize
            , force_horizontal   = self.args.force_horizontal
            , force_forward      = self.args.force_forward
            , force_backward     = self.args.force_backward
            , force_same_theta   = self.args.force_same_theta
            , relax_swr          = self.args.relax_swr
            , maxswr             = self.args.max_swr
            , copper_loading     = self.args.copper_loading
            , multiobjective     = self.args.multiobjective
            , nsga_iii           = self.args.nsga_iii
            , stagnation_max     = self.args.stagnation_max
            , DE_variant         = self.args.DE_variant
            , DE_crossover_prob  = self.args.DE_crossover_prob
            , DE_jitter          = self.args.DE_jitter
            , DE_dither          = self.args.DE_dither
            , epsilon_generation = self.args.epsilon_generation
            , min_gain           = self.args.min_gain
            , min_fb             = self.args.min_fb
            , use_mid            = self.args.use_mid
            , frq_min            = self.args.frq_min
            , frq_max            = self.args.frq_max
            )
        return d
    # end def default_optimization_args

    @property
    def default_antenna_args (self):
        frq_step_max = self.args.frq_step_max
        if self.args.action in ('frgain', 'necout'):
            frq_step_max = 3
        d = dict \
            ( avg_gain       = self.args.average_gain
            , frq_step_max   = frq_step_max
            , frq_step_nec   = self.args.frq_step_max
            , wire_radius    = self.args.wire_radius
            , copper_loading = self.args.copper_loading
            , frq_min        = self.args.frq_min
            , frq_max        = self.args.frq_max
            )
        return d
    # end def default_antenna_args

    def add_argument (self, *args, **kw):
        if 'help' in kw and kw.get ('type', None) == float and 'default' in kw:
            kw ['help'] = kw ['help'] + ' default=%(default)g'
        self.cmd.add_argument (*args, **kw)
    # end def add_argument

    def parse_args (self, *args, **kw):
        self.args = self.cmd.parse_args (*args, **kw)
        return self.args
    # end def parse_args

# end class Arg_Handler

def antenna_actions (cmd, args, antenna):
    if args.action == 'necout':
        print (antenna.as_nec ())
    elif args.action not in cmd.actions:
        cmd.print_usage ()
    else:
        if antenna.avg_gain:
            antenna.compute (avgain = True)
        antenna.compute ()
    if args.action == 'swr':
        antenna.swr_plot ()
    elif args.action == 'gain':
        for frq_idx in range (len (antenna.frq_ranges)):
            antenna.plot (frq_idx)
    elif args.action == 'frgain':
        for frq_idx in range (len (antenna.frq_ranges)):
            print ('\n'.join (antenna.show_gains ()))
# end def antenna_actions
