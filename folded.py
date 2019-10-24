from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

import math
import PyNEC

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

    def __init__ (self, refl_dist = 0.010) :
        self.refl_dist = refl_dist
        self.tag       = 1
        self.nec       = PyNEC.nec_context ()
        self.nec.set_extended_thin_wire_kernel (True)
        self.geo       = self.nec.get_geometry ()
        self.ex        = None
        a = ((-90, 90), (90, 270))
        for n, z in enumerate ((1, -1)) :
            a1, a2 = a [n]
            self.geo.arc \
                ( self.tag
                , self.segs_arc
                , self.dipole_radius
                , a1
                , a2
                , self.wire_radius
                )
            self.geo.move (0, 0, 0, z * self.lambda_4, 0, 0, self.tag, 0, 0)
            self.tag += 1
        for x in (-self.lambda_4, self.lambda_4) :
            for z in (self.dipole_radius, -self.dipole_radius) :
                self.geo.wire \
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
        self.geo.wire \
            ( self.tag
            , self.segs_boom
            , 0, 0,  self.dipole_radius
            , 0, 0, -self.dipole_radius
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # second part of boom from dipole to reflector
        self.geo.wire \
            ( self.tag
            , self.segs_boom
            , 0, 0, -self.dipole_radius
            , 0, 0, -(self.refl_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Reflector
        self.geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , self.reflector, 0, -(self.refl_dist + self.dipole_radius)
            , 0,              0, -(self.refl_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        self.geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , -self.reflector, 0, -(self.refl_dist + self.dipole_radius)
            , 0,               0, -(self.refl_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Turn around Y by 270 deg, move everything up by reflector length
        self.geo.move (0, 270, 0, 0, self.reflector, 0, 0, 0, 0)
        self.nec.geometry_complete (0)
        self.nec.ld_card (5, 0, 0, 0, 37735849, 0, 0)
        self.nec.ex_card (0, self.ex, 1, 0, 1, 0, 0, 0, 0, 0)
        self.nec.fr_card (0, self.frqidxmax, self.frqstart, self.frqinc)
        self.nec.rp_card (0, 37, 73, 0, 0, 0, 0, 0, 0, 5, 5, 0, 0)
    # end def __init__

    @property
    def frqidxrange (self) :
        return range (self.frqidxmax)
    # end def frqidxrange

    def plot (self, frqidx = 100) :
        self.rp = self.nec.get_radiation_pattern (frqidx)

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
        return (1. + gamma) / (1. - gamma)
    # end def vswr

# end class Folded_Dipole


f = Folded_Dipole ()
f.swr_plot ()
