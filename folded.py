from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

import math
import PyNEC

class Folded_Dipole (object) :

    wire_radius   = 1.5e-3
    dipole_radius = 2.2225e-2
    lambda_4      = 1.117475e-1
    segs_dipole   = 19
    segs_arc      = 17
    segs_boom     =  5
    reflector     = 0.3

    def __init__ (self, refl_dist = 7.77875e-2) :
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
            self.geo.move (0, 0, 0, z * self.lambda_4, 0, 0, 0, 0, self.tag)
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
                if self.ex is None :
                    self.ex = self.tag
                self.tag += 1
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
            , 0, 0, -self.refl_dist
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Reflector
        self.geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , self.reflector, 0, -self.refl_dist
            , 0,              0, -self.refl_dist
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        self.geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , -self.reflector, 0, -self.refl_dist
            , 0,               0, -self.refl_dist
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Turn around Y by 270 deg, move everything up by reflector length
        self.geo.move (0, 270, 0, 0, self.reflector, 0, 0, 0, 0)
        self.nec.geometry_complete (0)
        self.nec.ex_card (0, self.ex, 1, 0, 1, 0, 0, 0, 0, 0)
        self.nec.fr_card (0, 301, 4.550E+02, 0.1)
        self.nec.rp_card (0, 37, 73, 0, 0, 0, 0, 0, 0, 5, 5, 0, 0)
        #self.nec.rp_card (0, 7, 13, 0, 0, 0, 0, 0, 0, 30, 30, 0, 0)
        #self.nec.rp_card (0, 13, 13, 0, 0, 0, 0, 0, 0, 15, 30, 0, 0)
        #self.nec.rp_card (0, 7, 7, 0, 0, 0, 0, 0, 0, 30, 60, 0, 0)
        self.rp = self.nec.get_radiation_pattern (150)
    # end def __init__

# end class Folded_Dipole

f = Folded_Dipole ()
gains  = f.rp.get_gain ()
gains  = 10.0 ** (gains / 10.0)
thetas = f.rp.get_theta_angles () * 3.1415 / 180.0
phis   = f.rp.get_phi_angles ()   * 3.1415 / 180.0

P, T = np.meshgrid (phis, thetas)

X = np.cos (P) * np.sin (T) * gains
Y = np.sin (P) * np.sin (T) * gains
Z = np.cos (T) * gains

fig = plt.figure ()
ax  = fig.gca (projection='3d')

# Create cubic bounding box to simulate equal aspect ratio
max_range = np.array \
    ([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max() / 2.0

mid_x = (X.max () + X.min ()) * 0.5
mid_y = (Y.max () + Y.min ()) * 0.5
mid_z = (Z.max () + Z.min ()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

ax.plot_wireframe (X, Y, Z, color = 'r')
plt.show ()
