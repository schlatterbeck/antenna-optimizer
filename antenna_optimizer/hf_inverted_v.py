#!/usr/bin/python3
from __future__ import print_function

import numpy as np
from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions

class Inverted_V (Antenna_Model):

    wire_radius   = 1.5e-3 / 2.0
    dipole_half   =  5.3
    dipole_height =   5
    dipole_angle  =  90
    seg_len       = 0.05
    frq_ranges    = [(13.000, 14.350)]
    theta_range   = 90

    def __init__ \
        ( self
        , dipole_half   = dipole_half
        , dipole_height = dipole_height
        , dipole_angle  = dipole_angle
        , wire_radius   = wire_radius
        , ** kw
        ):
        self.dipole_half   = dipole_half
        self.dipole_height = dipole_height
        self.dipole_angle  = dipole_angle
        self.wire_radius   = wire_radius
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return \
            ( "--dipole-half %(dipole_half)1.4f "
              "--dipole-height %(dipole_height)1.4f "
              "--dipole-angle %(dipole_angle)1.4f"
            % self.__dict__
            )
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        self._geometry (geo)
        # Move everything up
        geo.move (0, 0, 0, 0, 0, self.dipole_height, 0, 0, 0)
    # end def geometry

    def _geometry (self, geo):
        # Use a horizontal 3-element wire in the middle for excitation
        geo.wire \
                ( self.tag
                , 3
                , -self.seg_len * 3 / 2, 0, 0
                ,  self.seg_len * 3 / 2, 0, 0
                , self.wire_radius
                , 1, 1
                )
        self.ex = Excitation (1, 1)
        self.tag += 1
        for sgn in (-1, 1):
            a = self.dipole_angle / 180 * np.pi
            s = sgn * self.seg_len * 3 / 2
            z = np.cos (a) * self.dipole_half
            x = sgn * np.sin (a) * self.dipole_half
            segs = int (self.dipole_half // self.seg_len)
            geo.wire \
                ( self.tag
                , segs
                , s, 0, 0
                , x, 0, -z
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1
    # end def _geometry

    def ground (self, nec = None):
        """ We currently asume a medium ground with sommerfield-norton
            ground model
        """
        if nec is None:
            nec = self.nec
        # ground-type 2: Sommerfield-Norton
        nec.gn_card (2, 0, 13, 0.005, 0, 0, 0, 0)
    # end def ground

# end class Inverted_V

def main ():
    cmd = Arg_Handler ()
    cmd.add_argument \
        ( '-d', '--dipole-half'
        , type = float
        , help = "Length of one dipole branch"
        , default = Inverted_V.dipole_half
        )
    cmd.add_argument \
        ( '-l', '--dipole-angle'
        , type    = float
        , help = "Angle of one dipole branch against vertical"
        , default = Inverted_V.dipole_angle
        )
    cmd.add_argument \
        ( '--dipole-height'
        , type    = float
        , help    = "Feedpoint height above ground"
        , default = Inverted_V.dipole_height
        )
    args = cmd.parse_args ()
    fd = Inverted_V \
        ( dipole_half   = args.dipole_half
        , dipole_angle  = args.dipole_angle
        , dipole_height = args.dipole_height
        , ** cmd.default_antenna_args
        )
    antenna_actions (cmd, args, fd)
# end def main

if __name__ == '__main__':
    main ()
