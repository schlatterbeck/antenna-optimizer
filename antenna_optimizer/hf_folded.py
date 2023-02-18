#!/usr/bin/python3
from __future__ import print_function

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions

class Folded_Dipole (Antenna_Model):

    wire_radius   = 1.5e-3 / 2.0
    dipole_dist   =   0.25
    dipole_len    =  10.06
    dipole_height =   5
    segs_dipole   = 121
    segs_end      =   5
    frq_ranges    = [(13.000, 14.350)]
    theta_range   = 90

    def __init__ \
        ( self
        , dipole_dist   = dipole_dist
        , dipole_len    = dipole_len
        , dipole_height = dipole_height
        , wire_radius   = wire_radius
        , ** kw
        ):
        self.dipole_dist   = dipole_dist
        self.dipole_len    = dipole_len
        self.dipole_height = dipole_height
        self.wire_radius   = wire_radius
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return \
            ( "-d %(dipole_dist)1.4f -l %(dipole_len)1.4f "
              "-h %(dipole_height)1.4f"
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
        for y in (0, self.dipole_dist):
            geo.wire \
                ( self.tag
                , self.segs_dipole
                , -self.dipole_len / 2, 0, y
                ,  self.dipole_len / 2, 0, y
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1
        self.ex = Excitation (1, self.segs_dipole // 2)
        for x in (-self.dipole_len / 2, self.dipole_len / 2):
            geo.wire \
                ( self.tag
                , self.segs_end
                , x, 0, 0
                , x, 0, self.dipole_dist
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

# end class Folded_Dipole

def main ():
    cmd = Arg_Handler ()
    cmd.add_argument \
        ( '-d', '--dipole-dist'
        , type = float
        , help = "Distance of folded dipole upper and lower wire"
        , default = Folded_Dipole.dipole_dist
        )
    cmd.add_argument \
        ( '-l', '--dipole-len'
        , type    = float
        , help = "Length of the folded dipole"
        , default = Folded_Dipole.dipole_len
        )
    cmd.add_argument \
        ( '--dipole-height'
        , type    = float
        , help    = "Feedpoint height above ground"
        , default = Folded_Dipole.dipole_height
        )
    args = cmd.parse_args ()
    fd = Folded_Dipole \
        ( dipole_dist   = args.dipole_dist
        , dipole_len    = args.dipole_len
        , dipole_height = args.dipole_height
        , ** cmd.default_antenna_args
        )
    antenna_actions (cmd, args, fd)
# end def main

if __name__ == '__main__':
    main ()
