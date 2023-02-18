#!/usr/bin/python3
from __future__ import print_function

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions

class Folded_Dipole (Antenna_Model):

    wire_radius   = 1.5e-3 / 2.0
    dipole_radius = 0.010
    lambda_4      = .146
    segs_dipole   = 19
    segs_arc      = 17
    segs_boom     =  5
    reflector     = 0.2

    def __init__ \
        ( self
        , refl_dist     = .01
        , dipole_radius = dipole_radius
        , lambda_4      = lambda_4
        , reflector     = reflector
        , wire_radius   = wire_radius
        , ** kw
        ):
        self.refl_dist     = refl_dist
        self.dipole_radius = dipole_radius
        self.lambda_4      = lambda_4
        self.reflector     = reflector
        self.wire_radius   = wire_radius
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return \
            ("-r %(dipole_radius)1.4f -d %(refl_dist)1.4f " \
             "-l %(reflector)1.4f -4 %(lambda_4)1.4f" % self.__dict__
            )
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        self._geometry (geo)
        # Turn around Y by 270 deg, move everything up
        geo.move (0, 270, 0, 0, 0, self.up, 0, 0, 0)
    # end def geometry

    def _geometry (self, geo):
        a = ((-90, 90), (90, 270))
        for n, z in enumerate ((1, -1)):
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
        for x in (-self.lambda_4, self.lambda_4):
            for z in (self.dipole_radius, -self.dipole_radius):
                geo.wire \
                    ( self.tag
                    , self.segs_dipole
                    , 0, 0, z
                    , x, 0, z
                    , self.wire_radius
                    , 1, 1
                    )
                self.tag += 1
        self.ex = Excitation (self.tag - 1, 1)
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
    # end def _geometry

    @property
    def up (self):
        """ move everything up by max (reflector length, lambda_4 + r)
        """
        return max (self.reflector, self.lambda_4 + self.dipole_radius)
    # end def up

# end class Folded_Dipole

class Folded_Dipole_Optimizer (Antenna_Optimizer):
    """ Optimize given folded dipole
        Length are encoded as integers with a resolution of .5mm
        We use:
        * 8mm  <= dipole_radius <=  5cm
        * 8mm  <= refl_dist     <= 10cm
        * 10cm <= reflector     <= 40cm
        * 10cm <= lambda_4      <= 20cm
        If we force a reflector, we use
        reflector = dipole_radius + lambda_4 + x
        with 5mm <= x <= 10cm
    """
    # The antenna class
    ant_cls  = Folded_Dipole

    def __init__ (self, force_reflector = False, **kw):
        self.force_reflector = force_reflector
        self.minmax = [(8e-3, 0.05), (8e-3, 0.1), (0.1, 0.4), (0.1, 0.2)]
        if self.force_reflector:
            self.minmax [2] = (5e-3, 0.1)
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop):
        dipole_radius = self.get_parameter (p, pop, 0)
        refl_dist     = self.get_parameter (p, pop, 1)
        reflector     = self.get_parameter (p, pop, 2)
        lambda_4      = self.get_parameter (p, pop, 3)
        if self.force_reflector:
            reflector += dipole_radius + lambda_4
        fd = self.ant_cls \
            ( dipole_radius = dipole_radius
            , refl_dist     = refl_dist
            , reflector     = reflector
            , lambda_4      = lambda_4
            , **self.antenna_args
            )
        return fd
    # end def compute_antenna

# end class Folded_Dipole_Optimizer

def main ():
    cmd = Arg_Handler ()
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
        ( '--force-reflector'
        , help    = "When optimizing force reflector "
                    " >= 5mm + lambda-len + dipole-radius"
        , action  = "store_true"
        )
    cmd.add_argument \
        ( '-l', '--reflector-length'
        , type = float
        , help = "(Half) Length of the reflector"
        , default = 0.2
        )
    cmd.add_argument \
        ( '-r', '--dipole-radius'
        , type    = float
        , help    = "Radius of the rounded corner of the folded dipole"
        , default = 0.01
        )
    args = cmd.parse_args ()
    if args.action == 'optimize':
        do = Folded_Dipole_Optimizer \
            ( force_reflector = args.force_reflector
            , ** cmd.default_optimization_args
            )
        do.run ()
    else:
        fd = Folded_Dipole \
            ( dipole_radius = args.dipole_radius
            , refl_dist     = args.reflector_distance
            , reflector     = args.reflector_length
            , lambda_4      = args.lambda_len
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, fd)
# end def main

if __name__ == '__main__':
    main ()
