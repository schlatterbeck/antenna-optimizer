#!/usr/bin/python3
from __future__ import print_function
import sys

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions

class Folded_Dipole (Antenna_Model):
    """ This is a folded dipole with a balcony rail as a large reflector.
        Note that most conductors modelled in the following are not
        actually round. The upper railing modelled as a 20.5mm diameter
        round conductor *is* the only conductor that really is round.
        The upper and lower horizontal conductors are 8.5x30.5mm. The
        left and right vertical conductors are 9x30.5mm. The inner
        vertical conductors are 8.5x25.5mm.
    """

    wire_radius   = 1.5e-3 / 2.0
    dipole_radius = 0.010
    lambda_4      = .146
    segs_dipole   = 19
    segs_arc      = 17
    segs_boom     =  5
    director      = 0.2
    refl_d4       = (53 + 8.5 / 2.0) * 1e-3

    def __init__ \
        ( self
        , refl_dist     = .01
        , dir_dist      = .01
        , ant_h         = ((690 - 8.5) / 2.0) * 1e-3
        , dipole_radius = dipole_radius
        , lambda_4      = lambda_4
        , director      = director
        , wire_radius   = wire_radius
        , refl_radius1  = (8.5  / 2.0) * 1e-3
        , refl_radius2  = (9.0  / 2.0) * 1e-3
        , refl_radius3  = (20.5 / 2.0) * 1e-3
        , refl_d1       = 61 * 1e-3
        , refl_d2       = 150 * 1e-3
        , refl_d3       = (690 - 8.5) * 1e-3
        , refl_d4       = refl_d4
        , refl_d5       = (96.5 + 8.5 / 2.0) * 1e-3
        , ** kw
        ):
        self.refl_dist     = refl_dist
        self.dir_dist      = dir_dist
        self.ant_h         = ant_h
        self.dipole_radius = dipole_radius
        self.lambda_4      = lambda_4
        self.director      = director
        self.wire_radius   = wire_radius
        self.refl_radius1  = refl_radius1
        self.refl_radius2  = refl_radius2
        self.refl_radius3  = refl_radius3
        self.refl_d1       = refl_d1
        self.refl_d2       = refl_d2
        self.refl_d3       = refl_d3
        self.refl_d4       = refl_d4
        self.refl_d5       = refl_d5
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return \
            ( "-r %(dipole_radius)1.4f -d %(refl_dist)1.4f "
              "-l %(director)1.4f -4 %(lambda_4)1.4f "
              "-D %(dir_dist)1.4f -H %(ant_h)1.4f" % self.__dict__
            )
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        self._geometry (geo)
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
        self.ex = Excitation (self.tag - 2, 1)
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
        # second part of boom from dipole to director
        geo.wire \
            ( self.tag
            , self.segs_boom
            , 0, 0, -self.dipole_radius
            , 0, 0, -(self.dir_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Director
        geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , self.director, 0, -(self.dir_dist + self.dipole_radius)
            , 0,             0, -(self.dir_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , -self.director, 0, -(self.dir_dist + self.dipole_radius)
            , 0,              0, -(self.dir_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Turn around Y by 270 deg, move everything in place
        geo.move \
            ( 0, 270, 0
            , self.dipole_radius + self.refl_dist, 0, self.refl_d4 + self.ant_h
            , 0, 0, 0
            )
        for y in (-1, 1):
            z = 0
            geo.wire \
                ( self.tag
                , 3
                , 0, y * 4 * self.refl_d5, z
                , 0, y * 4 * self.refl_d5, z + self.refl_d4
                , self.refl_radius2
                , 1, 1
                )
            self.tag += 1
            z += self.refl_d4
            geo.wire \
                ( self.tag
                , self.segs_dipole
                , 0, y * 4 * self.refl_d5, z
                , 0, y * 4 * self.refl_d5, z + self.refl_d3
                , self.refl_radius2
                , 1, 1
                )
            self.tag += 1
            z += self.refl_d3
            geo.wire \
                ( self.tag
                , 3
                , 0, y * 4 * self.refl_d5, z
                , 0, y * 4 * self.refl_d5, z + self.refl_d2
                , self.refl_radius2
                , 1, 1
                )
            self.tag += 1
            z += self.refl_d2
            geo.wire \
                ( self.tag
                , 3
                , 0, y * 4 * self.refl_d5, z
                , 0, y * 4 * self.refl_d5, z + self.refl_d1
                , self.refl_radius2
                , 1, 1
                )
            self.tag += 1
        geo.wire \
            ( self.tag
            , self.segs_dipole
            , 0, -4 * self.refl_d5, z
            , 0,  4 * self.refl_d5, z
            , self.refl_radius3
            , 1, 1
            )
        self.tag += 1
        for z in [self.refl_d4, self.refl_d4 + self.refl_d3]:
            for y in range (-4, 4):
                geo.wire \
                    ( self.tag
                    , 3
                    , 0, y * self.refl_d5,                z
                    , 0, y * self.refl_d5 + self.refl_d5, z
                    , self.refl_radius2
                    , 1, 1
                    )
                self.tag += 1
        for y in range (-3, 4):
            geo.wire \
                ( self.tag
                , self.segs_dipole
                , 0, y * self.refl_d5, self.refl_d4
                , 0, y * self.refl_d5, self.refl_d4 + self.refl_d3
                , self.refl_radius1
                , 1, 1
                )
            self.tag += 1
    # end def _geometry

    @property
    def up (self):
        """ move everything up by max (director length, lambda_4 + r)
        """
        return max (self.director, self.lambda_4 + self.dipole_radius)
    # end def up

# end class Folded_Dipole

class Folded_Dipole_Optimizer (Antenna_Optimizer):
    """ Optimize given folded dipole
        Length are encoded as integers with a resolution of .5mm
        We use:
        * 8mm  <= dipole_radius <=  5cm
        * 8mm  <= refl_dist     <= 15cm
        * 8mm  <= dir_dist      <= 15cm
        * 10cm <= director      <= 40cm
        * 10cm <= lambda_4      <= 20cm
        * 0mm  <= ant_h         <= ((690 - 8.5) / 2.0)mm
        Or with the large_refldist option we use
        * 2cm  <= refl_dist     <= 25cm
    """
    ant_cls = Folded_Dipole

    def __init__ (self, large_refldist = False, **kw):
        self.ant_height = (690 - 8.5) * 1e-3
        self.large_refldist = large_refldist
        self.minmax = \
            [ (8e-3, 0.05)
            , (8e-3, 0.15)
            , (8e-3, 0.15)
            , (0.1,  0.4)
            , (0.1,  0.2)
            , (0,    self.ant_height)
            ]
        if self.large_refldist:
            self.minmax [1] = (0.02, 0.25)
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop):
        dipole_radius = self.get_parameter (p, pop, 0)
        refl_dist     = self.get_parameter (p, pop, 1)
        dir_dist      = self.get_parameter (p, pop, 2)
        director      = self.get_parameter (p, pop, 3)
        lambda_4      = self.get_parameter (p, pop, 4)
        ant_h         = self.get_parameter (p, pop, 5) + self.ant_cls.refl_d4
        fd = self.ant_cls \
            ( dipole_radius = dipole_radius
            , dir_dist      = dir_dist
            , refl_dist     = refl_dist
            , director      = director
            , lambda_4      = lambda_4
            , ant_h         = ant_h
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
        ( '-D', '--director-distance'
        , type    = float
        , help    = "Distance of the director from nearest dipole part"
        , default = 0.01
        )
    cmd.add_argument \
        ( '--force-reflector'
        , help    = "When optimizing force reflector "
                    " >= 5mm + lambda-len + dipole-radius"
        , action  = "store_true"
        )
    cmd.add_argument \
        ( '-H', '--antenna-height'
        , type    = float
        , help    = "Height of the antenna relative to middle of reflector"
        , default = 0.03
        )
    cmd.add_argument \
        ( '-l', '--director-length'
        , type = float
        , help = "(Half) Length of the director"
        , default = 0.2
        )
    cmd.add_argument \
        ( '--large-refldist'
        , help    = "Use a larger (min and max) reflector distance in optimizer"
        , action  = "store_true"
        )
    cmd.add_argument \
        ( '-r', '--dipole-radius'
        , type    = float
        , help    = "Radius of the rounded corner of the folded dipole"
        , default = 0.01
        )
    cmd.add_argument \
        ( '--refl-radius1'
        , type    = float
        , help    = "Radius of smaller rectangular vertical elements"
        , default = (8.5  / 2.0) * 1e-3
        )
    cmd.add_argument \
        ( '--refl-radius2'
        , type    = float
        , help    = "Radius of larger rectangular vertical elements"
        , default = (9.0  / 2.0) * 1e-3
        )
    cmd.add_argument \
        ( '--refl-radius3'
        , type    = float
        , help    = "Radius of round horizontal element"
        , default = (20.5 / 2.0) * 1e-3
        )
    args = cmd.parse_args ()
    if args.action == 'optimize':
        do = Folded_Dipole_Optimizer \
            ( large_refldist = args.large_refldist
            , ** cmd.default_optimization_args
            )
        do.run ()
    else:
        fd = Folded_Dipole \
            ( dipole_radius = args.dipole_radius
            , refl_dist     = args.reflector_distance
            , dir_dist      = args.director_distance
            , ant_h         = args.antenna_height
            , director      = args.director_length
            , lambda_4      = args.lambda_len
            , refl_radius1  = args.refl_radius1
            , refl_radius2  = args.refl_radius2
            , refl_radius3  = args.refl_radius3
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, fd)
# end def main

if __name__ == '__main__':
    main ()
