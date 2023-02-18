#!/usr/bin/python3
from __future__ import print_function

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions

class Folded_Dipole (Antenna_Model):
    """ Broadcast Reception Antenna
        For 87.5MHz - 108 MHz and 75 Ohm
        Note that with a single dipole this is too narrowband to work
        reliably.
    """

    wire_radius   = 1.5e-3 / 2.0
    dipole_radius = 0.010
    lambda_4      = .146
    segs_dipole   = 19
    segs_arc      = 17
    segs_boom     =  5
    reflector     = 0.2
    frq_ranges    = [(87.5, 108.0)]

    def __init__ \
        ( self
        , refl_dist     = .01
        , dipole_radius = dipole_radius
        , lambda_4      = lambda_4
        , wire_radius   = wire_radius
        , impedance     = 75.0
        , use_boom      = False
        , ** kw
        ):
        self.dipole_radius = dipole_radius
        self.lambda_4      = lambda_4
        if self.lambda_4 < 0.0005:
            self.lambda_4 = 0.0
        self.wire_radius   = wire_radius
        self.impedance     = impedance
        self.use_boom      = use_boom
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return ("-r %(dipole_radius)1.4f -4 %(lambda_4)1.4f" % self.__dict__)
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        self._geometry (geo)
        # Move everything up
        geo.move (0, 0, 90, 0, 0, self.up, 0, 0, 0)
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
            roundtag = self.tag
            self.tag += 1
        # Interpolate the number of segments for very small lambda_4
        if self.lambda_4 > .4:
            segs = self.segs_dipole
        elif self.lambda_4 <= 0.05:
            segs = 3
        else:
            segs = self.lambda_4 * (self.segs_dipole - 3) / 0.35 \
                 + (24 - self.segs_dipole) / 7.0
            segs = int (segs)
            if segs % 2 == 0:
                segs += 1
        if self.lambda_4 > 0.0:
            for z in (self.dipole_radius, -self.dipole_radius):
                if self.use_boom:
                    for x in (-self.lambda_4, self.lambda_4):
                        geo.wire \
                            ( self.tag
                            , segs
                            , 0, 0, z
                            , x, 0, z
                            , self.wire_radius
                            , 1, 1
                            )
                        self.tag += 1
                else:
                    geo.wire \
                        ( self.tag
                        , segs
                        , -self.lambda_4, 0, z
                        ,  self.lambda_4, 0, z
                        , self.wire_radius
                        , 1, 1
                        )
                    self.tag += 1
        if self.use_boom:
            self.ex = Excitation (self.tag - 1, 1)
        elif self.lambda_4 >= 0.002:
            self.ex = Excitation (self.tag - 1, segs // 2 + 1)
        else:
            # If straight piece is too small use last round segment
            self.ex = Excitation (roundtag, self.segs_arc)
        # boom across folded part
        if self.use_boom:
            geo.wire \
                ( self.tag
                , self.segs_boom
                , 0, 0,  self.dipole_radius
                , 0, 0, -self.dipole_radius
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1
    # end def _geometry

    @property
    def up (self):
        """ move everything up by max (reflector length, lambda_4 + r)
        """
        return self.lambda_4 + self.dipole_radius
    # end def up

# end class Folded_Dipole

class Folded_Dipole_Optimizer (Antenna_Optimizer):
    """ Optimize given folded dipole
        Length are encoded as integers with a resolution of .5mm
        We use:
        * 8mm  <= dipole_radius <= 40cm
        * 40cm <= lambda_4      <= 150cm
    """
    ant_cls = Folded_Dipole

    def __init__ \
        (self, impedance = 75.0, allow_loop = False, use_boom = False, **kw):
        self.allow_loop = allow_loop
        self.impedance  = impedance
        self.use_boom   = use_boom
        self.minmax = [(8e-3, 0.4), (0.4, 1.5)]
        if allow_loop:
            self.minmax = [(8e-3, 0.75), (0.00, 1.5)]
        self.__super.__init__ (nofb = True, **kw)
    # end def __init__

    def compute_antenna (self, p, pop):
        dipole_radius = self.get_parameter (p, pop, 0)
        lambda_4      = self.get_parameter (p, pop, 1)
        fd = self.ant_cls \
            ( dipole_radius = dipole_radius
            , lambda_4      = lambda_4
            , impedance     = self.impedance
            , use_boom      = self.use_boom
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
        ( '-r', '--dipole-radius'
        , type    = float
        , help    = "Radius of the rounded corner of the folded dipole"
        , default = 0.01
        )
    cmd.add_argument \
        ( '--allow-loop'
        , help    = "Allow the antenna to degenerate into a loop"
        , action  = 'store_true'
        )
    cmd.add_argument \
        ( '--impedance'
        , type    = float
        , help    = "Impedance of antenna (and ideally feed wire)"
        , default = 75.0
        )
    cmd.add_argument \
        ( '--use-boom'
        , help    = "Use a boom in the middle of the structure"
        , action  = 'store_true'
        )
    args = cmd.parse_args ()
    if args.action == 'optimize':
        do = Folded_Dipole_Optimizer \
            ( allow_loop = args.allow_loop
            , impedance  = args.impedance
            , use_boom   = args.use_boom
            , ** cmd.default_optimization_args
            )
        do.run ()
    else:
        fd = Folded_Dipole \
            ( dipole_radius = args.dipole_radius
            , lambda_4      = args.lambda_len
            , impedance     = args.impedance
            , use_boom      = args.use_boom
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, fd)
# end def main

if __name__ == '__main__':
    main ()
