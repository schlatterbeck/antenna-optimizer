#!/usr/bin/python3
from __future__ import print_function

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions

class Fuchs_Antenna (Antenna_Model):
    """ End fed antenna with a resonant L/C circuit for matching.
        First published by J. Fuchs in 1927 [1], [2].
        [1] Josef Fuchs. Sendeanordnung für drahtlose Telegraphie.
            Patent AT110357, Österreichisches Patentamt, June 1927.
        [2] Josef Fuchs. Tests on a method of voltage feeding the
            antenna. QST, pages 37, 42, July 1928.
        Note that we're modeling this in free space, so no capacity to
        ground can distort the results :-)
    """

    wire_radius   = 1.25e-3
    dipole_half   = 5.125
    coil_radius   = 0.05
    coil_pitch    = 0.02 # height of one winding
    primary_n     = 9
    secondary_n   = 2
    capacity      = 30e-12

    seg_len       = 0.1
    frq_ranges    = [(14.000, 14.350)]

    def __init__ \
        ( self
        , dipole_half   = dipole_half
        , coil_radius   = coil_radius
        , coil_pitch    = coil_pitch
        , primary_n     = primary_n
        , secondary_n   = secondary_n
        , capacity      = capacity
        , ** kw
        ):
        self.dipole_half = dipole_half
        self.coil_radius = coil_radius
        self.coil_pitch  = coil_pitch
        self.primary_n   = primary_n
        self.secondary_n = secondary_n
        self.capacity    = capacity
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return \
            ("--dipole-half %(dipole_half).4f "
             "--coil-radius %(coil_radius).4f "
             "--coil-pitch %(coil_pitch).4f "
             "--primary-n %(primary_n)d "
             "--secondary-n %(secondary_n)d "
             "--capacity %(capacity)6e "
            % self.__dict__
            )
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1

        # Dipole Element
        geo.wire \
            ( self.tag
            , round (self.dipole_half * 2 / self.seg_len)
            , 0, 0, 0
            , 2 * self.dipole_half, 0, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        h1 = self.coil_pitch * self.primary_n
        geo.helix \
            ( self.tag
            , 10 * self.primary_n
            , self.coil_pitch
            , h1
            , self.coil_radius
            , self.coil_radius
            , self.coil_radius
            , self.coil_radius
            , self.wire_radius
            )
        geo.move (0, 0, 0, -self.coil_radius, 0, 0, self.tag, 0, 0)
        self.tag += 1
        h2 = self.coil_pitch * self.secondary_n
        hm = round ((self.primary_n - self.secondary_n) / 2) * self.coil_pitch
        geo.helix \
            ( self.tag
            , 10 * self.secondary_n
            , self.coil_pitch
            , h2
            , self.coil_radius
            , self.coil_radius
            , self.coil_radius
            , self.coil_radius
            , self.wire_radius
            )
        geo.move (0, 0, 180, -self.coil_radius, 0, hm, self.tag, 0, 0)
        self.tag += 1
        geo.wire \
            ( self.tag
            , 1
            , 0, 0,  0
            , 0, 0, h1
            , self.wire_radius
            , 1, 1
            )
        self.cap_tag = self.tag
        self.tag += 1
        geo.wire \
            ( self.tag
            , 1
            , -2 * self.coil_radius, 0, hm
            , -2 * self.coil_radius, 0, hm + h2
            , self.wire_radius
            , 1, 1
            )
        self.ex = Excitation (self.tag, 1)
    # end def geometry

    def nec_params_compute (self, nec = None):
        if nec is None:
            nec = self.nec
        self.__super.nec_params_compute (nec)
        nec.ld_card (0, self.cap_tag, 1, 1, 0, 0, self.capacity)
    # end def nec_params_compute

# end class Fuchs_Antenna

class Fuchs_Optimizer (Antenna_Optimizer):
    """ Optimize given Fuchs antenna
        We use:
        * 35mm   <= coil_radius <= 100mm
        *  6mm   <= coil_pitch  <=  30mm
        *  5     <= primary_n   <= 30  
        *  1     <= secondary_n <= 10
        *  5pf   <= capacity    <= 100pf
    """
    ant_cls = Fuchs_Antenna

    def __init__ (self, **kw):
        self.minmax = \
            [(0.035, 0.100), (0.006, 0.030), (5, 30.9), (1, 10.9),  (5, 100)]
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop):
        coil_radius   = self.get_parameter (p, pop, 0)
        coil_pitch    = self.get_parameter (p, pop, 1)
        primary_n     = int (self.get_parameter (p, pop, 2))
        secondary_n   = int (self.get_parameter (p, pop, 3))
        capacity      = self.get_parameter (p, pop, 4)
        ant = self.ant_cls \
            ( coil_radius   = coil_radius
            , coil_pitch    = coil_pitch
            , primary_n     = primary_n
            , secondary_n   = secondary_n
            , capacity      = capacity * 1e-12
            , **self.antenna_args
            )
        return ant
    # end def compute_antenna

# end class Fuchs_Optimizer

def main ():
    cmd = Arg_Handler \
        ( wire_radius  = Fuchs_Antenna.wire_radius
        , frq_step_max = 11
        )
    cmd.add_argument \
        ( '--capacity'
        , type    = float
        , help    = "Capacity for matching circuit (in F), default=%(default)s"
        , default = Fuchs_Antenna.capacity
        )
    cmd.add_argument \
        ( '--coil-radius'
        , type    = float
        , help    = "Radius of coil for matching circuit, default=%(default)s"
        , default = Fuchs_Antenna.coil_radius
        )
    cmd.add_argument \
        ( '--coil-pitch'
        , type    = float
        , help    = "Height of one winding of coil, default=%(default)s"
        , default = Fuchs_Antenna.coil_pitch
        )
    cmd.add_argument \
        ( '--primary-n'
        , type    = int
        , help    = "Number of primary windings, default=%(default)s"
        , default = Fuchs_Antenna.primary_n
        )
    cmd.add_argument \
        ( '--secondary-n'
        , type    = int
        , help    = "Number of secondary windings, default=%(default)s"
        , default = Fuchs_Antenna.secondary_n
        )
    args = cmd.parse_args ()
    if args.action == 'optimize':
        ao = Fuchs_Optimizer (** cmd.default_optimization_args)
        ao.run ()
    else:
        ant = Fuchs_Antenna \
            ( coil_radius   = args.coil_radius
            , coil_pitch    = args.coil_pitch
            , primary_n     = args.primary_n
            , secondary_n   = args.secondary_n
            , capacity      = args.capacity
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, ant)
# end def main

if __name__ == '__main__':
    main ()
