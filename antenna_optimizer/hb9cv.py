#!/usr/bin/python3
from __future__ import print_function

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions
from .transmission  import transmission_line_z

class HB9CV (Antenna_Model):
    """ The HB9CV antenna is a two-element yagi-uda antenna.
        It uses transmission lines for off-center feeding.
        The feedpoint is between the middle of the reflector and two
        connected transmission lines. When viewed with the reflector at
        the front and the director at the back, the length l5 is the
        distance from the middle if the reflector to the feedpoint on
        the right of the reflector. The length l4 is the distance from
        the middle of the director to the left feedpoint on the
        director.
    """

    director      = 317.2e-3
    reflector     = 344.8e-3
    refl_dist     = 86.2e-3
    l4            = 43.1e-3
    l5            = 46.55e-3
    stub_height   = 5e-3
    segs_dipole   = 19
    segs_stub     =  1 # When using transmission line we use 1 here
    segs_boom     =  5
    segs_h        =  5
    geotypes      = 'transmission-1 transmission-2 parallel'.split ()
    vf            = 0.9

    def __init__ \
        ( self
        , director      = director
        , reflector     = reflector
        , refl_dist     = refl_dist
        , l4            = l4
        , l5            = l5
        , stub_height   = stub_height
        , geotype       = geotypes [0]
        , vf            = vf
        , ** kw
        ):
        self.director    = director
        self.reflector   = reflector
        self.refl_dist   = refl_dist
        self.l4          = l4
        self.l5          = l5
        self.stub_height = stub_height
        self.geotype     = geotype
        self.vf          = vf
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        return \
            ("-l %(director)1.4f -r %(reflector)1.4f " \
             "-d %(refl_dist)1.4f -4 %(l4)1.4f -5 %(l5)1.4f "
             "-H %(stub_height)1.4f" % self.__dict__
            )
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None

        # Director
        geo.wire \
            ( self.tag
            , self.segs_dipole
            , self.refl_dist,                    0, 0
            , self.refl_dist, -self.director / 2.0, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Second half of Director needs to be modelled in two parts,
        # because stub connection must be at end of wire(s)
        geo.wire \
            ( self.tag
            , self.segs_stub
            , self.refl_dist,       0, 0
            , self.refl_dist, self.l4, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        if self.director / 2.0 > self.l4:
            geo.wire \
                ( self.tag
                , self.segs_stub
                , self.refl_dist,             self.l4, 0
                , self.refl_dist, self.director / 2.0, 0
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1

        # Reflector
        geo.wire \
            ( self.tag
            , self.segs_dipole
            , 0,                    0, 0
            , 0, self.reflector / 2.0, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Second half of Reflector needs to be modelled in two parts,
        # because stub connection must be between end of wire(s).
        geo.wire \
            ( self.tag
            , self.segs_stub
            , 0,        0, 0
            , 0, -self.l5, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        if self.reflector / 2.0 > self.l5:
            geo.wire \
                ( self.tag
                , self.segs_stub
                , 0,              -self.l5, 0
                , 0, -self.reflector / 2.0, 0
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1
        # Boom
        geo.wire \
            ( self.tag
            , self.segs_boom
            , 0,              0, 0
            , self.refl_dist, 0, 0
            , self.boom_radius
            , 1, 1
            )
        self.tag += 1
        # Single segment wire from boom to phasing on director side
        # Used for feed
        geo.wire \
            ( self.tag
            , 1
            , self.refl_dist, 0, 0
            , self.refl_dist, 0, self.stub_height
            , self.wire_radius
            , 1, 1
            )
        self.ex = Excitation (self.tag, 1)
        self.tag  += 1
        # Connect Director Stub to Director
        geo.wire \
            ( self.tag
            , self.segs_stub
            , self.refl_dist, self.l4,                0
            , self.refl_dist, self.l4, self.stub_height
            , self.wire_radius
            , 1, 1
            )
        self.director_stub_conn_tag = self.tag
        self.tag  += 1
        # Connect Reflector Stub to Reflector
        geo.wire \
            ( self.tag
            , self.segs_stub
            , 0, -self.l5,                0
            , 0, -self.l5, self.stub_height
            , self.wire_radius
            , 1, 1
            )
        self.reflector_stub_conn_tag = self.tag
        self.tag  += 1

        if self.geotype == 'transmission-1':
            # Wire from lower part of feedpoint to upper part of
            # director stub: Used as transmission-line endpoint
            h = 0
            if self.geotype == 'parallel':
                h = self.stub_height
            geo.wire \
                ( self.tag
                , 1
                , self.refl_dist,        0, 0
                ,              0, -self.l5, self.stub_height
                , self.wire_radius
                , 1, 1
                )
            self.reflector_stub_tag = self.tag
            self.tag += 1

        if self.geotype in ('transmission-1', 'parallel'):
            # Wire from lower (or upper depending on geotype) part of
            # feedpoint to upper part of director stub: Used as
            # transmission-line endpoint or parallel wire endpoint
            h = 0
            if self.geotype == 'parallel':
                h = self.stub_height

            # Wire from lower (or upper depending on geotype) part of
            # feedpoint to upper part of reflector stub: Used as
            # transmission-line endpoint or parallel wire endpoint
            geo.wire \
                ( self.tag
                , 1
                , self.refl_dist,       0, h
                , self.refl_dist, self.l4, self.stub_height
                , self.wire_radius
                , 1, 1
                )
            self.director_stub_tag = self.tag
            self.tag += 1

        if self.geotype == 'parallel':
            geo.wire \
                ( self.tag
                , 1
                , self.refl_dist, 0, self.stub_height
                ,              0, 0, self.stub_height
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1
            geo.wire \
                ( self.tag
                , 1
                , 0,        0, self.stub_height
                , 0, -self.l5, self.stub_height
                , self.wire_radius
                , 1, 1
                )
            self.tag += 1
    # end def geometry


    def geometry_complete (self, nec = None):
        if nec is None:
            nec = self.nec
        nec.geometry_complete (0)
        impedance = transmission_line_z \
            (self.wire_radius * 2, self.stub_height)
        if self.geotype == 'transmission-1':
            nec.tl_card \
                ( self.ex.tag, self.ex.segment
                , self.reflector_stub_tag, 1
                , impedance
                , (self.refl_dist + self.l5) * self.vf
                , 0, 0, 0, 0
                )
            nec.tl_card \
                ( self.ex.tag, self.ex.segment
                , self.director_stub_tag, 1
                , impedance
                , self.l4 * self.vf
                , 0, 0, 0, 0
                )
        if self.geotype == 'transmission-2':
            nec.tl_card \
                ( self.ex.tag, self.ex.segment
                , self.reflector_stub_conn_tag, 1
                , impedance
                , (self.refl_dist + self.l5) * self.vf
                , 0, 0, 0, 0
                )
            nec.tl_card \
                ( self.ex.tag, self.ex.segment
                , self.director_stub_conn_tag, 1
                , impedance
                , self.l4 * self.vf
                , 0, 0, 0, 0
                )
    # end def geometry_complete

# end class HB9CV

class HB9CV_Optimizer (Antenna_Optimizer):
    """ Optimize given folded dipole
        Length are encoded with a resolution of .5mm
        We use:
        * 25cm   <= director    <= 35cm
        * 25cm   <= reflector   <= 35cm
        *  5cm   <= refl_dist   <= 15cm
        *  3cm   <= l4          <= 10cm
        *  3cm   <= l5          <= 10cm
        *  5mm   <= stub_height <= 1.5cm
    """
    ant_cls = HB9CV

    def __init__ (self, vf = 0.9, **kw):
        self.minmax = \
            [ (0.25, 0.35), (0.25, 0.35), (0.05, 0.15)
            , (0.03, 0.1),  (0.03, 0.1),  (0.005, 0.015)
            ]
        self.vf = vf
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop):
        director      = self.get_parameter (p, pop, 0)
        reflector     = self.get_parameter (p, pop, 1)
        refl_dist     = self.get_parameter (p, pop, 2)
        l4            = self.get_parameter (p, pop, 3)
        l5            = self.get_parameter (p, pop, 4)
        h             = self.get_parameter (p, pop, 5)
        fd = self.ant_cls \
            ( director      = director
            , reflector     = reflector
            , refl_dist     = refl_dist
            , l4            = l4
            , l5            = l5
            , stub_height   = h
            , vf            = self.vf
            , **self.antenna_args
            )
        return fd
    # end def compute_antenna

# end class HB9CV_Optimizer

def main ():
    cmd = Arg_Handler ()
    cmd.add_argument \
        ( '-4', '--l4'
        , type    = float
        , help    = "Length of l4 (stub on director side)"
        , default = 0.0431
        )
    cmd.add_argument \
        ( '-5', '--l5'
        , type    = float
        , help    = "Length of l5 (stub on reflector side)"
        , default = 0.04655
        )
    cmd.add_argument \
        ( '-d', '--reflector-distance'
        , type    = float
        , help    = "Distance of the reflector from nearest dipole part"
        , default = 0.0862
        )
    cmd.add_argument \
        ( '-g', '--geotype'
        , help    = "Geometry type, one of %s" % ','.join (HB9CV.geotypes)
        , default = HB9CV.geotypes [0]
        )
    cmd.add_argument \
        ( '-H', '--stub-height'
        , type    = float
        , help    = "Height of stubs above antenna"
        , default = 0.017
        )
    cmd.add_argument \
        ( '-l', '--director-length'
        , type    = float
        , help    = "Length of the director"
        , default = 0.3172
        )
    cmd.add_argument \
        ( '-r', '--reflector-length'
        , type    = float
        , help    = "Length of the reflector"
        , default = 0.3448
        )
    cmd.add_argument \
        ( '--vf'
        , type    = float
        , help    = "Velocity factor of transmission line"
        , default = 0.9
        )
    args = cmd.parse_args ()
    if args.action == 'optimize':
        do = HB9CV_Optimizer (vf = args.vf, ** cmd.default_optimization_args)
        do.run ()
    else:
        fd = HB9CV \
            ( director      = args.director_length
            , reflector     = args.reflector_length
            , refl_dist     = args.reflector_distance
            , l4            = args.l4
            , l5            = args.l5
            , stub_height   = args.stub_height
            , geotype       = args.geotype
            , vf            = args.vf
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, fd)
# end def main

if __name__ == '__main__':
    main ()
