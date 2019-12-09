#!/usr/bin/python3
from __future__ import print_function

from argparse import ArgumentParser
from antenna_model import Antenna_Model, Antenna_Optimizer

class HB9CV (Antenna_Model) :

    director      = 317.2e-3
    reflector     = 344.8e-3
    refl_dist     = 86.2e-3
    l4            = 43.1e-3
    l5            = 46.55e-3
    stub_height   = 5e-3
    segs_dipole   = 19
    segs_stub     = 11
    segs_boom     =  5
    segs_h        =  5

    def __init__ \
        ( self
        , director      = director
        , reflector     = reflector
        , refl_dist     = refl_dist
        , l4            = l4
        , l5            = l5
        , stub_height   = stub_height
        , ** kw
        ) :
        self.director    = director
        self.reflector   = reflector
        self.refl_dist   = refl_dist
        self.l4          = l4
        self.l5          = l5
        self.stub_height = stub_height
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self) :
        return \
            ("-l %(director)1.4f -r %(reflector)1.4f " \
             "-d %(refl_dist)1.4f -4 %(l4)1.4f -5 %(l5)1.4f "
             "-H %(stub_height)1.4f" % self.__dict__
            )
    # end def cmdline

    def geometry (self, geo = None) :
        self.tag = 1
        self.ex  = None
        if geo is None :
            geo = self.nec.get_geometry ()

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
        if self.director / 2.0 > self.l4 :
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
        # because stub connection must be at end of wire(s)
        geo.wire \
            ( self.tag
            , self.segs_stub
            , 0,        0, 0
            , 0, -self.l5, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        if self.reflector / 2.0 > self.l5 :
            geo.wire \
                ( self.tag
                , self.segs_stub
                , 0,             -self.l5, 0
                , 0, -self.director / 2.0, 0
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
        # Phasing
        # Phasing parallel to boom
        geo.wire \
            ( self.tag
            , self.segs_boom
            , 0,              0, self.stub_height
            , self.refl_dist, 0, self.stub_height
            , self.wire_radius
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
        self.ex    = self.tag
        self.exseg = 1
        self.tag  += 1
        # Stub on director
        geo.wire \
            ( self.tag
            , self.segs_stub
            , self.refl_dist,       0, self.stub_height
            , self.refl_dist, self.l4, self.stub_height
            , self.wire_radius
            , 1, 1
            )
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
        self.tag  += 1
        # Stub on reflector
        geo.wire \
            ( self.tag
            , self.segs_stub
            , 0,        0, self.stub_height
            , 0, -self.l5, self.stub_height
            , self.wire_radius
            , 1, 1
            )
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
        self.tag  += 1

    # end def geometry

# end class HB9CV

class HB9CV_Optimizer (Antenna_Optimizer) :
    """ Optimize given folded dipole
        Length are encoded with a resolution of .5mm
        We use:
        * 25cm   <= director    <= 35cm
        * 25cm   <= reflector   <= 35cm
        *  5cm   <= refl_dist   <= 15cm
        *  3cm   <= l4          <= 10cm
        *  3cm   <= l5          <= 10cm
        *  1.7cm <= stub_height <= 4cm
    """

    def __init__ (self, **kw) :
        self.minmax = \
            [ (0.25, 0.35), (0.25, 0.35), (0.05, 0.15)
            , (0.03, 0.1),  (0.03, 0.1),  (0.017, 0.04)
            ]
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop) :
        director      = self.get_parameter (p, pop, 0)
        reflector     = self.get_parameter (p, pop, 1)
        refl_dist     = self.get_parameter (p, pop, 2)
        l4            = self.get_parameter (p, pop, 3)
        l5            = self.get_parameter (p, pop, 4)
        h             = self.get_parameter (p, pop, 5)
        fd = HB9CV \
            ( director      = director
            , reflector     = reflector
            , refl_dist     = refl_dist
            , l4            = l4
            , l5            = l5
            , stub_height   = h
            , frqidxmax     = 3
            , wire_radius   = self.wire_radius
            )
        return fd
    # end def compute_antenna

# end class HB9CV_Optimizer

if __name__ == '__main__' :
    actions = ['optimize', 'necout', 'swr', 'gain', 'frgain']
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( 'action'
        , help = "Action to perform, one of %s" % ', '.join (actions)
        )
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
        ( '-H', '--stub-height'
        , type    = float
        , help    = "Height of stubs above antenna"
        , default = '0.017'
        )
    cmd.add_argument \
        ( '-l', '--director-length'
        , type    = float
        , help    = "Length of the director"
        , default = 0.3172
        )
    cmd.add_argument \
        ( '-R', '--random-seed'
        , type    = int
        , help    = "Random number seed for optimizer, default=%(default)s"
        , default = 42
        )
    cmd.add_argument \
        ( '-r', '--reflector-length'
        , type    = float
        , help    = "Length of the reflector"
        , default = 0.3448
        )
    cmd.add_argument \
        ( '-w', '--wire-radius'
        , type    = float
        , help    = "Radius of the wire"
        , default = 0.00075
        )
    cmd.add_argument \
        ( '-v', '--verbose'
        , help    = "Verbose reporting in every generation"
        , action  = 'store_true'
        )
    args = cmd.parse_args ()
    if args.action == 'optimize' :
        do = HB9CV_Optimizer \
            ( verbose     = args.verbose
            , random_seed = args.random_seed
            , wire_radius = args.wire_radius
            )
        do.run ()
    else :
        frqidxmax = 201
        if args.action in ('frgain', 'necout') :
            frqidxmax = 3
        fd = HB9CV \
            ( director      = args.director_length
            , reflector     = args.reflector_length
            , refl_dist     = args.reflector_distance
            , l4            = args.l4
            , l5            = args.l5
            , stub_height   = args.stub_height
            , wire_radius   = args.wire_radius
            , frqidxmax     = frqidxmax
            )
        if args.action == 'necout' :
            print (fd.as_nec ())
        elif args.action not in actions :
            cmd.print_usage ()
        else :
            fd.compute ()
        if args.action == 'swr' :
            fd.swr_plot ()
        elif args.action == 'gain' :
            fd.plot ()
        elif args.action == 'frgain' :
            print ('\n'.join (fd.show_gains ()))
