#!/usr/bin/python3
from __future__ import print_function

from antenna_model import Antenna_Model, Antenna_Optimizer, default_args
from folded        import Folded_Dipole

class Folded_Dipole_3el (Folded_Dipole) :

    wire_radius   = 1.5e-3 / 2.0
    dipole_radius = 0.010
    lambda_4      = .146
    segs_dipole   = 19
    segs_arc      = 17
    segs_boom     =  5
    reflector     = 0.2
    director      = 0.2

    def __init__ \
        ( self
        , refl_dist     = .01
        , dir_dist      = .01
        , dipole_radius = dipole_radius
        , lambda_4      = lambda_4
        , reflector     = reflector
        , director      = director
        , wire_radius   = wire_radius
        , ** kw
        ) :
        self.dir_dist      = dir_dist
        self.director      = director
        self.__super.__init__ \
            ( refl_dist     = refl_dist
            , dipole_radius = dipole_radius
            , lambda_4      = lambda_4
            , reflector     = reflector
            , wire_radius   = wire_radius
            , ** kw
            )
    # end def __init__

    def cmdline (self) :
        return \
            ("-r %(dipole_radius)1.4f -d %(refl_dist)1.4f "
             "-D %(dir_dist)1.4f -L %(director)1.4f "
             "-l %(reflector)1.4f -4 %(lambda_4)1.4f" % self.__dict__
            )
    # end def cmdline

    def _geometry (self, geo) :
        self.__super._geometry (geo)
        # third part of boom from dipole to director
        geo.wire \
            ( self.tag
            , self.segs_boom
            , 0, 0, self.dipole_radius
            , 0, 0, (self.dir_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # Director
        geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , self.director, 0, (self.dir_dist + self.dipole_radius)
            , 0,             0, (self.dir_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , self.segs_dipole * 3
            , -self.director, 0, (self.dir_dist + self.dipole_radius)
            , 0,              0, (self.dir_dist + self.dipole_radius)
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
    # end def _geometry

    @property
    def up (self) :
        return max \
            (self.reflector, self.director, self.lambda_4 + self.dipole_radius)
    # end def up

# end class Folded_Dipole_3el

class Folded_Dipole_Optimizer (Antenna_Optimizer) :
    """ Optimize given folded dipole
        Length are encoded as integers with a resolution of .5mm
        We use:
        * 8mm  <= dipole_radius <=  5cm
        * 35mm <= refl_dist     <= 30cm
        * 35mm <= dir_dist      <= 30cm
        * 10cm <= director      <= 40cm
        * 1mm  <= refldiff      <= 15cm
        * 10cm <= lambda_4      <= 20cm
    """

    def __init__ (self, **kw) :
        self.minmax = \
            [ (8e-3,  0.05)
            , (0.035, 0.3)
            , (0.035, 0.3)
            , ( 0.1,  0.4)
            , (1e-3,  0.15)
            , ( 0.1,  0.2)
            ]
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop) :
        dipole_radius = self.get_parameter (p, pop, 0)
        refl_dist     = self.get_parameter (p, pop, 1)
        dir_dist      = self.get_parameter (p, pop, 2)
        director      = self.get_parameter (p, pop, 3)
        refldiff      = self.get_parameter (p, pop, 4)
        lambda_4      = self.get_parameter (p, pop, 5)
        fd = Folded_Dipole_3el \
            ( dipole_radius = dipole_radius
            , refl_dist     = refl_dist
            , dir_dist      = dir_dist
            , director      = director
            , reflector     = director + refldiff
            , lambda_4      = lambda_4
            , frqidxmax     = 3
            , wire_radius   = self.wire_radius
            )
        return fd
    # end def compute_antenna

# end class Folded_Dipole_Optimizer

if __name__ == '__main__' :
    cmd = default_args ()
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
        ( '-L', '--director-length'
        , type = float
        , help = "(Half) Length of the director"
        , default = 0.2
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
    if args.action == 'optimize' :
        do = Folded_Dipole_Optimizer \
            ( verbose     = args.verbose
            , random_seed = args.random_seed
            , wire_radius = args.wire_radius
            , use_rtr     = args.use_rtr
            )
        do.run ()
    else :
        frqidxmax = args.frqidxmax
        if args.action in ('frgain', 'necout') :
            frqidxmax = 3
        fd = Folded_Dipole_3el \
            ( dipole_radius = args.dipole_radius
            , refl_dist     = args.reflector_distance
            , dir_dist      = args.director_distance
            , reflector     = args.reflector_length
            , director      = args.director_length
            , lambda_4      = args.lambda_len
            , wire_radius   = args.wire_radius
            , frqidxmax     = frqidxmax
            , frqidxnec     = args.frqidxmax
            , avg_gain      = args.average_gain
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
