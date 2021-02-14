#!/usr/bin/python3
from __future__ import print_function

from antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from antenna_model import default_args

class Folded_Dipole (Antenna_Model) :
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
    frqstart      = 87.5
    frqend        = 108.0
    impedance     = 75.0

    def __init__ \
        ( self
        , refl_dist     = .01
        , dipole_radius = dipole_radius
        , lambda_4      = lambda_4
        , wire_radius   = wire_radius
        , ** kw
        ) :
        self.dipole_radius = dipole_radius
        self.lambda_4      = lambda_4
        self.wire_radius   = wire_radius
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self) :
        return ("-r %(dipole_radius)1.4f -4 %(lambda_4)1.4f" % self.__dict__)
    # end def cmdline

    def geometry (self, nec = None) :
        if nec is None :
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        self._geometry (geo)
        # Move everything up
        geo.move (0, 0, 90, 0, 0, self.up, 0, 0, 0)
        nec.geometry_complete (0)
    # end def geometry

    def _geometry (self, geo) :
        a = ((-90, 90), (90, 270))
        for n, z in enumerate ((1, -1)) :
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
        for x in (-self.lambda_4, self.lambda_4) :
            for z in (self.dipole_radius, -self.dipole_radius) :
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
    # end def _geometry

    @property
    def up (self) :
        """ move everything up by max (reflector length, lambda_4 + r)
        """
        return self.lambda_4 + self.dipole_radius
    # end def up

# end class Folded_Dipole

class Folded_Dipole_Optimizer (Antenna_Optimizer) :
    """ Optimize given folded dipole
        Length are encoded as integers with a resolution of .5mm
        We use:
        * 8mm  <= dipole_radius <= 40cm
        * 40cm <= lambda_4      <= 150cm
    """

    def __init__ (self, **kw) :
        self.minmax = [(8e-3, 0.4), (0.4, 1.5)]
        self.__super.__init__ (nofb = True, **kw)
    # end def __init__

    def compute_antenna (self, p, pop) :
        dipole_radius = self.get_parameter (p, pop, 0)
        lambda_4      = self.get_parameter (p, pop, 1)
        fd = Folded_Dipole \
            ( dipole_radius = dipole_radius
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
        fd = Folded_Dipole \
            ( dipole_radius = args.dipole_radius
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
