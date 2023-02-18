#!/usr/bin/python3
from __future__ import print_function

import numpy as np
from .antenna_model import Antenna_Model, Antenna_Optimizer, Arg_Handler
from .antenna_model import Excitation, antenna_actions

class Multi_Dipole (Antenna_Model):
    """ Multiple coupled dipoles for 10m-15m, inspired by:
        Martin Erger. Endgespeiste Vertikalantenne für 15 m, 12 m und 10 m.
        Funkamateur, (5):384–387, May 2022.
    """

    name          = 'Multi Dipole 10m, 12m, 15m'
    wire_radius   = 1.5e-3 / 2.0
    radius_feed   = 3e-3
    lu_15m        = 3.26
    ll_15m        = 3.4
    lu_12m        = 3.05
    ll_12m        = 2.26
    lu_10m        = 2.45
    ll_10m        = 2.35
    d_12_15       = 50e-3
    d_10_15       = 33e-3
    seglen        = 0.15
    feedpoint_h   = 5
    # Use the data parts of each band
    frq_ranges    = [(28.070, 28.190), (24.915, 24.940), (21.070, 21.149)]
    # Theta must be < 90° when we have ground
    theta_range   = 90

    def __init__ \
        ( self
        , lu_15m        = lu_15m
        , ll_15m        = ll_15m
        , lu_12m        = lu_12m
        , ll_12m        = ll_12m
        , lu_10m        = lu_10m
        , ll_10m        = ll_10m
        , d_12_15       = d_12_15
        , d_10_15       = d_10_15
        , wire_radius   = wire_radius
        , radius_feed   = radius_feed
        , feedpoint_h   = feedpoint_h
        , ** kw
        ):
        self.lu_15m      = lu_15m
        self.ll_15m      = ll_15m
        self.lu_12m      = lu_12m
        self.ll_12m      = ll_12m
        self.lu_10m      = lu_10m
        self.ll_10m      = ll_10m
        self.d_12_15     = d_12_15
        self.d_10_15     = d_10_15
        self.radius_feed = radius_feed
        self.feedpoint_h = feedpoint_h
        self.__super.__init__ (** kw)
    # end def __init__

    def cmdline (self):
        return \
            ("--lu15 %(lu_15m)1.4f --ll15 %(ll_15m)1.4f "
             "--lu12 %(lu_12m)1.4f --ll12 %(ll_12m)1.4f "
             "--lu10 %(lu_10m)1.4f --ll10 %(ll_10m)1.4f "
             "--d12 %(d_12_15)1.4f --d10 %(d_10_15)1.4f"
            % self.__dict__
            )
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        # 15m is in the middle
        geo.wire \
            ( self.tag
            , round (self.lu_15m / self.seglen)
            , 0, 0, self.seglen * 1.5
            , 0, 0, self.lu_15m
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , round (self.ll_15m / self.seglen)
            , 0, 0, -self.seglen * 1.5
            , 0, 0, -self.ll_15m
            , self.radius_feed
            , 1, 1
            )
        self.tag += 1
        # A very small piece for the feedpoint
        geo.wire \
            ( self.tag
            , 3
            , 0, 0,  self.seglen * 1.5
            , 0, 0, -self.seglen * 1.5
            , self.wire_radius
            , 1, 1
            )
        self.ex = Excitation (self.tag, 1)
        self.tag += 1
        # 10m is to the left with distance -d_10_15
        geo.wire \
            ( self.tag
            , round (self.lu_10m / self.seglen)
            , -self.d_10_15, 0, 0
            , -self.d_10_15, 0, self.lu_10m
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , round (self.ll_10m / self.seglen)
            , -self.d_10_15, 0, 0
            , -self.d_10_15, 0, -self.ll_10m
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        # 12m is to the right with distance d_12_15
        geo.wire \
            ( self.tag
            , round (self.lu_12m / self.seglen)
            , self.d_12_15, 0, 0
            , self.d_12_15, 0, self.lu_12m
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , round (self.ll_12m / self.seglen)
            , self.d_12_15, 0, 0
            , self.d_12_15, 0, -self.ll_12m
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.move (0, 0, 0, 0, 0, self.feedpoint_h, 0, 0, 0)
    # end def geometry

    def ground (self, nec = None):
        """ We currently asume a medium ground with sommerfield-norton
            ground model
        """
        if nec is None:
            nec = self.nec
        # ground-type 2: Sommerfield-Norton
        # No radials for ground screen
        # Dieelectric-constant:
        # Conductivity:
        # Infinite groundplane: 0
        # Infinite groundplane: 0
        # Infinite groundplane: 0
        # Infinite groundplane: 0
        nec.gn_card (2, 0, 13, 0.005, 0, 0, 0, 0)
    # end def ground

# end class Multi_Dipole

class Multi_Dipole_Optimizer (Antenna_Optimizer):
    """ Optimize given multi-dipole for several bands
        Length are encoded as integers with a resolution of .5mm
        We use:
        2m <= l <= 4m
        for all lengths of half-elements (not quite half because the
        elements may have different length above and below the
        feedpoint) and
        3cm <= d <= 40cm
        for distances of 10/12m dipole from 15m
    """

    ant_cls = Multi_Dipole

    def __init__ \
        (self, radius_feed = None, feedpoint_h = 5, **kw):
        self.minmax = \
            [ (2,    4)
            , (2,    4)
            , (2,    4)
            , (2,    4)
            , (2,    4)
            , (2,    4)
            , (0.03, 0.4)
            , (0.03, 0.4)
            ]
        self.feedpoint_h = feedpoint_h
        # Force multiobjective
        if not kw.get ('multiobjective'):
            kw ['multiobjective'] = True
        self.radius_feed = radius_feed or self.ant_cls.radius_feed
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop):
        lu_15m  = self.get_parameter (p, pop, 0)
        ll_15m  = self.get_parameter (p, pop, 1)
        lu_12m  = self.get_parameter (p, pop, 2)
        ll_12m  = self.get_parameter (p, pop, 3)
        lu_10m  = self.get_parameter (p, pop, 4)
        ll_10m  = self.get_parameter (p, pop, 5)
        d_12_15 = self.get_parameter (p, pop, 6)
        d_10_15 = self.get_parameter (p, pop, 7)
        md = self.ant_cls \
            ( lu_15m  = lu_15m
            , ll_15m  = ll_15m
            , lu_12m  = lu_12m
            , ll_12m  = ll_12m
            , lu_10m  = lu_10m
            , ll_10m  = ll_10m
            , d_12_15 = d_12_15
            , d_10_15 = d_10_15
            , radius_feed   = self.radius_feed
            , avg_gain      = self.avg_gain
            , feedpoint_h   = feedpoint_h
            , **self.antenna_args
            )
        return md
    # end def compute_antenna

    def get_eval_and_constraints (self):
        """ We have three constraints on VSWR, one for each frequency.
            Then we have three constraints on the gain in opposite PHI
            direction (but equal THETA): These should not be off by more
            than 0.5dB. The real eval function is the gain on each
            frequency. Note that we ignore self.min_fb.
        """
        add = bool (self.min_gain) * 3
        return dict (num_eval = 9 + add, num_constraint = 6 + add)
    # end def get_eval_and_constraints

    def evaluate (self, p, pop):
        phenos = self.phenotype (p, pop)
        assert self.multiobjective
        retval = []
        for pheno in phenos:
            swr_max = max (pheno.vswrs)
            # We take absolute value of F/B ratio: We want about same
            # gain in all directions, note that the term
            # pheno.gmax - pheno.rmax may <0 if the gain is not the same
            # on all frequency points in one range: The gmax takes the
            # *minimum* while rmax takes the *maximum* over the
            # computed frequencies in one range.
            v = [ pheno.gmax
                , abs (pheno.gmax - pheno.rmax) - 1 # Less than 1 dB
                , swr_max - self.maxswr
                ]
            if self.min_gain:
                v.append (self.min_gain - pheno.gmax)
            retval.append (v)
        return tuple (np.array (retval).T.flatten ())
    # end def evaluate

    def get_gain_fw_maxswr (self, sub_eval):
        """ For printing the middle value is the F/B ration minus 1, see
            evaluate above, and for the last value we need to re-add
            maxswr.
        """
        sub_eval [-1] += self.maxswr
        sub_eval [-2] += 1
        return tuple (sub_eval)
    # end def get_gain_fw_maxswr

# end class Multi_Dipole_Optimizer

def main ():
    cmd = Arg_Handler (multiobjective = True)
    cmd.add_argument \
        ( '--lu15'
        , type = float
        , help = "Upper part of 15m (above feedpoint)"
        , default = Multi_Dipole.lu_15m
        )
    cmd.add_argument \
        ( '--ll15'
        , type = float
        , help = "Lower part of 15m (below feedpoint)"
        , default = Multi_Dipole.ll_15m
        )
    cmd.add_argument \
        ( '--lu12'
        , type = float
        , help = "Upper part of 12m (above feedpoint)"
        , default = Multi_Dipole.lu_12m
        )
    cmd.add_argument \
        ( '--ll12'
        , type = float
        , help = "Lower part of 12m (below feedpoint)"
        , default = Multi_Dipole.ll_12m
        )
    cmd.add_argument \
        ( '--lu10'
        , type = float
        , help = "Upper part of 10m (above feedpoint)"
        , default = Multi_Dipole.lu_10m
        )
    cmd.add_argument \
        ( '--ll10'
        , type = float
        , help = "Lower part of 10m (below feedpoint)"
        , default = Multi_Dipole.ll_10m
        )
    cmd.add_argument \
        ( '--d12'
        , type    = float
        , help    = "Distance of 12m dipole from 15m dipole"
        , default = Multi_Dipole.d_12_15
        )
    cmd.add_argument \
        ( '--d10'
        , type    = float
        , help    = "Distance of 10m dipole from 15m dipole"
        , default = Multi_Dipole.d_10_15
        )
    args = cmd.parse_args ()
    if args.action == 'optimize':
        mo = Multi_Dipole_Optimizer (** cmd.default_optimization_args)
        mo.run ()
    else:
        md = Multi_Dipole \
            ( lu_15m  = args.lu15
            , ll_15m  = args.ll15
            , lu_12m  = args.lu12
            , ll_12m  = args.ll12
            , lu_10m  = args.lu10
            , ll_10m  = args.ll10
            , d_12_15 = args.d12
            , d_10_15 = args.d10
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, md)
# end def main

if __name__ == '__main__':
    main ()
