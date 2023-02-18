#!/usr/bin/python3
from __future__ import print_function

import sys

from .antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from .antenna_model import Arg_Handler, antenna_actions
from .transmission import transmission_line_z_square

class Logperiodic (Antenna_Model):
    """ Log Periodic Antenna Model, Default data taken from
        http://adl303-archiv.oevsv.at/technikecke/dlogpant/index.html
        Given lengths are with boom radius according to website above in
        the picture it states: "Die Elementlaengen sind von Mitte des
        Booms bis zum Element-Ende errechnet !!" Contrary to this in the
        table with the lengths and distances it is stated:
        "Achtung: halbe Boomdicke kommt noch zur Elementlaenge hinzu! (7mm)"
        Other sources indicate they already include the boom radius, e.g.
        http://www.dl6ka.de/log.htm tells us
        "Achtung: Die Elementlaengen sind von Mitte des Booms bis zum
        Element-Ende berechnet !" while giving the same numbers as above.
        Simulation without added 7mm is better but we currently model
        the boom as a round wire which is incorrect.

        The design claims 8.5 dBd at 144 MHz and 10.9 dBd at 433 MHz.
        The current simulation results are (with 7.5mm added element length):
        144 MHz: 11.28 dBi =  9.13 dBd  SWR: 1.84 .. 2.28
        440 MHz: 12.30 dBi = 10.15 dBd  SWR: 3.36 .. 3.35
        Without 7.5mm added we get:
        144 MHz: 11.26 dBi =  9.11 dBd  SWR: 1.71 .. 1.84 (min 1.67)
        440 MHz: 11.50 dBi =  9.35 dBd  SWR: 3.40 .. 3.35

        The original PDF from DL9HCG Website and other versions floating
        around give a completely different design. All are using 4mm
        Aluminium "Schweisstab" and a 15mmx15mm square aluminium tube
        for the boom. The elements are in the inner half of the boom.

        One measures from element end to middle of boom:
        53.04, 50.67, 48.30, 45.93, 43.56, 41.19, 38.82, 36.45, 34.08
        with distances *between* elements ("Lichter Abstand") of:
        7.43, 7.41, 7.38, 7.36, 7.33, 7.31, 7.28, 7.26
        So we would have to add element width of 4mm to each of these
        distances to model with NEC.
        (this one indicates 10.5 dBd for 2m and 12.3 dBd for 70cm)

        The other has the element numbering reversed and lists the
        following but measuring from element end to *surface* of boom:
        length:
        52.00, 49.75, 47.50, 45.25, 43.00, 40.75, 38.50, 36.25, 34.00
        with distances middle-to-middle of:
        9.10, 8.59, 8.12, 7.67, 7.24, 6.84, 6.46, 6.10
        (this indicates 10.1 dBd for 2m and 11.6 dBd for 70cm)
    """

    frqstart      = 144
    frqend        = 146
    wire_radius   = 4e-3 / 2.0
    lengths       = [ 0.53,  0.509, 0.489, 0.469, 0.45
                    , 0.432, 0.415, 0.398, 0.382
                    ]
    dists         = [0.123, 0.121, 0.119, 0.117, 0.115, 0.114, 0.112, 0.11]
    boom_r        = 0.016  / 2.0
    segs_dipole   = 19
    segs_boom     =  5
    d_boom        = 0.005

    # This would add half of the boom length, element lengths claim to
    # be with boom length included, but see comments above
    # for n in range (len (lengths)):
    #     lengths [n] = lengths [n] + boom_r

    def __init__ \
        ( self
        , lengths     = lengths
        , dists       = dists
        , boom_r      = boom_r
        , d_boom      = d_boom
        , wire_radius = wire_radius
        , tline       = False
        , ** kw
        ):
        self.lengths     = lengths
        self.dists       = dists
        self.boom_r      = boom_r
        self.d_boom      = d_boom + 2 * boom_r
        self.wire_radius = wire_radius
        self.tline       = tline
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self):
        c = 1
        r = []
        for l in self.lengths:
            s = ' '
            if c % 4 == 0:
                s = '\n'
            r.append ('-l %1.4f%s' % (l, s))
            c += 1
        for d in self.dists:
            s = ' '
            if c % 4 == 0:
                s = '\n'
            r.append ('-d %1.4f%s' % (l, s))
            c += 1
        return "".join (r)
    # end def cmdline

    def geometry (self, nec = None):
        if nec is None:
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        dists = [0.0] + self.dists
        bpos  = 0.0
        tltag = []
        for n, (l, d) in enumerate (zip (self.lengths, dists)):
            dm = 1 - 2 * (n & 1)
            for nb, boomdist in enumerate ((0, self.d_boom)):
                boom_inv = 1 - 2 * (nb & 1)
                dir = dm * boom_inv
                geo.wire \
                    ( self.tag
                    , self.segs_dipole
                    , boomdist, d + bpos, 0
                    , boomdist, d + bpos, l * dir
                    , self.wire_radius
                    , 1, 1
                    )
                self.tag += 1
                if d and not self.tline:
                    geo.wire \
                        ( self.tag
                        , self.segs_boom
                        , boomdist, bpos,     0
                        , boomdist, bpos + d, 0
                        , self.boom_r
                        , 1, 1
                        )
                    self.tag += 1
            if self.tline:
                geo.wire \
                    ( self.tag
                    , 1
                    ,        0, d + bpos, 0
                    , boomdist, d + bpos, 0
                    , self.wire_radius
                    , 1, 1
                    )
                tltag.append (self.tag)
                self.tag += 1
            bpos += d

        if self.tline:
            self.ex = Excitation (tltag [-1], 1)
        else:
            geo.wire \
                ( self.tag
                , 1 # only one segment for feeding
                , 0,        bpos, 0
                , boomdist, bpos, 0
                , self.boom_r
                , 1, 1
                )
            self.ex = Excitation (self.tag, 1)
            self.tag += 1
        # Turn around Y by 270 deg, move everything up
        #geo.move (0, 270, 0, 0, 0, self.up, 0, 0, 0)
        nec.geometry_complete (0)
        if self.tline:
            impedance = transmission_line_z_square \
                (2 * self.boom_r, 2 * self.boom_r + self.d_boom)
            impedance = 50.0
            last = None
            for n, t in enumerate (tltag):
                if last:
                    d = self.dists [n - 1]
                    nec.tl_card (last, 1, t, 1, -impedance, d, 0, 0, 0, 0)
                last = t
    # end def geometry

    @property
    def up (self):
        """ move everything up by max (reflector length, lambda_4 + r)
        """
        return max (self.reflector, self.lambda_4 + self.dipole_radius)
    # end def up

# end class Logperiodic

def main ():
    cmd = Arg_Handler ()
    cmd.add_argument \
        ( '-b', '--boom-distance'
        , type    = float
        , help    = "Distance between booms"
        , default = 0.002
        )
    cmd.add_argument \
        ( '-d', '--distance'
        , type    = float
        , nargs   = '+'
        , help    = "Distance between elements, option should be repeated"
        )
    cmd.add_argument \
        ( '-l', '--length'
        , type    = float
        , nargs   = '+'
        , help    = "Length of one element, option should be repeated"
        )
    cmd.add_argument \
        ( '-t', '--use-transmission-line'
        , help    = "Use transmission line instead of solid boom"
        , action  = 'store_true'
        )
    args = cmd.parse_args ()
    if not args.distance:
        args.distance = Logperiodic.dists
    if not args.length:
        args.length = Logperiodic.lengths
    if args.action == 'optimize':
        do = Folded_Dipole_Optimizer (** cmd.default_optimization_args)
        do.run ()
    else:
        ant = Logperiodic \
            ( lengths     = args.length
            , dists       = args.distance
            , d_boom      = args.boom_distance
            , tline       = args.use_transmission_line
            , ** cmd.default_antenna_args
            )
        antenna_actions (cmd, args, ant)
# end def main

if __name__ == '__main__':
    main ()
