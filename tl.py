#!/usr/bin/python3
from __future__ import print_function

from antenna_model import Antenna_Model, Antenna_Optimizer, Excitation
from antenna_model import Arg_Handler

class Transmission_Line_Match (Antenna_Model) :
    wire_radius = 2e-3
    wire_len    = 0.005
    # Use wire_len for transmission lines: Since we're specifying the
    # transmission line length as a parameter to the TL card anyway we
    # can simplify the structure by not reserving so much space.
    use_wlen    = True
    feed_len    =  9.86
    z0          = 50.0
    frqstart    = 3.499
    frqend      = 3.501
    phi_inc     = 180
    theta_inc   = 90
    theta_max   = int (180 / theta_inc + 1)
    phi_max     = int (360 / phi_inc   + 1)

    def __init__ (self, stub_dist, stub_len, is_open = False, **kw) :
        self.stub_dist = stub_dist
        self.stub_len  = stub_len
        self.is_open   = is_open
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self) :
        r = []
        if self.is_open :
            r.append ('--open')
        r.append ('-d %(stub_dist)1.4f')
        r.append ('-l %(stub_len)1.4f')
        return ' '.join (r) % self.__dict__
    # end def cmdline

    def geometry (self, nec = None) :
        if nec is None :
            nec = self.nec
        geo = nec.get_geometry ()
        self.tag = 1
        self.ex  = None
        geo.wire \
            ( self.tag
            , 1
            , 0,             0, 0
            , self.wire_len, 0, 0
            , self.wire_radius
            , 1, 1
            )
        self.load_wire_tag = self.tag
        self.tag += 1
        sd = self.stub_dist
        if self.use_wlen :
            sd = self.wire_len
        geo.wire \
            ( self.tag
            , 1
            , 0,             sd, 0
            , self.wire_len, sd, 0
            , self.wire_radius
            , 1, 1
            )
        self.stub_point_tag = self.tag
        self.tag += 1
        geo.wire \
            ( self.tag
            , 1
            , 0, sd,                 0
            , 0, sd + self.wire_len, 0
            , self.wire_radius
            , 1, 1
            )
        self.tag += 1
        geo.wire \
            ( self.tag
            , 1
            , self.wire_len, sd,                 0
            , self.wire_len, sd + self.wire_len, 0
            , self.wire_radius
            , 1, 1
            )
        self.stub_start_tag = self.tag
        self.tag += 1
        sl = self.stub_len
        if self.use_wlen :
            sl = self.wire_len
        geo.wire \
            ( self.tag
            , 1
            , self.wire_len + sl, sd,                 0
            , self.wire_len + sl, sd + self.wire_len, 0
            , self.wire_radius
            , 1, 1
            )
        self.stub_end_tag = self.tag
        self.tag += 1
        geo.wire \
            ( self.tag
            , 1
            , 0,             sd + self.wire_len, 0
            , self.wire_len, sd + self.wire_len, 0
            , self.wire_radius
            , 1, 1
            )
        self.feed_end_tag = self.tag
        self.tag += 1
        fl = self.feed_len
        if self.use_wlen :
            fl = self.wire_len
        geo.wire \
            ( self.tag
            , 1
            , 0,             sd + self.wire_len + fl, 0
            , self.wire_len, sd + self.wire_len + fl, 0
            , self.wire_radius
            , 1, 1
            )
        self.ex = Excitation (self.tag, 1)
        nec.geometry_complete (0)
        nec.tl_card \
            ( self.stub_point_tag, 1
            , self.load_wire_tag,  1
            , self.z0
            , self.stub_dist
            , 0, 0, 1.0/150.0, 0 
            )
        stub_termination = 0.0
        if self.is_open :
            stub_termination = 1e30
        nec.tl_card \
            ( self.stub_start_tag, 1
            , self.stub_end_tag,   1
            , self.z0
            , self.stub_len
            , 0, 0, stub_termination, 0 
            )
        nec.tl_card \
            ( self.ex.tag,       self.ex.segment
            , self.feed_end_tag, 1
            , self.z0
            , self.feed_len
            , 0, 0, 0, 0 
            )
    # end def geometry
# end class Transmission_Line_Match

class Transmission_Line_Optimizer (Antenna_Optimizer) :

    def __init__ (self, **kw) :
        self.minmax = [(0, 40.0), (0, 40.0)]
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop) :
        stub_dist = self.get_parameter (p, pop, 0)
        stub_len  = self.get_parameter (p, pop, 1)
        tl = Transmission_Line_Match \
            ( stub_dist   = stub_dist
            , stub_len    = stub_len
            , frqidxmax   = 3
            )
        return tl
    # end def compute_antenna

    def evaluate (self, p, pop) :
        self.neval += 1
        ant, vswrs, gmax, rmax, swr_eval, swr_med = self.phenotype (p, pop)
        return 1.0 / swr_med
    # end def evaluate

# end class Transmission_Line_Optimizer

if __name__ == '__main__' :
    cmd = Arg_Handler (wire_radius = 0.002, frqidxmax = 3)
    cmd.add_argument \
        ( '-d', '--stub-distance'
        , type    = float
        , help    = "Distance of matching stub from load"
        , default = 7.137
        )
    cmd.add_argument \
        ( '-l', '--stub-length'
        , type    = float
        , help    = "Length of matching stub"
        , default = 30.39
        )
    args = cmd.parse_args ()
    if args.action == 'optimize' :
        tlo = Transmission_Line_Optimizer (** cmd.default_optimization_args)
        tlo.run ()
    else :
        tl = Transmission_Line_Match \
            ( stub_dist = args.stub_distance
            , stub_len  = args.stub_length
            , ** cmd.default_antenna_args
            )
        if args.action == 'necout' :
            print (tl.as_nec ())
        elif args.action not in actions :
            cmd.print_usage ()
        else :
            tl.compute ()
        if args.action == 'swr' :
            tl.swr_plot ()
        elif args.action == 'gain' :
            tl.plot ()
        elif args.action == 'frgain' :
            print ('\n'.join (tl.show_gains ()))
