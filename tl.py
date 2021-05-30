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

    def __init__ \
        ( self
        , stub_dist
        , stub_len
        , is_open   = False
        , is_series = False
        , frequency = None
        , **kw
        ) :
        self.stub_dist = stub_dist
        self.stub_len  = stub_len
        self.is_open   = is_open
        self.is_series = is_series
        if frequency is not None :
            self.frqstart  = frequency - 0.01
            self.frqend    = frequency + 0.01
            self.frequency = frequency
        self.__super.__init__ (**kw)
    # end def __init__

    def cmdline (self) :
        r = []
        if self.is_open :
            r.append ('--is-open')
        if self.is_series :
            r.append ('--is-series')
        r.append ('-d %(stub_dist)1.10f')
        r.append ('-l %(stub_len)1.10f')
        return ' '.join (r) % self.__dict__
    # end def cmdline

    def geometry (self, nec = None) :
        if nec is None :
            nec = self.nec
        geo = nec.get_geometry ()
        sd = self.stub_dist
        if self.use_wlen :
            sd = self.wire_len
        sl = self.stub_len
        if self.use_wlen :
            sl = self.wire_len
        fl = self.feed_len
        if self.use_wlen :
            fl = self.wire_len
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
        if self.is_series :
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
        else :
            self.stub_start_tag = self.tag - 1
            self.feed_end_tag   = self.tag - 1
            geo.wire \
                ( self.tag
                , 1
                , 0,             sd, -sl
                , self.wire_len, sd, -sl
                , self.wire_radius
                , 1, 1
                )
            self.stub_end_tag = self.tag
            self.tag += 1
        feedlen = sd + fl
        if self.is_series :
            feedlen += self.wire_len
        geo.wire \
            ( self.tag
            , 1
            , 0,             feedlen, 0
            , self.wire_len, feedlen, 0
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

    def __init__ \
        ( self
        , is_open      = False
        , is_series    = False
        , add_lambda_4 = False
        , frequency    = None
        , **kw
        ) :
        self.is_open      = is_open
        self.is_series    = is_series
        self.add_lambda_4 = add_lambda_4
        self.frequency    = frequency
        c  = 3e8
        tl = Transmission_Line_Match
        if frequency is None :
            self.frequency = (tl.frqstart + tl.frqend) / 2
        self.lambda_4 = c / 1e6 / self.frequency / 4
        self.minmax = [(0, self.lambda_4), (0, 2 * self.lambda_4)]
        if self.add_lambda_4 :
            self.minmax [0] = (self.lambda_4, 2 * self.lambda_4)
        self.__super.__init__ (**kw)
    # end def __init__

    def compute_antenna (self, p, pop) :
        stub_dist = self.get_parameter (p, pop, 0)
        stub_len  = self.get_parameter (p, pop, 1)
        tl = Transmission_Line_Match \
            ( stub_dist      = stub_dist
            , stub_len       = stub_len
            , is_open        = self.is_open
            , is_series      = self.is_series
            , frqidxmax      = 3
            , wire_radius    = self.wire_radius
            , copper_loading = self.copper_loading
            , frequency      = self.frequency
            )
        return tl
    # end def compute_antenna

    def evaluate (self, p, pop) :
        self.neval += 1
        ant, vswrs, gmax, rmax, swr_eval, swr_med = self.phenotype (p, pop)
        return 1.0 / swr_med
    # end def evaluate

#    def print_string (self, file, p, pop) :
#        antenna, vswrs, gmax, rmax, swr_eval, swr_med = self.phenotype (p, pop)
#        print ("vswrs: %s" % vswrs, file = file)
#        print (antenna.as_nec (), file = file)
#    # end def print_string

# end class Transmission_Line_Optimizer

if __name__ == '__main__' :
    cmd = Arg_Handler \
        ( wire_radius    = 0.002
        , frqidxmax      = 3
        , copper_loading = False
        )
    cmd.add_argument \
        ( '-d', '--stub-distance'
        , type    = float
        , help    = "Distance of matching stub from load"
        , default = 7.106592973007906
        )
    cmd.add_argument \
        ( '-f', '--frequency'
        , type    = float
        , help    = "Frequency to match transmission line (MHz)"
        , default = 3.5
        )
    cmd.add_argument \
        ( '-l', '--stub-length'
        , type    = float
        , help    = "Length of matching stub"
        , default = 32.193495674862291
        )
    cmd.add_argument \
        ( '--is-open'
        , help    = "Use open stub"
        , action  = "store_true"
        )
    cmd.add_argument \
        ( '--is-series'
        , help    = "Use stub in series with transmission line"
        , action  = "store_true"
        )
    cmd.add_argument \
        ( '--add-lambda-4'
        , help    = "Add lambda/4 to the matching point (optimizer)"
        , action  = "store_true"
        )
    args = cmd.parse_args ()
    if args.action == 'optimize' :
        tlo = Transmission_Line_Optimizer \
            ( is_open      = args.is_open
            , is_series    = args.is_series
            , add_lambda_4 = args.add_lambda_4
            , frequency    = args.frequency
            , ** cmd.default_optimization_args
            )
        tlo.run ()
    else :
        tl = Transmission_Line_Match \
            ( stub_dist = args.stub_distance
            , stub_len  = args.stub_length
            , is_open   = args.is_open
            , is_series = args.is_series
            , frequency = args.frequency
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
