from __future__ import print_function
import sys
import re
import math
from rsclib.stateparser import Parser

# Parse statistics and create summary

class Optimization_Parser (Parser) :
    re_tbl = re.compile (r'^R\s+Gain')
    re_num = re.compile (r'^([0-9 .]+)')

    #      State   Pattern  new State Action
    matrix = \
        [ ["init", re_tbl,  "tbl",    "start_table"]
        , ["init", None,    "init",   None]
        , ["tbl",  re_num,  "tbl",    "table_line"]
        , ["tbl",  None,    "init",   "end_table"]
        ]

    def __init__ (self) :
        self.head = None
        self.__super.__init__ ()
    # end def __init__

    def start_table (self, state, new_state, match) :
        self.head = self.line.strip ().split ()
        assert len (self.head) == 7
        self.cols = []
        for k in self.head :
            self.cols.append ([])
    # end def start_table

    def table_line (self, state, new_state, match) :
        l = match.group (1).strip ().split ()
        if len (l) < 7 :
            return
        for n, k in enumerate (l) :
            self.cols [n].append (float (k))
            if n == 6 :
                break
    # end def table_line

    def end_table (self, state, new_state, match) :
        assert self.head
        for n, c in enumerate (self.cols) :
            mean = sum (c) / len (c)
            sdev = 0
            for v in c :
                sdev += (mean - v) ** 2
            sdev = math.sqrt (sdev / len (c))
            print ("%11s: %2.2f %2.2f" % (self.head [n], mean, sdev))
        self.head = None
    # end def end_table

# end class Optimization_Parser

if __name__ == '__main__' :
    op = Optimization_Parser ()
    op.parse (open (sys.argv [1]))
    if op.head :
        op.end_table ('', '', '')
