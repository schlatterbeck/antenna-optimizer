#!/usr/bin/python3
from __future__ import print_function
import sys
import re
import math
from argparse import ArgumentParser
from rsclib.stateparser import Parser

# Parse statistics and create summary

class Optimization_Parser (Parser):
    encoding = None

    re_tbl = re.compile (r'^R\s+Gain')
    re_num = re.compile (r'^([0-9 .]+)')

    #      State   Pattern  new State Action
    matrix = \
        [ ["init", re_tbl,  "tbl",    "start_table"]
        , ["init", None,    "init",   None]
        , ["tbl",  re_num,  "tbl",    "table_line"]
        , ["tbl",  None,    "init",   "end_table"]
        ]

    def __init__ (self):
        self.head = None
        self.__super.__init__ ()
    # end def __init__

    def start_table (self, state, new_state, match):
        self.head = self.line.strip ().split ()
        assert 7 <= len (self.head) <= 8
        self.cols = []
        for k in self.head:
            self.cols.append ([])
    # end def start_table

    def table_line (self, state, new_state, match):
        l = match.group (1).strip ().split ()
        if len (l) < len (self.cols):
            return
        for n, k in enumerate (l):
            self.cols [n].append (float (k))
            if n == len (self.cols) - 1:
                break
    # end def table_line

    def end_table (self, state, new_state, match):
        assert self.head
        for n, c in enumerate (self.cols):
            mean = sum (c) / len (c)
            sdev = 0
            for v in c:
                sdev += (mean - v) ** 2
            sdev = math.sqrt (sdev / len (c))
            print ("%11s: %9.2f %9.2f" % (self.head [n], mean, sdev))
        self.head = None
    # end def end_table

# end class Optimization_Parser

class Optimization_Result:

    def __init__ (self, eval):
        self.best_eval = eval
    # end def __init__

    def __str__ (self):
        s = []
        s.append ("%02d"       % self.random_seed)
        s.append ("%5.2f"      % self.gmax)
        s.append ("%5.2f"      % self.fb)
        s.append ("%5.2f"      % self.vswr [0])
        s.append ("%5.2f"      % self.vswr [2])
        s.append ("%7.2f"      % self.best_eval)
        s.append ("%4d"        % self.generations)
        s.append ("       %7d" % self.n_eval)
        return ' '.join (s)
    # end def __str__
    __repr__ = __str__

# end class Optimization_Result

class Result_Parser (Parser):
    encoding = None

    re_best    = re.compile (r'^The Best Evaluation:\s+([0-9.eE+]+)[.]$')
    re_cmdline = re.compile (r'^((-[^ ]\s+[0-9.]+\s*)+)')
    re_swr     = re.compile (r'^VSWR:\s+\[((([0-9.]+)[, +]*){3})\]')
    re_max     = re.compile (r'^GMAX:\s+([-0-9.]+),\s+RMAX:\s+([-0-9.]+)')
    re_cache   = re.compile (r'Cache hits:\s+([0-9]+)/([0-9]+)\s+([0-9.]+)%')
    re_iter    = re.compile \
        (r'Iter:\s+([0-9]+)\s+Evals:\s+([0-9]+)\s+Stag:\s+([0-9]+)')
    re_gene    = re.compile (r'^#\s+(.*)')
    re_bingene = re.compile (r'^\[\s*([01]+)\s*\]$')
    re_iiter   = re.compile (r'^([0-9]+)\s+Best\s+[0-9.+eE]+$')

    #      State   Pattern  new State Action
    matrix = \
        [ ["init",   re_best,       "result",    "eval"]
        , ["init",   re_iiter,      "init",      "intermediate_iter"]
        , ["init",   None,          "init",      None]
        , ["result", re_cmdline,    "result",    "cmd"]
        , ["result", re_swr,        "result",    "swr"]
        , ["result", re_max,        "result",    "maxgain"]
        , ["result", re_cache,      "result",    "cache"]
        , ["result", re_iter,       "result",    "iterations"]
        , ["result", re_gene,       "init",      "gene"]
        , ["result", re_bingene,    "init",      "gene"]
        , ["result", None,          "result",    None]
        ]

    def __init__ (self):
        self.opt_results = []
        self.opt         = None
        self.randseed    = None
        self.iiter       = 0
        self.__super.__init__ ()
    # end def __init__

    def __iter__ (self):
        for opt in self.opt_results:
            yield (opt)
    # end def __iter__

    def set_random_seed (self, random_seed):
        self.randseed = random_seed
    # end def set_random_seed

    def result_iter (self):
        yield ("R  Gain  f/b    SWR   SWR  Eval    Generations Evaluations")
        for opt in self:
            yield (str (opt))
    # end def result_iter

    def intermediate_iter (self, state, new_state, match):
        self.iiter = int (match.group (1))
    # end def intermediate_iter

    def eval (self, state, new_state, match):
        evaluation = float (match.group (1))
        self.opt   = Optimization_Result (evaluation)
        self.opt.generations = self.iiter
        self.iiter = 0
        if self.randseed:
            self.opt.random_seed = self.randseed
            self.randseed   = None
        self.opt_results.append (self.opt)
    # end def eval

    def cmd (self, state, new_state, match):
        self.opt.cmdline = match.group (1)
    # end def cmd

    def swr (self, state, new_state, match):
        self.opt.vswr = [float (x) for x in match.group (1).split (', ')]
    # end def swr

    def maxgain (self, state, new_state, match):
        self.opt.gmax = float (match.group (1))
        self.opt.rmax = float (match.group (2))
        self.opt.fb   = self.opt.gmax - self.opt.rmax
    # end def maxgain

    def cache (self, state, new_state, match):
        self.opt.cache_hits    = float (match.group (1))
        self.opt.n_eval        = float (match.group (2))
        self.opt.cache_percent = float (match.group (3))
    # end def cache

    def iterations (self, state, new_state, match):
        self.opt.generations   = float (match.group (1))
        # Will overwrite previous result from cache, the new value also
        # has the first popsize evals in the 0th generation
        self.opt.n_eval        = float (match.group (2))
        self.opt.stag_count    = float (match.group (3))
    # end def iterations

    def gene (self, state, new_state, match):
        self.opt.gene = match.group (1)
        self.opt      = None
    # end def gene

# end class Result_Parser

if __name__ == '__main__':
    cmd = ArgumentParser ()
    cmd.add_argument \
        ( 'files'
        , help  = "Files to parse"
        , nargs = "+"
        )
    cmd.add_argument \
        ( '-o', '--optimization-lines'
        , action  = "store_true"
        , help    = "Parse optimization result lines, not optimizer run output"
        )
    args = cmd.parse_args ()
    if args.optimization_lines:
        if len (args.files) > 1:
            print ("Only one file for -o option", file = sys.stderr)
            sys.exit (1)
        op = Optimization_Parser ()
        with open (args.files [0], 'r') as f:
            op.parse (f)
        if op.head:
            op.end_table ('', '', '')
    else:
        digi = re.compile (r'^.*[^0-9]([0-9]+)[^/0-9]*$')
        rp   = Result_Parser ()
        for fn in args.files:
            d = digi.search (fn)
            rp.set_random_seed (int (d.group (1)))
            #print (fn, d.group (1))
            with open (fn, 'r') as f:
                rp.parse (f)
        for line in rp.result_iter ():
            print (line)
        print ("")
        op = Optimization_Parser ()
        op.parse (rp.result_iter ())
        if op.head:
            op.end_table ('', '', '')

