#!/usr/bin/env python3
# Copyright (C) 2019-22 Dr. Ralf Schlatterbeck Open Source Consulting.
# Reichergasse 131, A-3411 Weidling.
# Web: http://www.runtux.com Email: office@runtux.com
# ****************************************************************************
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ****************************************************************************

import sys
from setuptools import setup
sys.path.insert (1, '.')
from antenna_optimizer import __version__

description = []
with open ('README.rst') as f:
    description = f.read ()

license     = 'BSD License'
setup \
    ( name             = "antenna-optimizer"
    , version          = __version__
    , description      = "Optimize antennas using NEC and a genetic algorithm"
    , long_description = ''.join (description)
    , license          = license
    , author           = "Ralf Schlatterbeck"
    , author_email     = "rsc@runtux.com"
    , install_requires = [ 'matplotlib', 'numpy', 'scipy'
                         , 'pgapy',  'pynec', 'rsclib'
                         ]
    , packages         = ['antenna_optimizer']
    , platforms        = 'Any'
    , python_requires  = '>=3.7'
    , entry_points     = dict
        ( console_scripts =
            [ 'coaxmodel=antenna_optimizer.coaxmodel:main'
            , 'folded-3ele-antenna=antenna_optimizer.folded_3ele:main'
            , 'folded-antenna=antenna_optimizer.folded:main'
            , 'folded-bc-antenna=antenna_optimizer.folded_bc:main'
            , 'folded-bigrefl-antenna=antenna_optimizer.folded_bigrefl:main'
            , 'hb9cv-antenna=antenna_optimizer.hb9cv:main'
            , 'hf-folded-dipole=antenna_optimizer.hf_folded:main'
            , 'hf-fuchs=antenna_optimizer.hf_fuchs:main'
            , 'hf-inverted-v=antenna_optimizer.hf_inverted_v:main'
            , 'logper-antenna=antenna_optimizer.logper:main'
            , 'multi-dipole=antenna_optimizer.multi_dipole:main'
            , 'transmission-line=antenna_optimizer.tl:main'
            ]
        )
    , url              = 'https://github.com/schlatterbeck/antenna-optimizer'
    , classifiers      =
        [ 'Development Status :: 5 - Production/Stable'
        , 'License :: OSI Approved :: ' + license
        , 'Operating System :: OS Independent'
        , 'Programming Language :: Python'
        , 'Intended Audience :: Developers'
        , 'Intended Audience :: Science/Research'
        , 'Programming Language :: Python :: 2'
        , 'Programming Language :: Python :: 2.7'
        , 'Programming Language :: Python :: 3'
        , 'Programming Language :: Python :: 3.5'
        , 'Programming Language :: Python :: 3.6'
        , 'Programming Language :: Python :: 3.7'
        , 'Programming Language :: Python :: 3.8'
        , 'Programming Language :: Python :: 3.9'
        ]
    )
