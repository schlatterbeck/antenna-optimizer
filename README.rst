Antenna Optimizer
=================

.. |ohm| unicode:: U+02126 .. Omega
.. |_| unicode:: U+00A0 .. Non-breaking space
   :trim:

This project can optimize anntennas using genetic algorithms.  It uses
my pgapy_ Python wrapper for PGApack_, the parallel genetic algorithm
library, originally by David Levine at Argonne National Laboratory and
currently maintained by me. It also uses PyNEC_, the Python wrapper for
NEC2++, the C++ version of the `Numerical Electromagnetics Code`_.

Originally this started out with a low-gain two-element antenna where
the driven element is a folded dipole. One of the requirements for that
antenna was that it should have 50 |_| |ohm| impedance and at least some
forward gain. Optimizing this antenna by hand soon turned out to be
tedious and I started experimenting with optimization by a genetic
algorithm.

You can find the original experiments in ``folded.nec`` and
``folded2.nec``. These are input files to the command-line NEC programs.
You can either use the ``nec2c`` command-line program which produces an
output file that can be viewed with ``xnecview`` or you can use
graphical ``xnec2c`` program. All run under Linux and are included in
the Debian Linux distribution. The ``.nec`` files *may* also be usable
with other NEC versions on other operating systems but I have not tried.

For optimizing the two-element Yagi-Uda you can use the command-line
tool ``folded_antenna`` where the optimizer is implemented in
``folded.py``. Later a 3-element antenna with the same principles was
added in ``folded_3ele.py`` callable with the command-line tool
``folded_3ele_antenna``.

The ``antenna_model.py`` factors out the common parts of the various
antennas.

The ``hb9cv.py`` with the command-line tool ``hb9cv_antenna`` models the
well-know HB9CV antenna but only very crudely. NEC2 isn't really suited
for modelling the phasing stubs of that antenna because it doesn't like
parallel wires that are too close. I mainly did this for comparing some
of the antennas resulting from optimization with the well-known
performance of the HB9CV.

In the file ``logper.py`` with the command-line tool ``logper_antenna``
you can find a 9-element log-periodic antenna.  It can currently not be
optimized, the performance of the real antenna is better than the
results obtained with NEC, so I didn't implement an optimizer for it
yet.

All the antenna programs take an action as mandatory argument. The
action is typically either ``optimize`` for running the optimizer or
``necout`` for creating a ``.nec`` file which can then be fed to one of
the nec programs mentioned above. When running the optimizer it makes
sense to experiment with different random seeds, each random seed will
usually produce a different antenna. In addition there are some
experimental actions, ``frgain`` prints the forward and backward gains
(in dBi) for the lowest, the middle, and the highest frequencies and the
VSWR for those. The ``gain`` action visualizes the 3D antenna gain
pattern and the ``swr`` action visualizes the VSWR over the given
frequency range. Note that both, the ``gain`` and the ``swr`` action
compute the antenna data over the whole frequency range using NEC and
that may take some time.

The output of the optimizer is text (usually redirected to a file) that
prints the evaluation, the VSWR, maximum gain, and forward/backward
ratio of the best antenna for every 10th generation of the genetic
algorithm. In addition the command-line options to create that antenna
are printed. When the genetic algorithm doesn't make any more progress,
the search terminates and the data of the best evaluation are given. An
example of the last lines of such a text is as follows. The data is from
one of the best 2-element antennas was obtained with the random seed 26
of an earlier version of the program::

    The Best Evaluation: 2.886437e+02.
    The Best String:
    -r 0.0364 -d 0.0444 -l 0.1704 -4 0.1075
    VSWR: [1.7901433511443068, 1.1495780609982815, 1.7995760521232753]
    GMAX: 6.69913175227, RMAX: -3.03663376703
    Cache hits: 5670/9243 61.34%
    Eval: 288.64
    [ 101011001100101001111000000010011 ]

This tells us the evaluation (which is meaningful only to the
genetic algorithm), the genetic algorithm *maximizes* this value.
The command-line options after the line ``The Best String:`` can be used
to create a ``.nec`` file for that antenna. The antenna in the example
has a voltage standing wave ratio of < |_| 1.8 at the band ends and
around 1.15 in the middle of the band (the 70cm band from 430 to
440 |_| MHz in that case). The forward gain (in the middle of the band)
is 6.7 |_| dBi. The RMAX value is the (maximum) backward gain (in a 30
degree area in the back). So the F/B ratio of that antenna is::

 6.7 dB - -3.0 dB = 9.7 dB

The last line of the text output contains the genetic representation of
that antenna.
The ``.nec`` file for the antenna above which was optimized with an early
version of this package can be created with the command::

 folded_antenna -r 0.0364 -d 0.0444 -l 0.1704 -4 0.1075 necout > folded-opt.nec

The command-line options specify the radius of the folded dipole, the
distance of the reflector from the folded dipole, the (half) length of
the reflector, and the (half) length of the straight part of the folded
dipole, respectively.

According to NEC it has a standing wave ratio (VSWR) of < |_| 1.8 from
430-440 |_| MHz, a forward gain of > |_| 6.5 dBi over the whole
frequency range and a Forward/Back Ratio of 8-11 dB.

If you want to implement an optimizer for your own antenna, look at the
file ``folded.py``: You need to implement a class that defines the
geometry of the new antenna and an optimizer class that initializes the
gene ranges and implements a ``compute_antenna`` method that returns an
instance of your antenna class with the parameters obtained from the
given gene. All lengths in the models are metric (in meters) as is the
default in NEC.

A recent addition to this package involves modelling of coax cables.
This uses information from an old article by Frank Witt [1]_ to derive
everything necessary to model a transmission line *with* loss from the
manufacturer cable data. The command-line tool for using these coax
models is named ``coaxmodel``. Again this command has several
sub-commands: 

- ``loss``: This displays the fitted loss-curves from the manufacturer
  data against the curve-fit algorithm used, you can see how much
  difference in dB the fitted curve has to the loss at certain
  frequencies given by the manufacturer data.
- ``match`` computes the impedance at the load and input (depending on
  which was given as an input the other is computed), the matched and
  total loss (the matched loss is the loss in the cable if the load is
  prefectly matched, the total loss is the sum of the matched loss and
  the additional loss due to reflections), the SWR and data for various
  stub-matches to get to the cable impedance Z0.
- ``resonator`` computes the resistance and Q-factor of a coax resonator
  at the given frequency. A resonator is a piece of cable that either
  has a short-circuit or an open-circuit at the far end. The sub-command
  computes resonators for quarter and half wave at the given frequency.
- ``stub`` computes the length, Q-factor, and resulting impedance as
  well as the inductance or capacitance of a stub for a given impedance
  and frequency. The impedance by default is -100j Ohm and can be
  changed with the -x option. The stub with the shortest resulting
  length is chosen, so for a negative reactance an open stub is chosen
  while for a positive reactance a closed stub is chosen.

For all these sub-commands you can specify the frequency, length of
cable and impedance (either at the load or at the end of the cable) to
be used in computing the results. You can specify complex impedances as
a python complex number in the format a+bj, e.g. 50-500j.

Finally the ``transmission_line`` program can optimize the stub-matching
for a transmission line using NEC. By default a lossless line is asumed.
Also by default a closed stub at the closest possible position is
searched.

Transmission lines can be modelled by NEC with the ``TL`` (transmission
line) card.  But NEC can also model arbitrary (symmetric, passive)
networks with the ``NT`` (network) card. We use this (and the code in
``coaxmodel.py``) to model a real cable *with loss* for stub matching.
It is instructive to compare the values for stub-matching obtained
analytically from ``coaxmodel`` with the values obtained from an
optimization with the genetic algorithm by ``transmission_line``. Note
that the coaxmodel takes frequecies in Hz while transmission_line (which
uses NEC) accepts frequencies in MHz. So we match a complex impedance of
75+15j |Ohm| for example::

 coaxmodel -c sytronic_RG_58_CU -f 435e6 -z 75+15j match

This yields a stub of length 8.007 cm attached 7.888 cm from the load when
matching with a closed stub. When optimizing with ``transmission_line``::

 transmission_line -c sytronic_RG_58_CU -f 435 -z 75+15j optimize

we get 8.016cm for the stub length and 7.6cm for the distance of the
stub from the load. We can visualize this over a given frequency range
by either producing NEC output::

 transmission_line -c sytronic_RG_58_CU -f 435 -z 75+15j -i 50 \
    -l 0.0816 -d 0.076 --frqstart=430 --frqend=440 necout > x1.nec
 transmission_line -c sytronic_RG_58_CU -f 435 -z 75+15j -i 50 \
    -l 0.08007 -d 0.07888 --frqstart=430 --frqend=440 necout > x2.nec

And the compute the nec model and display with::

 nec2c -i x1.nec > x1.out
 nec2c -i x2.nec > x2.out
 xnecview x1.out
 xnecview x2.out

Or directly display the VSWR curves with::

 transmission_line -c sytronic_RG_58_CU -f 435 -z 75+15j -i 50 \
    -l 0.0816 -d 0.076 --frqstart=430 --frqend=440 swr
 transmission_line -c sytronic_RG_58_CU -f 435 -z 75+15j -i 50 \
    -l 0.08007 -d 0.07888 --frqstart=430 --frqend=440 swr

Both are close enough, the SWR is below 1.1 over the whole frequency
range given. Note that this can change drastically if load impedances
with a higher VSWR are matched.

Also note that the NEC files produced in the example above have a different
NEC network *for each frequency*. This is because NEC models networks
using an `admittance matrix`_ which is frequency dependent.

This means the sequence of two ``NT``
cards, a ``TL`` card, a ``FR`` card and a ``RP`` card are repeated for
each frequency. Here the two ``NT`` cards define the network of the
cable from the load to the stub and the stub itself while the ``TL``
card defines the length of a lossless transmission line from the stub to
the source. The ``FR`` card specifies a single frequency and the ``RP``
card defines a radiation pattern and triggers computation. This format
is perfectly valid NEC code, but certain programs (like the popular
``xnec2c``) cannot deal with this format and display only a single
frequency.

.. [1] Frank Witt. Transmission line properties from manufacturer’s
   data. In R. Dean Straw, editor, The ARRL Antenna Compendium, volume 6,
   pages 179–183. American Radio Relay League (ARRL), 1999.

.. _pgapy:  https://pypi.org/project/PGAPy/
.. _PGApack: https://github.com/schlatterbeck/pgapack
.. _PyNEC: https://pypi.org/project/PyNEC/
.. _`PyNEC source`: https://github.com/tmolteno/python-necpp
.. _`Numerical Electromagnetics Code`:
    https://en.wikipedia.org/wiki/Numerical_Electromagnetics_Code
.. _`admittance matrix`: https://en.wikipedia.org/wiki/Admittance_parameters

Changes
-------

Version 0.3: Multi-objective optimization

- Switch to pyproject.toml instead of setup.py
- Add multi-objective optimization
- Allow NSGA-III for multi-objective optimization
- Allow to model ground
- Multiple frequency ranges
- Allow to use average gain when optimizing. Needs a bug-fix in pynec
  https://github.com/tmolteno/necpp/pull/73
- Add epsilon constrained optimization
  This allows to better find areas with good gain even if constraining
  the solutions to low SWR

Version 0.2: More cable data

- Fix setup to correctly specify dependencies
- Add more cable data the following command will list supported cable
  types::

    coaxmodel --help

Version 0.1: Initial Release
