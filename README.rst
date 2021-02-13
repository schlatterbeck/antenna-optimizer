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

The project currently lacks an installer and I'm experimenting in the
working directory. It is still very much work in progress.

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

For optimizing the two-element Yagi-Uda you can use ``folded.py``.
Later a 3-element antenna with the same principles was added in
``folded_3ele.py``.

The ``antenna_model.py`` factors out the common parts of the various
antennas.

The ``hb9cv.py`` models the well-know HB9CV antenna but only very
crudely. NEC2 isn't really suited for modelling the phasing stubs of that
antenna because it doesn't like parallel wires that are too close. I
mainly did this for comparing some of the antennas resulting from
optimization with the well-known performance of the HB9CV.

In the file ``logper.py`` you can find a 9-element log-periodic antenna.
It can currently not be optimized, the performance of the real antenna
is better than the results obtained with NEC, so I didn't implement an
optimizer for it yet.

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
may take some time.

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

This tells us that the evaluation (which is meaningful only to the
genetic algorithm), the genetic algorithm *maximizes* this value.
The command-line options after the line ``The Best String:`` can be used
to create a ``.nec`` file for that antenna. The antenna in the example
has a voltage standing wave ratio of < |_| 1.8 at the band ends and
around 1.15 in the middle of the band (the 70cm band from 430 to
440 |_| MHz in that case). The forward gain (in the middle of the band)
is 6.7 |_| dBi. The RMAX value is the (maximum) backward gain (in a 30
degree area in the back). So the F/B ratio of that antenna is:

.. math::
 6.7 \mbox{dB} - -3.0 \mbox{dB} = 9.7 \mbox{dB}

The last line of the text output contains the genetic representation of
that antenna.
The ``.nec`` file for the antenna above can be created with the command::

 ./folded.py -r 0.0364 -d 0.0444 -l 0.1704 -4 0.1075 necout > folded-opt.nec

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

.. _pgapy:  https://pypi.org/project/PGAPy/
.. _PGApack: https://github.com/schlatterbeck/pgapack
.. _PyNEC: https://pypi.org/project/PyNEC/
.. _`PyNEC source`: https://github.com/tmolteno/python-necpp
.. _`Numerical Electromagnetics Code`:
    https://en.wikipedia.org/wiki/Numerical_Electromagnetics_Code

