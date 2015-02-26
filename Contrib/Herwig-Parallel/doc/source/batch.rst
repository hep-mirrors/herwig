.. _batch:

Use Herwig-Parallel in Batch Mode
=================================

Apart from the commands ``herwig-parallel-read`` and ``herwig-parallel-run`` there also is the command ``herwig-parallel``.
This program takes the name of a run definition file as the argument and subsequently executes ``herwig-parallel-read`` and ``herwig-parallel-run``
in order to generate a parallel run according to the specifications provided in the run definition file.
An example run definition file called ``testrun.conf`` is provided in the ``examples`` directory of Herwig-Parallel.

When started as a background process this mechanism can be used to script and run in an automated manner a large number of jobs.
To start runs in the background do::
   
   herwig-parallel run1.conf > run1.log 2>&1 &
   herwig-parallel run2.conf > run2.log 2>&1 &
   herwig-parallel run3.conf > run3.log 2>&1 &

This will redirect all output to the respective log files.

Leaving empty the settings ``buildscript``, ``integratescript`` and/or ``runscript`` in a run definition file will instruct ``herwig-parallel`` to
use the default scripts from the ``<Herwig-Parallel>`` root directory (just as the option ``-d`` had done for the ``herwig-parallel-read`` and
``herwig-parallel-run`` commands.)
