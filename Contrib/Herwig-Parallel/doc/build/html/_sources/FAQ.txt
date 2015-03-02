Frequently Asked Questions
==========================

.. contents::

Why does herwig-parallel-read start the build job alright but wait forever for its completion and never start the integrate job(s)?
-----------------------------------------------------------------------------------------------------------------------------------

The condition for the successful completion of the build step and the start of the integrate step is the existence of the runfile ``<generator name>.run``. The name of the generator is set in the Herwig++ infile by a directive such as ``saverun LHC-Matchbox LHCGenerator`` which in this case uses ``LHC-Matchbox`` as the name of the generator and subsequently names the runfile ``LHC-Matchbox.run``. In order for Herwig-Parallel to work properly you have to call ``herwig-parallel-read`` with the option ``-g LHC-Matchbox`` in this particular case. Please check your Herwig++ infile to find out the appropriate generator name.

In case you are using the command ``herwig-parallel`` you have to specify the correct generator name in the run definition file (such as ``MyRun.conf``, see section :ref:`batch` for further details).
