Caveats, Bugs and Known Issues
==============================

.. warning::
   Please make sure you only use the 'CellGridSampler' but no other sampler in conjunction with Herwig-Parallel. The cell grid sampler does all adapting in the read step and does not
   change in the run step. The continuous adaption during the run step of other samplers compromises the correct working of the parallelization.

.. note::
   At the moment unfortunately only Rivet 1.9.0 is supported. This is due to the fact that currently the development and use of patches for specific Rivet versions is necessary
   in order to have the required functionality. Hopefully, this procedure will become obsolete with future releases of Rivet.


Bugs and Known Issues
---------------------

None.
