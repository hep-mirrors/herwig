Overview
========

Herwig-Parallel provides the following commands and functionality:

  ============================= =============================================================================================
  Command                       Description
  ============================= =============================================================================================
  herwig-parallel-read          Prepare a run by doing a read step where the generator is set up and the Monte Carlo
                                integration is warmed up. The read step may be done in parallel.
  herwig-parallel-run           Start a parallel run consisting of an arbitrary number of individual jobs and events per job.
  herwig-parallel-status        Monitor the progress of the run and the status of the jobs it consists of.
                                Monitoring information such as queue and hostnames, status, progress and the estimated
                                total cross section and the error are displayed for each job and a preliminary combined
                                total cross section are displayed. Furthermore it is possible to configure consistency
                                checks.
  herwig-parallel-monitor       Continously monitor the status of a run by calling ``herwig-parallel-status`` over and
                                over again after a certain time interval.
  herwig-parallel-addjobs       Add an arbitrary number of jobs to an already existing run to increase statistics.
  herwig-parallel-abort         Abort either a subset or all of the jobs of a parallel run.
  herwig-parallel-restart       Restart either a subset or all of the jobs of a parallel run.
  herwig-parallel-combine       Combine the jobs of a parallel run to get the complete total cross section and error
                                estimate as well as the combined histograms. It is possible to merge either all the jobs of
                                a parallel run or only a subset of them by specifying a list of job numbers and/or job
                                number ranges to include in or exclude from the combination.
  herwig-parallel-compress      Compress all the job folders that contain all files related to the individual parallel jobs
                                into a ``.tar.gz`` archive. This can save a significant amount of disk space.
  herwig-parallel-uncompress    This is the counterpart of the ``Herwig-parallel-compress.py`` script and extracts the job
                                folders from the archive again.
  ============================= =============================================================================================

Additionally, the command ``herwig-parallel`` is provided which takes the name of a run definition file (see ``examples/testrun.conf``
for a working example) as an argument and combines the commands ``herwig-parallel-read`` and ``herwig-parallel-run`` for batch use
(see section :ref:`batch`).
