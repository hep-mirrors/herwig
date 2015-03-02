Configuration
=============

Herwig-Parallel comes with three configuration files which reside in the ``<Herwig-source-folder>/Contrib/Herwig-Parallel/config`` directory.

- ``herwig-parallel.conf``
- ``clusters.conf``
- ``queues.conf``

General Configuration File: ``herwig-parallel.conf``
----------------------------------------------------

::

    [tools]

    mergeGrids        = mergeGrids.py
    combineRuns       = combineRuns
    makeDistributions = makeDistributions
    rivet             = rivet
    flat2aida         = flat2aida


    [defaults]

    ## --------------------
    ## herwig-parallel-read
    ## --------------------
    generator       = Herwig
    buildscript     = herwig-parallel-build.sh
    buildqueue      = ITP-short
    integratescript = herwig-parallel-integrate.sh
    integratequeue  = ITP-short
    integratejobs   = 1

    ## -----------------------------------------
    ## herwig-parallel-run / -addjobs / -restart
    ## -----------------------------------------
    runscript = herwig-parallel-run.sh
    runqueue  = ITP-all
    runjobs   = 1
    events    = 100k
    seedmode  = ascending
    seed      = 1

    ## -----------------------
    ## herwig-parallel-combine
    ## -----------------------
    builtin-analysis     = 
    hepmc-rivet-analyses = 


    [consistency-checks]

    standard-deviations = 4.0
    error-factor        = 2.0
  
- ``[tools]``
  
  In this section the commands for merging the grid files of the sampling algorithm, combining the histograms from different runs and
  creating the distributions for plotting are specified. The corresponding programs and scripts are located in the ``bin``-folder of
  your Herwig++ installation.

  Additionally, the commands to call ``Rivet`` and ``flat2aida`` which is part of Rivet need to be provided.

  In case the commands are located in your ``PATH`` environment variable the above commands will be sufficient.
  Otherwise the full path has to be specified.

- ``[defaults]``

  This section contains the default values for some of the options or arguments that may be given to the various Herwig-Parallel commands.
  Altering the values in this section of the configuration file may save you from specifying the same option over and over again.

- ``[consistency-checks]``
  
  Configuration of the consistency checks that may be applied to all parallel jobs of a run
  
  - ``standard-deviations``
  
    Limit for the number of standard deviations that the cross section of each parallel job may deviate from the median cross section
    of all jobs in order to still be considered consistent.
    
  - ``error-factor``
  
    Limit for the factor by which the statistical uncertainty on the total cross section of each parallel job may deviate from the median
    statistical uncertainty in order to still be considered consistent.


Cluster Configuration File: ``clusters.conf``
---------------------------------------------

::
  
    [ITP]
    jobid = awk '{printf "%%s", $3}'
    joblist = qstat.py
    status = awk '{printf "%%s", $3}'
    statusqueued = qw
    statusrunning = r
    abort = qdel @JOBID@
    
For each cluster you want to configure add a section to the cluster configuration file. The label of the section will be used as the cluster name.

- ``jobid``
  
  Command to extract the job id from the response output of the cluster's batch queue system when a job is submitted.
  
- ``joblist``
  
  Command to retrieve the list of all jobs from the clusters's batch queue.
  
- ``status``
  
  Command to extract the status of a job from the joblist generated with the ``joblist`` command.
  
- ``statusqueued``
  
  Code for queued jobs used in the joblist generated with the ``joblist`` command.
  
- ``statusrunning``
  
  Code for running jobs used in the joblist generated with the ``joblist`` command.
  
- ``abort``

  Command for aborting a certain job identified by its jobid. The placeholder ``@JOBID@`` will be replaced by Herwig-Parallel when aborting jobs.


Batch queue configuration file: ``queues.conf``
-----------------------------------------------

::

    [ITP-all]
    cluster = ITP
    submit = qsub -cwd -notify -l 'all=1 mem_free=2.5G' ./@SCRIPT@

    [ITP-albatros]
    cluster = ITP
    submit = qsub -cwd -notify -l mem_free=2.5G ./@SCRIPT@

    [ITP-short]
    cluster = ITP
    submit = qsub -cwd -notify -l 'short mem_free=2.5G' ./@SCRIPT@
  
For each queue you want to configure add a section to the queue configuration file. The label of the section will be used as the queue name.

- ``cluster``
  
  Name of the cluster the queue is associated to.
  
- ``submit``
  
  Command for the submission of jobs to the respective queue. The placeholder ``@SCRIPT@`` will be replaced by Herwig-Parallel when starting jobs.

