How to Use
==========

For using Herwig-Parallel to create parallel Herwig++ runs four files are necessary:

- The **infile** for Herwig++, e.g. ``LHC-Matchbox.in``, that is needed for any normal Herwig++ run anyway.
  
  However, please ensure that the following settings apply:
     
  - ``set /Herwig/Samplers/Sampler:AddUpSamplers Off``
    
    This is necessary to ensure that the cross section and error are only calculated from 'Herwig run' step and not from both 'Herwig read' and 'Herwig run' as this would compromise the statistical independence of the parallel runs.
  
  - ::
    
       create Herwig::ParallelRunAnalysis /Herwig/Analysis/ParallelRunAnalysis HwAnalysis.so
       insert /Herwig/Generators/LHCGenerator:AnalysisHandlers 0 /Herwig/Analysis/ParallelRunAnalysis
    
    These lines load an analysis handler that writes status information concerning the Herwig run to a file. This information is necessary for later on combining the jobs of the parallel run.
  
  - ::
  
       set /Herwig/Samplers/Sampler:BinSampler /Herwig/Samplers/CellGridSampler
       set /Herwig/Samplers/CellGridSampler:UseAllIterations No
    
  .. warning::
     Please make sure you only use the 'CellGridSampler' but no other sampler in conjunction with Herwig-Parallel. The cell grid sampler does all adapting in the read step and does not
     change in the run step. The continuous adaption of other samplers during the run step compromises the correct working of the parallelization.
     
- A bash script, referred to as '**build script**', that conducts the Herwig++ 'build' step.

- A bash script, referred to as '**integrate script**', that conducts the Herwig++ 'integrate' step.

- A bash script, referred to as '**run script**', that conducts the Herwig++ 'run' step.
  
Working examples for the infile and the three scripts are given in the ``examples`` directory of Herwig-Parallel.

Note the keywords ``@INFILE@``, ``@JOBS@``, ``@RUNFILE@``, ``@JOBID@``, ``@SETUPFILE@``, ``@EVENTS@``, ``@SEED@`` and ``@CLEANUP@`` in the shell scripts which are replaced by Herwig-Parallel when starting parallel runs.
``@CLEANUP@`` removes unnecessary files and folders after completion of the job in order to save disk space. Since runs with hundred or more jobs easily take up large amounts of disk space
it is strongly recommended to include this command at the end of the run script. Furthermore, all of the scripts may contain other instructions, e.g. source ``.bashrc`` or set the Rivet environment variables.

To execute the Herwig-Parallel commands enter the folder where you want to create a parallel run in a terminal and run the desired command.

- First, run ``herwig-parallel-read`` to set up a parallel run. This calls the build script which triggers the generation of the matrix elements and prepares the integration jobs.
  For each parallel integration job the integrate script is called. Once all integration jobs are finished, the resulting grid files are merged and ``herwig-parallel-read`` indicates
  the successful completion of all tasks.
  
- Next, run ``herwig-parallel-run`` to start the desired number of parallel jobs.

- You may monitor the status of your calculation with ``herwig-parallel-monitor`` or ``herwig-parallel-status``, add additional jobs to the run with ``herwig-parallel-addjobs`` and
  abort and restart some or all of the jobs with ``herwig-parallel-abort`` and ``herwig-parallel-restart``.
  
- Once your jobs are finished run ``herwig-parallel-combine`` to combine the results of all parallel jobs.

- Lastly, you may compress the run with ``herwig-parallel-compress`` in order to save disk space. This compresses all the run folders but the combined output and the used infiles will still
  be accessible. If you decide to run more statistics you can uncompress the run again with ``herwig-parallel-uncompress`` and do ``herwig-parallel-addjobs``.

Note that you don't need to have the same build, integrate and run scripts in all of your run folders. By passing the option ``-d`` to ``herwig-parallel-read`` and
``herwig-parallel-run`` you can use the default scripts that are located in the ``misc`` folder of the Herwig-Parallel directory tree. The names of the respective scripts are configurable
via the file ``config/herwig-parallel.conf``.
