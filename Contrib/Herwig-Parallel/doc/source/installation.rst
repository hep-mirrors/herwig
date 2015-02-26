Installation and Setup
======================

1. Preparations
   
   - Place the folder ``<Herwig-source-folder>/Contrib/Herwig-Parallel/source`` in your ``$PATH`` environment variable in order to enable easier access to the ``Herwig-Parallel`` commands.
   
2. Prepare Rivet 1.9.0

   .. note::
      At the moment unfortunately only version 1.9.0 of Rivet is supported. This is due to the fact that currently the development and use of patches for specific Rivet versions is necessary
      in order to have the required functionality. Hopefully, this procedure will become obsolete with future releases of Rivet.
   
   - Copy the following files to the following destinations in your Rivet source location
   
      - ``patches/Rivet/AIHistogram1D.h`` to ``<your-Rivet-source-folder>/include/LWH``
      - ``patches/Rivet/Histogram1D.h`` to ``<your-Rivet-source-folder>/include/LWH``
      - ``patches/Rivet/AnalysisHandler.hh`` to ``<your-Rivet-source-folder>/include/Rivet``
      - ``patches/Rivet/Analysis.hh`` to ``<your-Rivet-source-folder>/include/Rivet``
      - ``patches/Rivet/AnalysisHandler.cc`` to ``<your-Rivet-source-folder>/src/Core``
      - ``patches/Rivet/Analysis.cc`` to ``<your-Rivet-source-folder>/src/Core``
     
     and do ``make -j4 install``.
   
   - Put the command ``source <your-Rivet-source-folder>/rivetenv.sh`` in your ``.bashrc``. In case of problems it may be necessary to comment out all LaTeX-related lines.

3. Prepare ThePEG

   - Copy the files
      
      - ``patches/ThePEG/NLORivetAnalysis.cc``
      - ``patches/ThePEG/RivetAnalysis.cc``
     
     to the ``Analysis`` folder of your ThePEG source directory and do ``make -j4 install``.

