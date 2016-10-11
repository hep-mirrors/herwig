// -*- C++ -*-
//
// HerwigUI.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2016 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_UI_H
#define SRC_HERWIG_UI_H

#include <vector>
#include <string>
#include <iosfwd>

namespace Herwig {

/**
 * Definition of the possible run modes of Herwig
 */
namespace RunMode {
  enum Mode { ERROR, INIT, READ, BUILD, INTEGRATE, MERGEGRIDS, RUN };
}

/**
 * HerwigUI is an interface to abstract the command line parameters.
 *
 * This allows any inheriting class to configure Herwig wihtout actually 
 * having to interact via main()
 */
 class HerwigUI {

 public:

   /// Requested Herwig run mode
   virtual RunMode::Mode runMode() const = 0;

   /// Repository name to operate on
   virtual std::string repository() const = 0;

   /// Name of the file to be read
   virtual std::string inputfile() const = 0;

   /// Name of the setup file to be read, to modify the repository
   virtual std::string setupfile() const = 0;

   /// Try to resume execution from an earlier interrupted run.
   virtual bool resume() const = 0;

   /// Require verbose progress markers
   virtual bool tics() const = 0;

   /// A user-defined tag to append to the run name.
   virtual std::string tag() const = 0;

   /// An identifier for the integration job to be handled
   virtual std::string integrationList() const = 0;

   /// Directories from which Herwig reads input files,  will be prepended to the search path.
   virtual const std::vector<std::string> & prependReadDirectories() const = 0;

   /// Directories from which Herwig reads input files,  will be appended to the search path.
   virtual const std::vector<std::string> & appendReadDirectories() const = 0;

   /// The number of events to generate
   virtual long N() const = 0;

   /// The seed to use
   virtual int seed() const = 0;

   /// The number of jobs to fork
   virtual int jobs() const = 0;

   /// The number of subprocesses to integrate per integratoin job
   virtual unsigned int jobSize() const = 0;

   /// The maximum number of integration jobs
   virtual unsigned int maxJobs() const = 0;

   /// Bail out and print usage information
   virtual void quitWithHelp() const = 0;

   /// Bail out and be quiet
   virtual void quit() const = 0;

   /// Destructor
   virtual ~HerwigUI() {}

   /// Return true, if this is an integration job
   bool integrationJob() const {
     return runMode() == RunMode::INTEGRATE;
   }

   /// Return the standard out stream to be used
   virtual std::ostream& outStream() const = 0;

   /// Return the standard err stream to be used
   virtual std::ostream& errStream() const = 0;

   /// Return the standard in stream to be used
   virtual std::istream& inStream() const = 0;

 };

}

#endif
