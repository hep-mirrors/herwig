// -*- C++ -*-
//
// Herwig.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2015 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_UI_H
#define SRC_HERWIG_UI_H

#include <vector>
#include <string>

namespace Herwig {

/**
 * Definition of the possible run modes of Herwig
 */
namespace RunMode {
  enum Mode { ERROR, INIT, READ, BUILD, INTEGRATE, RUN };
}

/**
 * HerwigUI is an interface to abstract the command line parameters.
 * This allows an inherited class to configure Herwig wihtout actually 
 * having to interact via main()
 */
 class HerwigUI {
  public:

    /// Requested Herwig run mode
    virtual RunMode::Mode runMode() const = 0;

    bool integrationJob() const {
        return runMode() == RunMode::INTEGRATE;
    }

    /// Repository name to operate on
    virtual std::string repository() const = 0;

    /// Name of the file to be read
    virtual std::string inputfile() const = 0;

    /// Name of the setup file to be read, to modify the repository
    virtual std::string setupfile() const = 0;

    //bool readInWasSuccessful() const = 0;
    virtual bool resume() const = 0;
    virtual bool tics() const = 0;
    virtual std::string tag() const = 0;

    virtual std::string integrationList() const = 0;

    /// Directories from which Herwig reads input files,  will be prepended to the search path.
    virtual const std::vector<std::string> & prependReadDirectories() const = 0;

    /// Directories from which Herwig reads input files,  will be appended to the search path.
    virtual const std::vector<std::string> & appendReadDirectories() const = 0;

    virtual long N() const = 0;
    virtual int seed() const = 0;
    virtual int jobs() const = 0;
    virtual unsigned int jobSize() const = 0;
    virtual unsigned int maxJobs() const = 0;

    virtual void quitWithError() const = 0;

    virtual ~HerwigUI() {}

 };

}

#endif
