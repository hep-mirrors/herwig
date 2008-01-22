// -*- C++ -*-
//
// Herwig.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_H
#define HERWIG_H
#include "ThePEG/Config/ThePEG.h"
#include <string>

namespace Herwig {
  struct HerwigVersion {

    /**
     * The version string
     */
    static std::string versionstring;

    /**
     *  the location of data directory
     */
    static std::string pkgdatadir;
  };
}

// #define HERWIG_NO_DEBUG  1

// Debugging in Herwig may be switched off completely by this compilation 
// switched, eliminating possible overhead in error checking.
#ifndef HERWIG_NO_DEBUG
#define HERWIG_DEBUG_LEVEL Herwig::HwDebug::level
#else
#define HERWIG_DEBUG_LEVEL 0
#endif

#endif /* HERWIG_H */

