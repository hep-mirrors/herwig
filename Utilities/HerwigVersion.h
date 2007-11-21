// -*- C++ -*-
//
// HerwigVersion.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
