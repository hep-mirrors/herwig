// -*- C++ -*-
//
// RunDirectories.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "RunDirectories.h"

#include "ThePEG/Utilities/Exception.h"
#include <boost/filesystem.hpp>

using namespace ThePEG;
using namespace Herwig;

void RunDirectories::prefix(string p) {
  if ( *p.rbegin() != '/' )
    p += "/";
  thePrefix() = p;
}

const string& RunDirectories::prefix() {
  return thePrefix();
}

string& RunDirectories::thePrefix() {
  static string p = "./Matchbox/";
  return p;
}

string& RunDirectories::theBuildStorage() {
  static string builds = "";
  return builds;
}

const string& RunDirectories::buildStorage() {
  if ( !theBuildStorage().empty() )
    return theBuildStorage();
  theBuildStorage() = prefix();
  if ( theBuildStorage().empty() )
    theBuildStorage() = "./Matchbox/";
  else if ( *theBuildStorage().rbegin() != '/' )
    theBuildStorage() += "/";
  theBuildStorage() += "Build/";
  if ( boost::filesystem::exists(theBuildStorage()) ) {
    if ( !boost::filesystem::is_directory(theBuildStorage()) )
      throw Exception() << "Matchbox build storage '"
			<< theBuildStorage() << "' existing but not a directory."
			<< Exception::abortnow;
  } else {
    boost::filesystem::create_directories(theBuildStorage());
  }
  return theBuildStorage();
}

bool RunDirectories::empty() {
  return theRunDirectories().empty();
}

list<string>& RunDirectories::theRunDirectories() {
  static list<string> rundirs;
  return rundirs;
}

const string& RunDirectories::runStorage() {
  if ( theRunDirectories().empty() )
    throw Exception() << "No run directory stack has been allocated."
		      << Exception::abortnow;
  if ( boost::filesystem::exists(theRunDirectories().front()) ) {
    if ( !boost::filesystem::is_directory(theRunDirectories().front()) )
      throw Exception() << "Matchbox run storage '"
			<< theRunDirectories().front() << "' existing but not a directory."
			<< Exception::abortnow;
  } else {
    boost::filesystem::create_directories(theRunDirectories().front());
  }
  return theRunDirectories().front();
}

void RunDirectories::pushRunId(string p) {
  if ( *p.rbegin() != '/' )
    p += "/";
  if ( theRunDirectories().empty() ) {
    theRunDirectories().push_front(prefix() + p);
    return;
  }
  theRunDirectories().push_front(theRunDirectories().front() + p);
}

RunDirectories::RunDirectories()
  : directoriesLeft(theRunDirectories()) {}

string RunDirectories::nextRunStorage() {
  if ( directoriesLeft.empty() )
    throw Exception() << "No more run directories left to try."
		      << Exception::abortnow;
  string res = directoriesLeft.front();
  directoriesLeft.pop_front();
  return res;
}
