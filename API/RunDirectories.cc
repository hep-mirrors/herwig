// -*- C++ -*-
//
// RunDirectories.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "RunDirectories.h"

#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Handlers/SamplerBase.h"
#include "Filesystem.h"

using namespace ThePEG;
using namespace Herwig;

using std::string;
using std::list;

void RunDirectories::prefix(string p) {
  if ( *p.rbegin() != '/' )
    p += "/";
  thePrefix() = p;
}

const string& RunDirectories::prefix() {
  return thePrefix();
}

string& RunDirectories::thePrefix() {
  static string p = "./Herwig-cache/";
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
    theBuildStorage() = "./Herwig-cache/";
  else if ( *theBuildStorage().rbegin() != '/' )
    theBuildStorage() += "/";
  theBuildStorage() += "Build/";
  if ( filesystem::exists(theBuildStorage()) ) {
    if ( !filesystem::is_directory(theBuildStorage()) )
      throw Exception() << "Herwig build storage '"
			<< theBuildStorage() << "' exists but not a directory."
			<< Exception::abortnow;
  } else {
    filesystem::create_directory(theBuildStorage());
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
  if ( filesystem::exists(theRunDirectories().front()) ) {
    if ( !filesystem::is_directory(theRunDirectories().front()) )
      throw Exception() << "Herwig run storage '"
			<< theRunDirectories().front() << "' exists but not a directory."
			<< Exception::abortnow;
  } else {
    filesystem::create_directory(theRunDirectories().front());
  }
  return theRunDirectories().front();
}

const string& RunDirectories::interfaceStorage() {
  if ( SamplerBase::runLevel() == SamplerBase::IntegrationMode ||
       SamplerBase::runLevel() == SamplerBase::RunMode )
    return runStorage();
  return buildStorage();
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

