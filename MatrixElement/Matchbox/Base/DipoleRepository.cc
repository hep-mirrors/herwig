// -*- C++ -*-
//
// DipoleRepository.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "DipoleRepository.h"

using namespace Herwig;


vector<Ptr<SubtractionDipole>::ptr>& DipoleRepository::theDipoles() {
  static vector<Ptr<SubtractionDipole>::ptr> theDipoles_;
  return theDipoles_;
}

vector<Ptr<MatchboxInsertionOperator>::ptr>& DipoleRepository::theInsertionOperators() {
  static vector<Ptr<MatchboxInsertionOperator>::ptr> theInsertionOperators_;
  return theInsertionOperators_;
}


bool& DipoleRepository::initialized() {
  static bool theInitialized = false;
  return theInitialized;
}

void DipoleRepository::setup() {

  if ( initialized() )
    return;

  try {
    BaseRepository::CheckDirectory(HERWIG_MatchboxDipoles);
  } catch (RepositoryNoDirectory& d) {
    d.handle();
    BaseRepository::CreateDirectory(HERWIG_MatchboxDipoles);
  }

  try {
    BaseRepository::CheckDirectory(HERWIG_MatchboxInsertionOperators);
  } catch (RepositoryNoDirectory& d) {
    d.handle();
    BaseRepository::CreateDirectory(HERWIG_MatchboxInsertionOperators);
  }

  try {
    BaseRepository::CheckDirectory(HERWIG_MatchboxTildes);
  } catch (RepositoryNoDirectory& d) {
    d.handle();
    BaseRepository::CreateDirectory(HERWIG_MatchboxTildes);
  }

  initialized() = true;

}
