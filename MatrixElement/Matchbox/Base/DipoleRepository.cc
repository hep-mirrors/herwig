// -*- C++ -*-
//
// DipoleRepository.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "DipoleRepository.h"

using namespace Herwig;


vector<Ptr<SubtractionDipole>::ptr>& DipoleRepository::theDipoles(int id) {
  static map<int,vector<Ptr<SubtractionDipole>::ptr> > theDipoles_;
  return theDipoles_[id];
}

vector<Ptr<MatchboxInsertionOperator>::ptr>& DipoleRepository::theInsertionIOperators(int id) {
  static map<int,vector<Ptr<MatchboxInsertionOperator>::ptr> > theInsertionIOperators_;
  return theInsertionIOperators_[id];
}

vector<Ptr<MatchboxInsertionOperator>::ptr>& DipoleRepository::theInsertionPKOperators(int id) {
  static map<int,vector<Ptr<MatchboxInsertionOperator>::ptr> > theInsertionPKOperators_;
  return theInsertionPKOperators_[id];
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
    BaseRepository::CheckDirectory(HERWIG_MatchboxInsertionIOperators);
  } catch (RepositoryNoDirectory& d) {
    d.handle();
    BaseRepository::CreateDirectory(HERWIG_MatchboxInsertionIOperators);
  }

  try {
    BaseRepository::CheckDirectory(HERWIG_MatchboxInsertionPKOperators);
  } catch (RepositoryNoDirectory& d) {
    d.handle();
    BaseRepository::CreateDirectory(HERWIG_MatchboxInsertionPKOperators);
  }

  try {
    BaseRepository::CheckDirectory(HERWIG_MatchboxTildes);
  } catch (RepositoryNoDirectory& d) {
    d.handle();
    BaseRepository::CreateDirectory(HERWIG_MatchboxTildes);
  }

  initialized() = true;

}
