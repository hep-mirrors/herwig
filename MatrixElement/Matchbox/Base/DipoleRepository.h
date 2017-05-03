// -*- C++ -*-
//
// DipoleRepository.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleRepository_H
#define HERWIG_DipoleRepository_H

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Pointer/Ptr.h"

namespace Herwig {

using namespace ThePEG;

#define HERWIG_MatchboxDipoles "/Herwig/MatrixElements/Matchbox/Dipoles/"
#define HERWIG_MatchboxTildes "/Herwig/MatrixElements/Matchbox/TildeKinematics/"
#define HERWIG_MatchboxInsertionIOperators "/Herwig/MatrixElements/Matchbox/InsertionIOperators/"
#define HERWIG_MatchboxInsertionPKOperators "/Herwig/MatrixElements/Matchbox/InsertionPKOperators/"

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Repository of known subtraction dipoles.
 */
class DipoleRepository {

public:

  /**
   * Return the known dipoles
   */
  static const vector<Ptr<SubtractionDipole>::ptr>& dipoles(int id) {
    return theDipoles(id);
  }

  /**
   * Return the known I insertion operators
   */
  static const vector<Ptr<MatchboxInsertionOperator>::ptr>& insertionIOperators(int id) {
    return theInsertionIOperators(id);
  }

  /**
   * Return the known PK insertion operators
   */
  static const vector<Ptr<MatchboxInsertionOperator>::ptr>& insertionPKOperators(int id) {
    return theInsertionPKOperators(id);
  }

public:

  /**
   * Register a dipole with associated tilde kinematics
   */
  template<int id, class DipoleT, class TildeKinematicsT, class InvertedTildeKinematicsT>
  static void registerDipole(string name, string tildeName, string invertedTildeName) {

    setup();

    BaseRepository::PushDirectory(HERWIG_MatchboxTildes);

    typename Ptr<TildeKinematicsT>::ptr tilde;
    if ( !BaseRepository::GetPointer(HERWIG_MatchboxTildes + tildeName) ) {
      tilde = new_ptr(TildeKinematicsT());
      BaseRepository::Register(tilde,tildeName);
    } else {
      tilde = 
	dynamic_ptr_cast<typename Ptr<TildeKinematicsT>::ptr>
	(BaseRepository::GetPointer(HERWIG_MatchboxTildes + tildeName));
    }

    typename Ptr<InvertedTildeKinematicsT>::ptr itilde;
    if ( !BaseRepository::GetPointer(HERWIG_MatchboxTildes + invertedTildeName) ) {
      itilde = new_ptr(InvertedTildeKinematicsT());
      BaseRepository::Register(itilde,invertedTildeName);
    } else {
      itilde = 
	dynamic_ptr_cast<typename Ptr<InvertedTildeKinematicsT>::ptr>
	(BaseRepository::GetPointer(HERWIG_MatchboxTildes + invertedTildeName));
    }

    BaseRepository::PopDirectory();
    BaseRepository::PushDirectory(HERWIG_MatchboxDipoles);

    typename Ptr<DipoleT>::ptr dip = new_ptr(DipoleT());
    dip->tildeKinematics(tilde);
    dip->invertedTildeKinematics(itilde);
    BaseRepository::Register(dip,name);

    theDipoles(id).push_back(dip);

    BaseRepository::PopDirectory();

  }

  /**
   * Register an I insertion operator
   */
  template<int id, class InsertionOperatorT>
  static void registerInsertionIOperator(string name) {

    setup();

    BaseRepository::PushDirectory(HERWIG_MatchboxInsertionIOperators);
    typename Ptr<InsertionOperatorT>::ptr iop = new_ptr(InsertionOperatorT());
    BaseRepository::Register(iop,name);

    theInsertionIOperators(id).push_back(iop);

    BaseRepository::PopDirectory();

  }

  /**
   * Register an PK insertion operator
   */
  template<int id, class InsertionOperatorT>
  static void registerInsertionPKOperator(string name) {

    setup();

    BaseRepository::PushDirectory(HERWIG_MatchboxInsertionPKOperators);
    typename Ptr<InsertionOperatorT>::ptr iop = new_ptr(InsertionOperatorT());
    BaseRepository::Register(iop,name);

    theInsertionPKOperators(id).push_back(iop);

    BaseRepository::PopDirectory();

  }

private:

  /**
   * The known dipoles
   */
  static vector<Ptr<SubtractionDipole>::ptr>& theDipoles(int id);

  /**
   * The known I insertion operators
   */
  static vector<Ptr<MatchboxInsertionOperator>::ptr>& theInsertionIOperators(int id);

  /**
   * The known PK insertion operators
   */
  static vector<Ptr<MatchboxInsertionOperator>::ptr>& theInsertionPKOperators(int id);

  /**
   * True, if initialized.
   */
  static bool& initialized();

  /**
   * Setup directories in the repository
   */
  static void setup();

};

}

#endif // HERWIG_DipoleRepository_H
