// -*- C++ -*-
//
// DipoleRepository.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleRepository_H
#define HERWIG_DipoleRepository_H

#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig++/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Pointer/Ptr.h"

namespace Herwig {

using namespace ThePEG;

#define HERWIG_MatchboxDipoles "/Herwig/MatrixElements/Matchbox/Dipoles/"
#define HERWIG_MatchboxTildes "/Herwig/MatrixElements/Matchbox/TildeKinematics/"
#define HERWIG_MatchboxInsertionOperators "/Herwig/MatrixElements/Matchbox/InsertionOperators/"

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
  static const vector<Ptr<SubtractionDipole>::ptr>& dipoles() {
    return theDipoles();
  }

  /**
   * Return the known insertion operators
   */
  static const vector<Ptr<MatchboxInsertionOperator>::ptr>& insertionOperators() {
    return theInsertionOperators();
  }

public:

  /**
   * Register a dipole with associated tilde kinematics
   */
  template<class DipoleT, class TildeKinematicsT, class InvertedTildeKinematicsT>
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

    theDipoles().push_back(dip);

    BaseRepository::PopDirectory();

  }

  /**
   * Register an insertion operator
   */
  template<class InsertionOperatorT>
  static void registerInsertionOperator(string name) {

    setup();

    BaseRepository::PushDirectory(HERWIG_MatchboxInsertionOperators);
    typename Ptr<InsertionOperatorT>::ptr iop = new_ptr(InsertionOperatorT());
    BaseRepository::Register(iop,name);

    theInsertionOperators().push_back(iop);

    BaseRepository::PopDirectory();

  }

private:

  /**
   * The known dipoles
   */
  static vector<Ptr<SubtractionDipole>::ptr>& theDipoles();

  /**
   * The known insertion operators
   */
  static vector<Ptr<MatchboxInsertionOperator>::ptr>& theInsertionOperators();

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
