// -*- C++ -*-
//
// ShowerApproximation.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ShowerApproximation_H
#define Herwig_ShowerApproximation_H
//
// This is the declaration of the ShowerApproximation class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ShowerApproximation describes the shower emission to be used
 * in NLO matching.
 *
 */
class ShowerApproximation: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ShowerApproximation();

  /**
   * The destructor.
   */
  virtual ~ShowerApproximation();
  //@}

public:

  /**
   * Set the XComb object describing the Born process
   */
  void setBornXComb(tStdXCombPtr xc) { theBornXComb = xc; }

  /**
   * Return the XComb object describing the Born process
   */
  tStdXCombPtr bornXComb() const { return theBornXComb; }

  /**
   * Set the XComb object describing the real emission process
   */
  void setRealXComb(tStdXCombPtr xc) { theRealXComb = xc; }

  /**
   * Return the XComb object describing the real emission process
   */
  tStdXCombPtr realXComb() const { return theRealXComb; }

  /**
   * Set the dipole in charge for the emission
   */
  void setDipole(Ptr<SubtractionDipole>::tcptr);

  /**
   * Return the dipole in charge for the emission
   */
  Ptr<SubtractionDipole>::tcptr dipole() const;

public:

  /**
   * Return true if one of the recently encountered configutations was
   * below the infrared cutoff.
   */
  bool belowCutoff() const { return theBelowCutoff; }

  /**
   * Indicate that one of the recently encountered configutations was
   * below the infrared cutoff.
   */
  void wasBelowCutoff() { theBelowCutoff = true; }

  /**
   * Reset the below cutoff flag.
   */
  void resetBelowCutoff() { theBelowCutoff = false; }

public:

  /**
   * Return true, if the shower was able to generate an emission
   * leading from the given Born to the given real emission process.
   */
  virtual bool isInShowerPhasespace() const = 0;

  /**
   * Return true, if the shower emission leading from the given Born
   * to the given real emission process would have been generated
   * above the shower's infrared cutoff.
   */
  virtual bool isAboveCutoff() const = 0;

  /**
   * Return the shower approximation to the real emission cross
   * section for the given pair of Bron and real emission
   * configurations.
   */
  virtual CrossSection dSigHatDR() const = 0;

  /**
   * Return the shower approximated real emission cross section
   * divided by the Born cross section and the radiation phasespace
   * jacobian for the given pair of Born and real emission
   * configurations.
   */
  virtual double me2() const = 0;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The XComb object describing the Born process
   */
  tStdXCombPtr theBornXComb;

  /**
   * The XComb object describing the real emission process
   */
  tStdXCombPtr theRealXComb;

  /**
   * The dipole in charge for the emission
   */
  Ptr<SubtractionDipole>::tcptr theDipole;

  /**
   * True if one of the recently encountered configutations was below
   * the infrared cutoff.
   */
  bool theBelowCutoff;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerApproximation & operator=(const ShowerApproximation &);

};

}

#endif /* Herwig_ShowerApproximation_H */
