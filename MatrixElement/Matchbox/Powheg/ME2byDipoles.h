// -*- C++ -*-
//
// ME2byDipoles.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ME2byDipoles_H
#define HERWIG_ME2byDipoles_H
//
// This is the declaration of the ME2byDipoles class.
//

#include "ThePEG/Handlers/LastXCombInfo.h"

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxReweightBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig++/MatrixElement/Matchbox/Base/SubtractedME.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ME2byDipoles is the base class for all
 * quantities of type |M_R|^2 / \sum D
 *
 */
class ME2byDipoles: public MatchboxReweightBase, public LastXCombInfo<StandardXComb> {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ME2byDipoles();

  /**
   * The destructor.
   */
  virtual ~ME2byDipoles();
  //@}

public:

  /**
   * Set the XComb object steering the real emission.
   */
  virtual void setXComb(tStdXCombPtr real);

  /**
   * Return true, if applies to the process in the xcomb.
   */
  virtual bool apply() const { return projectionDipole()->apply(); }

  /**
   * Evaluate the ratio.
   */
  virtual double evaluate() const {
    double dummy; return evaluate(dummy);
  }

  /**
   * Evaluate the ratio.
   */
  double evaluate(double&) const;

public:

  /**
   * Return a scaled Born screening cross section
   */
  double scaledBornScreen() const;

  /**
   * Return a scaled Born cross section
   */
  double scaledBorn(Energy2 factorizationScale = ZERO) const;

protected:

  /**
   * Setup dependent xcombs for a yet unknown
   * real emission xcomb.
   */
  void getXCombs(tStdXCombPtr);

public:

  /**
   * Return the real emission matrix element
   */
  Ptr<MatchboxMEBase>::tcptr realME() const { return theRealME; }

  /**
   * Return the real emission matrix element
   */
  Ptr<MatchboxMEBase>::tptr realME() { return theRealME; }

  /**
   * Set the real emission matrix element
   */
  void realME(Ptr<MatchboxMEBase>::tptr me) { theRealME = me; }

  /**
   * Return the subtraction dipoles
   */
  const vector<Ptr<SubtractionDipole>::ptr>& dipoles() const { return theDipoles; }

  /**
   * Access the subtraction dipoles
   */
  vector<Ptr<SubtractionDipole>::ptr>& dipoles() { return theDipoles; }

  /**
   * Setup from a subtracted matrix element, calculating
   * |M_R|^2 / Dipoles
   */
  void setup(Ptr<SubtractedME>::tptr);

  /**
   * Setup from a dipole and subtracted matrix element and
   * dipole calculating dip / Dipoles
   */
  void setup(Ptr<SubtractionDipole>::tptr dip,
	     Ptr<SubtractedME>::tptr);

  /**
   * Set the dipole to obtain Born cross sections for
   * further reweighting.
   */
  void projectionDipole(Ptr<SubtractionDipole>::tptr dip) { theProjectionDipole = dip; }

  /**
   * Return the dipole to obtain Born cross sections for
   * further reweighting.
   */
  Ptr<SubtractionDipole>::tcptr projectionDipole() const { return theProjectionDipole; }

  /**
   * Return the dipole to obtain Born cross sections for
   * further reweighting.
   */
  Ptr<SubtractionDipole>::tptr projectionDipole() { return theProjectionDipole; }

public:

  /**
   * Dump the setup to an ostream
   */
  virtual void print(ostream&) const;

  /**
   * Print debug information on the last event
   */
  virtual void printLastEvent(ostream&) const;

  /**
   * Dump xcomb hierarchies.
   */
  void dumpInfo(const string& prefix = "") const;

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The real emission matrix element.
   */
  Ptr<MatchboxMEBase>::ptr theRealME;

  /**
   * A dipole to obtain Born cross sections for
   * further reweighting.
   */
  Ptr<SubtractionDipole>::ptr theProjectionDipole;

  /**
   * The dipoles to be considered.
   */
  vector<Ptr<SubtractionDipole>::ptr> theDipoles;

  /**
   * A map of real xcombs to Born xcombs to
   * be used by the dipoles.
   */
  map<StdXCombPtr,vector<StdDependentXCombPtr> > theXCombMap;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ME2byDipoles & operator=(const ME2byDipoles &);

};

}

#endif /* HERWIG_ME2byDipoles_H */
