// -*- C++ -*-
//
// MatchboxRambo.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxRambo_H
#define Herwig_MatchboxRambo_H
//
// This is the declaration of the MatchboxRambo class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxRambo implements RAMBO phase space generation.
 *
 */
class MatchboxRambo: public MatchboxPhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxRambo();

  /**
   * The destructor.
   */
  virtual ~MatchboxRambo();
  //@}

public:

  /**
   * Prepare a phase space generator for the given xcomb object.
   */
  virtual void setXComb(tStdXCombPtr);

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateTwoToNKinematics(const double*,
					  vector<Lorentz5Momentum>& momenta);

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDimPhasespace(int nFinal) const {
    if ( nFinal == 1 )
      return 1;
    return 4*nFinal;
  }

protected:

  /**
   * The function object defining the equation
   * to be solved.
   */
  struct ReshuffleEquation {

    typedef double ArgType;
    typedef Energy ValType;

    static double aUnit() { return 1.; }
    static Energy vUnit() { return 1.*GeV; }

    Energy operator() (double xi) const;

    Energy w;
    cPDVector::const_iterator dataBegin;
    cPDVector::const_iterator dataEnd;
    vector<Lorentz5Momentum>::const_iterator momentaBegin;
    vector<Lorentz5Momentum>::const_iterator momentaEnd;

    ReshuffleEquation(Energy q,
		      cPDVector::const_iterator dBegin,
		      cPDVector::const_iterator dEnd,
		      vector<Lorentz5Momentum>::const_iterator mBegin,
		      vector<Lorentz5Momentum>::const_iterator mEnd)
      : w(q),
	dataBegin(dBegin), dataEnd(dEnd),
	momentaBegin(mBegin),
	momentaEnd(mEnd) {}

  };

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxRambo & operator=(const MatchboxRambo &);

  /**
   * Whether or not we need to reshuffle.
   */
  bool needToReshuffle;

  /**
   * True, if a reference sample of phasespace points should be
   * generated.
   */
  bool theMakeReferenceSample;

  /**
   * Map processes to streams for reference samples
   */
  map<cPDVector,ofstream*> referenceSamples;

  /**
   * The stream to fill for the reference sample
   */
  ofstream* referenceSample;

  /**
   * Write the generated point to the reference sample
   */
  void dumpReference(const vector<Lorentz5Momentum>&, double) const;

};

}

#endif /* Herwig_MatchboxRambo_H */
