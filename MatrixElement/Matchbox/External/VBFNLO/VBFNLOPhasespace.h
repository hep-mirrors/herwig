// -*- C++ -*-
//
// VBFNLOPhasespace.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_VBFNLOPhasespace_H
#define Herwig_VBFNLOPhasespace_H
//
// This is the declaration of the VBFNLOPhasespace class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Michael Rauch
 *
 * \brief VBFNLOPhasespace is an interface to the internal phasespace generator
 * of VBFNLO. It uses the information passed via the BLHA interface to obtain
 * information on the required channels.
 *
 * @see \ref VBFNLOPhasespaceInterfaces "The interfaces"
 * defined for VBFNLOPhasespace.
 */
class VBFNLOPhasespace: public MatchboxPhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOPhasespace();

  /**
   * The destructor.
   */
  virtual ~VBFNLOPhasespace();
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
   * Generate a phase space point and return its weight.
   */
  virtual double generateKinematics(const double* random,
				    vector<Lorentz5Momentum>& momenta) {
    return generateTwoToNKinematics(random,momenta);
  }

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDimPhasespace(int nFinal) const;

  /**
   * Return true, if this phasespace generator will generate incoming
   * partons itself.
   */
  virtual bool haveX1X2() const { return true; }

  /**
   * Return true, if this phase space generator expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS() const { return false; }

  /**
   * Return true, if this phase space generator is invertible
   */
  virtual bool isInvertible() const { return false; }


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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}


private:

  /**
   * Last hadronic centre-of-mass energy, to update VBFNLO value only
   * when this has changed.
   */
  Energy lastSqrtS;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOPhasespace & operator=(const VBFNLOPhasespace &);

  /**
   * Whether or not we need to reshuffle.
   */
  bool needToReshuffle;

protected:

  /**
   * Location of the VBFNLO library
   */
  string VBFNLOlib_;

  /**
   *  load the VBFNLO library
   */
  void loadVBFNLO();

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

};

}

#endif /* Herwig_VBFNLOPhasespace_H */
