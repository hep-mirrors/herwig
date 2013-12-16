// -*- C++ -*-
//
// NLOJetPhasespace.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_NLOJetPhasespace_H
#define HERWIG_NLOJetPhasespace_H
//
// This is the declaration of the NLOJetPhasespace class.
//

#include "Herwig++/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

#include "NLOJetRandomWrapper.h"

#include "nlo++/bits/hep-lorentzvector.h"
#include "nlo++/bits/psg-phasespace.h"
#include "nlo++/bits/nlo-phasespace.h"
#include "nlo++/bits/nlo-event.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Convert nlo::lorentzvector<double> to Lorentz5Vector<Energy>,
 * assuming nlo::lorentzvector<double> is given in units og GeV
 */
struct NLOMomentumConverter {
  Lorentz5Momentum operator()(const nlo::lorentzvector<double>& v) const {
    return Lorentz5Momentum(v.X()*GeV,v.Y()*GeV,v.Z()*GeV,v.T()*GeV);
  }
  nlo::lorentzvector<double> operator()(const Lorentz5Momentum& v) const {
    return nlo::lorentzvector<double>(v.x()/GeV,v.y()/GeV,v.z()/GeV,v.t()/GeV);
  }
};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Keys for XComb meta information
 */
struct NLOMetaKeys {

  enum Keys {
    HadronicEvent
  };

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief NLOJetPhasespace provides an interface to 
 * nlojet's phasespace generator classes.
 *
 * @see \ref NLOJetPhasespaceInterfaces "The interfaces"
 * defined for NLOJetPhasespace.
 */
template<unsigned int N, unsigned int I, unsigned int F>
class NLOJetPhasespace: public Herwig::MatchboxPhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NLOJetPhasespace();

  /**
   * The destructor.
   */
  virtual ~NLOJetPhasespace();
  //@}

public:

  /**
   * Define the event type.
   */
  typedef typename nlo::hadronic_event<nlo::lorentzvector<double>,
				       nlo::hadronic_event_traits<N,I,F> > NLOEvent;

  /**
   * Define the type of the nlojet phasespace generator.
   */
  typedef typename nlo::basic_phasespace<NLOEvent> NLOPhasespace;

  /**
   * Get the phasepsace generator to wrap around.
   */
  NLOPhasespace* nloPhasespace() const { return theNLOPhasespace; }

  /**
   * Prepare a phase space generator for the given xcomb object.
   */
  virtual void prepare(tStdXCombPtr xc, bool verbose = false);

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateTwoToNKinematics(const double* r,
					  vector<Lorentz5Momentum>& momenta);

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDim(int nFinal) const;

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
   * Fill a diagram selector for the last phase space point.
   */
  virtual Selector<MEBase::DiagramIndex> selectDiagrams(const MEBase::DiagramVector&) const {
    return Selector<MEBase::DiagramIndex>();
  }

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
   * The number of outgoing partons from which on dipole generation is
   * used.
   */
  unsigned int nHard;

  /**
   * The phasepsace generator to wrap around.
   */
  NLOPhasespace* theNLOPhasespace;

  /**
   * The random number wrapper
   */
  NLOJetRandomWrapper nloRnd;

  /**
   * The event to be filled
   */
  NLOEvent* nloEvent;

  /**
   * The hard event to start with
   */
  NLOEvent* nloHardEvent;

  /**
   * The XCombPtr to nloEvent pointer map
   */
  map <  XCPtr, NLOEvent * > theXCtoEvents;
  /**
   * The XCombPtr to nloHardEvent pointer map
   */
  map < XCPtr, NLOEvent * > theXCtoHardEvents;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NLOJetPhasespace & operator=(const NLOJetPhasespace &);

};

}

#include "NLOJetPhasespace.tcc"

#endif /* HERWIG_NLOJetPhasespace_H */
