// -*- C++ -*-
//
// PowhegFactory.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PowhegFactory_H
#define HERWIG_PowhegFactory_H
//
// This is the declaration of the PowhegFactory class.
//

#include "ThePEG/Handlers/SubProcessHandler.h"

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxNLOME.h"
#include "Herwig++/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig++/MatrixElement/Matchbox/Powheg/PowhegInclusiveME.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief PowhegFactory automatically sets up
 * a POWHEG matching
 *
 * @see \ref PowhegFactoryInterfaces "The interfaces"
 * defined for PowhegFactory.
 */
class PowhegFactory: public SubProcessHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PowhegFactory();

  /**
   * The destructor.
   */
  virtual ~PowhegFactory();
  //@}

public:

  /**
   * Return the Born/virtual pieces to consider
   */
  const vector<Ptr<MatchboxNLOME>::ptr>& bornVirtuals() const {
    return theBornVirtuals;
  }

  /**
   * Access the Born/virtual pieces to consider
   */
  vector<Ptr<MatchboxNLOME>::ptr>& bornVirtuals() {
    return theBornVirtuals;
  }

  /**
   * Return the subtracted matrix elements to consider
   */
  const vector<Ptr<SubtractedME>::ptr>& subtractedMEs() const {
    return theSubtractedMEs;
  }

  /**
   * Access the subtracted matrix elements to consider
   */
  vector<Ptr<SubtractedME>::ptr>& subtractedMEs() {
    return theSubtractedMEs;
  }

  /**
   * Return true, if 'Born screening' should be done
   */
  bool bornScreening() const { return theBornScreening; }

  /**
   * Switch on 'Born screening'
   */
  void doBornScreening() { theBornScreening = true; }

  /**
   * Switch off 'Born screening'
   */
  void noBornScreening() { theBornScreening = false; }

  /**
   * Return the inclusive (BBar) ME's generated
   */
  const vector<Ptr<PowhegInclusiveME>::ptr>& inclusiveMEs() const {
    return theInclusiveMEs;
  }

  /**
   * Return the finite real emission contributions
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& realMEs() const {
    return theRealMEs;
  }

public:

  /**
   * Return true if this object needs to be initialized before all
   * other objects (except those for which this function also returns
   * true).  This default version always returns false, but subclasses
   * may override it to return true.
   */
  virtual bool preInitialize() const { return true; }

  /**
   * Setup everything.
   */
  void setup();

  /**
   * Dump the setup
   */
  void print(ostream&) const;

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
  //@}


private:

  /**
   * The Born/virtual pieces to consider
   */
  vector<Ptr<MatchboxNLOME>::ptr> theBornVirtuals;

  /**
   * The subtracted matrix elements to consider
   */
  vector<Ptr<SubtractedME>::ptr> theSubtractedMEs;

  /**
   * An optional NLO factory to pick matrix elements from
   */
  Ptr<MatchboxFactory>::ptr theMatchboxFactory;

  /**
   * Wether or not 'Born screening' should be performed
   */
  bool theBornScreening;

  /**
   * The inclusive (BBar) ME's generated
   */
  vector<Ptr<PowhegInclusiveME>::ptr> theInclusiveMEs;

  /**
   * The finite real emission contributions
   */
  vector<Ptr<MatchboxMEBase>::ptr> theRealMEs;

  /**
   * Switch on or off verbosity
   */
  bool theVerbose;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegFactory & operator=(const PowhegFactory &);

};

}

#endif /* HERWIG_PowhegFactory_H */
