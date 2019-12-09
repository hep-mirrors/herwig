// -*- C++ -*-
//
// MatchboxReweightBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxReweightBase_H
#define HERWIG_MatchboxReweightBase_H
//
// This is the declaration of the MatchboxReweightBase class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxReweightBase is the base class
 * for reweighting MatchboxMEBase matrix elements
 * as |M|^2 ( w_1 + ... + w_n )
 *
 */
class MatchboxReweightBase: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxReweightBase();

  /**
   * The destructor.
   */
  virtual ~MatchboxReweightBase();
  //@}

public:

  /**
   * Clone this reweight.
   */
  Ptr<MatchboxReweightBase>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxReweightBase>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
  virtual void cloneDependencies(const std::string& prefix = "");

  /**
   * Set the XComb object.
   */
  virtual void setXComb(tStdXCombPtr) = 0;

  /**
   * Return true, if applies to the process in the xcomb.
   */
  virtual bool apply() const = 0;

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches() = 0;

  /**
   * Evaluate the reweight.
   */
  virtual double evaluate() const = 0;

  /**
   * Set veto scales on the particles at the given
   * SubProcess which has been generated using this
   * matrix element.
   */
  virtual void setVetoScales(tSubProPtr) const {}

public:

  /**
   * Dump the setup to an ostream
   */
  virtual void print(ostream&) const {}

  /**
   * Print debug information on the last event
   */
  virtual void printLastEvent(ostream&) const {}

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxReweightBase & operator=(const MatchboxReweightBase &) = delete;

};

}


#endif /* HERWIG_MatchboxReweightBase_H */
