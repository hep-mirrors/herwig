// -*- C++ -*-
//
// PhasespaceCouplings.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_PhasespaceCouplings_H
#define Herwig_PhasespaceCouplings_H
//
// This is the declaration of the PhasespaceCouplings class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace Herwig {

using namespace ThePEG;

typedef boost::tuple<long,long,long> LTriple;

inline PersistentOStream& operator<<(PersistentOStream& os, const LTriple& t) {
  os << t.get<0>() << t.get<1>() << t.get<2>();
  return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is, LTriple& t) {
  is >> t.get<0>() >> t.get<1>() >> t.get<2>();
  return is;
}

/**
 *
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Store couplings for the phasespace generator.
 *
 * @see \ref PhasespaceCouplingsInterfaces "The interfaces"
 * defined for PhasespaceCouplings.
 */
class PhasespaceCouplings: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PhasespaceCouplings();

  /**
   * The destructor.
   */
  virtual ~PhasespaceCouplings();
  //@}

public:

  /**
   * Couplings to be used in diagram weighting
   */
  map<LTriple,double>& couplings() { return theCouplings; }

  /**
   * Couplings to be used in diagram weighting
   */
  const map<LTriple,double>& couplings() const { return theCouplings; }

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
   * Couplings to be used in diagram weighting
   */
  map<LTriple,double> theCouplings;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PhasespaceCouplings & operator=(const PhasespaceCouplings &);

};

}

#endif /* Herwig_PhasespaceCouplings_H */
