// -*- C++ -*-
//
// CKKWVeto.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licence under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_CKKWVeto_H
#define HERWIG_CKKWVeto_H
//
// This is the declaration of the CKKWVeto class.
//

#include "ShowerVeto.h"
#include "Herwig++/Shower/ShowerHandler.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 * 
 * Veto shower emissions according to transverse
 * momentum of a branching.
 *
 * Emissions may be vetoed, if the pt is above and/or
 * below given values.
 *
 * @see \ref CKKWVetoInterfaces "The interfaces"
 * defined for CKKWVeto.
 */
class CKKWVeto: public ShowerVeto {

  /**
   *  The ShowerHandler is a friend to set some parameters at initialisation
   */
  friend class ShowerHandler;

public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CKKWVeto() : ShowerVeto(ShowerVeto::Emission), pTVetoDefinition_(1),
	       reversepTVeto_(false), pTVeto_(ZERO), 
	       highestMult_(false), dynamicSuds_(false) {}
  //@}

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

public:

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoTimeLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			     const Branching&);

  /**
   * Return true, if the selected emission off the given
   * particle and progenitor is vetoed.
   */
  virtual bool vetoSpaceLike (tcShowerProgenitorPtr, tcShowerParticlePtr,
			      const Branching&);

  /**
   *  Set whether showering highest multiplicity configuration
   */ 
  void setHighest( bool isHighest ){
    highestMult_ = isHighest;
  }

  /**
   * Access to the veto scale for CKKW merging
   */
  Energy getVeto() {
    return pTVeto_;
  }

  /**
   * Set dynamic sudakovs on or off
   */
  void setDynamicSuds( bool dynamicSuds ) {
    dynamicSuds_ = dynamicSuds;
  }

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CKKWVeto> initCKKWVeto;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CKKWVeto & operator=(const CKKWVeto &);

private:

  /**
   * The transverse momentum definition to be used
   */
  unsigned int pTVetoDefinition_;
  
  /**
   * Whether to reverse the veto
   */
  bool reversepTVeto_; 

  /**
   * The scale at which the veto should be applied
   */
  Energy pTVeto_;

  /**
   * Whether we are showering highest multiplicity channe;
   */
  bool highestMult_;
  
  /**
   * Whether to generate the Sudakov weights dynamically by event vetoes
   */
  bool dynamicSuds_;

  /**
   * Apply the veto to timelike showering
   */
  bool vetoTimeLike_;

  /**
   * Apply the veto to spacelike showering
   */
  bool vetoSpaceLike_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CKKWVeto. */
template <>
struct BaseClassTrait<Herwig::CKKWVeto,1> {
  /** Typedef of the first base class of CKKWVeto. */
  typedef Herwig::ShowerVeto NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CKKWVeto class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CKKWVeto>
  : public ClassTraitsBase<Herwig::CKKWVeto> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CKKWVeto"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CKKWVeto is implemented. It may also include several, space-separated,
   * libraries if the class CKKWVeto depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_CKKWVeto_H */
