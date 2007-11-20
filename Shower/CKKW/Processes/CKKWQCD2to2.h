// -*- C++ -*-
//
// CKKWQCD2to2.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_CKKWQCD2to2_H
#define HERWIG_CKKWQCD2to2_H
//
// This is the declaration of the CKKWQCD2to2 class.
//

#include "Herwig++/Shower/CKKW/Clustering/CKKWHardProcess.h"
#include "CKKWQCD2to2.fh"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * QCD 2 to 2 scatterings as hard proceswses in ME/PS merging.
 *
 *@author Simon Plaetzer
 *
 * @see \ref CKKWQCD2to2Interfaces "The interfaces"
 * defined for CKKWQCD2to2.
 */
class CKKWQCD2to2: public CKKWHardProcess {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline CKKWQCD2to2();

  /**
   * The destructor.
   */
  virtual ~CKKWQCD2to2();
  //@}

public:

  /**
   * Return true, if the configuration given corresponds
   * to a valid signal process.
   */
  virtual bool reachedHard (const vector<ClusteringParticleData>&) const;

  /**
   * Return true, if the given configuration is not to be
   * used as the hard process.
   */
  virtual inline bool veto (const vector <tClusteringParticlePtr>&) const;

  /**
   * Wether or not to use colour information
   */
  inline bool useColour () const;

  /**
   * Return the running coupling weight for the hard process.
   */
  virtual double hardCouplings (const vector <tClusteringParticlePtr>&,double MEAlpha);

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * Return true, if colour conservation is fullfilled.
   */
  bool colourConservation (const vector<ClusteringParticleData>&) const;

  /**
   * Wether or not to use colour information
   */
  bool _useColour;

  /**
   * The running alpha to be used.
   */
  ShowerAlphaPtr _alpha;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CKKWQCD2to2> initCKKWQCD2to2;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CKKWQCD2to2 & operator=(const CKKWQCD2to2 &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CKKWQCD2to2. */
template <>
struct BaseClassTrait<Herwig::CKKWQCD2to2,1> {
  /** Typedef of the first base class of CKKWQCD2to2. */
  typedef Herwig::CKKWHardProcess NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CKKWQCD2to2 class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CKKWQCD2to2>
  : public ClassTraitsBase<Herwig::CKKWQCD2to2> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CKKWQCD2to2"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CKKWQCD2to2 is implemented. It may also include several, space-separated,
   * libraries if the class CKKWQCD2to2 depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSHower.so"; }
};

/** @endcond */

}

#include "CKKWQCD2to2.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWQCD2to2.tcc"
#endif

#endif /* HERWIG_CKKWQCD2to2_H */
