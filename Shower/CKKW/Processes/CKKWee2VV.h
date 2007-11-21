// -*- C++ -*-
//
// CKKWee2VV.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_CKKWee2VV_H
#define HERWIG_CKKWee2VV_H
//
// This is the declaration of the CKKWee2VV class.
//

#include "Herwig++/Shower/CKKW/Clustering/CKKWHardProcess.h"
#include "CKKWee2VV.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * The e+e- to V V hard process.
 *
 *@author Simon Plaetzer
 *
 * @see \ref CKKWee2VVInterfaces "The interfaces"
 * defined for CKKWee2VV.
 */
class CKKWee2VV: public CKKWHardProcess {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~CKKWee2VV();
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CKKWee2VV> initCKKWee2VV;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CKKWee2VV & operator=(const CKKWee2VV &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CKKWee2VV. */
template <>
struct BaseClassTrait<Herwig::CKKWee2VV,1> {
  /** Typedef of the first base class of CKKWee2VV. */
  typedef Herwig::CKKWHardProcess NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CKKWee2VV class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CKKWee2VV>
  : public ClassTraitsBase<Herwig::CKKWee2VV> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CKKWee2VV"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CKKWee2VV is implemented. It may also include several, space-separated,
   * libraries if the class CKKWee2VV depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "CKKWee2VV.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWee2VV.tcc"
#endif

#endif /* HERWIG_CKKWee2VV_H */
