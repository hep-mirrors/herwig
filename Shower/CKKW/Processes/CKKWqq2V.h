// -*- C++ -*-
#ifndef HERWIG_CKKWqq2V_H
#define HERWIG_CKKWqq2V_H
//
// This is the declaration of the CKKWqq2V class.
//

#include "Herwig++/Shower/CKKW/Clustering/CKKWHardProcess.h"
#include "CKKWqq2V.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * The q q to V hard process.
 *
 *@author Simon Plaetzer
 *
 * @see \ref CKKWqq2VInterfaces "The interfaces"
 * defined for CKKWqq2V.
 */
class CKKWqq2V: public CKKWHardProcess {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline CKKWqq2V();

  /**
   * The destructor.
   */
  virtual ~CKKWqq2V();
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
   * Wether or not to use colour information
   */
  bool _useColour;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CKKWqq2V> initCKKWqq2V;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CKKWqq2V & operator=(const CKKWqq2V &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CKKWqq2V. */
template <>
struct BaseClassTrait<Herwig::CKKWqq2V,1> {
  /** Typedef of the first base class of CKKWqq2V. */
  typedef Herwig::CKKWHardProcess NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CKKWqq2V class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CKKWqq2V>
  : public ClassTraitsBase<Herwig::CKKWqq2V> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CKKWqq2V"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CKKWqq2V is implemented. It may also include several, space-separated,
   * libraries if the class CKKWqq2V depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "CKKWqq2V.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWqq2V.tcc"
#endif

#endif /* HERWIG_CKKWqq2V_H */
