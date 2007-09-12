// -*- C++ -*-
#ifndef HERWIG_CKKWee2qq_H
#define HERWIG_CKKWee2qq_H
//
// This is the declaration of the CKKWee2qq class.
//

#include "Herwig++/Shower/CKKW/Clustering/CKKWHardProcess.h"
#include "CKKWee2qq.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * The e+e- to q qbar hard process.
 *
 *@author Simon Plaetzer
 *
 * @see \ref CKKWee2qqInterfaces "The interfaces"
 * defined for CKKWee2qq.
 */
class CKKWee2qq: public CKKWHardProcess {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline CKKWee2qq();

  /**
   * The destructor.
   */
  virtual ~CKKWee2qq();
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
  static ClassDescription<CKKWee2qq> initCKKWee2qq;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CKKWee2qq & operator=(const CKKWee2qq &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CKKWee2qq. */
template <>
struct BaseClassTrait<Herwig::CKKWee2qq,1> {
  /** Typedef of the first base class of CKKWee2qq. */
  typedef Herwig::CKKWHardProcess NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CKKWee2qq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CKKWee2qq>
  : public ClassTraitsBase<Herwig::CKKWee2qq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CKKWee2qq"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CKKWee2qq is implemented. It may also include several, space-separated,
   * libraries if the class CKKWee2qq depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "CKKWee2qq.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWee2qq.tcc"
#endif

#endif /* HERWIG_CKKWee2qq_H */
