// -*- C++ -*-
#ifndef HERWIG_CKKWHardProcess_H
#define HERWIG_CKKWHardProcess_H
//
// This is the declaration of the CKKWHardProcess class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "CKKWHardProcess.fh"

#include "ClusteringParticle.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * Defines an interface for any hard process of interest
 * to be reached at the end of a cascade reconstruction.
 *
 *@author Simon Plaetzer
 *
 * @see \ref CKKWHardProcessInterfaces "The interfaces"
 * defined for CKKWHardProcess.
 */
class CKKWHardProcess: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~CKKWHardProcess();
  //@}

public:

  /**
   * Return true, if the configuration given corresponds
   * to a valid signal process.
   */
  virtual bool reachedHard (const vector<ClusteringParticleData>&) const = 0;

  /**
   * Return true, if the given configuration is not to be
   * used as the hard process.
   */
  virtual bool veto (const vector <tClusteringParticlePtr>&) const = 0;

  /**
   * Return the running coupling weight for the hard process.
   * The default returns one.
   */
  virtual inline double hardCouplings (const vector <tClusteringParticlePtr>&,double MEAlpha);

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
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<CKKWHardProcess> initCKKWHardProcess;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CKKWHardProcess & operator=(const CKKWHardProcess &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CKKWHardProcess. */
template <>
struct BaseClassTrait<Herwig::CKKWHardProcess,1> {
  /** Typedef of the first base class of CKKWHardProcess. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CKKWHardProcess class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CKKWHardProcess>
  : public ClassTraitsBase<Herwig::CKKWHardProcess> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CKKWHardProcess"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CKKWHardProcess is implemented. It may also include several, space-separated,
   * libraries if the class CKKWHardProcess depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "CKKWHardProcess.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CKKWHardProcess.tcc"
#endif

#endif /* HERWIG_CKKWHardProcess_H */
