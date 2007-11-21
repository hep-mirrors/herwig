// -*- C++ -*-
//
// EWKQtildeClusterer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EWKQtildeClusterer_H
#define HERWIG_EWKQtildeClusterer_H
//
// This is the declaration of the EWKQtildeClusterer class.
//

#include "Herwig++/Shower/CKKW/Vertices/EWKEmitterSpectatorClusterer.h"
#include "EWKQtildeClusterer.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * EWK clusterings in the qtilde formalism
 *
 *@author Simon Plaetzer
 *
 * @see \ref EWKQtildeClustererInterfaces "The interfaces"
 * defined for EWKQtildeClusterer.
 */
class EWKQtildeClusterer: public EWKEmitterSpectatorClusterer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~EWKQtildeClusterer();
  //@}

public:

  /**@name Methods used for performing an actual clustering. */
  //@{
  /**
   * Set the scale, weight and possible other parameters
   * associated with the clustering resulting from
   * given children, parents and configuration. Return the
   * clustering, if done. Return null, if couldn't be handled.
   */
  virtual ClusteringPtr doScale (const vector<tClusteringParticlePtr>&,
				 const vector<ClusteringParticlePtr>&,
				 const tClusteringConfigurationPtr&);

  /**
   * The clustering was chosen to be performed:
   * Do all the kinematics and complete the clustering
   */
  virtual void doKinematics (const tClusteringPtr&);

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
  static ClassDescription<EWKQtildeClusterer> initEWKQtildeClusterer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EWKQtildeClusterer & operator=(const EWKQtildeClusterer &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EWKQtildeClusterer. */
template <>
struct BaseClassTrait<Herwig::EWKQtildeClusterer,1> {
  /** Typedef of the first base class of EWKQtildeClusterer. */
  typedef Herwig::EWKEmitterSpectatorClusterer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EWKQtildeClusterer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EWKQtildeClusterer>
  : public ClassTraitsBase<Herwig::EWKQtildeClusterer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EWKQtildeClusterer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EWKQtildeClusterer is implemented. It may also include several, space-separated,
   * libraries if the class EWKQtildeClusterer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "EWKQtildeClusterer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EWKQtildeClusterer.tcc"
#endif

#endif /* HERWIG_EWKQtildeClusterer_H */
