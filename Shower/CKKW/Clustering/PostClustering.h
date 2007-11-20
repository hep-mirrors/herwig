// -*- C++ -*-
//
// PostClustering.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PostClustering_H
#define HERWIG_PostClustering_H
//
// This is the declaration of the PostClustering class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "PostClustering.fh"

#include "ClusteringParticle.fh"
#include "Clustering.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * The interface for any transformation to be applied after
 * a clustering has been performed (such as boosting to the
 * CMS of the new incoming particles).
 *
 *@author Simon Plaetzer
 *
 * @see \ref PostClusteringInterfaces "The interfaces"
 * defined for PostClustering.
 */
class PostClustering: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~PostClustering();
  //@}

public:

  /**
   * Initialize the PostClustering object given the
   * particles after clustering has been performed.
   */
  virtual void initialize (const vector<tClusteringParticlePtr>&) = 0;

  /**
   * Apply the transformation
   */
  virtual void doTransform (tClusteringParticlePtr) = 0;

public:

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
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<PostClustering> initPostClustering;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PostClustering & operator=(const PostClustering &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PostClustering. */
template <>
struct BaseClassTrait<Herwig::PostClustering,1> {
  /** Typedef of the first base class of PostClustering. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PostClustering class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PostClustering>
  : public ClassTraitsBase<Herwig::PostClustering> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PostClustering"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PostClustering is implemented. It may also include several, space-separated,
   * libraries if the class PostClustering depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

//#include "PostClustering.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PostClustering.tcc"
#endif

#endif /* HERWIG_PostClustering_H */
