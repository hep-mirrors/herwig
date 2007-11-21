// -*- C++ -*-
//
// Clusterer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Clusterer_H
#define HERWIG_Clusterer_H
//
// This is the declaration of the Clusterer class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Clusterer.fh"

#include "ClusteringConfiguration.h"
#include "ClusteringParticle.h"
#include "Clustering.h"
#include "PostClustering.h"
#include "Herwig++/Shower/CKKW/Reweighting/JetMeasure.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup CKKW
 *
 * Clusterer is the base class for all objects which
 * perform clusterings during a parton shower history
 * reconstruction.
 *
 *@author Simon Plaetzer
 *
 * @see \ref ClustererInterfaces "The interfaces"
 * defined for Clusterer.
 */
class Clusterer: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~Clusterer();
  //@}

public:
  /**@name Methods used to determine the ClusteringGuide */
  //@{

  /**
   * Return the number of particles this clusterer clusters.
   */
  virtual unsigned int toClusterMultiplicity () const = 0;

  /**
   * Return the number of particles emerging from a clustering
   * being performed by this clusterer.
   */
  virtual unsigned int fromClusterMultiplicity () const = 0;

  /**
   * Return the possible clustering configurations resulting
   * from the given configuration. Return empty vector,
   * if configuration cannot be handled by this clusterer.
   */
  virtual vector<ClusteringConfigurationPtr> configurations (const vector<ClusteringParticleData>&) = 0;

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
				 const tClusteringConfigurationPtr&) = 0;

  /**
   * The clustering was chosen to be performed:
   * Do all the kinematics and complete the clustering
   */
  virtual void doKinematics (const tClusteringPtr&) = 0;

  //@}

  /**
   * Return the jet measure used
   */
  inline tJetMeasurePtr jetMeasure () const;

  /**
   * Set the jet measure used
   */
  inline void jetMeasure (JetMeasurePtr);

  /**
   * Reset this clusterer
   */
  virtual void reset () = 0;

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



protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}


private:

  /**
   * The jet measure used
   */ 
  JetMeasurePtr _jetMeasure;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<Clusterer> initClusterer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Clusterer & operator=(const Clusterer &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Clusterer. */
template <>
struct BaseClassTrait<Herwig::Clusterer,1> {
  /** Typedef of the first base class of Clusterer. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Clusterer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Clusterer>
  : public ClassTraitsBase<Herwig::Clusterer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Clusterer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Clusterer is implemented. It may also include several, space-separated,
   * libraries if the class Clusterer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "Clusterer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Clusterer.tcc"
#endif

#endif /* HERWIG_Clusterer_H */
