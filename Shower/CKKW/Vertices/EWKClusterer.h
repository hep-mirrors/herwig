// -*- C++ -*-
//
// EWKClusterer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EWKClusterer_H
#define HERWIG_EWKClusterer_H
//
// This is the declaration of the EWKClusterer class.
//

#include "Herwig++/Shower/CKKW/Clustering/Clusterer.h"
#include "EWKClusterer.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * Electroweak vertices in ME/PS merging
 *
 *@author Simon Plaetzer
 *
 * @see \ref EWKClustererInterfaces "The interfaces"
 * defined for EWKClusterer.
 */
class EWKClusterer: public Clusterer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline EWKClusterer();

  /**
   * The destructor.
   */
  virtual ~EWKClusterer();
  //@}

public:

  /**
   * Reset this clusterer
   */
  virtual inline void reset ();

  /**@name Methods used to determine the ClusteringGuide */
  //@{

  /**
   * Return the number of particles this clusterer clusters.
   */
  virtual inline unsigned int toClusterMultiplicity () const;

  /**
   * Return the number of particles emerging from a clustering
   * being performed by this clusterer.
   */
  virtual inline unsigned int fromClusterMultiplicity () const;

  /**
   * Return the possible clustering configurations resulting
   * from the given configuration. Return empty vector,
   * if configuration cannot be handled by this clusterer.
   */
  virtual vector<ClusteringConfigurationPtr> configurations (const vector<ClusteringParticleData>&);

  //@}

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

protected:

  /**
   * Return true, if colour connections are
   * to be taken into account.
   */
  inline bool useColour () const;

  /**
   * Return true, if hadronic clusterings should be used
   */
  inline bool useHadronic () const;

  /**
   * Return the internal line emerging from the
   * given particles according to a electroweak vertex.
   * The CKM matrix is assumed diagonal.
   */
  inline ClusteringParticleData emergingLine (ClusteringParticleData, ClusteringParticleData, bool& isEWK);

  /**
   * Return the electroweak order to be considered
   */
  inline unsigned int EWKOrder () const;

  /**
   * Indicate that a clustering is to be performed.
   * Return true, if this is allowed according to
   * the electroweak order considered.
   */
  inline bool performEWK ();

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
   * Return the internal line emerging from the
   * given particles according to a EWK vertex.
   */
  ClusteringParticleData doEmergingLine (ClusteringParticleData, ClusteringParticleData, bool& isEWK) const;

  /**
   * True, if colour information should be considered
   * when clustering particles.
   */
  bool _useColour;

  /**
   * Wether or not to use hadronic clusterings
   */
  bool _useHadronic;

  /**
   * The electroweak order to be considered
   */
  unsigned int _EWKOrder;

  /**
   * The current electroweak order
   */
  unsigned int _currentEWKOrder;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static ClassDescription<EWKClusterer> initEWKClusterer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EWKClusterer & operator=(const EWKClusterer &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EWKClusterer. */
template <>
struct BaseClassTrait<Herwig::EWKClusterer,1> {
  /** Typedef of the first base class of EWKClusterer. */
  typedef Herwig::Clusterer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EWKClusterer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EWKClusterer>
  : public ClassTraitsBase<Herwig::EWKClusterer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EWKClusterer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EWKClusterer is implemented. It may also include several, space-separated,
   * libraries if the class EWKClusterer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "EWKClusterer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EWKClusterer.tcc"
#endif

#endif /* HERWIG_EWKClusterer_H */
