// -*- C++ -*-
//
// JetMeasure.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_JetMeasure_H
#define HERWIG_JetMeasure_H
//
// This is the declaration of the JetMeasure class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "JetMeasure.fh"

#include "Herwig++/Shower/CKKW/Clustering/ClusteringParticle.h"
#include "Herwig++/Shower/CKKW/Clustering/ClusteringConfiguration.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * JetMeasure is the base class for doing cuts on
 * matrix elements, integrating splitting probabilities
 * for reweighting and vetoing shower emissions
 * in an CKKW-like ME/PS merging approach.
 *
 *@author Simon Plaetzer
 *
 * @see \ref JetMeasureInterfaces "The interfaces"
 * defined for JetMeasure.
 */
class JetMeasure: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline JetMeasure();

  /**
   * The destructor.
   */
  virtual ~JetMeasure();
  //@}

public:

  /**
   * Return true, if all particles are
   * resolvable. This is used for cuts on matrix elements.
   * The second member of the pair indicates final state (false)
   * or initial state (true).
   */
  virtual bool resolvable (const vector<pair<PPtr,bool> >&) = 0;

  /**
   * Return the clustering scale for the given configuration.
   */
  virtual Energy2 scale (const vector<tClusteringParticlePtr>&,
			 const vector<ClusteringParticlePtr>&,
			 const tClusteringConfigurationPtr&) = 0;

  /**
   * Return the momentum fraction for the given clustering.
   */
  virtual double z (const vector<tClusteringParticlePtr>&,
		    const vector<ClusteringParticlePtr>&,
		    const tClusteringConfigurationPtr&) = 0;

  /**
   * Return the minimum resolvable clustering scale,
   * This is set as the splitting scale of external legs.
   */
  virtual Energy2 minResolvableScale (long PDGId, bool initial = false) = 0;

  /**
   * Get the resolution parameter
   */
  inline Energy2 resolution () const;

  /**
   * Return the minimum number of partons
   * to be considered for evaluating cuts.
   */
  inline unsigned int minPartons () const;

  /**
   * Return the maximum number of partons
   * to be considered for evaluating cuts.
   */
  inline unsigned int maxPartons () const;

  /**
   * Set the production scales for the hard process
   * particles.
   */
  virtual bool hardScales (const vector <tClusteringParticlePtr>&) = 0;

protected:

  /**
   * Construct giving the maximum and minimum numbers
   * of partons to evaluate cuts for.
   */
  inline JetMeasure (unsigned int min, unsigned int max);

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
   * The resolution scale.
   */
  Energy2 _resolution;

  /**
   * The minimum number of partons to be
   * considered for cuts.
   */
  unsigned int _min;

  /**
   * The maximum number of partons to be
   * considered for cuts.
   */
  unsigned int _max;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<JetMeasure> initJetMeasure;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  JetMeasure & operator=(const JetMeasure &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of JetMeasure. */
template <>
struct BaseClassTrait<Herwig::JetMeasure,1> {
  /** Typedef of the first base class of JetMeasure. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the JetMeasure class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::JetMeasure>
  : public ClassTraitsBase<Herwig::JetMeasure> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::JetMeasure"; }
  /**
   * The name of a file containing the dynamic library where the class
   * JetMeasure is implemented. It may also include several, space-separated,
   * libraries if the class JetMeasure depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "JetMeasure.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "JetMeasure.tcc"
#endif

#endif /* HERWIG_JetMeasure_H */
