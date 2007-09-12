// -*- C++ -*-
#ifndef HERWIG_KtTildeMeasure_H
#define HERWIG_KtTildeMeasure_H
//
// This is the declaration of the KtTildeMeasure class.
//

#include "Herwig++/Shower/CKKW/Reweighting/DefaultJetMeasure.h"
#include "KtTildeMeasure.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 * This class implements the \tilde{k}_\perp
 * jet measure, which generalizes the Durham measure
 * to the Herwig++ evolution variable.
 *
 * @see \ref KtTildeMeasureInterfaces "The interfaces"
 * defined for KtTildeMeasure.
 */
class KtTildeMeasure: public DefaultJetMeasure {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline KtTildeMeasure();

  /**
   * The destructor.
   */
  virtual ~KtTildeMeasure();
  //@}

public:

  /**
   * Return true, if all particles are
   * resolvable. This is used for cuts on matrix elements.
   * The second member of the pair indicates final state (false)
   * or initial state (true).
   */
  virtual bool resolvable (const vector<pair<PPtr,bool> >&);

  /**
   * Return the clustering scale for the given configuration.
   */
  virtual Energy2 scale (const vector<tClusteringParticlePtr>&,
			 const vector<ClusteringParticlePtr>&,
			 const tClusteringConfigurationPtr&);

  /**
   * Return the momentum fraction for the given clustering.
   */
  virtual double z (const vector<tClusteringParticlePtr>&,
		    const vector<ClusteringParticlePtr>&,
		    const tClusteringConfigurationPtr&);

  /**
   * Return the minimum resolvable clustering scale,
   * This is set as the splitting scale of external legs.
   */
  virtual Energy2 minResolvableScale (long PDGId, bool initial = false);

  /**
   * Set the production scales for the hard process
   * particles.
   */
  virtual bool hardScales (const vector <tClusteringParticlePtr>&);

  /**
   * Return true, if the jet resolution can
   * handle the given configuration.
   */
  virtual bool canHandle (const IdList&, bool initial = false);

  /**
   * Return true, if the branching configuration
   * together with the scale and momentum fraction
   * is resolvable. This is used for integrating
   * Sudakov exponents.
   */
  virtual bool resolvable (const IdList&, Energy2, Energy2, double, bool initial = false);

  /**
   * Return the Jacobian for going from the parton shower
   * variables to the clustering variables as a function
   * of the clustering variables.
   */
  virtual double showerJacobian (const IdList&, Energy2, double, bool initial = false);

  /**
   * Return the shower scale and momentum fraction
   * as a function of the clustering variables.
   */
  virtual pair<Energy2,double> invertClustering (const IdList&, Energy2, double, bool initial = false);


  /**
   * Return true, if the given branching is resolvable.
   * This is used for vetoing shower emissions.
   */
  virtual bool resolvable (tcShowerParticlePtr, const Branching&, bool initial = false);

  /**
   * Return the bounds on the clustering momentum
   * fractions as a function of scale and the hard
   * scale.
   */
  virtual pair<double,double> zLimits (Energy2 q, Energy2 Q, const IdList&, bool initial = false);

protected:

  /**
   * Calculate z+/-
   */
  double zPM (Energy2 mi2, Energy2 mj2, Energy2 mij2);

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
   * Get the emitter's mass for QCD branchings
   */
  Energy2 emitterMass (long, long, bool&);

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<KtTildeMeasure> initKtTildeMeasure;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KtTildeMeasure & operator=(const KtTildeMeasure &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of KtTildeMeasure. */
template <>
struct BaseClassTrait<Herwig::KtTildeMeasure,1> {
  /** Typedef of the first base class of KtTildeMeasure. */
  typedef Herwig::DefaultJetMeasure NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KtTildeMeasure class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::KtTildeMeasure>
  : public ClassTraitsBase<Herwig::KtTildeMeasure> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::KtTildeMeasure"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KtTildeMeasure is implemented. It may also include several, space-separated,
   * libraries if the class KtTildeMeasure depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "KtTildeMeasure.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KtTildeMeasure.tcc"
#endif

#endif /* HERWIG_KtTildeMeasure_H */
