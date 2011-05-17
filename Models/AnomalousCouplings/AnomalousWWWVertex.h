// -*- C++ -*-
#ifndef HERWIG_AnomalousWWWVertex_H
#define HERWIG_AnomalousWWWVertex_H
//
// This is the declaration of the AnomalousWWWVertex class.
//

#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The AnomalousWWWVertex class implements the interaction
 * of three electroweak gauge bosons including
 * anomalous terms.
 *
 * @see \ref AnomalousWWWVertexInterfaces "The interfaces"
 * defined for AnomalousWWWVertex.
 */
class AnomalousWWWVertex: public AbstractVVVVertex {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  AnomalousWWWVertex();
  //@}
  
public:
  
  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   */
  Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2, const VectorWaveFunction & vec3);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const VectorWaveFunction & vec3,
			      complex<Energy> mass=-GeV,
			      complex<Energy> width=-GeV);
  //@}

protected:
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3,
			   double & g, double & kappa,
			   unsigned int & order);
  
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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<AnomalousWWWVertex> initAnomalousWWWVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AnomalousWWWVertex & operator=(const AnomalousWWWVertex &);

private:

  /**
   *  Couplings
   */
  //@{
  /**
   *  \f$g^Z\f$
   */
  double gZ_;

  /**
   *  \f$g^\gamma\f$
   */
  double gGamma_;

  /**
   *  \f$\kappa^Z\f$
   */
  double kappaZ_;

  /**
   *  \f$\kappa^\gamma\f$
   */
  double kappaGamma_;

  /**
   *  \f$\lambda^\gamma\f$
   */
  double lambda_;
  //@}

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The factor for the \f$Z\f$ vertex.
   */
  double _zfact;

  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AnomalousWWWVertex. */
template <>
struct BaseClassTrait<Herwig::AnomalousWWWVertex,1> {
  /** Typedef of the first base class of AnomalousWWWVertex. */
  typedef Helicity::AbstractVVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AnomalousWWWVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::AnomalousWWWVertex>
  : public ClassTraitsBase<Herwig::AnomalousWWWVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::AnomalousWWWVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * AnomalousWWWVertex is implemented. It may also include several, space-separated,
   * libraries if the class AnomalousWWWVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwAnomalousCouplings.so"; }
};

/** @endcond */

}

#endif /* HERWIG_AnomalousWWWVertex_H */
