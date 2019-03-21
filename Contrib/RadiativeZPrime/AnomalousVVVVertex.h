// -*- C++ -*-
#ifndef RADIATIVEZPRIME_AnomalousVVVVertex_H
#define RADIATIVEZPRIME_AnomalousVVVVertex_H
//
// This is the declaration of the AnomalousVVVVertex class.
//

#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"

namespace RadiativeZPrime {

using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The AnomalousVVVVertex class implements the anomalous Vector-Vector-Vector vertex.
 *
 *  Only the member which evaluates the matrix element is implemented. The vertex
 *  is defined to be 
 *  \[ i\epsilon_{\mu\nu\alpha\beta}
 *      \varepsilon^\mu_1 \varepsilon^\nu_2 \varepsilon^\alpha_2 p_1^\beta\]
 *
 *
 * @see \ref AnomalousVVVVertexInterfaces "The interfaces"
 * defined for AnomalousVVVVertex.
 */
class AnomalousVVVVertex: public AbstractVVVVertex {

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
  virtual Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const VectorWaveFunction & vec3);

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
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec2,
				      const VectorWaveFunction & vec3,
				      complex<Energy> mass=-GeV,
				      complex<Energy> width=-GeV);
  //@}

  /**
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3)=0;

  /**
   * Dummy setCouplings for a four point interaction 
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   */
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr) {
    assert(false);
  }

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<AnomalousVVVVertex> initAnomalousVVVVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AnomalousVVVVertex & operator=(const AnomalousVVVVertex &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AnomalousVVVVertex. */
template <>
struct BaseClassTrait<RadiativeZPrime::AnomalousVVVVertex,1> {
  /** Typedef of the first base class of AnomalousVVVVertex. */
  typedef Helicity::AbstractVVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AnomalousVVVVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<RadiativeZPrime::AnomalousVVVVertex>
  : public ClassTraitsBase<RadiativeZPrime::AnomalousVVVVertex> {
  /** Return a platform-independent class name */
  static string className() { return "RadiativeZPrime::AnomalousVVVVertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * AnomalousVVVVertex is implemented. It may also include several, space-separated,
   * libraries if the class AnomalousVVVVertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "RadiativeZPrime.so"; }
};

/** @endcond */

}

#endif /* RADIATIVEZPRIME_AnomalousVVVVertex_H */
