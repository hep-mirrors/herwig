// -*- C++ -*-
#ifndef HERWIG_SSSSVertex_H
#define HERWIG_SSSSVertex_H
//
// This is the declaration of the SSSSVertex class.
//
#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "SSSSVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 * 
 *  The SSSSVertex class is the implementation of the interaction of
 *  four scalars. It inherits from the VertexBase class for the storage 
 *  of the particles interacting at the vertex and implements the 
 *  helicity calculations.
 *
 *  Any classes implementating the vertex should inherit from it and implement
 *  the virtual set Coupling member.
 *
 *  The form of the vertex is
 * \f[ic\phi_1\phi_2\phi_3\phi_4\f]
 *
 *  @see VertexBase
 */
class SSSSVertex: public VertexBase {
  
public:
  
  /**
   * Default constructor.
   */
  inline SSSSVertex();
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{  
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param sca3 The wavefunction for the third  scalar.
   * @param sca4 The wavefunction for the fourth scalar.
   */
  Complex evaluate(Energy2 q2, const ScalarWaveFunction & sca1,
		   const ScalarWaveFunction & sca2, const ScalarWaveFunction & sca3,
		   const ScalarWaveFunction & sca4);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param sca3 The wavefunction for the third  scalar.
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const ScalarWaveFunction & sca1,
			      const ScalarWaveFunction & sca2,
			      const ScalarWaveFunction & sca3);
  //@}

  /**
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3,tcPDPtr part4)=0;

private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<SSSSVertex> initSSSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SSSSVertex & operator=(const SSSSVertex &);
  
};
}
}

#include "SSSSVertex.icc"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SSSSVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::SSSSVertex,1> {
  /** Typedef of the base class of SSSSVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::SSSSVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSSSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::SSSSVertex"; }
};

/** @endcond */

}


#endif /* HERWIG_SSSSVertex_H */
