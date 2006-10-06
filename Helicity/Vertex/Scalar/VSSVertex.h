// -*- C++ -*-
#ifndef HERWIG_VSSVertex_H
#define HERWIG_VSSVertex_H
//
// This is the declaration of the VSSVertex class.
//
#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The VSSVertex class is the implementation of the vector-scalar-scalar
 *  vertex. It inherits from the VertexBase class for storage of the particles 
 *  and implements the helicity calculations.
 *
 *  All such vertices should inherit from this class and implement the virtual
 *  setCoupling member
 *
 *  The form of the vertex is
 * \f[-ic\left(p_2-p_3\right)\cdot\epsilon_1\phi_2\phi_3\f]
 *
 *  @see VertexBase
 */
class VSSVertex: public VertexBase {
    
public:
  
  /**
   * Default constructor.
   */
  inline VSSVertex();
  
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
   * @param vec1 The wavefunction for the vector.
   * @param sca2 The wavefunction for the first  scalar.
   * @param sca3 The wavefunction for the second scalar.
   */
  Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
		   const ScalarWaveFunction & sca2, const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sca2 The wavefunction for the first  scalar.
   * @param sca3 The wavefunction for the second scalar.
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const ScalarWaveFunction & sca2,
			      const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the vector.
   * @param sca2 The wavefunction for the scalar.
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const ScalarWaveFunction & sca2);
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
  
private:
  
  /**
   * Describe an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<VSSVertex> initVSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VSSVertex & operator=(const VSSVertex &);
  
};

}
}

#include "VSSVertex.icc"

namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of VSSVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::VSSVertex,1> {
  /** Typedef of the base class of VSSVertex. */
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::VSSVertex>
    : public ClassTraitsBase<Herwig::Helicity::VSSVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::VSSVertex"; }
  };
  
}

#endif /* HERWIG_VSSVertex_H */
