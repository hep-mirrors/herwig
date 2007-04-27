// -*- C++ -*-
#ifndef HERWIG_VVVVertex_H
#define HERWIG_VVVVertex_H
//
// This is the declaration of the VVVVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVVVertex.fh"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;
  
/** \ingroup Helicity
 *
 *  The VVVVertex class is the base class for tripe vectro vertices
 *  using the perturbative form in Herwig++. 
 *  It inherits from the VertexBase class for the storage of the 
 *  particles allowed at the vertex.
 *
 *  Classes which implement a specific vertex should inherit from this and
 *  implement the virtual setCoupling member.
 *
 *  The form of the vertex is
 *  \f[ig\left[  (p_1-p_2)^\gamma g^{\alpha\beta }
 *              +(p_2-p_3)^\alpha g^{\beta \gamma}
 *              +(p_3-p_1)^\beta  g^{\alpha\gamma}
 *   \right]\epsilon_{1\alpha}\epsilon_{2\beta}\epsilon_{3\gamma}\f]
 *
 *  @see VertexBase
 */
class VVVVertex: public VertexBase {
    
public:
  
  /**
   * Default constructor.
   */
  inline VVVVertex();

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
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const VectorWaveFunction & vec3);
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
   * Describe an abstract base class with no persistent data.
   */
  static AbstractNoPIOClassDescription<VVVVertex> initVVVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVVertex & operator=(const VVVVertex &);
  
};

}
}

#include "VVVVertex.icc"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of VVVVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VVVVertex,1> {
  /** Typedef of the base class of VVVVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VVVVertex>
  : public ClassTraitsBase<Herwig::Helicity::VVVVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::VVVVertex"; }
};

/** @endcond */

}


#endif /* HERWIG_VVVVertex_H */
