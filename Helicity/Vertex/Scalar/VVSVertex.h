// -*- C++ -*-
#ifndef HERWIG_VVSVertex_H
#define HERWIG_VVSVertex_H
//
// This is the declaration of the VVSVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVSVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 * The VVSVertex class is the implementation of the vector-vector-scalar.
 * It inherits from the VertexBase class for the storage of the particles and
 * implements the helicity calculations.
 *
 * All interactions of this type should inherit from it and implement the virtual
 * setCoupling member.
 *
 *  The form of the vertex is
 *  \f[icg^{\mu\nu}\epsilon_{1\mu}\epsilon_{2\nu}\f]
 *
 * @see VertexBase
 */
class VVSVertex: public VertexBase {
  
public:
  
  /**
   * Default constructor.
   */
  inline VVSVertex();

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
   * @param sca3 The wavefunction for the scalar.
   */
  Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2, const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the vector.
   * @param sca3 The wavefunction for the scalar.
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const VectorWaveFunction & vec2);
  //@}

  /**
   * Calculate the couplings. This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  scalar.
   * @param part2 The ParticleData pointer for the second scalar.
   * @param part3 The ParticleData pointer for the third  scalar.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3)=0;
  
private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<VVSVertex> initVVSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVSVertex & operator=(const VVSVertex &);
  
};

}
}
#include "VVSVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of VVSVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::VVSVertex,1> {
  /** Typedef of the base class of VVSVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::VVSVertex>
  : public ClassTraitsBase<Herwig::Helicity::VVSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::VVSVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSVertex.so"; }

};

}


#endif /* HERWIG_VVSVertex_H */
