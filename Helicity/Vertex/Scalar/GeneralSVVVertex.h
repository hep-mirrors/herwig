// -*- C++ -*-
#ifndef HERWIG_GeneralSVVVertex_H
#define HERWIG_GeneralSVVVertex_H
//
// This is the declaration of the GeneralSVVVertex class.
//

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "GeneralSVVVertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
/**
 * The GeneralSVVVertex class implements a
 * general Scalar-Vector-Vector vertex allowing for decay modes 
 * that only enter at the one-loop level
 * 
 * The loop integral is calculated by Passarino-Veltman reduction 
 * and the coefficients are stored here. They must be calculated
 * in the inheriting class along with implementation of the
 * setCoupling member.
 */
class GeneralSVVVertex: public VertexBase {
  
public:

  /**
   * The default constructor.
   */
  inline GeneralSVVVertex();
  
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
  
  /** @name Member functions to calculate helicity amplitudes for vertices */
  //@{
  /**
   * Evaluate the vertex
   * @param q2 Scale at which to evaluate the coupling
   * @param sca Scalar wavefunction 
   * @param vec1 Wavefunction of first vector particle
   * @param vec2 Wavefunction of second vector particle
   */
  Complex evaluate(Energy2 q2, const ScalarWaveFunction & sca,
		   const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2);
  
  /**
   * Calculate coupling.
   *@param q2 Scale at which to evaluate couplings
   *@param part1 ParticleDataPointer to first particle 
   *@param part2 ParticleDataPointer to second particle
   *@param part3 ParticleDataPointer to third particle 
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3)=0;
  //@}

protected:
  
  /**@name Tensor Coefficients access and setting functions. */
  //@{
  /**
   * Set tensor coefficients 
   */
  inline void a00(const Complex & val);
  
  inline void a11(const Complex & val);
  
  inline void a12(const Complex & val);
  
  inline void a21(const Complex & val);
  
  inline void a22(const Complex & val);
  
  inline void aEp(const Complex & val);
  
  /**
   *Access to tensor coefficients
   */
  inline Complex a00() const;
  
  inline Complex a11() const;
  
  inline Complex a12() const;
  
  inline Complex a21() const;
  
  inline Complex a22() const;
  
  inline Complex aEp() const;
  //@}
  
private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<GeneralSVVVertex> initGeneralSVVVertex;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralSVVVertex & operator=(const GeneralSVVVertex &);
  
  /**@name Store tensor coefficients. */
  //@{
  /**
   *  The tensor coefficients
   */
  Complex _a00,_a11,_a12,_a21,_a22,_aEp;
  //@}
};
}
}

#include "GeneralSVVVertex.icc"

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralSVVVertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::GeneralSVVVertex,1> {
  /** Typedef of the first base class of GeneralSVVVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralSVVVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::GeneralSVVVertex>
  : public ClassTraitsBase<Herwig::Helicity::GeneralSVVVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::GeneralSVVVertex"; }
};

/** @endcond */

}


#endif /* HERWIG_GeneralSVVVertex_H */
