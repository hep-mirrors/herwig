// -*- C++ -*-
#ifndef HERWIG_FFSVertex_H
#define HERWIG_FFSVertex_H
//
// This is the declaration of the FFSVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "FFSVertex.fh"

namespace Herwig {
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;

namespace Helicity{
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The FFSVertex class is the implementation of the interact of a
 *  scalar boson and a fermion-antifermion pair. It inherits from the VertexBase
 *  class for storage of the particles interacting at the vertex and implements
 *  the helicity calculations.
 *
 *  Implementations of specific interactions should inherit from this and implement
 *  the virtual setCoupling member.
 *
 *  The form of the vertex is
 *  \f[ic\bar{f_2}a^\lambda P_\lambda f_1\phi_3\f]
 *  where \f$a^\pm\f$ are the right and left couplings and \f$P_\pm=(1\pm\gamma_5)\f$
 *  are the chirality projection operators.
 *
 *  @see VertexBase
 */
class FFSVertex: public VertexBase {
      
public:

  /**
   * Default constructor.
   */
  inline FFSVertex();
  
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
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param sca3  The wavefunction for the scalar.
   */
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2,
		   const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param sca3  The wavefunction for the scalar.
   * @param drep The Dirac matrix representation
   */
  SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1, 
			      const ScalarWaveFunction & sca3,
			      DiracRep drep=defaultDRep);

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param sca3  The wavefunction for the scalar.
   * @param drep The Dirac matrix representation
   */
  SpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				 const SpinorBarWaveFunction & sbar2,
				 const ScalarWaveFunction & sca3,
				 DiracRep drep=defaultDRep);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1, 
			      const SpinorBarWaveFunction & sbar2);
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

protected:

  /**
   *  Set and get the couplings
   */
  //@{
  /**
   * Set the left coupling.
   */
  inline void setLeft(Complex);

  /**
   * Set the right coupling.
   */
  inline void setRight(Complex);

  /**
   * Get the left coupling.
   */
  inline Complex getLeft();

  /**
   * Get the right coupling.
   */
  inline Complex getRight();
  //@}
  
private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<FFSVertex> initFFSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  FFSVertex & operator=(const FFSVertex &);
  
private:

  /**
   * Storage of the left coupling.
   */
  Complex _left;

  /**
   * Storage of the right coupling.
   */
  Complex _right;

};
}
}

#include "FFSVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of FFSVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::FFSVertex,1> {
  /** Typedef of the base class of FFSVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/** 
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::FFSVertex>
  : public ClassTraitsBase<Herwig::Helicity::FFSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::FFSVertex"; }
};

}


#endif /* HERWIG_FFSVertex_H */
