// -*- C++ -*-
#ifndef HERWIG_FFVVertex_H
#define HERWIG_FFVVertex_H
//
// This is the declaration of the FFVVertex class.

#include <Herwig++/Helicity/Vertex/VertexBase.h>
#include <Herwig++/Helicity/WaveFunction/VectorWaveFunction.h>
#include <Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h>
#include <Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h>
#include "FFVVertex.fh"

namespace Herwig {
 
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;

namespace Helicity{

using namespace ThePEG; 

/** \ingroup Helicity
 *
 *  The FFVVertex class is the base class for all helicity amplitude
 *  vertices which used he renormalisable form for the 
 *  fermion-fermion-vector vertex. 
 *
 *  Any such vertices should inherit from this class and implement the virtual
 *  setcoupling member function. The base VertexBase class is used to store the
 *  particles allowed to interact at the vertex.
 *
 *  The form of the vertex is
 *  \f[ic\bar{f_2}\gamma^\mu a^\lambda P_\lambda f_1\epsilon_{3\mu}\f]
 *
 *  @see VertexBase
 */
class FFVVertex: public VertexBase {
  
public:
  
  /**
   * Default constructor.
   */
  inline FFVVertex();
  
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
   * @param vec3  The wavefunction for the vector.
   */
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2,const VectorWaveFunction & vec3);

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param drep The Dirac matrix representation
   */
  SpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				 const SpinorBarWaveFunction & sbar2,
				 const VectorWaveFunction & vec3,
				 DiracRep drep=defaultDRep);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const SpinorBarWaveFunction & sbar2);

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param vec3  The wavefunction for the vector.
   * @param drep The Dirac matrix representation
   */
  SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const VectorWaveFunction & vec3,DiracRep drep=defaultDRep);
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
  inline void setLeft(const Complex &);

  /**
   * Set the right coupling.
   */
  inline void setRight(const Complex &);

  /**
   * Get the left coupling.
   */
  inline const Complex & getLeft();

  /**
   * Get the right coupling.
   */
  inline const Complex & getRight();
  //@}

private:
  
  /**
   * Describe an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<FFVVertex> initFFVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  FFVVertex & operator=(const FFVVertex &);
  
private:

  /**
   * Left coupling.
   */
  Complex _left;

  /**
   * Right coupling.
   */
  Complex _right;
  
};
}
}

#include "FFVVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of FFVVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::FFVVertex,1> {
  /** Typedef of the base class of FFVVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::FFVVertex>
  : public ClassTraitsBase<Herwig::Helicity::FFVVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::FFVVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVVertex.so"; }

};

}


#endif /* HERWIG_FFVVertex_H */
