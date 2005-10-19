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
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline FFVVertex();

  /**
   * Copy-constructor.
   */
  inline FFVVertex(const FFVVertex &);

  /**
   * Destructor.
   */
  virtual ~FFVVertex();
  //@}
  
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

protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}
  
private:
  
  /**
   * Describe an abstract class with persistent data.
   */
  static AbstractClassDescription<FFVVertex> initFFVVertex;
  
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
  static string className() { return "Herwig++::Helicity::FFVVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwVVertex.so"; }

};

}


#endif /* HERWIG_FFVVertex_H */
