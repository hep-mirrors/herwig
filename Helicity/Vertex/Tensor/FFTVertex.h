// -*- C++ -*-
#ifndef HERWIG_FFTVertex_H
#define HERWIG_FFTVertex_H
//
// This is the declaration of the FFTVertex class.
//
#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/TensorWaveFunction.h"

namespace Herwig {
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::defaultDRep;
namespace Helicity {

using namespace ThePEG;

/** \ingroup Helicity
 *
 *  The FFTVertex class is the implementation of the fermion-fermion-tensor
 *  vertex. It inherits from the VertexBase class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  The vertex has the form
 * \f[-\frac{i\kappa}8\bar{f_2}\left[
 *  \gamma_\mu(k_1-k_2)_\nu+\gamma_\nu(k_1-k_2)_\mu
 * -2g_{\mu\nu}(k\!\!\!\!\!\not\,\,\,_1-k\!\!\!\!\!\not\,\,\,_2)+4g_{\mu\nu}m_{f}
 * \right]f_1\epsilon^{\mu\nu}_3\f] 
 *
 *  @see VertexBase
 */
class FFTVertex: public VertexBase {
      
public:
    
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline FFTVertex();

  /**
   * Copy-constructor.
   */
  inline FFTVertex(const FFTVertex &);

  /**
   * Destructor.
   */
  virtual ~FFTVertex();
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
   * @param ten3  The wavefunction for the tensor.
   */
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2, 
		   const TensorWaveFunction & ten3);

  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   */
  TensorWaveFunction evaluate(Energy2 q2, int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const SpinorBarWaveFunction & sbar2);

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param ten3  The wavefunction for the tensor.
   * @param drep The Dirac matrix representation
   */
  SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const TensorWaveFunction & ten3,DiracRep drep=defaultDRep);

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param ten3  The wavefunction for the tensor.
   * @param drep The Dirac matrix representation
   */
  SpinorBarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				 const SpinorBarWaveFunction & sbar2,
				 const TensorWaveFunction& ten3,
				 DiracRep drep=defaultDRep);
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
  static AbstractClassDescription<FFTVertex> initFFTVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  FFTVertex & operator=(const FFTVertex &);
  
};

}
}

#include "FFTVertex.icc"

namespace ThePEG {

  /**
   * The following template specialization informs ThePEG about the
   * base class of FFTVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::FFTVertex,1> {
  /** Typedef of the base class of FFTVertex. */
    typedef Herwig::Helicity::VertexBase NthBase;
  };

  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::FFTVertex>
    : public ClassTraitsBase<Herwig::Helicity::FFTVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::Helicity::FFTVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwTVertex.so"; }

  };
  
}

#endif /* HERWIG_FFTVertex_H */
