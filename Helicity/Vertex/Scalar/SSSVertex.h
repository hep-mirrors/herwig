// -*- C++ -*-
#ifndef HERWIG_SSSVertex_H
#define HERWIG_SSSVertex_H
//
// This is the declaration of the SSSVertex class.
//
#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *  
 *  The SSSVertex class is the implementation of the interaction of
 *  three scalars. It inherits from the VertexBase class for the storage of the
 *  particles interacting at the vertex and implements the helicity calculations.
 *
 *  Any classes implementating the vertex should inherit from it and implement
 *  the virtual set Coupling member.
 *
 *  The form of the vertex is
 * \f[ic\phi_1\phi_2\phi_3\f]
 *
 *  @see VertexBase
 */
class SSSVertex: public VertexBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline SSSVertex();

  /**
   * Copy-constructor.
   */
  inline SSSVertex(const SSSVertex &);

  /**
   * Destructor.
   */
  virtual ~SSSVertex();
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
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param sca3 The wavefunction for the third  scalar.
   */
  Complex evaluate(Energy2 q2,const ScalarWaveFunction & sca1,
		   const ScalarWaveFunction & sca2,const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out, 
			      const ScalarWaveFunction & sca1,
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
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<SSSVertex> initSSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SSSVertex & operator=(const SSSVertex &);
  
};

}
}

#include "SSSVertex.icc"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of SSSVertex.
 */
template <>
struct BaseClassTrait<Herwig::Helicity::SSSVertex,1> {
  /** Typedef of the base class of SSSVertex. */
  typedef Herwig::Helicity::VertexBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::Helicity::SSSVertex>
  : public ClassTraitsBase<Herwig::Helicity::SSSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::SSSVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSVertex.so"; }

};

}


#endif /* HERWIG_SSSVertex_H */
