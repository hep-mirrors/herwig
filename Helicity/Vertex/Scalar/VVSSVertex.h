// -*- C++ -*-
#ifndef HERWIG_VVSSVertex_H
#define HERWIG_VVSSVertex_H
//
// This is the declaration of the VVSSVertex class.
//
#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
namespace Helicity {

/** \ingroup Helicity
 *
 *  The VVSSVertex class is the implementation of the coupling of two
 *  vectors and two scalars. It inherits from the VertexBase class for the 
 *  storage of the particles and implements the helicity calculations.
 *
 *  All classes implementing the vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  The form of the vertex is \f[icg^{\mu\nu}\epsilon_{1\mu}\epsilon_{2\nu}\f]
 *
 *  @see VertexBase
 */
class VVSSVertex: public VertexBase {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline VVSSVertex();

  /**
   * Copy-constructor.
   */
  inline VVSSVertex(const VVSSVertex &);

  /**
   * Destructor.
   */
  virtual ~VVSSVertex();
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
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param sca3 The wavefunction for the first  scalar.
   * @param sca4 The wavefunction for the second scalar.
   */
  Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2, const ScalarWaveFunction & sca3,
		   const ScalarWaveFunction & sca4);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the second vector.
   * @param sca3 The wavefunction for the first  scalar.
   * @param sca4 The wavefunction for the second scalar.
   */
  VectorWaveFunction evaluate(Energy2 q2, int iopt,tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const ScalarWaveFunction & sca3,
			      const ScalarWaveFunction & sca4);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param sca3 The wavefunction for the second scalar.
   */
  ScalarWaveFunction evaluate(Energy2 q2, int iopt,tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const VectorWaveFunction & vec2,
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
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4)=0;
  
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
  static AbstractClassDescription<VVSSVertex> initVVSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVSSVertex & operator=(const VVSSVertex &);
  
};

}
}

#include "VVSSVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of VVSSVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVSSVertex,1> {
  /** Typedef of the base class of VVSSVertex. */
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::VVSSVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVSSVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::VVSSVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwSVertex.so"; }

  };

}


#endif /* HERWIG_VVSSVertex_H */
