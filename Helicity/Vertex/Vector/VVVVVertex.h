// -*- C++ -*-
#ifndef HERWIG_VVVVVertex_H
#define HERWIG_VVVVVertex_H
//
// This is the declaration of the VVVVVertex class.

#include "Herwig++/Helicity/Vertex/VertexBase.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVVVVertex.fh"

namespace Herwig {
namespace Helicity{
using namespace ThePEG;

/** \ingroup Helicity
 *
 * This is the implementation of the four vector vertex. 
 * It is based on the VertexBase class for the storage of particles 
 * which are allowed to interact at the vertex.
 * Classes implementation a specific vertex should inherit from this 
 * one and implement the virtual setCoupling member.
 *
 * The form of the vertex is
 * \f[ic^2\left[
 *   2\epsilon_1\cdot\epsilon_2\epsilon_3\cdot\epsilon_4
 *  - \epsilon_1\cdot\epsilon_3\epsilon_2\cdot\epsilon_4
 *  - \epsilon_1\cdot\epsilon_4\epsilon_2\cdot\epsilon_3
 * \right]\f]
 *  optional the additional diagrams from the three point vertices can be included.
 *
 * @see VertexBase
 */
class VVVVVertex: public VertexBase {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline VVVVVertex();

  /**
   * Copy-constructor.
   */
  inline VVVVVertex(const VVVVVertex &);

  /**
   * Destructor.
   */
  virtual ~VVVVVertex();
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
   * @param iopt Evaluation option, 0 just evaluate the four point vertex, 1
   * include all the three point diagrams as well.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param vec4 The wavefunction for the fourth vector.
   */
  Complex evaluate(Energy2 q2, int iopt,
		   const VectorWaveFunction & vec1, const VectorWaveFunction & vec2,
		   const VectorWaveFunction & vec3, const VectorWaveFunction & vec4);
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
  
protected:
  
  /**
   * Set the order of the particles.
   * @param id1 The PDG code of the first  particle.
   * @param id2 The PDG code of the second particle.
   * @param id3 The PDG code of the third  particle.
   * @param id4 The PDG code of the fourth particle.
   */
  inline void setOrder(int id1,int id2,int id3,int id4);

  /**
   * Set the type of the vertex.
   * @param itype The type of vertex (QCD=1 or electroweak=2).
   */
  inline void setType(int itype);

  /**
   * Set the intermediate particles if including s/u/t channel terms.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param c1 The coupling for the first  particle.
   * @param c2 The coupling for the second particle.
   */
  inline void setIntermediate(tcPDPtr part1,tcPDPtr part2,Complex c1,Complex c2);
  
private:
  
  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<VVVVVertex> initVVVVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVVVertex & operator=(const VVVVVertex &);
  
private:

  /**
   * Type of vertex 1=QCD 2=EW.
   */
  int _itype;

  /**  
   * Order of the particles.
   */
  int _iorder[4];

  /**
   *  Intermediate particles
   */
  tcPDPtr _inter[2];

  /**
   * Couplings of the intermediate particles.
   */
  Complex _coup[2];

};
}
}

#include "VVVVVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of VVVVVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::VVVVVertex,1> {
  /** Typedef of the base class of VVVVVertex. */
    typedef Herwig::Helicity::VertexBase NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::VVVVVertex>
    : public ClassTraitsBase<Herwig::Helicity::VVVVVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::Helicity::VVVVVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwVVertex.so"; }

  };
  
}

#endif /* HERWIG_VVVVVertex_H */
