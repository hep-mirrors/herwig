// -*- C++ -*-
#ifndef HERWIG_SMFFWVertex_H
#define HERWIG_SMFFWVertex_H
//
// This is the declaration of the SMFFWVertex class.

#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/StandardModel/CKMBase.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"
namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the Standard model coupling 
 *  of the W to the fermions.
 *
 *  @see FFVVertex
 *  @see VertexBase
 */
class SMFFWVertex: public FFVVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline SMFFWVertex();

  /**
   * Copy-constructor.
   */
  inline SMFFWVertex(const SMFFWVertex &);

  /**
   * Destructor.
   */
  virtual ~SMFFWVertex();
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
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
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
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SMFFWVertex> initSMFFWVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMFFWVertex & operator=(const SMFFWVertex &);

private:

  /**
   * Pointer to the Standard Model object.
   */
  SMPtr _theSM;

  /**
   * Pointer to the CKM object.
   */
  Ptr<CKMBase>::pointer _theCKM;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The elements of the CKM matrix.
   */
  Complex _ckm[3][3];

  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;
  //@}
};

}
}
#include "SMFFWVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SMFFWVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMFFWVertex,1> {
    /** Typedef of the base class of SMFFWVertex. */
    typedef Herwig::Helicity::FFVVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMFFWVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMFFWVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::Helicity::SMFFWVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */
    static string library() { return "libHwSMVertex.so"; }

  };
  
}

#endif /* HERWIG_SMFFWVertex_H */
