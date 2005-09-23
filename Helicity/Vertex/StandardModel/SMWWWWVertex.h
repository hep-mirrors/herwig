// -*- C++ -*-
#ifndef HERWIG_SMWWWWVertex_H
#define HERWIG_SMWWWWVertex_H
//
// This is the declaration of the SMWWWWVertex class.
//
#include "Herwig++/Helicity/Vertex/Vector/VVVVVertex.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Rebinder.h"


namespace Herwig {
namespace Helicity {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The SMWWWWVertex class is the implementation of the
 *  coupling of four electroweak gauge bosons in the SM. 
 *  It inherits from VVVVVertex nad implements the setCoupling member.
 *
 *  @see VVVVVVertex
 *  @see VertexBase
 */
class SMWWWWVertex: public VVVVVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline SMWWWWVertex();

  /**
   * Copy-constructor.
   */
  inline SMWWWWVertex(const SMWWWWVertex &);

  /**
   * Destructor.
   */
  virtual ~SMWWWWVertex();
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
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4);
  
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
  static ClassDescription<SMWWWWVertex> initSMWWWWVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SMWWWWVertex & operator=(const SMWWWWVertex &);

  /**
   * Pointer to the SM object.
   */
  SMPtr _theSM;

  /**
   *  Intermediate particles
   */
  //@{
  /**
   * The ParticleData pointer for the photon
   */
  tcPDPtr _gamma;

  /**
   * The ParticleData pointer for the \f$Z^0\f$ boson.
   */
  tcPDPtr _Z0;

  /**
   * The ParticleData pointer for the \f$W^+\f$ boson.
   */
  tcPDPtr _wplus;

  /**
   * The ParticleData pointer for the \f$W^-\f$ boson.
   */
  tcPDPtr _wminus;
  //@}

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the electroweak coupling calculated.
   */
  Complex _couplast;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 _q2last;

  /**
   *  The factors for the different bosons.
   */
  double _vfact[4];

  /**
   *  \f$\sin^2\theta_W\f$.
   */
  double _sw2;

  /**
   *  \f$\cos^2\theta_W\f$.
   */
  double _cw2;
  //@}
};
 
}
}

#include "SMWWWWVertex.icc"

namespace ThePEG {
  
  /**
   * The following template specialization informs ThePEG about the
   * base class of SMWWWWVertex.
   */
  template <>
  struct BaseClassTrait<Herwig::Helicity::SMWWWWVertex,1> {
    /** Typedef of the base class of SMWWWWVertex. */
    typedef Herwig::Helicity::VVVVVertex NthBase;
  };
  
  /**
   * The following template specialization informs ThePEG about the
   * name of this class and the shared object where it is defined.
   */
  template <>
  struct ClassTraits<Herwig::Helicity::SMWWWWVertex>
    : public ClassTraitsBase<Herwig::Helicity::SMWWWWVertex> {

    /**
     * Return the class name.
     */
    static string className() { return "Herwig++::Helicity::SMWWWWVertex"; }

    /**
     * Return the name of the shared library to be loaded to get
     * access to this class and every other class it uses
     * (except the base class).
     */ 
    static string library() { return "HwSMVertex.so"; }

  };
  
}


#endif /* HERWIG_SMWWWWVertex_H */
