// -*- C++ -*-
#ifndef HERWIG_RSModelVVVGRVertex_H
#define HERWIG_RSModelVVVGRVertex_H
//
// This is the declaration of the RSModelVVVGRVertex class.

#include "ThePEG/Helicity/Vertex/Tensor/VVVTVertex.h"
#include "Herwig++/Models/RSModel/RSModel.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 *
 *  The RSModelVVVGRVertex class is the implementation of the 
 *  triple vector graviton couling in the RS model. 
 *  It inherits from VVVTVertex and implements the setCoupling member.
 *
 *  @see VVVTVertex
 *  @see VertexBase
 */
class RSModelVVVGRVertex: public VVVTVertex {
  
public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline RSModelVVVGRVertex();
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
   * @param part4 The ParticleData pointer for the foruth particle.
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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<RSModelVVVGRVertex> initRSModelVVVGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelVVVGRVertex & operator=(const RSModelVVVGRVertex &);

private:

  /**
   * Pointer to the model object.
   */
  tcSMPtr _theModel;

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   * The graviton coupling.
   */
  InvEnergy _theKappa;

  /**
   *  The last value of the coupling/
   */
  vector<Complex> _couplast;

  /**
   *  The last value of the scale, \f$q^2\f$.
   */
  vector<Energy2> _q2last;

  /**
   *  The prefactor for the \f$Z\f$ vertex.
   */
  double _zfact;
  //@}
};
}


#include "RSModelVVVGRVertex.icc"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModelVVVGRVertex. 
 */
template <>
struct BaseClassTrait<Herwig::RSModelVVVGRVertex,1> {
    /** Typedef of the base class of RSModelVVVGRVertex. */
  typedef ThePEG::Helicity::VVVTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::RSModelVVVGRVertex>
  : public ClassTraitsBase<Herwig::RSModelVVVGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig++::RSModelVVVGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwRSModel.so"; }

};

/** @endcond */

}


#endif /* HERWIG_RSModelVVVGRVertex_H */
