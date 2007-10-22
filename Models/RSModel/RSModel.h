// -*- C++ -*-
#ifndef HERWIG_RSModel_H
#define HERWIG_RSModel_H
// This is the declaration of the RSModel class.

#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Helicity/Vertex/Tensor/FFTVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/VVTVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/SSTVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/FFVTVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/VVVTVertex.h"
#include "RSModel.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Models
 *
 *  This is the class to be used instead of the Standard Model class for
 *  the Randell Sundrum model.
 *
 * @see \ref RSModelInterfaces "The interfaces"
 * defined for RSModel.
 * @see StandardModel
 * @see StandardModelBase
 * 
 */
class RSModel: public StandardModel {

public:

  /**
   * The default constructor 
   */
  RSModel();
  
  /**
   * Return the gravition coupling
   */
  inline Energy lambda_pi() const;


  /** @name Vertices */
  //@{
  /**
   * Pointer to the object handling the \f$G\to f\bar{f}\f$ vertex.
   */
  inline tFFTVertexPtr   vertexFFGR() const;

  /**
   * Pointer to the object handling the \f$G\to VV\f$ vertex.
   */
  inline tVVTVertexPtr   vertexVVGR() const;

  /**
   * Pointer to the object handling the \f$G\to SS\f$ vertex.
   */
  inline tSSTVertexPtr   vertexSSGR() const;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}V\f$ vertex.
   */
  inline tFFVTVertexPtr  vertexFFVGR() const;

  /**
   * Pointer to the object handling the \f$G\to VVV\f$ vertex.
   */
  inline tVVVTVertexPtr  vertexVVVGR() const;
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}
  
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

private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<RSModel> initRSModel;
  
    /**
     * Private and non-existent assignment operator.
     */
  RSModel & operator=(const RSModel &);
  
  /**
   * Coupling of the graviton
   */
  Energy _theLambda_pi;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}\f$ vertex.
   */
  FFTVertexPtr  _theFFGRVertex;

  /**
   * Pointer to the object handling the \f$G\to VV\f$ vertex.
   */
  VVTVertexPtr  _theVVGRVertex;

  /**
   * Pointer to the object handling the \f$G\to SS\f$ vertex.
   */
  SSTVertexPtr  _theSSGRVertex;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}V\f$ vertex.
   */
  FFVTVertexPtr _theFFVGRVertex;

  /**
   * Pointer to the object handling the \f$G\to VVV\f$ vertex.
   */
  VVVTVertexPtr _theVVVGRVertex;
  
};
}
#include "RSModel.icc"


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModel.
 */
template <>
struct BaseClassTrait<Herwig::RSModel,1> {
  /** Typedef of the base class of RSModel. */
  typedef Herwig::StandardModel NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::RSModel>
  : public ClassTraitsBase<Herwig::RSModel> {
  /** Return the class name.*/
  static string className() { return "Herwig::RSModel"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwRSModel.so"; }
  
};

/** @endcond */
  
}


#endif /* HERWIG_RSModel_H */
