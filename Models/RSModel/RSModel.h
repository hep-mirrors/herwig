// -*- C++ -*-
//
// RSModel.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RSModel_H
#define HERWIG_RSModel_H
// This is the declaration of the RSModel class.

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVTVertex.h"
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
class RSModel: public BSMModel {

public:

  /**
   * The default constructor 
   */
  RSModel(): Lambda_pi_(10000*GeV) {
    useMe();
  }
  
  /**
   * Return the gravition coupling
   */
  Energy lambda_pi() const {return Lambda_pi_;}


  /** @name Vertices */
  //@{
  /**
   * Pointer to the object handling the \f$G\to f\bar{f}\f$ vertex.
   */
  tAbstractFFTVertexPtr   vertexFFGR() const {return FFGRVertex_;}

  /**
   * Pointer to the object handling the \f$G\to VV\f$ vertex.
   */
  tAbstractVVTVertexPtr   vertexVVGR() const {return VVGRVertex_;}

  /**
   * Pointer to the object handling the \f$G\to SS\f$ vertex.
   */
  tAbstractSSTVertexPtr   vertexSSGR() const {return SSGRVertex_;}

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}g\f$ vertex.
   */
  tAbstractFFVTVertexPtr  vertexFFGGR() const {return FFGGRVertex_;}

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}W^\pm/Z^0/\gamma\f$ vertex.
   */
  tAbstractFFVTVertexPtr  vertexFFWGR() const {return FFWGRVertex_;}

  /**
   * Pointer to the object handling the \f$G\to W^+W^-Z^0/\gamma\f$ vertex.
   */
  tAbstractVVVTVertexPtr  vertexWWWGR() const {return WWWGRVertex_;}

  /**
   * Pointer to the object handling the \f$G\to ggg\f$ vertex.
   */
  tAbstractVVVTVertexPtr  vertexGGGGR() const {return GGGGRVertex_;}
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
  virtual void doinit();
  //@}
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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

private:
  
  /**
   * Coupling of the graviton
   */
  Energy Lambda_pi_;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}\f$ vertex.
   */
  AbstractFFTVertexPtr  FFGRVertex_;

  /**
   * Pointer to the object handling the \f$G\to VV\f$ vertex.
   */
  AbstractVVTVertexPtr  VVGRVertex_;

  /**
   * Pointer to the object handling the \f$G\to SS\f$ vertex.
   */
  AbstractSSTVertexPtr  SSGRVertex_;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}g\f$ vertex.
   */
  AbstractFFVTVertexPtr FFGGRVertex_;

  /**
   * Pointer to the object handling the \f$G\to f\bar{f}W/Z^0\gamma\f$ vertex.
   */
  AbstractFFVTVertexPtr FFWGRVertex_;

  /**
   * Pointer to the object handling the \f$G\to W^+W^-Z^0\gamma\f$ vertex.
   */
  AbstractVVVTVertexPtr WWWGRVertex_;

  /**
   * Pointer to the object handling the \f$G\to ggg\f$ vertex.
   */
  AbstractVVVTVertexPtr GGGGRVertex_;
  
};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModel.
 */
template <>
struct BaseClassTrait<Herwig::RSModel,1> {
  /** Typedef of the base class of RSModel. */
  typedef Herwig::BSMModel NthBase;
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
