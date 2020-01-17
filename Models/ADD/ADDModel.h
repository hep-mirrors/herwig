// -*- C++ -*-
//
// ADDModel.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ADDModel_H
#define HERWIG_ADDModel_H
// This is the declaration of the ADDModel class.

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVTVertex.h"
#include "ADDModel.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Models
 *
 *  This is the class to be used instead of the Standard Model class for
 *  the ADD model.
 *
 * @see \ref ADDModelInterfaces "The interfaces"
 * defined for ADDModel.
 * @see StandardModel
 * @see StandardModelBase
 * 
 */
class ADDModel: public BSMModel {

public:

  /**
   * The default constructor 
   */
  ADDModel() : delta_(2), mPlanckBar_(2.4e18*GeV),
	       md_(1000.*GeV), lambdaT_(1000.*GeV) {
    useMe();
  }
  
  /**
   * Number of extrac dimensions
   */
  unsigned int delta() const {return delta_;}

  /**
   *  The reduced Planck mass in 4d
   */
  Energy MPlanckBar() const {return mPlanckBar_;}

  /**
   *  The d-dimension Planck mass
   */
  Energy MD() const {return md_;}

  /**
   *  The cut-off for virtual gravition processes
   */
  Energy LambdaT() const {return lambdaT_;}

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
     * Private and non-existent assignment operator.
     */
  ADDModel & operator=(const ADDModel &) = delete;

private:
  
  /**
   * Number of extrac dimensions
   */
  unsigned int delta_;

  /**
   *  The reduced Planck mass in 4d
   */
  Energy mPlanckBar_;

  /**
   *  The d-dimension Planck mass
   */
  Energy md_;

  /**
   *  Cut-off parameter for virtual gravitons
   */
  Energy lambdaT_;

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

#endif /* HERWIG_ADDModel_H */
