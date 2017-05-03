// -*- C++ -*-
//
// MEvv2ss.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEvv2ss_H
#define HERWIG_MEvv2ss_H
//
// This is the declaration of the MEvv2ss class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;

/**
 * This is the implementation of the matrix element for the process
 * vector-vector to scalar-scalar. It inherits from GeneralHardME and
 * implements the required virtual functions.
 *
 * @see \ref MEff2ffInterfaces "The Interfaces"
 * defined for MEff2ff.
 * @see GeneralHardME
 */
class MEvv2ss: public GeneralHardME {

public:

  /** A vector of VectorWaveFunction objects*/
  typedef vector<VectorWaveFunction> VBVector;

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
  //@}

  /**
   * Set the Hardvertex for the spin correlations
   * @param sub
   */
  virtual void constructVertex(tSubProPtr sub);

private:

  /**
   * Calculate the matrix element.
   * @param v1 A vector of VectorWaveFunction objects for the first boson
   * @param v2 A vector of VectorWaveFunction objects for the second boson
   * @param sca1 A ScalarWaveFunction for the first outgoing
   * @param sca2 A ScalarWaveFunction for the second outgoing
   * @param me2 The value of the spin-summed matrix element squared
   * (to be calculated)
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   */
  ProductionMatrixElement vv2ssME(const VBVector & v1, const VBVector & v2,
				  const ScalarWaveFunction & sca1, 
				  const ScalarWaveFunction & sca2, 
				  double & me2, bool first) const;
protected:
  
  /**
   * A debugging function to test the value of me2 against an
   * analytic function.
   * @param me2 The value of the \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  virtual void debug(double me2) const;

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEvv2ss> initMEvv2ss;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2ss & operator=(const MEvv2ss &);

private:

  /** @name The dynamically casted vertices. */
  //@{
  /**
   * Intermediate s-channel scalar
   */
  vector<pair<AbstractVVSVertexPtr, AbstractSSSVertexPtr> > scalar1_;

  /**
   * Intermediate t-channel scalar
   */
  vector<pair<AbstractVSSVertexPtr, AbstractVSSVertexPtr> > scalar2_;

  /**
   * Intermediate t-channel scalar
   */
  vector<pair<AbstractVVSVertexPtr, AbstractVVSVertexPtr> > scalar3_;

  /**
   * Intermediate s-channel vector
   */
  vector<pair<AbstractVVVVertexPtr, AbstractVSSVertexPtr> > vector_;

  /**
   * Intermediate s-channel tensor
   */
  vector<pair<AbstractVVTVertexPtr, AbstractSSTVertexPtr> > tensor_;
  
  /**
   * The contact vertex 
   */
  AbstractVVSSVertexPtr contact_;
  //@}
  
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2ss. */
template <>
struct BaseClassTrait<Herwig::MEvv2ss,1> {
  /** Typedef of the first base class of MEvv2ss. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2ss class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2ss>
  : public ClassTraitsBase<Herwig::MEvv2ss> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEvv2ss"; }
};

/** @endcond */

}

#endif /* HERWIG_MEvv2ss_H */
