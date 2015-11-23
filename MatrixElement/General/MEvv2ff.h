// -*- C++ -*-
//
// MEvv2ff.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEvv2ff_H
#define HERWIG_MEvv2ff_H
//
// This is the declaration of the MEvv2ff class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;

/**
 * This class is designed to implement the matrix element for the 
 * \f$2 \rightarrow 2\f$ process vector-vector to fermion-antifermion pair. It
 * inherits from GeneralHardME and implements the me2() virtual function.
 *
 * @see \ref MEvv2ffInterfaces "The Interfaces"
 * defined for MEvv2ff.
 * @see GeneralHardME
 * 
 */
class MEvv2ff: public GeneralHardME {

public:
  
  /** A Vector of VectorWaveFunction objects. */
  typedef vector<VectorWaveFunction> VBVector;

  /** A vector of SpinorBarWaveFunction objects. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** A vector of SpinorBarWaveFunction objects. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

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
   * Construct the vertex information for the spin correlations
   * @param sub Pointer to the relevent SubProcess
   */
  virtual void constructVertex(tSubProPtr sub);

private:

  /**
   * Calculate the value of the matrix element 
   */
  ProductionMatrixElement vv2ffME(const VBVector & v1, const VBVector & v2,
				  const SpinorBarVector & sbar,
				  const SpinorVector & sp, 
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
  static ClassDescription<MEvv2ff> initMEvv2ff;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2ff & operator=(const MEvv2ff &);

private:
  
  /** @name Dynamically casted vertices. */
  //@{
  /**
   *  Intermediate scalar
   */
  vector<pair<AbstractVVSVertexPtr, AbstractFFSVertexPtr > > scalar_;
  /**
   * Intermediate fermion 
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > fermion_;

  /**
   * Intermediate vector
   */
  vector<pair<AbstractVVVVertexPtr, AbstractFFVVertexPtr> > vector_;
  
  /**
   * Intermediate tensor
   */
  vector<pair<AbstractVVTVertexPtr, AbstractFFTVertexPtr> > tensor_;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2ff. */
template <>
struct BaseClassTrait<Herwig::MEvv2ff,1> {
  /** Typedef of the first base class of MEvv2ff. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2ff class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2ff>
  : public ClassTraitsBase<Herwig::MEvv2ff> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEvv2ff"; }
};

/** @endcond */

}

#endif /* HERWIG_MEvv2ff_H */
