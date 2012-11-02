// -*- C++ -*-
//
// MEff2ss.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEff2ss_H
#define HERWIG_MEff2ss_H
//
// This is the declaration of the MEff2ss class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::ScalarWaveFunction;

/**
 * The MEff2ss class is designed to implement the matrix element for a
 * fermion-antifermion to scalar-scalar hard process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions for this 
 * specific spin combination.
 *
 * @see \ref MEff2ssInterfaces "The interfaces"
 * defined for MEff2ss.
 * @see GeneralHardME
 */
class MEff2ss: public GeneralHardME {

public:

  /** Vector of SpinorWaveFunctions objects */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunction objects. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

public:

  /**
   * The default constructor.
   */
  MEff2ss() : fermion_(0), vector_(0), tensor_(0) {}

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
  static ClassDescription<MEff2ss> initMEff2ss;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2ss & operator=(const MEff2ss &);

private:

  /**
   * Calculate the matrix element
   * @param sp A vector of SpinorWaveFunction objects
   * @param sbar A vector of SpinorBarWaveFunction objects
   * @param sca1 A ScalarWaveFunction for an outgoing scalar
   * @param sca2 A ScalarWaveFunction for the other outgoing scalar
   * @param me2 The spin averaged matrix element
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   */
  ProductionMatrixElement ff2ssME(const SpinorVector & sp, 
				  const SpinorBarVector & sbar, 
				  const ScalarWaveFunction & sca1,
				  const ScalarWaveFunction & sca2,
				  double & me2, bool first) const;

private:

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * fermion
   */
  vector<pair<AbstractFFSVertexPtr, AbstractFFSVertexPtr> > fermion_;

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * vector
   */
  vector<pair<AbstractFFSVertexPtr, AbstractSSSVertexPtr> > scalar_;

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * vector
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVSSVertexPtr> > vector_;
  
  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * tensor
   */
  vector<pair<AbstractFFTVertexPtr, AbstractSSTVertexPtr> > tensor_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEff2ss. */
template <>
struct BaseClassTrait<Herwig::MEff2ss,1> {
  /** Typedef of the first base class of MEff2ss. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2ss class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2ss>
  : public ClassTraitsBase<Herwig::MEff2ss> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEff2ss"; }
};

/** @endcond */

}

#endif /* HERWIG_MEff2ss_H */
