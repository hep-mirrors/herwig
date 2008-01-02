// -*- C++ -*-
//
// MEff2vs.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEff2vs_H
#define HERWIG_MEff2vs_H
//
// This is the declaration of the MEff2vs class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ProductionMatrixElement.h"
#include "MEff2vs.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::ScalarWaveFunction;

/**
 * The MEff2vs class is designed to implement the matrix element for a
 * fermion-antifermion to vector-scalar hard process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions for this 
 * specific spin combination.
 *
 * @see \ref MEff2vsInterfaces "The interfaces"
 * defined for MEff2vs.
 * @see GeneralHardME
 */
class MEff2vs: public GeneralHardME {

public:

  /** @name Typedefs */
  //@{
  /**
   * A vector of SpinorWaveFunctions 
   */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /**
   * A vector of SpinorWaveBarFunctions 
   */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

  /**
   * A vector of VectorWaveFunctions 
   */
  typedef vector<VectorWaveFunction> VBVector;
  //@}

public:

  /**
   * The default constructor.
   */
  inline MEff2vs();

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

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}

  /**
   * Construct the vertex information for the spin correlations
   * @param sub Pointer to the relevent SubProcess
   */
  virtual void constructVertex(tSubProPtr sub);

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
  virtual void doinit() throw(InitException);
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEff2vs> initMEff2vs;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2vs & operator=(const MEff2vs &);

private:

  /** @name Functions to compute the ProductionMatrixElement. */
  //@{
  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\Psi\bar{\Psi}\f$
   * @param sp Spinors for first incoming particle
   * @param spbar SpinorBar Wavefunctions for second incoming particle
   * @param vec VectorWaveFunctions for outgoing vector
   * @param sca Outgoing ScalarWaveFunction
   * @param me2 colour averaged, spin summed ME
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ffb2vsHeME(SpinorVector & sp, SpinorBarVector & spbar,
	     VBVector & vec, ScalarWaveFunction & sca, 
	     double & me2) const;
  //@}


private:

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * scalar
   */
  vector<pair<AbstractFFSVertexPtr, AbstractVSSVertexPtr> > theSca;

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * vector
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVVSVertexPtr> > theVec;
  
  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * fermion
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFSVertexPtr> > theFerm;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEff2vs. */
template <>
struct BaseClassTrait<Herwig::MEff2vs,1> {
  /** Typedef of the first base class of MEff2vs. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2vs class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2vs>
  : public ClassTraitsBase<Herwig::MEff2vs> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEff2vs"; }
};

/** @endcond */

}

#include "MEff2vs.icc"

#endif /* HERWIG_MEff2vs_H */
