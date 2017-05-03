// -*- C++ -*-
//
// MEff2sv.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEff2sv_H
#define HERWIG_MEff2sv_H
//
// This is the declaration of the MEff2sv class.
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
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::ScalarWaveFunction;

/**
 * The MEff2sv class is designed to implement the matrix element for a
 * fermion-antifermion to vector-scalar hard process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions for this 
 * specific spin combination.
 *
 * @see \ref MEff2svInterfaces "The interfaces"
 * defined for MEff2sv.
 * @see GeneralHardME
 */
class MEff2sv: public GeneralHardME {

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
  MEff2sv() : scalar_(0), vector_(0), fermion_(0) {}

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
  static ClassDescription<MEff2sv> initMEff2sv;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2sv & operator=(const MEff2sv &);

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
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ffb2svHeME(SpinorVector & sp, SpinorBarVector & spbar,
	     ScalarWaveFunction & sca, VBVector & vec, 
	     double & me2,bool first) const;
  //@}


private:

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * scalar
   */
  vector<pair<AbstractFFSVertexPtr, AbstractVSSVertexPtr> > scalar_;

  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * vector
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVVSVertexPtr> > vector_;
  
  /**
   * Storage for dynamically cast vertices for a diagram with intermediate
   * fermion
   */
  vector<pair<AbstractFFSVertexPtr, AbstractFFVVertexPtr> > fermion_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEff2sv. */
template <>
struct BaseClassTrait<Herwig::MEff2sv,1> {
  /** Typedef of the first base class of MEff2sv. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2sv class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2sv>
  : public ClassTraitsBase<Herwig::MEff2sv> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEff2sv"; }
};

/** @endcond */

}

#endif /* HERWIG_MEff2sv_H */
