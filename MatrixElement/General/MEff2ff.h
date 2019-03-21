// -*- C++ -*-
//
// MEff2ff.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEff2ff_H
#define HERWIG_MEff2ff_H
//
// This is the declaration of the MEff2ff class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;

/**
 * This is the implementation of the \f$ 2\to 2\f$ matrix element for
 * a \f$ \Psi \Psi \to \Psi \Psi\f$ process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions.
 *
 * @see \ref MEff2ffInterfaces "The Interfaces"
 * defined for MEff2ff.
 * @see GeneralHardME
 */
class MEff2ff: public GeneralHardME {

public:
  
  /** Vector of SpinorWaveFunctions. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunctions. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

public:

  /**
   * The default constructor.
   */
  MEff2ff() : scalar_(0), vector_(0), tensor_(0), spin_(4), sbar_(4) 
  {}

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

private:
  
  /** @name Functions to compute the ProductionMatrixElement. */
  //@{
  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\Psi\bar{\Psi}\f$
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ffb2ffbHeME(double & me2, bool first) const;

  /**
   * Compute the matrix element for \f$\Psi\Psi\to\Psi\Psi\f$
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement ff2ffHeME(double & me2, bool first) const;
  
  /**
   * Compute the matrix element for 
   * \f$\bar{\Psi}\bar{\Psi}\to\bar{\Psi}\bar{\Psi}\f$
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement fbfb2fbfbHeME(double & me2, bool first) const;

  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\lambda\lambda\f$
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement 
  ffb2mfmfHeME(double & me2, bool first) const;
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2ff & operator=(const MEff2ff &) = delete;

private:
  
  /**
   * Store the vector of FFSVertex pairs
   */
  vector<pair<AbstractFFSVertexPtr, AbstractFFSVertexPtr> > scalar_;

  /**
   * Store the vector of FFVVertex pairs
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > vector_;

  /**
   * Store the vector of FFTVertex pairs
   */
  vector<pair<AbstractFFTVertexPtr, AbstractFFTVertexPtr> > tensor_;

  /**
   *  Spinors
   */
  mutable vector<vector<SpinorWaveFunction> > spin_;

  /**
   *  Barred spinors
   */
  mutable vector<vector<SpinorBarWaveFunction> > sbar_;
};

}

#endif /* HERWIG_MEff2ff_H */
