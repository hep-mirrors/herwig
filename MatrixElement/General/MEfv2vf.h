// -*- C++ -*-
//
// MEfv2vf.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEfv2vf_H
#define HERWIG_MEfv2vf_H
//
// This is the declaration of the MEfv2vf class.
//

#include "GeneralHardME.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVVertex.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the matrix element for a fermion and a vector
 * boson to a fermion and a vector boson. It inherits from GeneralHardME
 * and implements the appropriate virtual functions. 
 * 
 * @see GeneralHardME
 *
 */
class MEfv2vf: public GeneralHardME {

public:

  /** A vector of SpinorWaveFunctions. */
  typedef vector<Helicity::SpinorWaveFunction> SpinorVector;

  /** A vector of SpinorBarWaveFunctions. */
  typedef vector<Helicity::SpinorBarWaveFunction> SpinorBarVector;

  /** A vector of VectorWaveFunctions. */
  typedef vector<Helicity::VectorWaveFunction> VBVector;

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

  /** @name Functions to calculate the Helicity MatrixElement.*/
  //@{
  /**
   * Calculate the matrix element for an incoming fermion
   * @param spIn A vector of spinors for the incoming fermion
   * @param vecIn A vector of VectorWaveFunctions for the incoming boson
   * @param spbOut A vector of SpinorBarWaveFunctions for the outgoing fermion
   * @param vecOut A vector of VectorWaveFunctions for the outgoing boson
   * @param mc If the outgoing vector is massless or not
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param mesq The matrix element squared
  */
  ProductionMatrixElement
  fv2vfHeME(const SpinorVector & spIn,  const VBVector & vecIn, 
	    const VBVector & vecOut, bool mc,
	    const SpinorBarVector & spbOut, 
	    double & mesq, bool first) const;

  /**
   * Calculate the matrix element for an incoming anti-fermion
   * @param spbIn A vector of SpinorBarWaveFunctions for the incoming anti-fermion
   * @param vecIn A vector of VectorWaveFunctions for the incoming boson
   * @param spOut A vector of Spinors for the outgoing antifermion
   * @param vecOut A vector of VectorWaveFunctions for the outgoing boson
   * @param mc If the outgoing vector is massless or not
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param mesq The matrix element squared
  */
  ProductionMatrixElement
  fbv2vfbHeME(const SpinorBarVector & spbIn,  const VBVector & vecIn, 
	      const VBVector & vecOut, bool mc,
	      const SpinorVector & spOut, 
	      double & mesq, bool first) const;
  //@}

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
  MEfv2vf & operator=(const MEfv2vf &) = delete;

private:
  
  /** @name Store dynamically casted vertices. */
  //@{
  /**
   * A pair off FFVVertex pointers 
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > fermion_;

  /**
   * A pair of FFVVertex, VVVertex pointers 
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVVVVertexPtr> > vector_;

  /**
   *  Four point vertices
   */
  vector<AbstractFFVVVertexPtr> four_;
  //@}

};

}

#endif /* HERWIG_MEfv2vf_H */
