// -*- C++ -*-
#ifndef HERWIG_MEff2tv_H
#define HERWIG_MEff2tv_H
//
// This is the declaration of the MEff2tv class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVTVertex.h"

namespace Herwig {

using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::TensorWaveFunction;

/**
 * Here is the documentation of the MEff2tv class.
 *
 * @see \ref MEff2tvInterfaces "The interfaces"
 * defined for MEff2tv.
 */
class MEff2tv: public GeneralHardME {

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

  /**
   * A vector of VectorWaveFunctions 
   */
  typedef vector<TensorWaveFunction> TBVector;
  //@}

public:

  /**
   * The default constructor.
   */
  MEff2tv() : fermion_(0), vector_(0), fourPoint_(0) {}

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2tv & operator=(const MEff2tv &) = delete;

private:

  /** @name Functions to compute the ProductionMatrixElement. */
  //@{
  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\Psi\bar{\Psi}\f$
   * @param sp Spinors for first incoming particle
   * @param spbar SpinorBar Wavefunctions for second incoming particle
   * @param vec VectorWaveFunctions for outgoing vector
   * @param ten Outgoing TensorWaveFunction
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ffb2tvHeME(SpinorVector & sp, SpinorBarVector & spbar,
	     TBVector & ten, VBVector & vec,
	     double & me2,bool first) const;
  //@}

  /**
   * A debugging function to test the value of me2 against an
   * analytic function.
   * @param me2 The value of the \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  virtual void debug(double me2) const;

private:

  /**
   * Store a pair of  FFTVertex and FFVVertex pointers  
   */
  vector<pair<AbstractFFTVertexPtr, AbstractFFVVertexPtr> > fermion_;

  /**
   *  Store a pair of FFTVertex and VVTVertex pointers
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVVTVertexPtr> > vector_;

  /**
   *  The four point vertex
   */
  vector<AbstractFFVTVertexPtr> fourPoint_;

};

}

#endif /* HERWIG_MEff2tv_H */
