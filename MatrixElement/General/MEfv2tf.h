// -*- C++ -*-
#ifndef HERWIG_MEfv2tf_H
#define HERWIG_MEfv2tf_H
//
// This is the declaration of the MEfv2tf class.
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
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;

/**
 * Here is the documentation of the MEfv2tf class.
 *
 * @see \ref MEfv2tfInterfaces "The interfaces"
 * defined for MEfv2tf.
 */
class MEfv2tf: public GeneralHardME {

  /** Vector of SpinorWaveFunctions. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunctions. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

  /** Vector of VectorWaveFunctions. */
  typedef vector<VectorWaveFunction> VBVector;

  /** Vector of TensorWaveFunctions. */
  typedef vector<TensorWaveFunction> TBVector;

public:

  /**
   * The default constructor.
   */
  MEfv2tf() : fermion_(0), vector_(0), fourPoint_(0) {}

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

  /** @name Functions to calculate production matrix elements and me2. */
  //@{
  /**
   * Calculate me2 and the production matrix element for the normal mode.
   * @param spIn Vector of SpinorWaveFunction for the incoming fermion
   * @param vecIn Vector of VectorWaveFunction for incoming boson
   * @param spbOut Vector of SpinorBarWaveFunction for outgoing fermion
   * @param tenOut TensorWaveFunction for outgoing tensor.
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement fv2tfHeME(const SpinorVector & spIn, 
				    const VBVector & vecIn,
				    const TBVector & tenOut,
				    const SpinorBarVector & spbOut,
				    double & full_me, bool first) const;

  /**
   * Calculate me2 and the production matrix element for the cc mode.
   * @param spbIn Vector of SpinorBarWaveFunction for the incoming fermion
   * @param vecIn Vector of VectorWaveFunction for incoming boson
   * @param spOut Vector of SpinorWaveFunction for outgoing fermion
   * @param tenOut TensorWaveFunction for outgoing tensor.
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement fbv2tfbHeME(const SpinorBarVector & spbIn, 
				      const VBVector & vecIn,
				      const TBVector & tenOut,
				      const SpinorVector & spOut,
				      double & full_me, bool first) const;
  //@}

  /**
   * A debugging function to test the value of me2 against an
   * analytic function.
   * @param me2 The value of the \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  virtual void debug(double me2) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEfv2tf & operator=(const MEfv2tf &) = delete;

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

#endif /* HERWIG_MEfv2tf_H */
