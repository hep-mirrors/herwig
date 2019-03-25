// -*- C++ -*-
#ifndef Herwig_MEfv2rs_H
#define Herwig_MEfv2rs_H
//
// This is the declaration of the MEfv2rs class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * This class is designed to implement the matrix element for 
 * fermion-vector to RS fermion scalar. It inherits from GeneralHardME 
 * and implements the required virtual functions.
 *
 * @see GeneralHardME
 */
class MEfv2rs: public GeneralHardME {

  /** Vector of SpinorWaveFunctions. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunctions. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

  /** Vector of RSSpinorWaveFunctions. */
  typedef vector<RSSpinorWaveFunction> RSSpinorVector;

  /** Vector of RSSpinorBarWaveFunctions. */
  typedef vector<RSSpinorBarWaveFunction> RSSpinorBarVector;

  /** Vector of VectorWaveFunctions. */
  typedef vector<VectorWaveFunction> VecWFVector;

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
   * @param subp Pointer to the relevent SubProcess
   */
  virtual void constructVertex(tSubProPtr subp);

private:

  /** @name Functions to calculate production matrix elements and me2. */
  //@{
  /**
   * Calculate me2 and the production matrix element for the normal mode.
   * @param spIn Vector of SpinorWaveFunction for the incoming fermion
   * @param vecIn Vector of VectorWaveFunction for incoming boson
   * @param spbOut Vector of SpinorBarWaveFunction for outgoing fermion
   * @param scaOut ScalarWaveFunction for outgoing scalar.
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement fv2rbsHeME(const SpinorVector & spIn, 
				     const VecWFVector & vecIn,
				     const RSSpinorBarVector & spbOut,
				     const ScalarWaveFunction & scaOut,
				     double & full_me, bool first) const;
  
  /**
   * Calculate me2 and the production matrix element for the cc mode.
   * @param spbIn Vector of SpinorBarWaveFunction for the incoming fermion
   * @param vecIn Vector of VectorWaveFunction for incoming boson
   * @param spOut Vector of SpinorWaveFunction for outgoing fermion
   * @param scaOut ScalarWaveFunction for outgoing scalar.
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement fbv2rsHeME(const SpinorBarVector & spbIn, 
				     const VecWFVector & vecIn,
				     const RSSpinorVector & spOut,
				     const ScalarWaveFunction & scaOut,
				     double & full_me, bool first) const;
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
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEfv2rs & operator=(const MEfv2rs &) = delete;

private:

  /**
   * Store a pair of  FFSVertex and VSSVertex pointers  
   */
  vector<pair<AbstractRFSVertexPtr, AbstractVSSVertexPtr> > scalar_;

  /**
   * Store a pair of  FFSVertex and FFVVertex pointers  
   */
  vector<pair<AbstractFFVVertexPtr, AbstractRFSVertexPtr> > fermion1_;

  /**
   * Store a pair of  FFSVertex and FFVVertex pointers  
   */
  vector<pair<AbstractFFSVertexPtr, AbstractRFVVertexPtr> > fermion2_;

  /**
   * Store a pair of  VVSVertex and FFVVertex pointers  
   */
  vector<pair<AbstractRFVVertexPtr,AbstractVVSVertexPtr> > vector_;

  /**
   *  Store the 4-point vertices
   */
  vector<AbstractRFVSVertexPtr> four_;
};

}

#endif /* Herwig_MEfv2rs_H */
