// -*- C++ -*-
#ifndef Herwig_MEff2rf_H
#define Herwig_MEff2rf_H
//
// This is the declaration of the MEff2rf class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVVertex.h"

namespace Herwig {

using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using Helicity::RSSpinorWaveFunction;
using Helicity::RSSpinorBarWaveFunction;

/**
 * This is the implementation of the \f$ 2\to 2\f$ matrix element for
 * a fermion-fermion to fermion RS fermion process. It inherits from 
 * GeneralHardME and implements the appropriate virtual functions.
 *
 * @see GeneralHardME
 */
class MEff2rf: public GeneralHardME {

public:
  
  /** Vector of SpinorWaveFunctions. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunctions. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;
  
  /** Vector of RSSpinorWaveFunctions. */
  typedef vector<RSSpinorWaveFunction> RSSpinorVector;

  /** Vector of RSSpinorBarWaveFunctions. */
  typedef vector<RSSpinorBarWaveFunction> RSSpinorBarVector;

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
  ffb2rfbHeME(double & me2, bool first) const;
  
  /**
   * Compute the matrix element for \f$\Psi\bar{\Psi}\to\Psi\bar{\Psi}\f$
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement
  ffb2rbfHeME(double & me2, bool first) const;

  /**
   * Compute the matrix element for \f$\Psi\Psi\to\Psi\Psi\f$
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement ff2rfHeME(double & me2, bool first) const;
  
  /**
   * Compute the matrix element for 
   * \f$\bar{\Psi}\bar{\Psi}\to\bar{\Psi}\bar{\Psi}\f$
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement fbfb2rbfbHeME(double & me2, bool first) const;

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
  MEff2rf & operator=(const MEff2rf &) = delete;

private:
  
  /**
   * Store the vector of FFSVertex pairs
   */
  vector<pair<AbstractFFSVertexPtr, AbstractRFSVertexPtr> > scalar_;

  /**
   * Store the vector of FFVVertex pairs
   */
  vector<pair<AbstractFFVVertexPtr, AbstractRFVVertexPtr> > vector_;

  /**
   *  Spinors
   */
  mutable std::array<SpinorVector,3> spin_;

  /**
   *  Barred spinors
   */
  mutable std::array<SpinorBarVector,3> sbar_;

  mutable RSSpinorVector rs_;
  mutable RSSpinorBarVector rsbar_;
};

}

#endif /* Herwig_MEff2rf_H */
