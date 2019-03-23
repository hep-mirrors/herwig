// -*- C++ -*-
#ifndef Herwig_MEfv2rv_H
#define Herwig_MEfv2rv_H
//
// This is the declaration of the MEfv2rv class.
//

#include "GeneralHardME.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * This class implements the matrix element for a fermion and a vector
 * boson to a vector boson and a RS fermion. It inherits from GeneralHardME
 * and implements the appropriate virtual functions. 
 * 
 * @see GeneralHardME
 *
 */
class MEfv2rv: public GeneralHardME {

public:

  /** A vector of SpinorWaveFunctions. */
  typedef vector<Helicity::SpinorWaveFunction> SpinorVector;

  /** A vector of SpinorBarWaveFunctions. */
  typedef vector<Helicity::SpinorBarWaveFunction> SpinorBarVector;

  /** A vector of SpinorWaveFunctions. */
  typedef vector<Helicity::RSSpinorWaveFunction> RSSpinorVector;

  /** A vector of SpinorBarWaveFunctions. */
  typedef vector<Helicity::RSSpinorBarWaveFunction> RSSpinorBarVector;

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
  fv2rvHeME(const SpinorVector & spIn,  const VBVector & vecIn, 
	    const RSSpinorBarVector & spbOut, 
	    const VBVector & vecOut, bool mc,
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
  fbv2rbvHeME(const SpinorBarVector & spbIn,  const VBVector & vecIn,
	      const RSSpinorVector & spOut,
	      const VBVector & vecOut, bool mc, 
	      double & mesq, bool first) const;
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
  MEfv2rv & operator=(const MEfv2rv &) = delete;

private:
  
  /** @name Store dynamically casted vertices. */
  //@{
  /**
   * A pair off FFVVertex pointers 
   */
  vector<pair<AbstractFFVVertexPtr, AbstractRFVVertexPtr> > fermion_;

  /**
   * A pair of FFVVertex, VVVertex pointers 
   */
  vector<pair<AbstractRFVVertexPtr, AbstractVVVVertexPtr> > vector_;

  /**
   *  Four point vertices
   */
  vector<AbstractRFVVVertexPtr> four_;
  //@}

};

}

#endif /* Herwig_MEfv2rv_H */
