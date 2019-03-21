// -*- C++ -*-
#ifndef Herwig_MEvv2rf_H
#define Herwig_MEvv2rf_H
//
// This is the declaration of the MEvv2rf class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/Vertex/AbstractRFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractRFVVVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * This class is designed to implement the matrix element for the 
 * \f$2 \rightarrow 2\f$ process vector-vector to fermion- RS fermion. It
 * inherits from GeneralHardME and implements the me2() virtual function.
 *
 * @see GeneralHardME
 * 
 */
class MEvv2rf: public GeneralHardME {

public:
  
  /** A Vector of VectorWaveFunction objects. */
  typedef vector<VectorWaveFunction> VBVector;

  /** A vector of SpinorBarWaveFunction objects. */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** A vector of SpinorBarWaveFunction objects. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;
  
  /** A vector of SpinorBarWaveFunction objects. */
  typedef vector<RSSpinorWaveFunction> RSSpinorVector;

  /** A vector of SpinorBarWaveFunction objects. */
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
   * Calculate the value of the matrix element 
   */
  ProductionMatrixElement vv2rfME(const VBVector & v1, const VBVector & v2,
				  const RSSpinorBarVector & sbar,
				  const SpinorVector & sp, 
				  double & me2, bool first) const;
  
  /**
   * Calculate the value of the matrix element 
   */
  ProductionMatrixElement vv2frME(const VBVector & v1, const VBVector & v2,
				  const SpinorBarVector & sbar,
				  const RSSpinorVector & sp, 
				  double & me2, bool first) const;
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2rf & operator=(const MEvv2rf &) = delete;
private:
  
  /** @name Dynamically casted vertices. */
  //@{
  /**
   *  Intermediate scalar
   */
  vector<pair<AbstractVVSVertexPtr, AbstractRFSVertexPtr > > scalar_;
  
  /**
   * Intermediate fermion 
   */
  vector<pair<AbstractRFVVertexPtr, AbstractFFVVertexPtr> > fermion_;

  /**
   * Intermediate vector
   */
  vector<pair<AbstractVVVVertexPtr, AbstractRFVVertexPtr> > vector_;

  /**
   *  Four point vertices
   */
  vector<AbstractRFVVVertexPtr> four_;
  //@}

};

}

#endif /* Herwig_MEvv2rf_H */
