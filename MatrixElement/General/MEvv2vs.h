// -*- C++ -*-
#ifndef HERWIG_MEvv2vs_H
#define HERWIG_MEvv2vs_H
//
// This is the declaration of the MEvv2vs class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVSVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;
using Helicity::VectorWaveFunction;
using Helicity::ScalarWaveFunction;

/**
 * This is the implementation of the matrix element for 
 * \f$2\to 2\f$ massless vector-boson pair to a vector and scalar boson.
 * It inherits from GeneralHardME and implements the appropriate virtual
 * member functions.
 *
 * @see \ref MEvv2vsInterfaces "The interfaces"
 * defined for MEvv2vs.
 */
class MEvv2vs: public GeneralHardME {

public:

  /**
   *  Typedef for VectorWaveFunction
   */
  typedef vector<VectorWaveFunction> VBVector;

public:

  /** @name Virtual functions required by the GeneralHardME class. */
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

  /**
   * Compute the matrix element for \f$V\, V\to V\, V\f$
   * @param vin1 VectorWaveFunctions for first incoming particle
   * @param vin2 VectorWaveFunctions for second incoming particle
   * @param vout1 VectorWaveFunctions for first outgoing particle
   * @param mc Whether vout1 is massless or not
   * @param sout2  ScalarWaveFunction for outgoing particle
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement 
  vv2vsHeME(VBVector & vin1, VBVector & vin2, 
	    VBVector & vout1, bool mc, ScalarWaveFunction & sout2,
	    double & me2, bool first ) const;

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
  MEvv2vs & operator=(const MEvv2vs &) = delete;

private:

  /**
   * Store the dynamically casted VVSVertex and VSSVertex pointers
   */
  vector<pair<AbstractVVSVertexPtr, AbstractVSSVertexPtr> > scalar_;

  /**
   * Store the dynamically casted VVVVertex and VVSVertex pointers
   */
  vector<pair<AbstractVVVVertexPtr, AbstractVVSVertexPtr> > vector_;

  /**
   * Store the dynamically casted VVVSVertex pointer
   */
  vector<AbstractVVVSVertexPtr> four_;

};

}

#endif /* HERWIG_MEvv2vs_H */
