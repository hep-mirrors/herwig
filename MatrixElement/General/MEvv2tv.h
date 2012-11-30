// -*- C++ -*-
#ifndef HERWIG_MEvv2tv_H
#define HERWIG_MEvv2tv_H
//
// This is the declaration of the MEvv2tv class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVTVertex.h"

namespace Herwig {

using namespace ThePEG;
using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::TensorWaveFunction;

/**
 * Here is the documentation of the MEvv2tv class.
 *
 * @see \ref MEvv2tvInterfaces "The interfaces"
 * defined for MEvv2tv.
 */
class MEvv2tv: public GeneralHardME {

  /** Vector of VectorWaveFunctions. */
  typedef vector<VectorWaveFunction> VBVector;

  /** Vector of TensorWaveFunctions. */
  typedef vector<TensorWaveFunction> TBVector;

public:

  /**
   * The default constructor.
   */
  MEvv2tv() : vector_(0), fourPoint_(0) {}

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
   * @param vec1 Vector of VectorWaveFunction for the 1st incoming boson
   * @param vec2 Vector of VectorWaveFunction for the 2nd incoming boson
   * @param ten TensorWaveFunction for outgoing tensor.
   * @param vec3 Vector of VectorWaveFunction for the outgoing boson
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param full_me The value of me2 calculation
   */
  ProductionMatrixElement vv2tvHeME(const VBVector & vec1,
				    const VBVector & vec2,
				    const TBVector & ten,
				    const VBVector & vec3,
				    double & full_me, bool first) const;

  /**
   * A debugging function to test the value of me2 against an
   * analytic function.
   * @param me2 The value of the \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  virtual void debug(double me2) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEvv2tv> initMEvv2tv;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2tv & operator=(const MEvv2tv &);

private:

  /**
   * Store a pair of  FFTVertex and FFVVertex pointers  
   */
  vector<pair<AbstractVVVVertexPtr, AbstractVVTVertexPtr> > vector_;

  /**
   *  The four point vertex
   */
  vector<AbstractVVVTVertexPtr> fourPoint_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2tv. */
template <>
struct BaseClassTrait<Herwig::MEvv2tv,1> {
  /** Typedef of the first base class of MEvv2tv. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2tv class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2tv>
  : public ClassTraitsBase<Herwig::MEvv2tv> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEvv2tv"; }
};

/** @endcond */

}

#endif /* HERWIG_MEvv2tv_H */
