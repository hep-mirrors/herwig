// -*- C++ -*-
//
// MEff2vv.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEff2vv_H
#define HERWIG_MEff2vv_H
//
// This is the declaration of the MEff2vv class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::VectorWaveFunction;

/**
 * This class implements the matrix element calculation for a generic
 * \f$\Psi \Psi \rightarrow V^{\mu} V^{\nu}\f$ process. 
 *
 * @see \ref MEff2vvInterfaces "The interfaces"
 * defined for MEff2vv.
 */
class MEff2vv: public GeneralHardME {
public:
  
  /** Vector of SpinorWaveFunctions objects */
  typedef vector<SpinorWaveFunction> SpinorVector;

  /** Vector of SpinorBarWaveFunction objects. */
  typedef vector<SpinorBarWaveFunction> SpinorBarVector;

  /** Vector of VectorWaveFunction objects. */
  typedef vector<VectorWaveFunction> VBVector;

public:

  /**
   * The default constructor.
   */
  MEff2vv() : vector_(0), tensor_(0), scalar_(0) {}

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


protected:
  
  /**
   * A debugging function to test the value of me2 against an
   * analytic function.
   * @param me2 \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  virtual void debug(double me2) const;

private: 
  
  /**
   * Compute the production matrix element.
   * @param sp Spinors for first incoming fermion
   * @param sbar SpinorBar Wavefunctions for incoming anti-fermion
   * @param v1 A vector of VectorWaveFunction objects for the first vector
   * @param m1 Whether v1 is massless or not
   * @param v2 A vector of VectorWaveFunction objects for the second vector
   * @param m2 Whether v2 is massless or not
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @param me2 The value of the \f$ |\bar{\mathcal{M}}|^2 \f$
   */
  ProductionMatrixElement 
  ff2vvME(const SpinorVector & sp, const SpinorBarVector sbar, 
	  const VBVector & v1, bool m1, const VBVector & v2, bool m2,
	  double & me2, bool first) const;

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEff2vv> initMEff2vv;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEff2vv & operator=(const MEff2vv &);

private:

  /**
   * Storage for a dynamically cast vertices for a tchannel vector
   * intermediate
   */
  vector<pair<AbstractFFVVertexPtr, AbstractFFVVertexPtr> > fermion_;

  /**
   * Storage for a dynamically cast vertices for a schannel vector
   * intermediate
   */
  vector<pair<AbstractFFVVertexPtr, AbstractVVVVertexPtr> > vector_;

  /**
   * Storage for a dynamically cast vertices for a schannel scalar
   * intermediate
   */
  vector<pair<AbstractFFTVertexPtr, AbstractVVTVertexPtr> > tensor_;

  /**
   * Storage for a dynamically cast vertices for a schannel scalar
   * intermediate for massless external vector bosons
   */
  vector<pair<AbstractFFSVertexPtr, AbstractVVSVertexPtr> > scalar_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEff2vv. */
template <>
struct BaseClassTrait<Herwig::MEff2vv,1> {
  /** Typedef of the first base class of MEff2vv. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEff2vv class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEff2vv>
  : public ClassTraitsBase<Herwig::MEff2vv> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEff2vv"; }
};

/** @endcond */

}

#endif /* HERWIG_MEff2vv_H */
