// -*- C++ -*-
//
// MEvv2vv.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEvv2vv_H
#define HERWIG_MEvv2vv_H
//
// This is the declaration of the MEvv2vv class.
//

#include "GeneralHardME.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VectorWaveFunction;

/**
 * This is the implementation of the matrix element for 
 * \f$2\to 2\f$ massless vector-boson pair to vector-boson pair. It inherits from
 * GeneralHardME and implements the appropriate virtual member functions.
 *
 * @see \ref MEvv2vvInterfaces "The interfaces"
 * defined for MEvv2vv.
 */
class MEvv2vv: public GeneralHardME {

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
   * @param vout2  VectorWaveFunctions for outgoing particle
   * @param md Whether vout2 is massless or not
   * @param me2 colour averaged, spin summed ME
   * @param first Whether or not first call to decide if colour decomposition etc
   * should be calculated
   * @return ProductionMatrixElement containing results of 
   * helicity calculations
   */
  ProductionMatrixElement 
  vv2vvHeME(VBVector & vin1, VBVector & vin2, 
	    VBVector & vout1, bool mc, VBVector & vout2, bool md,
	    double & me2, bool first ) const;

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
  static ClassDescription<MEvv2vv> initMEvv2vv;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEvv2vv & operator=(const MEvv2vv &);

private:

  /**
   * Store the dynamically casted VVSVertex pointers
   */
  vector<pair<AbstractVVSVertexPtr, AbstractVVSVertexPtr> > scalar_;

  /**
   * Store the dynamically casted VVVVertex pointers
   */
  vector<pair<AbstractVVVVertexPtr, AbstractVVVVertexPtr> > vector_;

  /**
   * Store the dynamically casted VVTVertex pointers
   */
  vector<pair<AbstractVVTVertexPtr, AbstractVVTVertexPtr> > tensor_;

  /**
   * Store the dynamically casted VVVVVertex pointer
   */
  AbstractVVVVVertexPtr fourPointVertex_;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEvv2vv. */
template <>
struct BaseClassTrait<Herwig::MEvv2vv,1> {
  /** Typedef of the first base class of MEvv2vv. */
  typedef Herwig::GeneralHardME NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEvv2vv class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEvv2vv>
  : public ClassTraitsBase<Herwig::MEvv2vv> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEvv2vv"; }
};

/** @endcond */

}

#endif /* HERWIG_MEvv2vv_H */
