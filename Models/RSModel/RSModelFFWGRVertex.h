// -*- C++ -*-
//
// RSModelFFWGRVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RSModelFFWGRVertex_H
#define HERWIG_RSModelFFWGRVertex_H
//
// This is the declaration of the RSModelFFWGRVertex class.

#include "ThePEG/Helicity/Vertex/Tensor/FFVTVertex.h"
#include "RSModel.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Helicity
 *
 *  This is the implementation of the fermion-antifermion-vector-graviton
 *  vertex for the Randell-Sundrum model
 *
 *  @see FFVTVertex
 *  @see VertexBase
 */
class RSModelFFWGRVertex: public FFVTVertex {
  
public:
  
  /**
   * Default constructor.
   */
  RSModelFFWGRVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the foruth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
			   tcPDPtr part3, tcPDPtr part4);
  
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
   * Standard Init function used to initialize the interfaces.
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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
    
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<RSModelFFWGRVertex> initRSModelFFWGRVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RSModelFFWGRVertex & operator=(const RSModelFFWGRVertex &);

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The charges of the Standard Model fermions.
   */
  vector<double> charge_;

  /**
   *  The left couplings of the Standard Model fermions.
   */
  vector<double> gl_;

  /**
   *  The right couplings of the Standard Model fermions.
   */
  vector<double> gr_;

  /**
   *  The elements of the CKM matrix.
   */
  vector<vector<Complex> > ckm_;

  /**
   *  The last value of the coupling/
   */
  Complex couplast_;

  /**
   *  The last value of the scale, \f$q^2\f$.
   */
  Energy2 q2last_;

  /**
   * The graviton coupling.
   */
  InvEnergy kappa_;
  //@}
};
}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of RSModelFFWGRVertex.
 */
template <>
struct BaseClassTrait<Herwig::RSModelFFWGRVertex,1> {
    /** Typedef of the base class of RSModelFFWGRVertex. */
  typedef ThePEG::Helicity::FFVTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::RSModelFFWGRVertex>
  : public ClassTraitsBase<Herwig::RSModelFFWGRVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "Herwig::RSModelFFWGRVertex"; }

  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwRSModel.so"; }

};

/** @endcond */

}


#endif /* HERWIG_RSModelFFWGRVertex_H */
