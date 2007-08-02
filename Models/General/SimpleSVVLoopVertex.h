// -*- C++ -*-
#ifndef HERWIG_SimpleSVVLoopVertex_H
#define HERWIG_SimpleSVVLoopVertex_H
//
// This is the declaration of the SimpleSVVLoopVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/GeneralSVVVertex.h"
#include "ThePEG/PDT/PDT.h"
#include "SimpleSVVLoopVertex.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The <code>SimpleSVVLoopVertex</code> is designed to
 * calculate the coefficents for the terms in the
 * Passarino-Veltman tensor reduction scheme. A vertex
 * class should inherit from this and implement it's own 
 * setCoupling member from which the SimpleSVVLoopVertex  
 * setCoupling member is called.
 */
class SimpleSVVLoopVertex: public GeneralSVVVertex {

public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SimpleSVVLoopVertex();
  
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
  
  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2, tcPDPtr part3);
  
protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;
  
  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}
  
  /**
   *  Member functions to calculate the loop functions
   */
  //@{
  /**
   *  The \f$A_1(s)\f$ function of NPB297 (1988) 221-243.
   * @param s   The invariant
   * @param mf2 The fermion mass squared
   */
  Complex A1(Energy2 s,Energy2 mf2) const;

  /**
   *  The \f$W_2(s)\f$ function of NPB297 (1988) 221-243.
   * @param s   The invariant
   * @param mf2 The fermion mass squared
   */
  Complex W2(Energy2 s,Energy2 mf2) const;
  //@}
  Complex W3(Energy mh, Energy mf) const;

protected:
  
  /**
   * Vector of loop masses
   */
  vector<Energy> masses;

  
  /**
   * Vector of loop types
   */
  vector<PDT::Spin> type;
  
  /**
   * Vector of left couplings
   */
  vector<double> left;
  
  /**
   * Vector of right couplings
   */
  vector<double> right;
  
  /**
   * The number of particles in the loop 
   */
  unsigned int theNpart;

private:
  
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<SimpleSVVLoopVertex> initSimpleSVVLoopVertex;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleSVVLoopVertex & operator=(const SimpleSVVLoopVertex &);
  
};

}

#include "SimpleSVVLoopVertex.icc"

// CLASSDOC OFF

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SimpleSVVLoopVertex. */
template <>
struct BaseClassTrait<Herwig::SimpleSVVLoopVertex,1> {
  /** Typedef of the first base class of SVVLoopVertex. */
  typedef Helicity::GeneralSVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SVVLoopVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SimpleSVVLoopVertex>
  : public ClassTraitsBase<Herwig::SimpleSVVLoopVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SimpleSVVLoopVertex"; }
};

}
 
#endif /* HERWIG_SimpleSVVLoopVertex_H */
