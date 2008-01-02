// -*- C++ -*-
#ifndef HERWIG_SimpleVVSLoopVertex_H
#define HERWIG_SimpleVVSLoopVertex_H
//
// This is the declaration of the SimpleVVSLoopVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/GeneralVVSVertex.h"
#include "SimpleVVSLoopVertex.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The <code>SimpleVVSLoopVertex</code> is designed to
 * calculate the coefficents for the terms in the
 * Passarino-Veltman tensor reduction scheme. A vertex
 * class should inherit from this and implement it's own 
 * setCoupling member from which the SimpleVVSLoopVertex  
 * setCoupling member is called.
 */
class SimpleVVSLoopVertex: public GeneralVVSVertex {

public:

  /**
   * The default constructor.
   */
  inline SimpleVVSLoopVertex();

  /** 
   * Calculate couplings
   *@param q2 Scale at which to evaluate coupling
   *@param part1 ParticleData pointer to first particle
   *@param part2 ParticleData pointer to first particle
   *@param part3 ParticleData pointer to first particle
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1, tcPDPtr part2, tcPDPtr part3);

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:
  
  /**
   *  Member functions to calculate the loop functions
   */
  //@{
  /**
   *  The \f$A_1(s)\f$ function of NPB297 (1988) 221-243.
   * @param s   The invariant
   * @param mf2 The fermion mass squared
   */
  Complex A1(Energy2 s, Energy2 mf2) const;

  /**
   *  The \f$W_2(s)\f$ function of NPB297 (1988) 221-243.
   * @param s   The invariant
   * @param mf2 The fermion mass squared
   */
  Complex W2(Energy2 s, Energy2 mf2) const;

  /**
   *  The \f$I_q\f$ function of V.D.Berger and R.N.Phillips Collider Physics, p .434
   * N.B. is not used in the code...
   * @param mh The Higgs mass (invariant)
   * @param mf The fermion mass
   */
   Complex W3(Energy mh, Energy mf) const;
  //@}

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
   * The left and right couplings for a fermion loop 
   */
  vector<pair<Complex, Complex> > couplings;

  /**
   * The number of particles in the loop 
   */
  unsigned int theNpart;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<SimpleVVSLoopVertex> initSimpleVVSLoopVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleVVSLoopVertex & operator=(const SimpleVVSLoopVertex &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SimpleVVSLoopVertex. */
template <>
struct BaseClassTrait<Herwig::SimpleVVSLoopVertex,1> {
  /** Typedef of the first base class of SimpleVVSLoopVertex. */
  typedef Helicity::GeneralVVSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SimpleVVSLoopVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SimpleVVSLoopVertex>
  : public ClassTraitsBase<Herwig::SimpleVVSLoopVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SimpleVVSLoopVertex"; }
};

/** @endcond */

}

#include "SimpleVVSLoopVertex.icc"

#endif /* HERWIG_SimpleVVSLoopVertex_H */
