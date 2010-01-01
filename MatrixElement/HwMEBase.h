// -*- C++ -*-
#ifndef HERWIG_HwMEBase_H
#define HERWIG_HwMEBase_H
//
// This is the declaration of the HwMEBase class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig++/Shower/Base/ShowerKinematics.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/HardTree.h"
#include "HwMEBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HwMEBase class.
 *
 * @see \ref HwMEBaseInterfaces "The interfaces"
 * defined for HwMEBase.
 */
class HwMEBase: public MEBase {

public:

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual bool hasPOWHEGCorrection() {return false;}

  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return false;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(ShowerTreePtr , double & ,
				      double & ) {}

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr) {}

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param parent The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr,
				     ShowerParticlePtr,Branching) {
    return false;
  }

  /**
   *  Apply the POWHEG style correction
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr) {
    return HardTreePtr();
  }
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<HwMEBase> initHwMEBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HwMEBase & operator=(const HwMEBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HwMEBase. */
template <>
struct BaseClassTrait<Herwig::HwMEBase,1> {
  /** Typedef of the first base class of HwMEBase. */
  typedef MEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HwMEBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HwMEBase>
  : public ClassTraitsBase<Herwig::HwMEBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HwMEBase"; }
};

/** @endcond */

}

#endif /* HERWIG_HwMEBase_H */
