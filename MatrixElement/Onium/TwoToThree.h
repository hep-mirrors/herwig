// -*- C++ -*-
#ifndef Herwig_TwoToThree_H
#define Herwig_TwoToThree_H
//
// This is the declaration of the TwoToThree class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The TwoToThree class implements simple \f$2\to3\f$ phase space
 * for the three-body \f$B_c\f$ and diquark prodcution processes
 *
 * @see \ref TwoToThreeInterfaces "The interfaces"
 * defined for TwoToThree.
 */
class TwoToThree: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  TwoToThree() : mOpt_(1)
  {}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the scale associated with the last set phase space point.
   */
  Energy2 scale() const {
    return sHat();
  }

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  void setKinematics() {
    HwMEBase::setKinematics();
  }

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  int nDim() const {
    if(mOpt_==1 && massGen_) return 5;
    else return 4;
  }

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  CrossSection dSigHatDR() const;
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

  /**
   * Set the mass generatror
   */
  void setMassGenerator(GenericMassGeneratorPtr in) {
    massGen_=in;
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TwoToThree & operator=(const TwoToThree &) = delete;

private:
  
  /**
   *  Option for the mass handling
   */
  unsigned int mOpt_;

  /**
   *  Mass generator for the state
   */
  GenericMassGeneratorPtr massGen_;

};

}

#endif /* Herwig_TwoToThree_H */
