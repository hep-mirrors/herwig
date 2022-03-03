// -*- C++ -*-
#ifndef Herwig_MassiveIncoming_H
#define Herwig_MassiveIncoming_H
//
// This is the declaration of the MassiveIncoming class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MassiveIncoming class provides the implementation fo the kinematics for massive incoming particles.
 *
 * @see \ref MassiveIncomingInterfaces "The interfaces"
 * defined for MassiveIncoming.
 */
class MassiveIncoming: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MassiveIncoming() : mOpt_(1)
  {}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;
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
  MassiveIncoming & operator=(const MassiveIncoming &) = delete;

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

#endif /* Herwig_MassiveIncoming_H */
