// -*- C++ -*-
#ifndef Herwig_MEPP2OniumPowheg_H
#define Herwig_MEPP2OniumPowheg_H
//
// This is the declaration of the MEPP2OniumPowheg class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "OniumParameters.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2OniumPowheg class is designed to allow easy implementation of the matrix elements for quarkonium production inthe POWHEG scheme
 *
 * @see \ref MEPP2OniumPowhegInterfaces "The interfaces"
 * defined for MEPP2OniumPowheg.
 */
class MEPP2OniumPowheg: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2OniumPowheg();

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const {
    return 2;
  }

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const {
    return 0;
  }

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

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
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const {
    return 2;
  }

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;
  
  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}

public:
 
  // /** @name Member functions for the generation of hard QCD radiation */
  // //@{
  // /**
  //  *  Has a POWHEG style correction
  //  */
  // virtual POWHEGType hasPOWHEGCorrection() {return ISR;}

  // /**
  //  *  Apply the POWHEG style correction
  //  */
  // virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr,
  // 						 ShowerInteraction);
  // //@}

protected:

  /**
   *  The leading-order matrix element
   */
  virtual Energy2 leadingOrderME2() const = 0;

  /**
   *  set the state
   */
  void setParticleData(long ioff) {
    unsigned int iq = 4+state_;
    onium_ = getParticleData(long(iq*110+ioff + (n_-1)*100000));
  }

  /**
   *  Access to the quarkonium parameters
   */
  const OniumParametersPtr parameters() const {
    return params_;
  }

  /**
   *  Type of quarkonium state
   */
  OniumState state() const {
    return state_;
  }

  /**
   *  Principal quantum number
   */
  unsigned int principalQuantumNumber() const {
    return n_;
  }
  
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
  MEPP2OniumPowheg & operator=(const MEPP2OniumPowheg &) = delete;

private:

  /**
   *   Parameters for the form-factors
   */
  //@{
  /**
   *  Access to the parameters for the quarkonium states
   */
  OniumParametersPtr params_;

  /**
   *  Type of state
   */
  OniumState state_;

  /**
   *  Principal quantum number
   */
  unsigned int n_;

  /**
   *  The particle data object for the state being produced
   */
  PDPtr onium_;

  /**
   *  The mass generator for the onium state
   */
  GenericMassGeneratorPtr massGen_;
  //@}

private:

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Scales etc
   */
  //@{
  /**
   * Selects a dynamic (sHat) or fixed factorization scale
   */
  unsigned int scaleopt_;

  /**
   * The factorization  scale 
   */
  Energy mu_F_;

  /**
   * The renormalization scale
   */
  Energy mu_UV_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;
  //@}

  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$\tilde{x}\f$ variable
   */
  double xt_;

  /**
   *  The \f$y\f$ angular variable
   */
  double y_;
  //@}
};

}

#endif /* Herwig_MEPP2OniumPowheg_H */
