// -*- C++ -*-
#ifndef HERWIG_HwMEBase_H
#define HERWIG_HwMEBase_H
//
// This is the declaration of the HwMEBase class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig/Shower/Base/ShowerParticle.fh"
#include "Herwig/Shower/Base/ShowerProgenitor.fh"
#include "Herwig/Shower/Base/ShowerTree.fh"
#include "Herwig/Shower/Base/HardTree.fh"
#include "Herwig/Shower/ShowerConfig.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "HwMEBase.fh"

namespace Herwig {

struct Branching;

using namespace ThePEG;

typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

/**
 * The HwMEBase class serves a number of purposes
 * - it implements the phase space for \f$2\to2\f$ scattering processes
 * - it provides virtual members for the implementation of hard radiation
 * - it gives us greater control over the masses of the outgoing
 *   particles so that they can be
 *   - set massless where required by gauge invariance
 *   - have their off-shell masses generated using the sophisticated approaches
 *     available in Herwig.
 *
 * @see \ref HwMEBaseInterfaces "The interfaces"
 * defined for HwMEBase.
 */
class HwMEBase: public MEBase {

public:

  /**
   * Default constructor.
   */
  HwMEBase() : lastTHat_(ZERO), lastUHat_(ZERO),
	       lastPhi_(0.0), rescaleOption_(1)
  {}

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given 'nDim()' uniform
   * random numbers in the interval ]0,1[. To help the phase space
   * generator, the 'dSigHatDR()' should be a smooth function of these
   * numbers, although this is not strictly necessary. Return
   * false if the chosen points failed the kinematical cuts.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object.
   */
  virtual void setKinematics();
  //@}

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{

  /**
   * Type of POWHEG correction
   */
  enum POWHEGType {No, ISR, FSR, Both};

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType  hasPOWHEGCorrection() {return No;}

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
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,
				     Branching br);

  /**
   *  Apply the POWHEG style correction
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr,vector<ShowerInteraction::Type>);
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

  /** @name Access cached values in of the last set phase space point. */
  //@{
  /**
   * Return the \f$\hat{t}\f$ of the last set phase space point.
   */
  Energy2 tHat() const { return lastTHat_; }

  /**
   * Return the \f$\hat{u}\f$ of the last set phase space point.
   */
  Energy2 uHat() const { return lastUHat_; }

  /**
   * Return the azimuth angle of the last set phase space point.
   */
  double phi() const { return lastPhi_; }
  //@}

  /** @name Set the cached values in of the last set phase space point. */
  //@{
  /**
   * Set the \f$\hat{t}\f$ of the last set phase space point.
   */
  void tHat(Energy2 e2) { lastTHat_ = e2; }

  /**
   * Set the \f$\hat{u}\f$ of the last set phase space point.
   */
  void uHat(Energy2 e2) { lastUHat_ = e2; }

  /**
   * Set the azimuth angle of the last set phase space point.
   */
  void phi(double phi) { lastPhi_ = phi; }
  //@}

  /**
   *  Set the treatment of the outgoing masses
   * @param iopt The option for the treatment of the mass
   */
  void massOption(vector<unsigned int> iopt) {
    massOption_ = iopt;
  }

  /**
   *  Rescaled momenta for the helicity ME
   */
  //@{
  /**
   * Set the treatment of the rescaling of the momenta for 
   * the matrix element calculation
   * @param iopt The rescaling option
   */
  void rescalingOption(unsigned int iopt) {
    rescaleOption_=iopt;
  }

  /**
   *  rescale the momenta for the computation of the helicity matrix element
   */
  bool rescaleMomenta(const vector<Lorentz5Momentum> &,
		      const cPDVector &);

  /**
   *  Access to the rescaled momenta
   */
  const vector<Lorentz5Momentum> & rescaledMomenta() const {
    return rescaledMomenta_;
  }
  //@}

  /**
   *  Generate the masses of the particles
   */ 
  bool generateMasses(vector<Energy> & masses, double & mjac,
		      const double *r);

  /**
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  virtual double getCosTheta(double cthmin, double cthmax, const double r);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractClassDescription<HwMEBase> initHwMEBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HwMEBase & operator=(const HwMEBase &);

private:

  /**
   *  Option for the treatment of the particle masses
   */
  vector<unsigned int> massOption_;
  
  /**
   * The \f$\hat{t}\f$ of the last set phase space point.
   */
  Energy2 lastTHat_;
  
  /**
   * The \f$\hat{u}\f$ of the last set phase space point.
   */
  Energy2 lastUHat_;

  /**
   * The azimuth angle of the last set phase space point.
   */
  double lastPhi_;

  /**
   *  Produced to produce rescaled momenta
   */
  unsigned int rescaleOption_;

  /**
   *  Rescaled momenta for use in ME calculations
   */
  vector<Lorentz5Momentum> rescaledMomenta_;

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
