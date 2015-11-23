// -*- C++ -*-
#ifndef HERWIG_DrellYanBase_H
#define HERWIG_DrellYanBase_H
//
// This is the declaration of the DrellYanBase class.
//

#include "HwMEBase.h"
#include "Herwig/Shower/ShowerConfig.h"
#include "Herwig/Shower/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DrellYanBase class class provides a base class for the implemented
 * of Drell-Yan type processes and provides the matrix element and POWHEG
 * style hard corrections
 *
 * @see \ref DrellYanBaseInterfaces "The interfaces"
 * defined for DrellYanBase.
 */
class DrellYanBase: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  DrellYanBase();

  /**
   *  Has a POWHEG style correction
   */
  //virtual bool hasPOWHEGCorrection() {return _alpha;}

  virtual POWHEGType hasPOWHEGCorrection() {return ISR;}


  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return _alpha;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(ShowerTreePtr, double & initial,
				      double & final) {
    final   = 1.;
    initial = 1.;
  }

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

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
  virtual HardTreePtr generateHardest(ShowerTreePtr,
				      vector<ShowerInteraction::Type>);

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object.
   */
  virtual void setKinematics() {
    HwMEBase::setKinematics();
    mb2_ = sHat();
  }

protected:

  /**
   *  Return the momenta and type of hard matrix element correction
   * @param quarks The original incoming particles.
   * @param beams The BeamParticleData objects
   * @param boson The momentum of the original outgoing gauge boson
   * @param iemit Whether the first (0) or second (1) particle emitted
   * the radiation
   * @param itype The type of radiated particle (0 is gluon, 1 is quark 
   *              and 2 is antiquark)
   * @param pnew The momenta of the new particles
   * @param trans The LorentzRotation from the boson rest frame to the new lab
   * @param xnew The new values of the momentuym fractions
   * @return Whether or not the matrix element correction needs to be applied
   */
  bool applyHard(ShowerParticleVector quarks,
		 vector<tcBeamPtr> beams,
		 Lorentz5Momentum boson,unsigned int & iemit,
		 unsigned int & itype,vector<Lorentz5Momentum> & pnew,
		 LorentzRotation & trans, pair<double,double> & xnew);

  /**
   * Returns the matrix element for a given type of process,
   * rapidity of the jet \f$y_j\f$ and transverse momentum \f$p_T\f$
   * @param emis_type the type of emission,
   * (0 is \f$q\bar{q}\to Vg\f$, 1 is \f$qg\to Vq\f$ and 2 is \f$g\bar{q}\to V\bar{q}\f$)
   * @param pt The transverse momentum of the jet
   * @param yj The rapidity of the jet
   */
  double getResult(int emis_type, Energy pt, double yj);
 
  /**
   *  generates the hardest emission (yj,p)
   * @param pnew The momenta of the new particles
   * @param emissiontype The type of emission, as for getResult
   * @return Whether not an emission was generated
   */
  bool getEvent(vector<Lorentz5Momentum> & pnew,int & emissiontype);
  
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

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DrellYanBase> initDrellYanBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DrellYanBase & operator=(const DrellYanBase &);

private:

  /**
   *  Mass squared of the vector boson
   */
  Energy2 mb2_;

  /**
   *  Parameters for the old-style ME correction
   */
  //@{
  /**
   *  Relative weight for the \f$q\bar{q}\f$ and \f$q/\bar{q}g\f$ channels
   */
  double _channelwgtA;

  /**
   *  Relative weight for the \f$qg\f$ and \f$\bar{q}g\f$ channels 
   */
  double _channelwgtB;

  /**
   *  Weights for the channels as a vector
   */
  vector<double> _channelweights;
  
  /**
   *  Number of weights greater than 1
   */
  unsigned int _nover;

  /**
   *  Maximum weight
   */
  double _maxwgt;
  //@}

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   * The power, \f$n\f$, for the sampling
   */
  double _power;

  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double _preqqbar;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double _preqg;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double _pregqbar;

  /**
   *  The prefactors as a vector for easy use
   */
  vector<double> _prefactor;
  //@}

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> _beams;
  
  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> _partons;
  //@}

  /**
   *  Properties of the boson and jets
   */
  //@{
  /**
   *  The rapidity of the gauge boson
   */
  double _yb;

  /**
   *  The mass of the gauge boson
   */
  Energy _mass;

  /**
   *  Whether the quark is in the + or - z direction
   */
  bool _quarkplus;

  /**
   *  the rapidity of the jet
   */
  double _yj;

  /**
   *  The transverse momentum of the jet
   */
  Energy _pt;
  //@}

  /**
   *  The transverse momentum of the jet
   */
  Energy _min_pt;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alpha;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DrellYanBase. */
template <>
struct BaseClassTrait<Herwig::DrellYanBase,1> {
  /** Typedef of the first base class of DrellYanBase. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DrellYanBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DrellYanBase>
  : public ClassTraitsBase<Herwig::DrellYanBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DrellYanBase"; }
};

/** @endcond */

}

#endif /* HERWIG_DrellYanBase_H */
