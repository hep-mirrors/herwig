// -*- C++ -*-
//
// MSSM.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MSSM_H
#define HERWIG_MSSM_H
//
// This is the declaration of the MSSM class.
//

#include "SusyBase.h"
#include "MSSM.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The MSSM class provides the main model class to replace the Standard Model 
 * when using the Minimal Supersymmetric Standard Model.
 *
 * @see \ref MSSMInterfaces "The interfaces"
 * defined for MSSM.
 */
class MSSM: public SusyBase {

public:

  /**
   *  Default constructor
   */
  MSSM() : createDiagonalMixing_(false) {}

public:

  /**
   * Value of Higgs mixing angle \f$\alpha\f$.
   */
  double higgsMixingAngle() const {return theAlpha;}

  /**
   * Value of up-type trilinear couplings
   */
  const complex<Energy> & topTrilinear() const {return theAtop;}

  /**
   * Value of down-type trilinear couplings
   */
  const complex<Energy> & bottomTrilinear() const {return theAbottom;}

  /**
   * Value of lepton trilinear couplings
   */
  const complex<Energy> & tauTrilinear() const {return theAtau;}

  /**
   * The stop mixing matrix
   */
  const MixingMatrixPtr & stopMix() const {return theStopMix;}

  /**
   * The sbottom chargino mixing matrix
   */
  const MixingMatrixPtr & sbottomMix() const {return theSbotMix;}

  /**
   * The stau mixing matrix
   */
  const MixingMatrixPtr & stauMix() const {return theStauMix;}

  /**
   * Mixing matrix for the neutral CP-even Higgs bosons
   */
  const MixingMatrixPtr & CPevenHiggsMix() const {return theHiggsMix;}

  /**
   * Mixing matrix for the neutral CP-odd Higgs bosons
   */
  const MixingMatrixPtr & CPoddHiggsMix() const {
    return HiggsAMix_;
  }

  /**
   * Mixing matrix for the charged Higgs bosons
   */
  const MixingMatrixPtr & ChargedHiggsMix() const {
    return HiggsPMix_;
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
   *  Create the mixing matrices for the model
   */
  virtual void createMixingMatrices();

  /**
   *  Extract the parameters from the input blocks
   */
  virtual void extractParameters(bool checkModel=true);

  /**
   * Adjust row of Mixing Matrix if a negative mass occurs in LHA file
   * @param id The PDG code of the particle with a negative mass
   */
  virtual void adjustMixingMatrix(long id);

  /**
   * Mixing matrix for the neutral CP-even Higgs bosons
   */
  void CPevenHiggsMix(MixingMatrixPtr in) {theHiggsMix = in;}

  /**
   * Mixing matrix for the neutral CP-odd Higgs bosons
   */
  void CPoddHiggsMix(MixingMatrixPtr in) {HiggsAMix_ = in;}

  /**
   * Mixing matrix for the charged Higgs bosons
   */
  void ChargedHiggsMix(MixingMatrixPtr in) {HiggsPMix_ = in;}

  /**
   * Set the stop mixing matrix
   */
  void stopMix(MixingMatrixPtr in) {theStopMix =in;}

  /**
   * The sbottom chargino mixing matrix
   */
  void sbottomMix(MixingMatrixPtr in) {theSbotMix = in;}

  /**
   * The stau mixing matrix
   */
  void stauMix(MixingMatrixPtr in) {theStauMix = in;}

  /**
   * Value of Higgs mixing angle \f$\alpha\f$.
   */
  void higgsMixingAngle(double in) {theAlpha = in;}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MSSM & operator=(const MSSM &) = delete;

private:

  /**
   *  Third generation squark and slepton mixing matrices
   */
  //@{
  /**
   *  The \f$\tilde{t}\f$ mixing matrix
   */
  MixingMatrixPtr theStopMix;

  /**
   *  The \f$\tilde{b}\f$ mixing matrix
   */
  MixingMatrixPtr theSbotMix;

  /**
   *  The \f$\tilde{\tau}\f$ mixing matrix
   */
  MixingMatrixPtr theStauMix;
  //@}

  /**
   * Trilinear couplings stored as vector of complex numbers to make use
   * of routine already available to read complex matrices
   */
  //@{
  /**
   *  For the up type squarks
   */
  complex<Energy> theAtop;

  /**
   *  For the down type squarks
   */
  complex<Energy> theAbottom;

  /**
   *  For the charged sleptons
   */
  complex<Energy> theAtau;
  //@}

  /**
   * Value of higgs mixing angle.
   */
  double theAlpha;

  /**
   *  Higgs boson mixing matrices
   */
  //@{
  /**
   *  Scalar Higgs mixing matrix
   */
  MixingMatrixPtr theHiggsMix;
  /**
   *  Pseudoscalar Higgs mixing
   */
  MixingMatrixPtr HiggsAMix_;

  /**
   *  Charged Higgs mixing
   */
  MixingMatrixPtr HiggsPMix_;
  //@}

  /**
   *  Create diagonal the stop, bottom and stau mixing matrices if needed
   */
  bool createDiagonalMixing_;

};

}

#endif /* HERWIG_MSSM_H */
