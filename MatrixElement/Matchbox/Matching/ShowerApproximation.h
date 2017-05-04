// -*- C++ -*-
//
// ShowerApproximation.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ShowerApproximation_H
#define Herwig_ShowerApproximation_H
//
// This is the declaration of the ShowerApproximation class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/MatrixElement/Matchbox/Utility/ColourBasis.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.fh"
#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.fh"
#include "HardScaleProfile.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ShowerApproximation describes the shower emission to be used
 * in NLO matching.
 *
 */
class ShowerApproximation: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ShowerApproximation();

  /**
   * The destructor.
   */
  virtual ~ShowerApproximation();
  //@}

public:

  /**
   * Return true, if this shower approximation will require a
   * splitting generator
   */
  virtual bool needsSplittingGenerator() const { return false; }

  /**
   * Return true, if this shower approximation will require
   * H events
   */
  virtual bool hasHEvents() const { return true; }

  /**
   * Return true, if this shower approximation will require tilde
   * XCombs for the real phase space point generated
   */
  virtual bool needsTildeXCombs() const { return false; }

  /**
   * Return true, if this shower approximation will require 
   * a truncated parton shower
   */
  virtual bool needsTruncatedShower() const { return false; }

  /**
   * Return the tilde kinematics object returning the shower
   * kinematics parametrization if different from the nominal dipole
   * mappings.
   */
  virtual Ptr<TildeKinematics>::tptr showerTildeKinematics() const;

  /**
   * Return the tilde kinematics object returning the shower
   * kinematics parametrization if different from the nominal dipole
   * mappings.
   */
  virtual Ptr<InvertedTildeKinematics>::tptr showerInvertedTildeKinematics() const;

public:

  /**
   * Set the XComb object describing the Born process
   */
  void setBornXComb(tStdXCombPtr xc) { theBornXComb = xc; }

  /**
   * Return the XComb object describing the Born process
   */
  tStdXCombPtr bornXComb() const { return theBornXComb; }

  /**
   * Return the XComb object describing the Born process
   */
  tcStdXCombPtr bornCXComb() const { return theBornXComb; }

  /**
   * Set the XComb object describing the real emission process
   */
  void setRealXComb(tStdXCombPtr xc) { theRealXComb = xc; }

  /**
   * Return the XComb object describing the real emission process
   */
  tStdXCombPtr realXComb() const { return theRealXComb; }

  /**
   * Return the XComb object describing the real emission process
   */
  tcStdXCombPtr realCXComb() const { return theRealXComb; }

  /**
   * Set the tilde xcomb objects associated to the real xcomb
   */
  void setTildeXCombs(const vector<StdXCombPtr>& xc) { theTildeXCombs = xc; }

  /**
   * Return the tilde xcomb objects associated to the real xcomb
   */
  const vector<StdXCombPtr>& tildeXCombs() const { return theTildeXCombs; }

  /**
   * Set the dipole in charge for the emission
   */
  void setDipole(Ptr<SubtractionDipole>::tptr);

  /**
   * Return the dipole in charge for the emission
   */
  Ptr<SubtractionDipole>::tptr dipole() const;

  /**
   * Return true, if this matching is capable of spin correlations.
   */
  virtual bool hasSpinCorrelations() const { return false; }

public:

  /**
   * Return true if one of the recently encountered configutations was
   * below the infrared cutoff.
   */
  bool belowCutoff() const { return theBelowCutoff; }

  /**
   * Indicate that one of the recently encountered configutations was
   * below the infrared cutoff.
   */
  void wasBelowCutoff() { theBelowCutoff = true; }

  /**
   * Reset the below cutoff flag.
   */
  void resetBelowCutoff() { theBelowCutoff = false; }

  /**
   * Return the pt cut to be applied for final-final dipoles.
   */
  Energy ffPtCut() const { return theFFPtCut; }

  /**
   * Return the pt cut to be applied for final-initial dipoles.
   */
  Energy fiPtCut() const { return theFIPtCut; }

  /**
   * Return the pt cut to be applied for initial-initial dipoles.
   */
  Energy iiPtCut() const { return theIIPtCut; }

  /**
   * Return the pt cut to be applied for initial-initial dipoles.
   */
  Energy safeCut() const { return theSafeCut;}

  /**
   * Return the screening scale to be applied for final-final dipoles.
   */
  Energy ffScreeningScale() const { return theFFScreeningScale; }

  /**
   * Return the screening scale to be applied for final-initial dipoles.
   */
  Energy fiScreeningScale() const { return theFIScreeningScale; }

  /**
   * Return the screening scale to be applied for initial-initial dipoles.
   */
  Energy iiScreeningScale() const { return theIIScreeningScale; }

  /**
   * Return the shower renormalization scale
   */
  virtual Energy2 showerEmissionScale() const;

  /**
   * Return the shower renormalization scale
   */
  Energy2 showerRenormalizationScale() const {
    return sqr(renormalizationScaleFactor())*showerEmissionScale();
  }

  /**
   * Return the shower factorization scale
   */
  Energy2 showerFactorizationScale() const {
    return sqr(factorizationScaleFactor())*showerEmissionScale();
  }

  /**
   * Return the Born renormalization scale
   */
  Energy2 bornRenormalizationScale() const;

  /**
   * Return the Born factorization scale
   */
  Energy2 bornFactorizationScale() const;

  /**
   * Return the real emission renormalization scale
   */
  Energy2 realRenormalizationScale() const;

  /**
   * Return the real emission factorization scale
   */
  Energy2 realFactorizationScale() const;

  /**
   * Enumerate possible scale choices
   */
  enum ScaleChoices {

    bornScale = 0,
    /** Use the born scales */

    realScale = 1,
    /** Use the real scales */

    showerScale = 2
    /** Use the shower scales */

  };

  /**
   * Return the scale choice in the real emission cross section to be
   * used in the matching subtraction.
   */
  int realEmissionScaleInSubtraction() const { return theRealEmissionScaleInSubtraction; }

  /**
   * Return the scale choice in the born cross section to be
   * used in the matching subtraction.
   */
  int bornScaleInSubtraction() const { return theBornScaleInSubtraction; }

  /**
   * Return the scale choice in the emission contribution to be
   * used in the matching subtraction.
   */
  int emissionScaleInSubtraction() const { return theEmissionScaleInSubtraction; }

  /**
   * Return the scale choice in the real emission cross section to be
   * used in the splitting.
   */
  int realEmissionScaleInSplitting() const { return theRealEmissionScaleInSplitting; }

  /**
   * Return the scale choice in the born cross section to be
   * used in the splitting.
   */
  int bornScaleInSplitting() const { return theBornScaleInSplitting; }

  /**
   * Return the scale choice in the emission contribution to be
   * used in the splitting.
   */
  int emissionScaleInSplitting() const { return theEmissionScaleInSplitting; }

  /**
   * Return the scale weight
   */
  double scaleWeight(int rScale, int bScale, int eScale) const;

  /**
   * Return the scale weight for the matching subtraction
   */
  double subtractionScaleWeight() const {
    return scaleWeight(realEmissionScaleInSubtraction(),
		       bornScaleInSubtraction(),
		       emissionScaleInSubtraction());
  }

  /**
   * Return the scale weight for the splitting
   */
  double splittingScaleWeight() const {
    return scaleWeight(realEmissionScaleInSplitting(),
		       bornScaleInSplitting(),
		       emissionScaleInSplitting());
  }

public:

  /**
   * Return true, if the phase space restrictions of the dipole shower should
   * be applied.
   */
  bool restrictPhasespace() const { return theRestrictPhasespace; }

  /**
   * Indicate that the phase space restrictions of the dipole shower should
   * be applied.
   */
  void restrictPhasespace(bool yes) { theRestrictPhasespace = yes; }

  /**
   * Return profile scales
   */
  Ptr<HardScaleProfile>::tptr profileScales() const { return theHardScaleProfile; }

  /**
   * Set profile scales
   */
  void profileScales(Ptr<HardScaleProfile>::ptr prof) { theHardScaleProfile = prof; }

  /**
   * Return true if maximum pt should be deduced from the factorization scale
   */
  bool hardScaleIsMuF() const { return maxPtIsMuF; }

  /**
   * Indicate that maximum pt should be deduced from the factorization scale
   */
  void hardScaleIsMuF(bool yes) { maxPtIsMuF = yes; }

  /**
   * Return the scale factor for the hard scale
   */
  double hardScaleFactor() const { return theHardScaleFactor; }

  /**
   * Set the scale factor for the hard scale
   */
  void hardScaleFactor(double f) { theHardScaleFactor = f; }

  /**
   * Get the factorization scale factor
   */
  double factorizationScaleFactor() const { return theFactorizationScaleFactor; }

  /**
   * Get the renormalization scale factor
   */
  double renormalizationScaleFactor() const { return theRenormalizationScaleFactor; }

  /**
   * Set the factorization scale factor
   */
  void factorizationScaleFactor(double f) { theFactorizationScaleFactor = f; }

  /**
   * Set the renormalization scale factor
   */
  void renormalizationScaleFactor(double f) { theRenormalizationScaleFactor = f; }

  /**
   * Determine if the configuration is below or above the cutoff.
   */
  virtual void checkCutoff();

  /**
   * Determine all kinematic variables which are not provided by the
   * dipole kinematics; store all shower variables in the respective
   * dipole object for later use.
   */
  virtual void getShowerVariables();

  /**
   * Return the shower approximation to the real emission cross
   * section for the given pair of Born and real emission
   * configurations.
   */
  virtual CrossSection dSigHatDR() const = 0;

  /**
   * Return the shower approximation splitting kernel for the given
   * pair of Born and real emission configurations in units of the
   * Born center of mass energy squared, and including a weight to
   * project onto the splitting given by the dipole used.
   */
  virtual double me2() const = 0;

  /**
   * Return the Born PDF weight
   */
  double bornPDFWeight(Energy2 muF) const;

  /**
   * Return the real emission PDF weight
   */
  double realPDFWeight(Energy2 muF) const;

protected:

  /**
   * Return true, if the shower was able to generate an emission
   * leading from the given Born to the given real emission process.
   */
  virtual bool isInShowerPhasespace() const;

  /**
   * Return true, if the shower emission leading from the given Born
   * to the given real emission process would have been generated
   * above the shower's infrared cutoff.
   */
  virtual bool isAboveCutoff() const;

  /**
   * Return the relevant hard scale
   */
  virtual Energy hardScale() const;

public:

  /**
   * Generate a weight for the given dipole channel
   */
  virtual double channelWeight(int emitter, int emission, 
			       int spectator, int bemitter) const;

  /**
   * Generate a normalized weight taking into account all channels
   */
  virtual double channelWeight() const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

public:

  /**
   * A large-N colour basis to be used when reproducing the shower
   * kernels.
   */
  Ptr<ColourBasis>::tptr largeNBasis() const { return theLargeNBasis; }

protected:

  /**
   * A large-N colour basis to be used when reproducing the shower
   * kernels.
   */
  Ptr<ColourBasis>::ptr theLargeNBasis;

  /**
   * Set the large-N basis
   */
  void setLargeNBasis();

  /**
   * The x value from which on we extrapolate PDFs for numerically stable ratios.
   */
  double theExtrapolationX;

private:

  /**
   * The XComb object describing the Born process
   */
  tStdXCombPtr theBornXComb;

  /**
   * The XComb object describing the real emission process
   */
  tStdXCombPtr theRealXComb;

  /**
   * The tilde xcomb objects associated to the real xcomb
   */
  vector<StdXCombPtr> theTildeXCombs;

  /**
   * The dipole in charge for the emission
   */
  Ptr<SubtractionDipole>::tptr theDipole;

  /**
   * True if one of the recently encountered configutations was below
   * the infrared cutoff.
   */
  bool theBelowCutoff;

  /**
   * The pt cut to be applied for final-final dipoles.
   */
  Energy theFFPtCut;

  /**
   * An optional screening scale for final-final dipoles; see
   * DipoleSplittingKernel
   */
  Energy theFFScreeningScale;

  /**
   * The pt cut to be applied for final-initial dipoles.
   */
  Energy theFIPtCut;

  /**
   * An optional screening scale for final-initial dipoles; see
   * DipoleSplittingKernel
   */
  Energy theFIScreeningScale;

  /**
   * The pt cut to be applied for initial-initial dipoles.
   */
  Energy theIIPtCut;

  /**
   * An optional screening scale for initial-initial dipoles; see
   * DipoleSplittingKernel
   */
  Energy theIIScreeningScale;

  /**
   * The cut to be applied as an enhanced shower cutoff.
   */
  Energy theSafeCut;

  /**
   * True, if the phase space restrictions of the dipole shower should
   * be applied.
   */
  bool theRestrictPhasespace;

  /**
   * The scale factor for the hard scale
   */
  double theHardScaleFactor;

  /**
   * The scale factor for the renormalization scale
   */
  double theRenormalizationScaleFactor;

  /**
   * The scale factor for the factorization scale
   */
  double theFactorizationScaleFactor;

  /**
   * The scale choice in the real emission cross section to be
   * used in the matching subtraction.
   */
  int theRealEmissionScaleInSubtraction;

  /**
   * The scale choice in the born cross section to be
   * used in the matching subtraction.
   */
  int theBornScaleInSubtraction;

  /**
   * The scale choice in the emission contribution to be
   * used in the matching subtraction.
   */
  int theEmissionScaleInSubtraction;

  /**
   * The scale choice in the real emission cross section to be
   * used in the splitting.
   */
  int theRealEmissionScaleInSplitting;

  /**
   * The scale choice in the born cross section to be
   * used in the splitting.
   */
  int theBornScaleInSplitting;

  /**
   * The scale choice in the emission contribution to be
   * used in the splitting.
   */
  int theEmissionScaleInSplitting;

  /**
   * A freezing value for the renormalization scale
   */
  Energy theRenormalizationScaleFreeze;

  /**
   * A freezing value for the factorization scale
   */
  Energy theFactorizationScaleFreeze;

  /**
   * True if maximum pt should be deduced from the factorization scale
   */
  bool maxPtIsMuF;

  /**
   * The profile scales
   */
  Ptr<HardScaleProfile>::ptr theHardScaleProfile;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerApproximation & operator=(const ShowerApproximation &);

};

}

#endif /* Herwig_ShowerApproximation_H */
