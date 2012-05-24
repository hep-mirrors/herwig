// -*- C++ -*-
//
// DipoleSplittingKernel.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleSplittingKernel_H
#define HERWIG_DipoleSplittingKernel_H
//
// This is the declaration of the DipoleSplittingKernel class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/StandardModel/AlphaSBase.h"
#include "ThePEG/PDF/PDF.h"

#include "Herwig++/DipoleShower/Utility/PDFRatio.h"
#include "Herwig++/DipoleShower/Base/DipoleSplittingInfo.h"
#include "Herwig++/DipoleShower/Kinematics/DipoleSplittingKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleSplittingKernel is the base class for all kernels
 * used within the dipole shower.
 *
 * @see \ref DipoleSplittingKernelInterfaces "The interfaces"
 * defined for DipoleSplittingKernel.
 */
class DipoleSplittingKernel: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleSplittingKernel();

  /**
   * The destructor.
   */
  virtual ~DipoleSplittingKernel();
  //@}

public:

  /**
   * Return the alpha_s to be used
   */
  Ptr<AlphaSBase>::tptr alphaS() const { return theAlphaS; }

  /**
   * Set the alpha_s to be used
   */
  void alphaS(Ptr<AlphaSBase>::tptr ap) { theAlphaS = ap; }

  /**
   * Return the splitting kinematics object
   */
  Ptr<DipoleSplittingKinematics>::tptr splittingKinematics() const { return theSplittingKinematics; }

  /**
   * Return the mc check object
   */
  Ptr<DipoleMCCheck>::ptr mcCheck() const { return theMCCheck; }

  /**
   * Set the splitting kinematics object
   */
  void splittingKinematics(Ptr<DipoleSplittingKinematics>::tptr sp) { theSplittingKinematics = sp; }

  /**
   * Return the PDFRatio object
   */
  Ptr<PDFRatio>::tptr pdfRatio() const { return thePDFRatio; }

  /**
   * Set the PDFRatio object
   */
  void pdfRatio(Ptr<PDFRatio>::tptr sp) { thePDFRatio = sp; }

  /**
   * Return the number of additional parameter
   * random variables needed to evaluate this kernel
   * except the momentum fractions of incoming partons.
   * These will be accessible through the 
   * lastSplittingParameters() container of the splitting
   * info object.
   */
  virtual int nDimAdditional() const { return 0; }

public:

  /**
   * Return true, if this splitting kernel
   * applies to the given dipole index.
   */
  virtual bool canHandle(const DipoleIndex&) const = 0;

  /**
   * Return true, if this splitting kernel is
   * the same for the given index a, as the given
   * splitting kernel for index b.
   */
  virtual bool canHandleEquivalent(const DipoleIndex& a,
				   const DipoleSplittingKernel& sk,
				   const DipoleIndex& b) const = 0;

  /**
   * Return the emitter data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr emitter(const DipoleIndex&) const = 0;

  /**
   * Return the emission data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr emission(const DipoleIndex&) const = 0;

  /**
   * Return the spectator data after splitting, given
   * a dipole index.
   */
  virtual tcPDPtr spectator(const DipoleIndex&) const = 0;

  /**
   * Return the flavour produced, if this cannot
   * be determined from the dipole.
   */
  PDPtr flavour() const { return theFlavour; }

  /**
   * Return true, if this splitting kernel is supposed to work in a
   * strict large-N limit, i.e. replacing C_F by C_A/2
   */
  bool strictLargeN() const { return theStrictLargeN; }

public:

  /**
   * Inform this splitting kernel, that it is being
   * presampled until a call to stopPresampling
   */
  virtual void startPresampling(const DipoleIndex&) {}

  /**
   * Inform this splitting kernel, that it is not being
   * presampled until a call to startPresampling
   */
  virtual void stopPresampling(const DipoleIndex&) {}

  /**
   * Return the number of points to presample this
   * splitting generator.
   */
  unsigned long presamplingPoints() const { return thePresamplingPoints; }

  /**
   * Return the maximum number of trials
   * to generate a splitting.
   */
  unsigned long maxtry() const { return theMaxtry; }

  /**
   * Evaluate this splitting kernel for the given
   * dipole splitting.
   */
  virtual double evaluate(const DipoleSplittingInfo&) const = 0;

  /**
   * Return true, if this kernel is capable of
   * delivering an overestimate to the kernel, and
   * of inverting the integral over the overestimate
   * w.r.t. the phasepsace provided by the given
   * DipoleSplittingInfo object.
   */
  virtual bool haveOverestimate(const DipoleSplittingInfo&) const { return false; }

  /**
   * Return the overestimate to this splitting kernel 
   * for the given dipole splitting.
   */
  virtual double overestimate(const DipoleSplittingInfo&) const { return -1.; }

  /**
   * Invert the integral over the overestimate 
   * w.r.t. the phasepsace provided by the given
   * DipoleSplittingInfo object to equal
   * the given value.
   */
  virtual double invertOverestimateIntegral(const DipoleSplittingInfo&, double) const { return -1.; }

public:

  /**
   * Get the factorization scale factor
   */
  double factorizationScaleFactor() const { return theFactorizationScaleFactor; }

  /**
   * Set the factorization scale factor
   */
  void factorizationScaleFactor(double f) { theFactorizationScaleFactor = f; }

  /**
   * Get the renormalization scale factor
   */
  double renormalizationScaleFactor() const { return theRenormalizationScaleFactor; }

  /**
   * Set the renormalization scale factor
   */
  void renormalizationScaleFactor(double f) { theRenormalizationScaleFactor = f; }

protected:

  /**
   * Return the common factor of (alphas/2pi)*(pdf ratio)
   */
  double alphaPDF(const DipoleSplittingInfo&) const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The alpha_s to be used.
   */
  Ptr<AlphaSBase>::ptr theAlphaS;

  /**
   * An optional 'colour screening' scale
   * for alternative intrinsic pt generation.
   */
  Energy theScreeningScale;

  /**
   * The splitting kinematics to be used.
   */
  Ptr<DipoleSplittingKinematics>::ptr theSplittingKinematics;

  /**
   * An optional PDF ratio object to be used
   * when evaluating this kernel.
   */
  Ptr<PDFRatio>::ptr thePDFRatio;

  /**
   * The number of points to presample this
   * splitting generator.
   */
  unsigned long thePresamplingPoints;

  /**
   * The maximum number of trials
   * to generate a splitting.
   */
  unsigned long theMaxtry;

  /**
   * The flavour produced, if this cannot
   * be determined from the dipole.
   */
  PDPtr theFlavour;

  /**
   * Pointer to a check histogram object
   */
  Ptr<DipoleMCCheck>::ptr theMCCheck;

  /**
   * True, if this splitting kernel is supposed to work in a
   * strict large-N limit, i.e. replacing C_F by C_A/2
   */
  bool theStrictLargeN;

  /**
   * The factorization scale factor.
   */
  double theFactorizationScaleFactor;

  /**
   * The renormalization scale factor.
   */
  double theRenormalizationScaleFactor;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DipoleSplittingKernel> initDipoleSplittingKernel;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleSplittingKernel & operator=(const DipoleSplittingKernel &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleSplittingKernel. */
template <>
struct BaseClassTrait<Herwig::DipoleSplittingKernel,1> {
  /** Typedef of the first base class of DipoleSplittingKernel. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleSplittingKernel class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleSplittingKernel>
  : public ClassTraitsBase<Herwig::DipoleSplittingKernel> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleSplittingKernel"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleSplittingKernel is implemented. It may also include several, space-separated,
   * libraries if the class DipoleSplittingKernel depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleSplittingKernel_H */
