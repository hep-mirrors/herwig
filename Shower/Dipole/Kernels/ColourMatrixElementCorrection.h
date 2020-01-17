// -*- C++ -*-
//
// ColourMatrixElementCorrection.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ColourMatrixElementCorrection_H
#define Herwig_ColourMatrixElementCorrection_H
//
// This is the declaration of the ColourMatrixElementCorrection class.
//

#include "Herwig/Shower/Dipole/Base/DipoleSplittingReweight.h"

#include <tuple>

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Johan Thoren, Simon Platzer
 *
 * \brief ColourMatrixElementCorrection is implementing colour matrix element
 * corrections through the weighted Sudakov algorithm
 *
 * @see \ref ColourMatrixElementCorrectionInterfaces "The interfaces"
 * defined for ColourMatrixElementCorrection.
 */
class ColourMatrixElementCorrection: public DipoleSplittingReweight {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ColourMatrixElementCorrection();

  /**
   * The destructor.
   */
  virtual ~ColourMatrixElementCorrection();
  //@}

public:

  /**
   * Calculate and cache the colour matrix element correction factor
   * for the given splitting type.
   */
  double cmec(const DipoleSplittingInfo&) const;

  /**
   * Return the reweighting factor for the given splitting type.
   */
  virtual double evaluate(const DipoleSplittingInfo& s) const {
    return cmec(s);
  }

  /**
   * Return the absolute value of the colour matrix element correction
   * as an enhancement hint for the sampling of the un-reweighted
   * splitting kernel.
   */
  virtual double hint(const DipoleSplittingInfo& s) const {
    if ( hintOnly(s) )
      return cmec(s);
    return abs(cmec(s))*lambda;
  }

  /**
   * Return true, if the reweight can be entirely absorbed into the hint. A
   * possible detuning will be switched off.
   */
  virtual bool hintOnly(const DipoleSplittingInfo& s) const {
    return cmec(s) > 0.;
  }

  /**
   * Set the factor in front of enhance used by the veto algorithm.
   */
  virtual void reweightFactor(const double c) {
    assert(c > 0.0);
    lambda = c;
  }

  /**
   * Set the factor in front of enhance used by the veto algorithm.
   */
  virtual void negativeScaling(const double c) {
    assert(c >= 0.0);
    negCMECScaling = c;
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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * Factor to shuffle the magnitude of the CMEC between the splitting kernel
   * and the weight in the reweighted veto algorithm.
   */
  double lambda;

  /**
   * Scaling factor multiplying all of the negative colour matrix element 
   * corrections. The physically sensible value is 1.0, but this factor can
   * be used to examine the effects of the negative contributions.
   */
  double negCMECScaling;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ColourMatrixElementCorrection & operator=(const ColourMatrixElementCorrection &);

};

}

#endif /* Herwig_ColourMatrixElementCorrection_H */
