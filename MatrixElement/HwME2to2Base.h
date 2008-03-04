// -*- C++ -*-
//
// HwME2to2Base.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ME2to2Base_H
#define HERWIG_ME2to2Base_H
// This is the declaration of the ME2to2Base class.

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Switch.fh"

namespace Herwig {
using namespace ThePEG;

/**
 * HwME2to2Base can be used as a base class for any matrix element class
 * implementing 2\f$\rightarrow\f$ 2 processes.
 *
 * It is heavily based on the ME2to2Base class of ThePEG but makes a number of
 * changes to give us greater control over the masses of the outgoing
 * particles so that they can be
 * - set massless where required by gauge invariance
 * - have their off-shell masses generated using the sophisticated approaches
 *   available in Herwig++.
 *
 * @see \ref HwME2to2BaseInterfaces "The interfaces"
 * defined for HwME2to2Base.
 * @see MEBase
 */
class HwME2to2Base: public ME2to2Base {

public:

  /**
   * Default constructor.
   */
  inline HwME2to2Base();

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
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object.
   */
  virtual void setKinematics();
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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

 protected:

  /**
   *  Set the treatment of the outgoing masses
   * @param first is the first outgoing particle
   * @param iopt The option for the treatment of the mass
   */
  inline void massOption(bool first, unsigned int iopt);

  /**
   * Set the treatment of the rescaling of the momenta for 
   * the matrix element calculation
   * @param iopt The rescaling option
   */
  inline void rescalingOption(unsigned int iopt);

  /**
   *  rescale the momenta for the computation of the helicity matrix element
   */
  bool rescaleMomenta(const vector<Lorentz5Momentum> &,
		      const cPDVector &);

  /**
   *  Access to the rescaled momenta
   */
  inline const vector<Lorentz5Momentum> & rescaledMomenta() const;

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<HwME2to2Base> initHwME2to2Base;

  /**
   *  Private and non-existent assignment operator.
   */
  HwME2to2Base & operator=(const HwME2to2Base &);

private:

  /**
   *  Treatment of the mass of the first particle
   */
  unsigned int _massopt1;

  /**
   *  Treatment of the mass of the second particle
   */
  unsigned int _massopt2;

  /**
   *  Produced to produce rescaled momenta
   */
  unsigned int _rescaleOption;

  /**
   *  Rescaled momenta for use in ME calculations
   */
  vector<Lorentz5Momentum> _rescaledMomenta;
};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs Herwig about the
 * base class of HwME2to2Base.
 */
template <>
struct BaseClassTrait<Herwig::HwME2to2Base,1>: public ClassTraitsType {
  /** Typedef of the base class of HwME2to2Base. */
  typedef ME2to2Base NthBase;
};

/**
 * This template specialization informs Herwig about the name of the
 * HwME2to2Base class.
 */
template <>
struct ClassTraits<Herwig::HwME2to2Base>
  : public ClassTraitsBase<Herwig::HwME2to2Base> {
  /** Return the class name. */
  static string className() { return "Herwig::HwME2to2Base"; }
};

/** @endcond */

}

#include "HwME2to2Base.icc"

#endif /* HERWIG_HwME2to2Base_H */
