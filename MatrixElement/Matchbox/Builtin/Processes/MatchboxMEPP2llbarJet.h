// -*- C++ -*-
//
// MatchboxMEPP2llbarJet.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEPP2llbarJet_H
#define HERWIG_MatchboxMEPP2llbarJet_H
//
// This is the declaration of the MatchboxMEPP2llbarJet class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Builtin/Processes/MatchboxMEllbarqqbarg.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMEPP2llbarJet implements charged lepton
 * pair plus jet production in hadron-hadron collisions.
 *
 * @see \ref MatchboxMEPP2llbarJetInterfaces "The interfaces"
 * defined for MatchboxMEPP2llbarJet.
 */
class MatchboxMEPP2llbarJet: public MatchboxMEBase, public MatchboxMEllbarqqbarg {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMEPP2llbarJet();

  /**
   * The destructor.
   */
  virtual ~MatchboxMEPP2llbarJet();
  //@}

public:

  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const { return phasespace() ? phasespace()->nDim(3) : 5; }

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return 1; }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return 2; }

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector. Derived classes
   * must call this method once internal degrees of freedom are setup
   * and finally return the result of this method.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Return the number of light flavours, this matrix
   * element is calculated for.
   */
  virtual unsigned int nLight() const { return theQuarkFlavours.size(); }

  /**
   * Return the colour correlated matrix element squared with
   * respect to the given two partons as appearing in mePartonData(),
   * suitably scaled by sHat() to give a dimension-less number.
   */
  virtual double colourCorrelatedME2(pair<int,int>) const;

  /**
   * Return the colour and spin correlated matrix element squared for
   * the gluon indexed by the first argument with and the
   * given correlation tensor build as p1^\mu p2^\nu
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& p2) const;

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

protected:

  /**
   * The lepton flavours to be considered.
   */
  PDVector theLeptonFlavours;

  /**
   * The quark flavours to be considered.
   */
  PDVector theQuarkFlavours;

  /**
   * The Z mass squared
   */
  Energy2 theZMass2;

  /**
   * The Z width squared
   */
  Energy2 theZWidth2;

private:

  /**
   * A user defined scale; if zero, the center of mass
   * energy is used.
   */
  Energy theUserScale;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<MatchboxMEPP2llbarJet> initMatchboxMEPP2llbarJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMEPP2llbarJet & operator=(const MatchboxMEPP2llbarJet &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MatchboxMEPP2llbarJet. */
template <>
struct BaseClassTrait<Herwig::MatchboxMEPP2llbarJet,1> {
  /** Typedef of the first base class of MatchboxMEPP2llbarJet. */
  typedef Herwig::MatchboxMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MatchboxMEPP2llbarJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MatchboxMEPP2llbarJet>
  : public ClassTraitsBase<Herwig::MatchboxMEPP2llbarJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MatchboxMEPP2llbarJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MatchboxMEPP2llbarJet is implemented. It may also include several, space-separated,
   * libraries if the class MatchboxMEPP2llbarJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MatchboxMEPP2llbarJet_H */
