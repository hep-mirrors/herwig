// -*- C++ -*-
//
// MatchboxMEllbar2qqbar.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEllbar2qqbar_H
#define HERWIG_MatchboxMEllbar2qqbar_H
//
// MatchboxMEllbar2qqbar
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Builtin/Processes/MatchboxMEllbarqqbar.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Martin Stoll
 *
 * \brief MatchboxMEllbar2qqbar implements the
 * matrix element for charged lepton annihilation
 * into a light q qbar pair.
 *
 * @see \ref MatchboxMEllbar2qqbarInterfaces "The interfaces"
 * defined for MatchboxMEllbar2qqbar.
 */
class MatchboxMEllbar2qqbar: public MatchboxMEBase, public MatchboxMEllbarqqbar {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMEllbar2qqbar();

  /**
   * The destructor.
   */
  virtual ~MatchboxMEllbar2qqbar();
  //@}

public:

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return 0; }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return 2; }

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;
  
  /**
   * Evaluate the finite virtual correction for the
   * variables supplied through the Born XComb object
   * and possible additional random numbers.
   */
  virtual double oneLoopInterference() const;
  
  /**
   * Return true, if this matrix element is capable of calculating
   * one-loop (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return true; };

  /**
   * Return true, if one loop corrections have been calculated in
   * dimensional reduction. Otherwise conventional dimensional
   * regularization is assumed. Note that renormalization is always
   * assumed to be MSbar.
   */
  virtual bool isDR() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of the integrated dipoles.
   */
  virtual bool isCS() const { return true; }

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
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { return false; }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const { return false; }

  /**
   * Return the colour correlated matrix element squared with
   * respect to the given two partons as appearing in mePartonData(),
   * suitably scaled by sHat() to give a dimension-less number.
   */
  virtual double colourCorrelatedME2(pair<int,int>) const;

  /**
   * Return the colour and spin correlated matrix element squared for
   * the gluon indexed by the first argument using the given
   * correlation tensor.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

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
   * The lepton flavours to be considered.
   */
  PDVector theLeptonFlavours;

  /**
   * The quark flavours to be considered.
   */
  PDVector theQuarkFlavours;

  /**
   * A user defined scale; if zero, the center of mass
   * energy is used.
   */
  Energy theUserScale;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MatchboxMEllbar2qqbar> initMatchboxMEllbar2qqbar;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMEllbar2qqbar & operator=(const MatchboxMEllbar2qqbar &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MatchboxMEllbar2qqbar. */
template <>
struct BaseClassTrait<Herwig::MatchboxMEllbar2qqbar,1> {
  /** Typedef of the first base class of MatchboxMEllbar2qqbar. */
  typedef Herwig::MatchboxMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MatchboxMEllbar2qqbar class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MatchboxMEllbar2qqbar>
  : public ClassTraitsBase<Herwig::MatchboxMEllbar2qqbar> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MatchboxMEllbar2qqbar"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MatchboxMEllbar2qqbar is implemented. It may also include several, space-separated,
   * libraries if the class MatchboxMEllbar2qqbar depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MatchboxMEllbar2qqbar_H */
