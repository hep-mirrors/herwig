// -*- C++ -*-
//
// VBFNLOMEqq2hqq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOMEqq2hqq_H
#define HERWIG_VBFNLOMEqq2hqq_H
//
// This is the declaration of the VBFNLOMEqq2hqq class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOMEVVJJNeutralBase.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOMEqq2hqq is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see \ref VBFNLOMEqq2hqqInterfaces "The interfaces"
 * defined for VBFNLOMEqq2hqq.
 */
class VBFNLOMEqq2hqq: public VBFNLOMEPP2hJetJet {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOMEqq2hqq();

  /**
   * The destructor.
   */
  virtual ~VBFNLOMEqq2hqq();
  //@}

public:


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

  /**
   * Temporary entry for rebuilding subprocess management
   */
  virtual bool hasSplitSubprocs() { return true; }

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VBFNLOMEqq2hqq> initVBFNLOMEqq2hqq;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOMEqq2hqq & operator=(const VBFNLOMEqq2hqq &);

  /*
   * Control if we include virtual corrections
   */
  int theNLOmode;

  /*
   * Control the decay channel simulated by VBFNLO
   */
  int theDecayChannel;

  /**
   * Should this ME calculate the Higgs boson decay 
   * in narrow width approximation?
   */
  bool theNarrowWidth;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VBFNLOMEqq2hqq. */
template <>
struct BaseClassTrait<Herwig::VBFNLOMEqq2hqq,1> {
  /** Typedef of the first base class of VBFNLOMEqq2hqq. */
  typedef Herwig::VBFNLOMEVVJJNeutralBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOMEqq2hqq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOMEqq2hqq>
  : public ClassTraitsBase<Herwig::VBFNLOMEqq2hqq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOMEqq2hqq"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOMEqq2hqq is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOMEqq2hqq depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOMEqq2hqq_H */
