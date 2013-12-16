// -*- C++ -*-
//
// DC_Hjj_NC_qBqB.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DC_Hjj_NC_qBqB_H
#define HERWIG_DC_Hjj_NC_qBqB_H
//
// This is the declaration of the DC_Hjj_NC_qBqB class.
//

#include "DiagramContainer.h"


namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * DC_Hjj_NC_qBqB is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see ef DC_Hjj_NC_qBqBInterfaces "The interfaces"
 * defined for DC_Hjj_NC_qBqB.
 */
class DC_Hjj_NC_qBqB: public DiagramContainer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DC_Hjj_NC_qBqB(const MEPtr me) : DiagramContainer(me) {};

  /**
   * The destructor.
   */
  virtual ~DC_Hjj_NC_qBqB() {};
  //@}

  virtual vector<DiagPtr> getDiagrams() const;

  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  void setQuarkFlavours ( PDVector input ) { theQuarkFlavours = input; };
  string className() const {return "DC_Hjj_NC_qBqB";};

  protected:
  /**
   * The quark flavours to be considered.
   */
  PDVector theQuarkFlavours;

  virtual IBPtr clone() const;

  virtual IBPtr fullclone() const;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DC_Hjj_NC_qBqB. */
template <>
struct BaseClassTrait<Herwig::DC_Hjj_NC_qBqB,1> {
  /** Typedef of the first base class of DC_Hjj_NC_qBqB. */
  typedef Herwig::DiagramContainer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DC_Hjj_NC_qBqB class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DC_Hjj_NC_qBqB>
  : public ClassTraitsBase<Herwig::DC_Hjj_NC_qBqB> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DC_Hjj_NC_qBqB"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DC_Hjj_NC_qBqB is implemented. It may also include several, space-separated,
   * libraries if the class DC_Hjj_NC_qBqB depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DC_Hjj_NC_qBqB_H */
