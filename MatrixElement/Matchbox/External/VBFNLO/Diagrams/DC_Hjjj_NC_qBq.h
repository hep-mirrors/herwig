// -*- C++ -*-
//
// DC_Hjjj_NC_qBq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DC_Hjjj_NC_qBq_H
#define HERWIG_DC_Hjjj_NC_qBq_H
//
// This is the declaration of the DC_Hjjj_NC_qBq class.
//

#include "DiagramContainer.h"


namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * DC_Hjjj_NC_qBq is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see ef DC_Hjjj_NC_qBqInterfaces "The interfaces"
 * defined for DC_Hjjj_NC_qBq.
 */
class DC_Hjjj_NC_qBq: public DiagramContainer {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DC_Hjjj_NC_qBq(const MEPtr me) : DiagramContainer(me) {};

  /**
   * The destructor.
   */
  virtual ~DC_Hjjj_NC_qBq() {};
  //@}

  virtual vector<DiagPtr> getDiagrams() const;

  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  void setQuarkFlavours ( PDVector input ) { theQuarkFlavours = input; };
  string className() const {return "DC_Hjjj_NC_qBq";};

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
 *  base classes of DC_Hjjj_NC_qBq. */
template <>
struct BaseClassTrait<Herwig::DC_Hjjj_NC_qBq,1> {
  /** Typedef of the first base class of DC_Hjjj_NC_qBq. */
  typedef Herwig::DiagramContainer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DC_Hjjj_NC_qBq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DC_Hjjj_NC_qBq>
  : public ClassTraitsBase<Herwig::DC_Hjjj_NC_qBq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DC_Hjjj_NC_qBq"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DC_Hjjj_NC_qBq is implemented. It may also include several, space-separated,
   * libraries if the class DC_Hjjj_NC_qBq depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DC_Hjjj_NC_qBq_H */
