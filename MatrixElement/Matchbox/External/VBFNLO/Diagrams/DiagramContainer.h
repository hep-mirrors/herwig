// -*- C++ -*-
//
// DiagramContainer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DiagramContainer_H
#define HERWIG_DiagramContainer_H
//
// This is the declaration of the DiagramContainer class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/MatrixElement/MEBase.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * DiagramContainer is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see \ref DiagramContainerInterfaces "The interfaces"
 * defined for DiagramContainer.
 */
  class DiagramContainer : public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DiagramContainer(const MEPtr me);

  /**
   * The destructor.
   */
  virtual ~DiagramContainer();
  //@}
    
  virtual vector<DiagPtr> getDiagrams() const {assert(false);};

    virtual string className() const { return "DiagramContainer";};
  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr) const {assert(false);};
 
  void setQuarkFlavours ( PDVector input ) { theQuarkFlavours = input; };

protected:
  /**
   * The quark flavours to be considered.
   */
  PDVector theQuarkFlavours;

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


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DiagramContainer> initDiagramContainer;

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

  const MEPtr theME;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DiagramContainer. */
template <>
struct BaseClassTrait<Herwig::DiagramContainer,1> {
  /** Typedef of the first base class of DiagramContainer. */
  typedef Herwig::Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DiagramContainer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DiagramContainer>
  : public ClassTraitsBase<Herwig::DiagramContainer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DiagramContainer"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DiagramContainer is implemented. It may also include several, space-separated,
   * libraries if the class DiagramContainer depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DiagramContainer_H */
