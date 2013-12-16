// -*- C++ -*-
//
// DiagramCollection.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DiagramCollection_H
#define HERWIG_DiagramCollection_H
//
// This is the declaration of the DiagramCollection class.
//
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"


namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * DiagramCollection is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see \ref DiagramCollectionInterfaces "The interfaces"
 * defined for DiagramCollection.
 */
  class DiagramCollection : public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DiagramCollection();

  /**
   * The destructor.
   */
  virtual ~DiagramCollection();
  //@}
    
  virtual vector<DiagPtr> getDiagrams() const {assert(false);};

    virtual string className() const { return "DiagramCollection";};
  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const {assert(false);};
 
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
  static ClassDescription<DiagramCollection> initDiagramCollection;

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

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DiagramCollection. */
template <>
struct BaseClassTrait<Herwig::DiagramCollection,1> {
  /** Typedef of the first base class of DiagramCollection. */
  typedef Herwig::Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DiagramCollection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DiagramCollection>
  : public ClassTraitsBase<Herwig::DiagramCollection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DiagramCollection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DiagramCollection is implemented. It may also include several, space-separated,
   * libraries if the class DiagramCollection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DiagramCollection_H */
