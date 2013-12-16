// -*- C++ -*-
//
// VBFNLOMEqq2Zqqg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOMEqq2Zqqg_H
#define HERWIG_VBFNLOMEqq2Zqqg_H
//
// This is the declaration of the VBFNLOMEqq2Zqqg class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOMEPP2ZJetJetJet.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOMEqq2Zqqg is the concrete implementation of an interface to
 * the q q -> Z0 q q g matrix element of VBFNLO
 *
 * @see \ref VBFNLOMEqq2ZqqgInterfaces "The interfaces"
 * defined for VBFNLOMEqq2Zqqg.
 */
class VBFNLOMEqq2Zqqg: public VBFNLOMEPP2ZJetJetJet {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOMEqq2Zqqg();

  /**
   * The destructor.
   */
  virtual ~VBFNLOMEqq2Zqqg();
  //@}

public:

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

  /*
   * Move the momenta from meMomenta() into an array that can be read by VBFNLO
   */
  virtual void prepareMomenta(double (&pbar)[14][4], double (&qbar)[5],int) const;

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static ClassDescription<VBFNLOMEqq2Zqqg> initVBFNLOMEqq2Zqqg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOMEqq2Zqqg & operator=(const VBFNLOMEqq2Zqqg &);

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
 *  base classes of VBFNLOMEqq2Zqqg. */
template <>
struct BaseClassTrait<Herwig::VBFNLOMEqq2Zqqg,1> {
  /** Typedef of the first base class of VBFNLOMEqq2Zqqg. */
  typedef Herwig::VBFNLOMEPP2ZJetJetJet NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOMEqq2Zqqg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOMEqq2Zqqg>
  : public ClassTraitsBase<Herwig::VBFNLOMEqq2Zqqg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOMEqq2Zqqg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOMEqq2Zqqg is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOMEqq2Zqqg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOMEqq2Zqqg_H */
