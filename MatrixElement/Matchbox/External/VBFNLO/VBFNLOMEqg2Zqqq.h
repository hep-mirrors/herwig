// -*- C++ -*-
//
// VBFNLOMEqg2Zqqq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOMEqg2Zqqq_H
#define HERWIG_VBFNLOMEqg2Zqqq_H
//
// This is the declaration of the VBFNLOMEqg2Zqqq class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOMEPP2ZJetJetJet.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOMEqg2Zqqq is the concrete implementation of an interface to
 * the q g -> Z0 q q q matrix element of VBFNLO with one of
 *
 * @see \ref VBFNLOMEqg2ZqqqInterfaces "The interfaces"
 * defined for VBFNLOMEqg2Zqqq.
 */
class VBFNLOMEqg2Zqqq: public VBFNLOMEPP2ZJetJetJet {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOMEqg2Zqqq();

  /**
   * The destructor.
   */
  virtual ~VBFNLOMEqg2Zqqq();
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
  static ClassDescription<VBFNLOMEqg2Zqqq> initVBFNLOMEqg2Zqqq;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOMEqg2Zqqq & operator=(const VBFNLOMEqg2Zqqq &);

  /**
   * Set which of the incoming partons
   * is the incoming gluon.
   */
  int theWhichGluon;

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
 *  base classes of VBFNLOMEqg2Zqqq. */
template <>
struct BaseClassTrait<Herwig::VBFNLOMEqg2Zqqq,1> {
  /** Typedef of the first base class of VBFNLOMEqg2Zqqq. */
  typedef Herwig::VBFNLOMEPP2ZJetJetJet NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOMEqg2Zqqq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOMEqg2Zqqq>
  : public ClassTraitsBase<Herwig::VBFNLOMEqg2Zqqq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOMEqg2Zqqq"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOMEqg2Zqqq is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOMEqg2Zqqq depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOMEqg2Zqqq_H */
