// -*- C++ -*-
//
// VBFNLOVirtualMEVVJJNeutral.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOVirtualMEVVJJNeutral_H
#define HERWIG_VBFNLOVirtualMEVVJJNeutral_H
//
// This is the declaration of the VBFNLOVirtualMEVVJJNeutral class.
//
#include "Herwig++/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "VBFNLOMEVVJJNeutralBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * \brief VBFNLOVirtualMEVVJJNeutral is the concrete class of virtual corrections
 * to matrix elements derieved from VBFNLOMEVVJJNeutralBase
 *
 * @see \ref VBFNLOVirtualMEVVJJNeutralInterfaces "The interfaces"
 * defined for VBFNLOVirtualMEVVJJNeutral.
 */
  class VBFNLOVirtualMEVVJJNeutral: public MatchboxInsertionOperator{

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOVirtualMEVVJJNeutral();

  /**
   * The destructor.
   */
  virtual ~VBFNLOVirtualMEVVJJNeutral();
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

public:

  /**
   * Return the Born matrix element this class represents 
   * virtual corrections to.
   */
    //  Ptr<VBFNLOMEVVJJNeutralBase>::tptr BornMEVVJJNeutral() { return theBornMEVVJJNeutral; }

    Ptr<VBFNLOMEVVJJNeutralBase>::tptr BornMEVVJJNeutral() { 
      return dynamic_ptr_cast<Ptr<VBFNLOMEVVJJNeutralBase>::tptr>(MatchboxInsertionOperator::lastBorn());
    }

  /**
   * Return the Born matrix element this class represents 
   * virtual corrections to.
   */  
    Ptr<VBFNLOMEVVJJNeutralBase>::tcptr BornMEVVJJNeutral() const { 
      return dynamic_ptr_cast<Ptr<VBFNLOMEVVJJNeutralBase>::tcptr>(MatchboxInsertionOperator::lastBorn());
    }
    //  Ptr<VBFNLOMEVVJJNeutralBase>::tcptr BornMEVVJJNeutral() const { return theBornMEVVJJNeutral; }

  /**
   * Evaluate the finite virtual correction for the
   * variables supplied through the Born XComb object
   * and possible additional random numbers.
   */
  virtual double me2() const;
    
  /**
   * Return true, if this virtual correction
   * applies to the given process.
   */
  virtual bool apply(const cPDVector&) const;

  /**
   * Return true, if this virtual correction
   * has been calculated using dimensional reduction.
   * CDR is assumed otherwise.
   */
  virtual bool isDR() const { return true; }

  /**
   * Return true, if the virtual correction has been calculated in the
   * dipole convention.
   */
  virtual bool isCS() const { return true; }

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

  //@}


private:

  /**
   * The Born matrix element this class represents 
   * virtual corrections to.
   */
  Ptr<MatchboxMEBase>::tptr theLastBorn;

  /**
   * The Born matrix element this class represents 
   * virtual corrections to. 
   */
  Ptr<VBFNLOMEVVJJNeutralBase>::ptr theBornMEVVJJNeutral;

  /**
   * The exchanged current
   */
  int theCurrent;

  /**
   * The possible exchanged currents
   */
  enum currents {
    neutral = 0,
    charged = 1
  };

  /**
   * Control if parton 1 is a particle or antiparticle
   */
  bool theIncoming1;

  /**
   * Control if parton 2 is a particle or antiparticle
   */
  bool theIncoming2;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static ClassDescription<VBFNLOVirtualMEVVJJNeutral> initVBFNLOVirtualMEVVJJNeutral;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOVirtualMEVVJJNeutral & operator=(const VBFNLOVirtualMEVVJJNeutral &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VBFNLOVirtualMEVVJJNeutral. */
template <>
struct BaseClassTrait<Herwig::VBFNLOVirtualMEVVJJNeutral,1> {
  /** Typedef of the first base class of VBFNLOVirtualMEVVJJNeutral. */
  typedef Herwig::MatchboxInsertionOperator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOVirtualMEVVJJNeutral class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOVirtualMEVVJJNeutral>
  : public ClassTraitsBase<Herwig::VBFNLOVirtualMEVVJJNeutral> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOVirtualMEVVJJNeutral"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOVirtualMEVVJJNeutral is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOVirtualMEVVJJNeutral depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOVirtualMEVVJJNeutral_H */
