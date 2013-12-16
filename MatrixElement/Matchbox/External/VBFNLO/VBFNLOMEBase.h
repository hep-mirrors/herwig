// -*- C++ -*-
//
// VBFNLOMEBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOMEBase_H
#define HERWIG_VBFNLOMEBase_H
//
// This is the declaration of the VBFNLOMEBase class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Diagrams/DiagramContainer.h"
#include <iostream>

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOMEBase is the base class for matrix elements
 * in the context of the matchbox NLO interface.
 *
 * @see \ref VBFNLOMEBaseInterfaces "The interfaces"
 * defined for VBFNLOMEBase.
 */
class VBFNLOMEBase: public MatchboxMEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOMEBase();

  /**
   * The destructor.
   */
  virtual ~VBFNLOMEBase();
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
  virtual unsigned int orderInAlphaEW() const { return 1; }

  /**
   * Return the symmetry factor for identical final state particles.
   */
  virtual double finalStateSymmetry() const { return 1;}

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return lastSHat(); }

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /*
   * Initialize the VBFNLO phasespace generator for this matrix element.
   */
  virtual void initPSGen() = 0;

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

public:

  /**
   * Convert momenta to VBFNLO compliant format
   */
  void L5MomToDouble(Lorentz5Momentum L5Mom, double * DoubleMom) const{
    *DoubleMom++ = L5Mom.t()/GeV;
    *DoubleMom++ = L5Mom.x()/GeV;
    *DoubleMom++ = L5Mom.y()/GeV;
    *DoubleMom++ = L5Mom.z()/GeV;
    return;
  };

protected:

  /**
   * Copy the defined coupling constants into the Fortran common blocks.
   */
  void initCouplings();
  
  /**
   * A derieved class may override this method to figure out an integer which
   * may then be used by VBFNLO to calculate a corresponding subprocess.
   */
  virtual int getSubprocessID(int) const = 0;

  /**
   * A derieved class may override this method to initialize values
   * in common blocks once per run.
   */
  virtual void initProcess(const int &) const = 0;

  /**
   * Initialize diagram containers
   */
  virtual void initDiagramContainers() const = 0;

protected:

  /**
   * The quark flavours to be considered.
   */
  PDVector theQuarkFlavours;

  /**
   * A user defined scale; if zero, the center of mass
   * energy is used.
   */
  Energy theUserScale;

  /**
   * A vector with all allowed diagrams for all possible 
   * configurations defined by the interfaces.
   */
  mutable vector<DiagPtr> allPossibleDiagrams;

  /**
   * A container from which we get the allowed diagrams
   */
  //DiagramContainer* theDiagramContainer;
  mutable RCPtr<DiagramContainer> theDiagramContainer;

  int theInteger; //debugging

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static AbstractClassDescription<VBFNLOMEBase> initVBFNLOMEBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOMEBase & operator=(const VBFNLOMEBase &);
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VBFNLOMEBase. */
template <>
struct BaseClassTrait<Herwig::VBFNLOMEBase,1> {
  /** Typedef of the first base class of VBFNLOMEBase. */
  typedef Herwig::MatchboxMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOMEBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOMEBase>
  : public ClassTraitsBase<Herwig::VBFNLOMEBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOMEBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOMEBase is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOMEBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOMEBase_H */
