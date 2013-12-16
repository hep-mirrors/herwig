// -*- C++ -*-
//
// VBFNLOAmplitudeVVJJNeutralBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOAmplitudeVVJJNeutralBase_H
#define HERWIG_VBFNLOAmplitudeVVJJNeutralBase_H
//
// This is the declaration of the VBFNLOAmplitudeVVJJNeutralBase class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOAmplitudeBase.h"


namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOAmplitudeVVJJNeutralBase is the base class for matrix elements
 * of neutral boson production via vector boson fusion interfaced
 * from VBFNLO
 *
 * @see \ref VBFNLOAmplitudeVVJJNeutralBaseInterfaces "The interfaces"
 * defined for VBFNLOAmplitudeVVJJNeutralBase.
 */
class VBFNLOAmplitudeVVJJNeutralBase: public VBFNLOAmplitudeBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOAmplitudeVVJJNeutralBase();

  /**
   * The destructor.
   */
  virtual ~VBFNLOAmplitudeVVJJNeutralBase();
  //@}

public:

  /**
   * Return the matrix element squared.
   */
  virtual double me2() const; 

  /**
   * The number of jets that this matrix element generates
   */
  virtual int Njets() const = 0;

  /**
   * Access on whether particle 1 is a particle or antiparticle
   */
  virtual bool Incoming1() const { return theIncoming1; };

  /**
   * Access on whether particle 2 is a particle or antiparticle
   */
  virtual bool Incoming2() const { return theIncoming2; };

  /**
   * Temporary entry for rebuilding subprocess management
   */
  virtual bool hasSplitSubprocs() { return false; }

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const;

  // /**
  //  * Return true, if the amplitude implements diagrams.
  //  */
  // virtual bool haveDiagrams() const { return false; }

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

  /*
   * Move the momenta from meMomenta() into an array that can be read by VBFNLO
   */
  virtual void prepareMomenta(double (&pbar)[14][4], double (&qbar)[5], int gluonIndex = -1) const = 0;

  /**
   * Check if this matrix element object allows the given parton as incoming parton on beam 1
   */
  virtual bool requestedAsIncoming1(PDPtr) const;

  /**
   * Check if this matrix element object allows the given parton as incoming parton on beam 2
   */
  virtual bool requestedAsIncoming2(PDPtr) const;

  /**
   * Check if the given diagram is allowed for the interfaced
   * flavour configuration and parameters
   */
  bool allowedDiagram(DiagPtr, int gluonIndex = -1) const;

  /**
   * Check if the given diagram is allowed for the interfaced
   * flavour configuration and parameters
   */
  virtual bool allowedProcess(const PDVector &) const = 0;

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
   * Indicates that this is a concrete class with persistent data.
   */
  static AbstractClassDescription<VBFNLOAmplitudeVVJJNeutralBase> initVBFNLOAmplitudeVVJJNeutralBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOAmplitudeVVJJNeutralBase & operator=(const VBFNLOAmplitudeVVJJNeutralBase &);

protected:

  /**
   * A derieved class may override this method to initialize values
   * in common blocks once per run.
   */
  virtual void initProcess(const int &) const = 0;

  /**
   * Get the subprocess ID determined by incoming and outgoing flavours
   */
  virtual int getSubprocessID(int gluonIndex = 2) const;

  /**
   * Get the subprocess ID determined by incoming and outgoing flavours
   */
  virtual int getSubprocessID(const PDVector& data, int gluonIndex = 2) const;

  /**
   * A derieved class should overwrite this method with a call of the
   * corresponding VBFNLO matrix element subroutine
   */
  virtual void VbfnloMe2(double pbar[14][4],int sign[14],double qbar[5],int,int,double&,double&,double&,double&,double&,double&) const = 0;
  
public:

  /**
   * Return true, if one loop corrections have been calculated in
   * dimensional reduction. Otherwise conventional dimensional
   * regularization is assumed. Note that renormalization is always
   * assumed to be MSbar.
   */
  virtual bool isDR() const { return true; }

  /**
   * Return true, if the virtual correction has been calculated in the
   * dipole convention.
   */
  virtual bool isCS() const { return true; }

  /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const;

protected:

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

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VBFNLOAmplitudeVVJJNeutralBase. */
template <>
struct BaseClassTrait<Herwig::VBFNLOAmplitudeVVJJNeutralBase,1> {
  /** Typedef of the first base class of VBFNLOAmplitudeVVJJNeutralBase. */
  typedef Herwig::VBFNLOAmplitudeBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOAmplitudeVVJJNeutralBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOAmplitudeVVJJNeutralBase>
  : public ClassTraitsBase<Herwig::VBFNLOAmplitudeVVJJNeutralBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOAmplitudeVVJJNeutralBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOAmplitudeVVJJNeutralBase is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOAmplitudeVVJJNeutralBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOAmplitudeVVJJNeutralBase_H */
