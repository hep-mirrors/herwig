// -*- C++ -*-
//
// VBFNLOMEVVJJChargedBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOMEVVJJChargedBase_H
#define HERWIG_VBFNLOMEVVJJChargedBase_H
//
// This is the declaration of the VBFNLOMEVVJJChargedBase class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOMEBase.h"

//#include "VBFNLOVirtualMEVVJJCharged.h"





namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOMEVVJJChargedBase is the base class for matrix elements
 * in the context of the matchbox NLO interface.
 *
 * @see \ref VBFNLOMEVVJJChargedBaseInterfaces "The interfaces"
 * defined for VBFNLOMEVVJJChargedBase.
 */
class VBFNLOMEVVJJChargedBase: public VBFNLOMEBase {

public:

  friend class VBFNLOVirtualMEVVJJCharged;

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOMEVVJJChargedBase();

  /**
   * The destructor.
   */
  virtual ~VBFNLOMEVVJJChargedBase();
  //@}

public:

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * The number of jets that this matrix element generates
   */
  virtual int Njets() const = 0;

  /*
   * Should the next call of this ME include virtual corrections?
   */
  int NLOmode() const { return theNLOmode; };
  
  /*
   * Set whether we are including virtual corrections
   */
  void setNLOmode(int n) { theNLOmode = n; }; 

  /**
   * Access on whether particle 1 is a particle or antiparticle
   */
  virtual bool Incoming1() const { return theIncoming1; };

  /**
   * Access on whether particle 2 is a particle or antiparticle
   */
  virtual bool Incoming2() const { return theIncoming2; };

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

  /*
   * Initialize the VBFNLO phasespace generator for this matrix element.
   */
  virtual void initPSGen() = 0;

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
   * Check if this matrix element object allows the given parton as incoming parton on beam 1
   */
  virtual bool requestedAsIncoming1tc(tcPDPtr) const;

  /**
   * Check if this matrix element object allows the given parton as incoming parton on beam 2
   */
  virtual bool requestedAsIncoming2tc(tcPDPtr) const;

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
  static AbstractClassDescription<VBFNLOMEVVJJChargedBase> initVBFNLOMEVVJJChargedBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOMEVVJJChargedBase & operator=(const VBFNLOMEVVJJChargedBase &);

  /**
   * Control if we include virtual corrections
   */
  int theNLOmode;

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
   * A derieved class should overwrite this method with a call of the
   * corresponding VBFNLO matrix element subroutine
   */
  virtual void VbfnloMe2(double pbar[14][4],int sign[6],double qbar[5],int,int,int,double&,double&,double&,double&) const = 0;

  /**
   * The quark flavours to be considered.
   */
  PDPtr theExchangedBoson;

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
 *  base classes of VBFNLOMEVVJJChargedBase. */
template <>
struct BaseClassTrait<Herwig::VBFNLOMEVVJJChargedBase,1> {
  /** Typedef of the first base class of VBFNLOMEVVJJChargedBase. */
  typedef Herwig::VBFNLOMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOMEVVJJChargedBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOMEVVJJChargedBase>
  : public ClassTraitsBase<Herwig::VBFNLOMEVVJJChargedBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOMEVVJJChargedBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOMEVVJJChargedBase is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOMEVVJJChargedBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOMEVVJJChargedBase_H */
