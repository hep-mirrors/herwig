// -*- C++ -*-
//
// VBFNLOAmplitudePP2hJetJetJet.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOAmplitudePP2hJetJetJet_H
#define HERWIG_VBFNLOAmplitudePP2hJetJetJet_H
//
// This is the declaration of the VBFNLOAmplitudePP2hJetJetJet class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOAmplitudeVVJJNeutralBase.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOAmplitudePP2hJetJetJet is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see \ref VBFNLOAmplitudePP2hJetJetJetInterfaces "The interfaces"
 * defined for VBFNLOAmplitudePP2hJetJetJet.
 */
class VBFNLOAmplitudePP2hJetJetJet: public VBFNLOAmplitudeVVJJNeutralBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOAmplitudePP2hJetJetJet();

  /**
   * The destructor.
   */
  virtual ~VBFNLOAmplitudePP2hJetJetJet();
  //@}

public:

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return 1; }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return 3; }

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const {return 1;}

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const {return 3;}

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { return true; }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const { return true; }

  /**
   * Return the colour correlated matrix element squared with
   * respect to the given two partons as appearing in mePartonData(),
   * suitably scaled by sHat() to give a dimension-less number.
   */
  virtual double colourCorrelatedME2(pair<int,int>) const;

  /**
   * Return the colour and spin correlated matrix element squared for
   * the gluon indexed by the first argument using the given
   * correlation tensor.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;

  /**
   * This function calls VBFNLO and returns the requested matrix element
   * squared.
   */
  virtual void VbfnloMe2(double pbar[14][4],int sign[14],double qbar[5],
			 int,int,double&,double&,double&,double&,double&,double&) const;
  
  /**
   * The number of jets that this matrix element generates
   */
  virtual int Njets() const;

  /**
   * Return the number of decay products for the selected decay channel
   */
  int NDecayProducts() const;

  /**
   * Return the branching ratio for the selected decay channel
   */
  double BranchingRatio() const;

protected:

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize VBFNLO process specific common blocks
   */
  virtual void initProcess(const int &) const;

protected:

  /**
   * Access the decay channel of ths matrix element
   */
  virtual int DecayChannel() const {return theDecayChannel;}; 


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
  static AbstractClassDescription<VBFNLOAmplitudePP2hJetJetJet> initVBFNLOAmplitudePP2hJetJetJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOAmplitudePP2hJetJetJet & operator=(const VBFNLOAmplitudePP2hJetJetJet &);

 protected:

  /*
   * Control the decay channel simulated by VBFNLO
   */
  int theDecayChannel;

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
 *  base classes of VBFNLOAmplitudePP2hJetJetJet. */
template <>
struct BaseClassTrait<Herwig::VBFNLOAmplitudePP2hJetJetJet,1> {
  /** Typedef of the first base class of VBFNLOAmplitudePP2hJetJetJet. */
  typedef Herwig::VBFNLOAmplitudeVVJJNeutralBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOAmplitudePP2hJetJetJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOAmplitudePP2hJetJetJet>
  : public ClassTraitsBase<Herwig::VBFNLOAmplitudePP2hJetJetJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOAmplitudePP2hJetJetJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOAmplitudePP2hJetJetJet is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOAmplitudePP2hJetJetJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOAmplitudePP2hJetJetJet_H */

