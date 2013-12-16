// -*- C++ -*-
//
// VBFNLOAmplitudePP2hJetJet.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOAmplitudePP2hJetJet_H
#define HERWIG_VBFNLOAmplitudePP2hJetJet_H
//
// This is the declaration of the VBFNLOAmplitudePP2hJetJet class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOAmplitudeVVJJNeutralBase.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOAmplitudePP2hJetJet is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see \ref VBFNLOAmplitudePP2hJetJetInterfaces "The interfaces"
 * defined for VBFNLOAmplitudePP2hJetJet.
 */
class VBFNLOAmplitudePP2hJetJet: public VBFNLOAmplitudeVVJJNeutralBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOAmplitudePP2hJetJet();

  /**
   * The destructor.
   */
  virtual ~VBFNLOAmplitudePP2hJetJet();
  //@}

public:

  /**
   * Return a ME instance appropriate for this amplitude and the given
   * subprocesses
   */
  virtual Ptr<MatchboxMEBase>::ptr makeME(const PDVector&) const;

  /**
   * Check if the given diagram is allowed for the interfaced
   * flavour configuration and parameters
   */
  virtual bool allowedProcess(const PDVector &) const;

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return 0; }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const { return 3; }

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const {return 0;}

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const {return 3;}

  // /**
  //  * Return true, if this matrix element will generate
  //  * momenta for the incoming partons itself. 
  //  * The matrix element is required to store the incoming 
  //  * parton momenta in meMomenta()[0,1]. No mapping
  //  * in tau and y is performed by the PartonExtractor object,
  //  * if a derived class returns true here.
  //  */
  // virtual bool haveX1X2() const { return (phasespace() ? phasespace()->haveX1X2() : true); }

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
   * If this matrix element is considered an underlying Born matrix
   * element in the context of a subtracted real emission, but
   * actually neglecting a subclass of the contributing diagrams,
   * return true if the given emitter-spectator configuration
   * should not be considered when setting up subtraction dipoles.
   */
  virtual bool noDipole(int,int) const;

  /**
   * This function calls VBFNLO and returns the requested matrix element
   * squared.
   */
  virtual void VbfnloMe2(double pbar[14][4],int sign[14],double qbar[5],
			 int,int,double&,double&,double&,double&,double&,double&) const;


  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return true; }

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

  // /*
  //  * Initialize the VBFNLO phasespace generator for this matrix element.
  //  */
  // virtual void initPSGen();

  /**
   * Initialize VBFNLO process specific common blocks
   */
  virtual void initProcess(const int &) const;

  // /**
  //  * Initialize diagram containers
  //  */
  // virtual void initDiagramContainers() const;

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
  virtual void prepareMomenta(double (&pbar)[14][4],double (&qbar)[5],int) const;

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
  static ClassDescription<VBFNLOAmplitudePP2hJetJet> initVBFNLOAmplitudePP2hJetJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOAmplitudePP2hJetJet & operator=(const VBFNLOAmplitudePP2hJetJet &);

  /*
   * Control if we include virtual corrections
   */
  int theNLOmode;

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
 *  base classes of VBFNLOAmplitudePP2hJetJet. */
template <>
struct BaseClassTrait<Herwig::VBFNLOAmplitudePP2hJetJet,1> {
  /** Typedef of the first base class of VBFNLOAmplitudePP2hJetJet. */
  typedef Herwig::VBFNLOAmplitudeVVJJNeutralBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOAmplitudePP2hJetJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOAmplitudePP2hJetJet>
  : public ClassTraitsBase<Herwig::VBFNLOAmplitudePP2hJetJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOAmplitudePP2hJetJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOAmplitudePP2hJetJet is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOAmplitudePP2hJetJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOAmplitudePP2hJetJet_H */
