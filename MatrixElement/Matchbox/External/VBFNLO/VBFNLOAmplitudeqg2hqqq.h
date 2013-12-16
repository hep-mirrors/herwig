// -*- C++ -*-
//
// VBFNLOAmplitudeqg2hqqq.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VBFNLOAmplitudeqg2hqqq_H
#define HERWIG_VBFNLOAmplitudeqg2hqqq_H
//
// This is the declaration of the VBFNLOAmplitudeqg2hqqq class.
//

#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOAmplitudeVVJJNeutralBase.h"
#include "Herwig++/MatrixElement/Matchbox/External/VBFNLO/VBFNLOAmplitudePP2hJetJetJet.h"

namespace Herwig {


using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold
 *
 * VBFNLOAmplitudeqg2hqqq is the concrete implementation of an interface to
 * the PP -> Higgs Jet Jet matrix element of VBFNLO.
 *
 * @see \ref VBFNLOAmplitudeqg2hqqqInterfaces "The interfaces"
 * defined for VBFNLOAmplitudeqg2hqqq.
 */
class VBFNLOAmplitudeqg2hqqq: public VBFNLOAmplitudePP2hJetJetJet {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  VBFNLOAmplitudeqg2hqqq();

  /**
   * The destructor.
   */
  virtual ~VBFNLOAmplitudeqg2hqqq();
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

  /* /\** */
  /*  * Return the order in \f$\alpha_S\f$ in which this matrix element */
  /*  * is given. */
  /*  *\/ */
  /* virtual unsigned int orderInAlphaS() const { return 0; } */

  /* /\** */
  /*  * Return the order in \f$\alpha_{EM}\f$ in which this matrix */
  /*  * element is given. Returns 0. */
  /*  *\/ */
  /* virtual unsigned int orderInAlphaEW() const { return 3; } */

  /* /\** */
  /*  * Return the (tree-level) order in \f$g_S\f$ in which this matrix */
  /*  * element is given. */
  /*  *\/ */
  /* virtual unsigned int orderInGs() const {return 0;} */

  /* /\** */
  /*  * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix */
  /*  * element is given. */
  /*  *\/ */
  /* virtual unsigned int orderInGem() const {return 3;} */

  // /**
  //  * Return true, if this matrix element will generate
  //  * momenta for the incoming partons itself. 
  //  * The matrix element is required to store the incoming 
  //  * parton momenta in meMomenta()[0,1]. No mapping
  //  * in tau and y is performed by the PartonExtractor object,
  //  * if a derived class returns true here.
  //  */
  // virtual bool haveX1X2() const { return (phasespace() ? phasespace()->haveX1X2() : true); }

  /* /\** */
  /*  * Return true, if this matrix element provides the PDF */
  /*  * weight for the first incoming parton itself. */
  /*  *\/ */
  /* virtual bool havePDFWeight1() const { return true; } */

  /* /\** */
  /*  * Return true, if this matrix element provides the PDF */
  /*  * weight for the second incoming parton itself. */
  /*  *\/ */
  /* virtual bool havePDFWeight2() const { return true; } */

  /* /\** */
  /*  * Return the colour correlated matrix element squared with */
  /*  * respect to the given two partons as appearing in mePartonData(), */
  /*  * suitably scaled by sHat() to give a dimension-less number. */
  /*  *\/ */
  /* virtual double colourCorrelatedME2(pair<int,int>) const; */

  /* /\** */
  /*  * Return the colour and spin correlated matrix element squared for */
  /*  * the gluon indexed by the first argument using the given */
  /*  * correlation tensor. */
  /*  *\/ */
  /* virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator, */
  /* 					 const SpinCorrelationTensor& c) const; */

  /* /\** */
  /*  * If this matrix element is considered an underlying Born matrix */
  /*  * element in the context of a subtracted real emission, but */
  /*  * actually neglecting a subclass of the contributing diagrams, */
  /*  * return true if the given emitter-spectator configuration */
  /*  * should not be considered when setting up subtraction dipoles. */
  /*  *\/ */
  /* virtual bool noDipole(int,int) const; */

  // /**
  //  * Add all possible diagrams with the add() function.
  //  */
  // virtual vector<Ptr<Tree2toNDiagram>::ptr> getDiagrams() const;

  // /**
  //  * With the information previously supplied with the
  //  * setKinematics(...) method, a derived class may optionally
  //  * override this method to weight the given diagrams with their
  //  * (although certainly not physical) relative probabilities.
  //  */
  // virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;

  // /**
  //  * Return a Selector with possible colour geometries for the selected
  //  * diagram weighted by their relative probabilities.
  //  */
  // virtual Selector<const ColourLines *>
  // colourGeometries(tcDiagPtr diag) const;

  // /**
  //  * Generate internal degrees of freedom given nDim() uniform random
  //  * numbers in the interval ]0,1[. To help the phase space generator,
  //  * the 'dSigHatDR' should be a smooth function of these numbers,
  //  * although this is not strictly necessary. The return value should
  //  * be true of the generation succeeded. If so the generated momenta
  //  * should be stored in the meMomenta() vector.
  //  */
  // virtual bool generateKinematics(const double * r);

  // /**
  //  * The number of internal degreed of freedom used in the matrix
  //  * element.
  //  */
  // virtual int nDim() const;

  // /**
  //  * This function calls VBFNLO and returns the requested matrix element
  //  * squared.
  //  */
  // virtual void VbfnloMe2(double pbar[14][4],int sign[14],double qbar[5],
  // 			 int,int,int,double&,double&,double&,double&,double&,double&) const;
  
  // /**
  //  * Return true, if this matrix element expects
  //  * the incoming partons in their center-of-mass system
  //  */
  // virtual bool wantCMS() const { return phasespace() ? phasespace()->wantCMS() : false; }

  /* /\** */
  /*  * The number of jets that this matrix element generates */
  /*  *\/ */
  /* virtual int Njets() const; */

  /* /\** */
  /*  * Return the number of decay products for the selected decay channel */
  /*  *\/ */
  /* int NDecayProducts() const; */

  /* /\** */
  /*  * Return the branching ratio for the selected decay channel */
  /*  *\/ */
  /* double BranchingRatio() const; */

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

  // /*
  //  * Initialize the VBFNLO phasespace generator for this matrix element.
  //  */
  // virtual void initPSGen();

  /* /\** */
  /*  * Initialize VBFNLO process specific common blocks */
  /*  *\/ */
  /* virtual void initProcess(const int &) const; */

  // /**
  //  * Initialize diagram containers
  //  */
  // virtual void initDiagramContainers() const;

protected:

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
  static ClassDescription<VBFNLOAmplitudeqg2hqqq> initVBFNLOAmplitudeqg2hqqq;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFNLOAmplitudeqg2hqqq & operator=(const VBFNLOAmplitudeqg2hqqq &);

  /**
   * Set which of the incoming partons
   * is the incoming gluon.
   */
  int theWhichGluon;

  /* /\* */
  /*  * Control if we include virtual corrections */
  /*  *\/ */
  /* int theNLOmode; */

  /* /\* */
  /*  * Control the decay channel simulated by VBFNLO */
  /*  *\/ */
  /* int theDecayChannel; */

  /* /\** */
  /*  * Should this ME calculate the Higgs boson decay  */
  /*  * in narrow width approximation? */
  /*  *\/ */
  /* bool theNarrowWidth; */
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VBFNLOAmplitudeqg2hqqq. */
template <>
struct BaseClassTrait<Herwig::VBFNLOAmplitudeqg2hqqq,1> {
  /** Typedef of the first base class of VBFNLOAmplitudeqg2hqqq. */
  typedef Herwig::VBFNLOAmplitudePP2hJetJetJet NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFNLOAmplitudeqg2hqqq class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFNLOAmplitudeqg2hqqq>
  : public ClassTraitsBase<Herwig::VBFNLOAmplitudeqg2hqqq> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFNLOAmplitudeqg2hqqq"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFNLOAmplitudeqg2hqqq is implemented. It may also include several, space-separated,
   * libraries if the class VBFNLOAmplitudeqg2hqqq depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMatchbox.so HwMatchboxVBFNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VBFNLOAmplitudeqg2hqqq_H */
