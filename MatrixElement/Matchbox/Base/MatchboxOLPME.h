// -*- C++ -*-
//
// MatchboxOLPME.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxOLPME_H
#define Herwig_MatchboxOLPME_H
//
// This is the declaration of the MatchboxOLPME class.
//

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxOLPME implements OLP interfaces.
 */
class MatchboxOLPME: public MatchboxAmplitude {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxOLPME();

  /**
   * The destructor.
   */
  virtual ~MatchboxOLPME();
  //@}

public:

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector& p,
			 Ptr<MatchboxFactory>::tptr,
			 bool) const;

  /**
   * Set the (tree-level) order in \f$g_S\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGs(unsigned int ogs) { theOrderInGs = ogs; }

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const { return theOrderInGs; }

  /**
   * Set the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGem(unsigned int oge) { theOrderInGem = oge; }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const { return theOrderInGem; }

  /**
   * Return true, if this amplitude is handled by a BLHA one-loop provider
   */
  virtual bool isOLPTree() const { return true; }

  /**
   * Return true, if this amplitude is handled by a BLHA one-loop provider
   */
  virtual bool isOLPLoop() const { return true; }

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { return false; }

  /**
   * Set the xcomb object.
   */
  virtual void setXComb(tStdXCombPtr xc);

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr) {}

  /**
   * Return the matrix element squared.
   */
  virtual double me2() const;

  /**
   * Return the colour correlated matrix element.
   */
  virtual double colourCorrelatedME2(pair<int,int> ij) const;

  /**
   * Return the large-N colour correlated matrix element.
   */
  virtual double largeNColourCorrelatedME2(pair<int,int>,
					   Ptr<ColourBasis>::tptr) const;

  /**
   * Return the largeN matrix element squared.
   */
  virtual double largeNME2(Ptr<ColourBasis>::tptr largeNBasis) const;

  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> ij,
					 const SpinCorrelationTensor& c) const;

  /**
   * Return the spin correlated matrix element.
   */
  virtual double spinCorrelatedME2(pair<int,int> ij,
				   const SpinCorrelationTensor& c) const;

  /**
   * Return true, if tree-level contributions will be evaluated at amplitude level.
   */
  virtual bool treeAmplitudes() const { return false; }

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return true; }

  /**
   * Return true, if this amplitude only provides
   * one-loop (QCD) corrections.
   */
  virtual bool onlyOneLoop() const { return false; }

  /**
   * Return true, if one-loop contributions will be evaluated at amplitude level.
   */
  virtual bool oneLoopAmplitudes() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of everything expanded.
   */
  virtual bool isExpanded() const { return true; }

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const;

  /**
   * Indicate that this amplitude is running alphas by itself.
   */
  virtual bool hasRunningAlphaS() const;

  /**
   * Indicate that this amplitude is running alphaew by itself.
   */
  virtual bool hasRunningAlphaEW() const;

  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const;

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const;

  /**
   * Calculate the one-loop amplitudes for the phasespace point
   * stored in lastXComb, if provided.
   */
  virtual void prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr) {}

  /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const;

public:

  /**
   * Call OLP_EvalSubProcess and fill in the results
   */
  virtual void evalSubProcess() const = 0;

  /**
   * Fill in results for the given colour correlator
   */
  virtual void evalColourCorrelator(pair<int,int> ij) const = 0;

  /**
   * Fill in results for the given colour/spin correlator
   */
  virtual void evalSpinColourCorrelator(pair<int,int> ij) const = 0;

  /**
   * Fill in results for the given spin correlator; may not be supported
   */
  virtual void evalSpinCorrelator(pair<int,int> ij) const;

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

  /**
   * Set an optional contract file name to be used
   */
  static string& optionalContractFile() {
    static string s = "";
    return s;
  }

  /**
   * Indicate that the OLP has been started
   */
  static bool& didStartOLP() {
    static bool f = false;
    return f;
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxOLPME & operator=(const MatchboxOLPME &) = delete;

  /**
   * The (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  unsigned int theOrderInGs;

  /**
   * The (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  unsigned int theOrderInGem;

  /**
   * Set the value of the dimensional regularization parameter 
   * to the value of the renormalization scale
   */
  bool theSetMuToMuR;

  /**
   * Use the running alpha_s instead of the reference alpha_s.
   * This also sets hasRunningAlphaS() to true.
   */
  bool theUseRunningAlphaS;

  /**
   * Use the running alpha_ew instead of the reference alpha_ew.
   * This also sets hasRunningAlphaEW() to true.
   */
  bool theUseRunningAlphaEW;

};

}

#endif /* Herwig_MatchboxOLPME_H */
