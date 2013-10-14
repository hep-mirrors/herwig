// -*- C++ -*-
//
// MatchboxHybridAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxHybridAmplitude_H
#define Herwig_MatchboxHybridAmplitude_H
//
// This is the declaration of the MatchboxHybridAmplitude class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxHybridAmplitude unifies two amplitude objects to
 * provide tree and one-loop matrix elements.
 *
 * @see \ref MatchboxHybridAmplitudeInterfaces "The interfaces"
 * defined for MatchboxHybridAmplitude.
 */
class MatchboxHybridAmplitude: public Herwig::MatchboxAmplitude {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxHybridAmplitude();

  /**
   * The destructor.
   */
  virtual ~MatchboxHybridAmplitude();
  //@}

public:

  /**
   * Return the amplitude object to provide tree-level amplitudes.
   */
  Ptr<MatchboxAmplitude>::tptr treeLevelAmplitude() const { return theTreeLevelAmplitude; }

  /**
   * Set the amplitude object to provide tree-level amplitudes.
   */
  void treeLevelAmplitude(Ptr<MatchboxAmplitude>::tptr amp) { theTreeLevelAmplitude = amp; }

  /**
   * Return the amplitude object to provide one-loop amplitudes.
   */
  Ptr<MatchboxAmplitude>::tptr oneLoopAmplitude() const { return theOneLoopAmplitude; }

  /**
   * Set the amplitude object to provide one-loop amplitudes.
   */
  void oneLoopAmplitude(Ptr<MatchboxAmplitude>::tptr amp) { theOneLoopAmplitude = amp; }

  /**
   * Return true, if the two amplitude objects can be used in a
   * consistent way.
   */
  bool isConsistent() const;

public:

  /** @name Subprocess information */
  //@{

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector& p,
			 Ptr<MatchboxFactory>::tptr f) const { 
    return 
      treeLevelAmplitude()->canHandle(p,f) &&
      oneLoopAmplitude()->canHandle(p,f) &&
      isConsistent();
  }

  /**
   * Return the number of random numbers required to evaluate this
   * amplitude at a fixed phase space point.
   */
  virtual int nDimAdditional() const { 
    return
      treeLevelAmplitude()->nDimAdditional() ?
      treeLevelAmplitude()->nDimAdditional() :
      oneLoopAmplitude()->nDimAdditional();
  }

  /**
   * Return a ME instance appropriate for this amplitude and the given
   * subprocesses
   */
  virtual Ptr<MatchboxMEBase>::ptr makeME(const PDVector& p) const {
    return treeLevelAmplitude()->makeME(p);
  }

  /**
   * Set the (tree-level) order in \f$g_S\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGs(unsigned int n) {
    treeLevelAmplitude()->orderInGs(n);
    oneLoopAmplitude()->orderInGs(n);
  }

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const {
    return treeLevelAmplitude()->orderInGs();
  }

  /**
   * Set the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGem(unsigned int n) {
    treeLevelAmplitude()->orderInGem(n);
    oneLoopAmplitude()->orderInGem(n);
  }

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const {
    return treeLevelAmplitude()->orderInGem();
  }

  /**
   * Tell whether the outgoing partons should be sorted when determining
   * allowed subprocesses. Otherwise, all permutations are counted as
   * separate subprocesses.
   */
  virtual bool sortOutgoing() {
    return treeLevelAmplitude()->sortOutgoing();
  }

  /**
   * Return true, if this amplitude already includes averaging over
   * incoming parton's quantum numbers.
   */
  virtual bool hasInitialAverage() const { 
    return treeLevelAmplitude()->hasInitialAverage();
  }

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const {
    return treeLevelAmplitude()->hasFinalStateSymmetry();
  }

  /**
   * Return true, if this amplitude is handled by a BLHA one-loop provider
   */
  virtual bool isOLPTree() const { return false; }

  /**
   * Return true, if this amplitude is handled by a BLHA one-loop provider
   */
  virtual bool isOLPLoop() const { 
    return oneLoopAmplitude()->isOLPLoop();
  }

  /**
   * Return true, if colour and spin correlated matrix elements should
   * be ordered from the OLP
   */
  virtual bool needsOLPCorrelators() const { return false; }

  /**
   * Start the one loop provider, if appropriate. This default
   * implementation writes an BLHA 2.0 order file and starts the OLP
   */
  virtual bool startOLP(const map<pair<Process,int>,int>& procs) {
    return oneLoopAmplitude()->startOLP(procs);
  }

  //@}

  /** @name Colour basis. */
  //@{

  /**
   * Return the colour basis.
   */
  virtual Ptr<ColourBasis>::tptr colourBasis() const { 
    return treeLevelAmplitude()->colourBasis();
  }

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { 
    return treeLevelAmplitude()->haveColourFlows();
  }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *> colourGeometries(tcDiagPtr diag) const {
    return treeLevelAmplitude()->colourGeometries(diag);
  }

  //@}

  /** @name Phasespace point, crossing and helicities */
  //@{

  /**
   * Set the xcomb object.
   */
  virtual void setXComb(tStdXCombPtr xc) {
    treeLevelAmplitude()->setXComb(xc);
    oneLoopAmplitude()->setXComb(xc);
  }

  /**
   * Perform a normal ordering of external legs and fill the
   * crossing information as. This default implementation sorts
   * lexicographically in (abs(colour)/spin/abs(charge)), putting pairs
   * of particles/anti-particles where possible.
   */
  virtual void fillCrossingMap(size_t shift = 0) {
    treeLevelAmplitude()->fillCrossingMap(shift);
  }

  //@}

  /** @name Tree-level amplitudes */
  //@{

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return the matrix element squared.
   */
  virtual double me2() const {
    return treeLevelAmplitude()->me2();
  }

  /**
   * Return the colour correlated matrix element.
   */
  virtual double colourCorrelatedME2(pair<int,int> ij) const {
    return treeLevelAmplitude()->colourCorrelatedME2(ij);
  }

  /**
   * Return the large-N colour correlated matrix element.
   */
  virtual double largeNColourCorrelatedME2(pair<int,int> ij,
					   Ptr<ColourBasis>::tptr largeNBasis) const {
    return treeLevelAmplitude()->largeNColourCorrelatedME2(ij,largeNBasis);
  }

  /**
   * Return a positive helicity polarization vector for a gluon of
   * momentum p (with reference vector n) to be used when evaluating
   * spin correlations.
   */
  virtual LorentzVector<Complex> plusPolarization(const Lorentz5Momentum& p,
						  const Lorentz5Momentum& n,
						  int id = -1) const {
    return treeLevelAmplitude()->plusPolarization(p,n,id);
  }

  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const {
    return treeLevelAmplitude()->spinColourCorrelatedME2(emitterSpectator,c);
  }


  /**
   * Return true, if tree-level contributions will be evaluated at amplitude level.
   */
  virtual bool treeAmplitudes() const { 
    return treeLevelAmplitude()->treeAmplitudes();
  }

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t a, const vector<int>& hel, Complex& largeN) { 
    return treeLevelAmplitude()->evaluate(a,hel,largeN);
  }

  //@}

  /** @name One-loop amplitudes */
  //@{

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
  virtual bool oneLoopAmplitudes() const { 
    return oneLoopAmplitude()->oneLoopAmplitudes();
  }

  /**
   * Return true, if one loop corrections have been calculated in
   * dimensional reduction. Otherwise conventional dimensional
   * regularization is assumed. Note that renormalization is always
   * assumed to be MSbar.
   */
  virtual bool isDR() const { 
    return oneLoopAmplitude()->isDR();
  }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of the integrated dipoles.
   */
  virtual bool isCS() const { 
    return oneLoopAmplitude()->isCS();
  }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */
  virtual bool isBDK() const { 
    return oneLoopAmplitude()->isBDK();
  }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of everything expanded.
   */
  virtual bool isExpanded() const { 
    return oneLoopAmplitude()->isExpanded();
  }

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { 
    return oneLoopAmplitude()->mu2();
  }

  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const { 
    return oneLoopAmplitude()->oneLoopDoublePole();
  }

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const { 
    return oneLoopAmplitude()->oneLoopSinglePole();
  }

  /**
   * Calculate the one-loop amplitudes for the phasespace point
   * stored in lastXComb, if provided.
   */
  virtual void prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const {
    return oneLoopAmplitude()->oneLoopInterference();
  }

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluateOneLoop(size_t a, const vector<int>& hel) { 
    return oneLoopAmplitude()->evaluateOneLoop(a,hel);
  }

  //@}

  /** @name Caching and helpers to setup amplitude objects. */
  //@{

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {
    treeLevelAmplitude()->flushCaches();
    oneLoopAmplitude()->flushCaches();
  }

  //@}

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The amplitude object to provide tree-level amplitudes.
   */
  Ptr<MatchboxAmplitude>::ptr theTreeLevelAmplitude;

  /**
   * The amplitude object to provide one-loop amplitudes.
   */
  Ptr<MatchboxAmplitude>::ptr theOneLoopAmplitude;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxHybridAmplitude & operator=(const MatchboxHybridAmplitude &);

};

}

#endif /* Herwig_MatchboxHybridAmplitude_H */
