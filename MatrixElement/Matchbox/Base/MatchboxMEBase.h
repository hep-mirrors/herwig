// -*- C++ -*-
//
// MatchboxMEBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEBase_H
#define HERWIG_MatchboxMEBase_H
//
// This is the declaration of the MatchboxMEBase class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/Tree2toNGenerator.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/MatchboxMECache.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxReweightBase.h"

namespace Herwig {

using namespace ThePEG;

class SubtractionDipole;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Keys for XComb meta information
 */
struct MatchboxMetaKeys {

  enum Keys {
    BornME,
    ColourCorrelatedMEs
  };

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxMEBase is the base class for matrix elements
 * in the context of the matchbox NLO interface.
 *
 * @see \ref MatchboxMEBaseInterfaces "The interfaces"
 * defined for MatchboxMEBase.
 */
class MatchboxMEBase: public MEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxMEBase();

  /**
   * The destructor.
   */
  virtual ~MatchboxMEBase();
  //@}

public:

  /** @name Subprocess and diagram information. */
  //@{

  /**
   * Return the subprocesses.
   */
  const vector<PDVector>& subProcesses() const { return theSubprocesses; }

  /**
   * Access the subprocesses.
   */
  vector<PDVector>& subProcesses() { return theSubprocesses; }

  /**
   * Return the diagram generator.
   */
  Ptr<Tree2toNGenerator>::tptr diagramGenerator() const { return theDiagramGenerator; }

  /**
   * Set the diagram generator.
   */
  void diagramGenerator(Ptr<Tree2toNGenerator>::ptr gen) { theDiagramGenerator = gen; }

  /**
   * Return true, if this matrix element does not want to
   * make use of mirroring processes; in this case all
   * possible partonic subprocesses with a fixed assignment
   * of incoming particles need to be provided through the diagrams
   * added with the add(...) method.
   */
  virtual bool noMirror () const { return true; }

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;
  using MEBase::getDiagrams;

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;
  using MEBase::diagrams;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const;
  using MEBase::orderInAlphaS;

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const;
  using MEBase::orderInAlphaEW;  

  /**
   * Return the number of light flavours, this matrix
   * element is calculated for.
   */
  virtual unsigned int nLight() const { return theNLight; }

  /**
   * Set the number of light flavours, this matrix
   * element is calculated for.
   */
  void nLight(unsigned int n) { theNLight = n; }

  //@}

  /** @name Phasespace generation */
  //@{

  /**
   * Return the phase space generator to be used.
   */
  Ptr<MatchboxPhasespace>::tptr phasespace() const { return thePhasespace; }

  /**
   * Set the phase space generator to be used.
   */
  void phasespace(Ptr<MatchboxPhasespace>::ptr ps) { thePhasespace = ps; }

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr xc);

  /**
   * Generate incoming parton momenta. This default
   * implementation performs the standard mapping
   * from x1,x2 -> tau,y making 1/tau flat; incoming
   * parton momenta are stored in meMomenta()[0,1],
   * only massless partons are supported so far;
   * return the Jacobian of the mapping
   */
  double generateIncomingPartons(const double* r1, const double* r2);

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector. Derived classes
   * must call this method once internal degrees of freedom are setup
   * and finally return the result of this method.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Return true, if this matrix element will generate momenta for the
   * incoming partons itself.  The matrix element is required to store
   * the incoming parton momenta in meMomenta()[0,1]. No mapping in
   * tau and y is performed by the PartonExtractor object, if a
   * derived class returns true here. The phase space jacobian is to
   * include a factor 1/(x1 x2).
   */
  virtual bool haveX1X2() const { return (phasespace() ? phasespace()->haveX1X2() : false); }

  /**
   * Return true, if this matrix element expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS() const { return phasespace() ? phasespace()->wantCMS() : true; }

  /**
   * Return true, if the XComb steering this matrix element
   * should keep track of the random numbers used to generate
   * the last phase space point
   */
  virtual bool keepRandomNumbers() const { return true; }

  /**
   * Return the meMomenta as generated at the last
   * phase space point.
   */
  const vector<Lorentz5Momentum>& lastMEMomenta() const { return meMomenta(); }

  /**
   * Access the meMomenta.
   */
  vector<Lorentz5Momentum>& lastMEMomenta() { return meMomenta(); }

  //@}

  /** @name Scale choices, couplings and PDFs */
  //@{

  /**
   * Set the scale choice object
   */
  void scaleChoice(Ptr<MatchboxScaleChoice>::ptr sc) { theScaleChoice = sc; }

  /**
   * Return the scale choice object
   */
  Ptr<MatchboxScaleChoice>::tptr scaleChoice() const { return theScaleChoice; }

  /**
   * Set scales and alphaS
   */
  void setScale() const;

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const { return lastScale(); }

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Get the factorization scale factor
   */
  virtual double factorizationScaleFactor() const { return theFactorizationScaleFactor; }

  /**
   * Set the factorization scale factor
   */
  void factorizationScaleFactor(double f) { theFactorizationScaleFactor = f; }

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Get the renormalization scale factor
   */
  virtual double renormalizationScaleFactor() const { return theRenormalizationScaleFactor; }

  /**
   * Set the renormalization scale factor
   */
  void renormalizationScaleFactor(double f) { theRenormalizationScaleFactor = f; }

  /**
   * Set veto scales on the particles at the given
   * SubProcess which has been generated using this
   * matrix element.
   */
  virtual void setVetoScales(tSubProPtr) const;

  /**
   * Return true, if fixed couplings are used.
   */
  bool fixedCouplings() const { return theFixedCouplings; }

  /**
   * Switch on fixed couplings.
   */
  void setFixedCouplings(bool on = true) { theFixedCouplings = on; }

  /**
   * Return the value of \f$\alpha_S\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaS(scale()).
   */
  virtual double alphaS() const { return lastAlphaS(); }

  /**
   * Return the value of \f$\alpha_EM\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaEM(scale()).
   */
  virtual double alphaEM() const { return lastAlphaEM(); }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { 
    return diagrams().front()->partons()[0]->coloured();
  }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const { 
    return diagrams().front()->partons()[1]->coloured();
  }

  /**
   * Set the PDF weight.
   */
  void getPDFWeight(Energy2 factorizationScale = ZERO) const;

  /**
   * Supply the PDF weight for the first incoming parton.
   */
  double pdf1(Energy2 factorizationScale = ZERO) const;

  /**
   * Supply the PDF weight for the second incoming parton.
   */
  double pdf2(Energy2 factorizationScale = ZERO) const;

  //@}

  /** @name Amplitude information and matrix element evaluation */
  //@{

  /**
   * Return the amplitude.
   */
  Ptr<MatchboxAmplitude>::tptr matchboxAmplitude() const { return theAmplitude; }

  /**
   * Set the amplitude.
   */
  void matchboxAmplitude(Ptr<MatchboxAmplitude>::ptr amp) { theAmplitude = amp; }

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

  /**
   * Return the symmetry factor for identical final state particles.
   */
  virtual double finalStateSymmetry() const;

  /**
   * Return the normalizing factor for the matrix element averaged
   * over quantum numbers and including running couplings.
   */
  double me2Norm(unsigned int addAlphaS = 0) const;

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  //@}

  /** @name One-loop corrections */
  //@{

  /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const;

  /**
   * Return true, if this matrix element is capable of calculating
   * one-loop (QCD) corrections.
   */
  virtual bool haveOneLoop() const;

  /**
   * Return true, if this matrix element only provides
   * one-loop (QCD) corrections.
   */
  virtual bool onlyOneLoop() const;

  /**
   * Return true, if one loop corrections have been calculated in
   * dimensional reduction. Otherwise conventional dimensional
   * regularization is assumed. Note that renormalization is always
   * assumed to be MSbar.
   */
  virtual bool isDR() const;

  /**
   * Return true, if one loop corrections are given in the conventions
   * of the integrated dipoles.
   */
  virtual bool isCS() const;

  //@}

  /** @name Dipole subtraction */
  //@{

  /**
   * If this matrix element is considered a real
   * emission matrix element, return all subtraction
   * dipoles needed given a set of subtraction terms
   * and underlying Born matrix elements to choose
   * from.
   */
  vector<Ptr<SubtractionDipole>::ptr> 
  getDipoles(const vector<Ptr<SubtractionDipole>::ptr>&,
	     const vector<Ptr<MatchboxMEBase>::ptr>&) const;

  /**
   * If this matrix element is considered a real emission matrix
   * element, but actually neglecting a subclass of the contributing
   * diagrams, return true if the given emitter-emission-spectator
   * configuration should not be considered when setting up
   * subtraction dipoles.
   */
  virtual bool noDipole(int,int,int) const { return false; }

  /**
   * If this matrix element is considered an underlying Born matrix
   * element in the context of a subtracted real emission, but
   * actually neglecting a subclass of the contributing diagrams,
   * return true if the given emitter-spectator configuration
   * should not be considered when setting up subtraction dipoles.
   */
  virtual bool noDipole(int,int) const { return false; }

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
   * Return true, if colour correlated matrix elements should be calculated
   * for later use
   */
  bool getColourCorrelatedMEs() const { return theGetColourCorrelatedMEs; }

  /**
   * Switch on calculation of colour correlated matrix elements for
   * later use
   */
  void doColourCorrelatedMEs() { theGetColourCorrelatedMEs = true; }

  /**
   * Calculate colour correlated matrix elements for later use
   */
  void storeColourCorrelatedMEs(double xme2 = -1.) const;

  //@}

  /** @name Caching and diagnostic information */
  //@{

  /**
   * Set the ME cache object
   */
  void cache(Ptr<MatchboxMECache>::ptr c) { theCache = c; }

  /**
   * Get the ME cache object
   */
  Ptr<MatchboxMECache>::tptr cache() const { return theCache; }

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches();

  /**
   * Return true, if the matrix element needs to be 
   * recalculated for the given phase space point.
   * If not, return the cached value in the given reference.
   */
  bool calculateME2(double& xme2,
		    const pair<int,int>& corr = make_pair(0,0)) const {
    if ( !cache() ) {
      xme2 = 0.;
      return true;
    }
    cache()->setXComb(lastXCombPtr());
    return cache()->calculateME2(xme2,corr);
  }

  /**
   * Cache a calculated matrix element
   * for the last phase space point.
   */
  void cacheME2(double xme2,
		const pair<int,int>& corr = make_pair(0,0)) const {
    if ( !cache() )
      return;
    cache()->cacheME2(xme2,corr);
  }

  /**
   * Return true, if verbose
   */
  bool verbose() const { return theVerbose; }

  /**
   * Switch on diagnostic information.
   */
  void setVerbose(bool on = true) { theVerbose = on; }

  /**
   * Print debug information on the last event
   */
  virtual void printLastEvent(ostream&) const;

  /**
   * Write out diagnostic information for
   * generateKinematics
   */
  void logGenerateKinematics(const double * r) const;

  /**
   * Write out diagnostic information for
   * setting scales
   */
  void logSetScale() const;

  /**
   * Write out diagnostic information for
   * pdf evaluation
   */
  void logPDFWeight() const;

  /**
   * Write out diagnostic information for
   * me2 evaluation
   */
  void logME2() const;

  /**
   * Write out diagnostic information
   * for dsigdr evaluation
   */
  void logDSigHatDR() const;

  /**
   * Dump xcomb hierarchies.
   */
  void dumpInfo(const string& prefix = "") const;

  //@}

  /** @name Reweight objects */
  //@{

  /**
   * Insert a reweight object
   */
  void addReweight(Ptr<MatchboxReweightBase>::ptr rw) { theReweights.push_back(rw); }

  /**
   * Return the reweights
   */
  const vector<Ptr<MatchboxReweightBase>::ptr>& reweights() const { return theReweights; }

  /**
   * Access the reweights
   */
  vector<Ptr<MatchboxReweightBase>::ptr>& reweights() { return theReweights; }

  //@}

  /** @name Methods used to setup MatchboxMEBase objects */
  //@{

  /**
   * Return true if this object needs to be initialized before all
   * other objects (except those for which this function also returns
   * true).  This default version always returns false, but subclasses
   * may override it to return true.
   */
  virtual bool preInitialize() const { return true; }

  /**
   * Clone this matrix element.
   */
  Ptr<MatchboxMEBase>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
  void cloneDependencies(const std::string& prefix = "");

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

  /**
   * The final state symmetry factors.
   */
  mutable map<tStdXCombPtr,double> symmetryFactors;

private:

  /**
   * A vector of reweight objects the sum of which
   * should be applied to reweight this matrix element
   */
  vector<Ptr<MatchboxReweightBase>::ptr> theReweights;

  /**
   * The phase space generator to be used.
   */
  Ptr<MatchboxPhasespace>::ptr thePhasespace;

  /**
   * The amplitude to be used
   */
  Ptr<MatchboxAmplitude>::ptr theAmplitude;

  /**
   * The diagram generator to be used.
   */
  Ptr<Tree2toNGenerator>::ptr theDiagramGenerator;

  /**
   * The scale choice object
   */
  Ptr<MatchboxScaleChoice>::ptr theScaleChoice;

  /**
   * The ME cache object
   */
  Ptr<MatchboxMECache>::ptr theCache;

  /**
   * The subprocesses to be considered, if a diagram generator is
   * present.
   */
  vector<PDVector> theSubprocesses;

  /**
   * The factorization scale factor.
   */
  double theFactorizationScaleFactor;

  /**
   * The renormalization scale factor.
   */
  double theRenormalizationScaleFactor;

  /**
   * Wether or not diagnostic information
   * should be written to the generator log
   */
  bool theVerbose;

  /**
   * Use non-running couplings.
   */
  bool theFixedCouplings;

  /**
   * The number of light flavours, this matrix
   * element is calculated for.
   */
  unsigned int theNLight;

  /**
   * True, if colour correlated matrix elements should be calculated
   * for later use
   */
  bool theGetColourCorrelatedMEs;

  /**
   * Map xcomb's to storage of Born matrix elements squared.
   */
  mutable map<tStdXCombPtr,double> bornMEs;

  /**
   * Map xcomb's to storage of colour correlated matrix elements.
   */
  mutable map<tStdXCombPtr,map<pair<int,int>,double> > colourCorrelatedMEs;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMEBase & operator=(const MatchboxMEBase &);

};

}

#endif /* HERWIG_MatchboxMEBase_H */
