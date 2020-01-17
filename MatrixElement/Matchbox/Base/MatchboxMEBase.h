// -*- C++ -*-
//
// MatchboxMEBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxMEBase_H
#define HERWIG_MatchboxMEBase_H
//
// This is the declaration of the MatchboxMEBase class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Utility/Tree2toNGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"
#include "Herwig/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxReweightBase.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig/MatrixElement/Matchbox/Base/MergerBase.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"
#include "Herwig/MatrixElement/Matchbox/Utility/LastMatchboxXCombInfo.h"
#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxXComb.h"

namespace Herwig {

using namespace ThePEG;

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
class MatchboxMEBase: 
    public MEBase, public LastMatchboxXCombInfo {

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

  /**
   * Return the factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr factory() const;

  /** @name Subprocess and diagram information. */
  //@{

  /**
   * Return the subprocess.
   */
  const Process& subProcess() const { return theSubprocess; }

  /**
   * Access the subprocess.
   */
  Process& subProcess() { return theSubprocess; }

  /**
   * Return the diagram generator.
   */
  Ptr<Tree2toNGenerator>::tptr diagramGenerator() const;

  /**
   * Return the process data.
   */
  Ptr<ProcessData>::tptr processData() const;

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
   * Return true, if this amplitude is capable of consistently filling
   * the rho matrices for the spin correllations
   */
  virtual bool canFillRhoMatrix() const { 
    if ( matchboxAmplitude() )
      return matchboxAmplitude()->canFillRhoMatrix();
    return false;
  }

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr) {}

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub, const ColourLines* cl);

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
   * Return true, if this amplitude already includes averaging over
   * incoming parton's quantum numbers.
   */
  virtual bool hasInitialAverage() const { 
    return matchboxAmplitude() ? matchboxAmplitude()->hasInitialAverage() : false;
  }

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const { 
    return matchboxAmplitude() ? matchboxAmplitude()->hasFinalStateSymmetry() : false; 
  }


  /**
   * Return the number of light flavours, this matrix
   * element is calculated for.
   */
  virtual unsigned int getNLight() const;

  /**
   * Return the vector that contains the PDG ids of 
   * the light flavours, which are contained in the
   * jet particle group.
   */
  virtual vector<long> getNLightJetVec() const;

  /**
   * Return the vector that contains the PDG ids of 
   * the heavy flavours, which are contained in the
   * jet particle group.
   */
  virtual vector<long> getNHeavyJetVec() const;

  /**
   * Return the vector that contains the PDG ids of 
   * the light flavours, which are contained in the
   * proton particle group.
   */
  virtual vector<long> getNLightProtonVec() const;

  /**
   * Return true, if this matrix element is handled by a BLHA one-loop provider
   */
  virtual bool isOLPTree() const { 
    return matchboxAmplitude() ? matchboxAmplitude()->isOLPTree() : false;
  }

  /**
   * Return true, if this matrix element is handled by a BLHA one-loop provider
   */
  virtual bool isOLPLoop() const { 
    return matchboxAmplitude() ? matchboxAmplitude()->isOLPLoop() : false;
  }

  /**
   * Return true, if colour and spin correlated matrix elements should
   * be ordered from the OLP
   */
  virtual bool needsOLPCorrelators() const { 
    return matchboxAmplitude() ? matchboxAmplitude()->needsOLPCorrelators() : true;
  }

  /**
   * Return the process index, if this is an OLP handled matrix element
   */
  const vector<int>& olpProcess() const { return theOLPProcess; }

  /**
   * Set the process index, if this is an OLP handled matrix element
   */
  void olpProcess(int pType, int id) { 
    if ( theOLPProcess.empty() )
      theOLPProcess.resize(5,0);
    theOLPProcess[pType] = id;
  }

  /**
   * Return true, if this is a real emission matrix element which does
   * not require colour correlators.
   */
  bool noCorrelations() const {
    return theNoCorrelations;
  }

  /**
   * Indicate that this is a real emission matrix element which does
   * not require colour correlators.
   */
  void needsNoCorrelations() {
    theNoCorrelations = true;
  }

  /**
   * Indicate that this is a virtual matrix element which does
   * require colour correlators.
   */
  void needsCorrelations() {
    theNoCorrelations = false;
  }

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
   * Return true, if the XComb steering this matrix element
   * should keep track of the random numbers used to generate
   * the last phase space point
   */
  virtual bool keepRandomNumbers() const { return true; }

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
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * Clear the information previously provided by a call to
   * setKinematics(...).
   */
  virtual void clearKinematics();

  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element for generating a Born phase space point
   */
  virtual int nDimBorn() const;

  /**
   * Return true, if this matrix element will generate momenta for the
   * incoming partons itself.  The matrix element is required to store
   * the incoming parton momenta in meMomenta()[0,1]. No mapping in
   * tau and y is performed by the PartonExtractor object, if a
   * derived class returns true here. The phase space jacobian is to
   * include a factor 1/(x1 x2).
   */
  virtual bool haveX1X2() const { 
    return 
      (phasespace() ? phasespace()->haveX1X2() : false) ||
      diagrams().front()->partons().size() == 3;
  }

  /**
   * Return true, if this matrix element expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS() const { 
    return 
      (phasespace() ? phasespace()->wantCMS() : true) &&
      diagrams().front()->partons().size() != 3; }

  /**
   * Return the meMomenta as generated at the last
   * phase space point.
   */
  const vector<Lorentz5Momentum>& lastMEMomenta() const { return meMomenta(); }

  /**
   * Access the meMomenta.
   */
  vector<Lorentz5Momentum>& lastMEMomenta() { return meMomenta(); }
  
  
  /**
   * leg size
   */
  
  int legsize() const {return int(meMomenta().size());}

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
  void setScale(Energy2 ren=ZERO,Energy2 fac=ZERO) const;

  /**
   * Indicate that this matrix element is running alphas by itself.
   */
  virtual bool hasRunningAlphaS() const { 
    if ( matchboxAmplitude() )
      return matchboxAmplitude()->hasRunningAlphaS();
    return false;
  }

  /**
   * Indicate that this matrix element  is running alphaew by itself.
   */
  virtual bool hasRunningAlphaEW() const {
    if ( matchboxAmplitude() )
      return matchboxAmplitude()->hasRunningAlphaEW();
    return false;
  }

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
  virtual double factorizationScaleFactor() const;
      
      
  /**
    * Get the factorization scale factor
    */
  virtual double facFac() const{return factorizationScaleFactor();}

  /**
   * Return the (QCD) renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 renormalizationScale() const;

  /**
   * Get the renormalization scale factor
   */
  virtual double renormalizationScaleFactor() const;
      
      
  /**
   * Get the renormalization scale factor
   */
  virtual double renFac() const{return renormalizationScaleFactor();}

  /**
   * Return the QED renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 renormalizationScaleQED() const;

  /**
   * Return the shower scale for the last generated phasespace point.
   */
  virtual Energy2 showerScale() const;

  /**
   * Set veto scales on the particles at the given
   * SubProcess which has been generated using this
   * matrix element.
   */
  virtual void setVetoScales(tSubProPtr) const;

  /**
   * Return true, if fixed couplings are used.
   */
  bool fixedCouplings() const;

  /**
   * Return true, if fixed couplings are used.
   */
  bool fixedQEDCouplings() const;

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
  virtual bool havePDFWeight1() const;

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const;

  /**
   * Set the PDF weight.
   */
  void getPDFWeight(Energy2 factorizationScale = ZERO) const;

  /**
   * Supply the PDF weight for the first incoming parton.
   */
  double pdf1(Energy2 factorizationScale = ZERO,
	      double xEx = 1., double xFactor = 1.) const;

  /**
   * Supply the PDF weight for the second incoming parton.
   */
  double pdf2(Energy2 factorizationScale = ZERO,
	      double xEx = 1., double xFactor = 1.) const;

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
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double largeNME2(Ptr<ColourBasis>::tptr largeNBasis) const;

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

  /**
   * Same prefactor for all dSigHat
   **/
  CrossSection prefactor()const;

  /**
   * Born part of the cross section
   **/
  CrossSection dSigHatDRB() const ;
      
  /**
   * Virtual corrections of the cross section
   **/
  CrossSection dSigHatDRV() const ;
      
  /**
   * Insertion operators of the cross section
   **/
  CrossSection dSigHatDRI() const ;
      
  /**
   * If diffAlpha is not 1 and the matrix element has insertion operators 
   * this routine adds the difference between the insertion operator calculated 
   * with an alpha-Parameter to the insertion operator without alpha-parameter.
   */ 
  CrossSection dSigHatDRAlphaDiff(double alpha) const ;
      

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
   * Return true, if the amplitude is DRbar renormalized, otherwise
   * MSbar is assumed.
   */
  virtual bool isDRbar() const;

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

  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */
  virtual bool isBDK() const;

  /**
   * Return true, if one loop corrections are given in the conventions
   * of everything expanded.
   */
  virtual bool isExpanded() const;

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const;

  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const;

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const;

  /**
   * Return true, if cancellationn of epsilon poles should be checked.
   */
  bool checkPoles() const;

  /**
   * Simple histogram for accuracy checks
   */
  struct AccuracyHistogram {

    /**
     * The lower bound
     */
    double lower;

    /**
     * The upper bound
     */
    double upper;

    /**
     * The bins, indexed by upper bound.
     */
    map<double,double> bins;

    /**
     * The number of points of same sign
     */
    unsigned long sameSign;

    /**
     * The number of points of opposite sign
     */
    unsigned long oppositeSign;

    /**
     * The number of points being nan or inf
     */
    unsigned long nans;

    /**
     * The overflow
     */
    unsigned long overflow;

    /**
     * The underflow
     */
    unsigned long underflow;

    /**
     * Constructor
     */
    AccuracyHistogram(double low = -40.,
		      double up = 0.,
		      unsigned int nbins = 80);

    /**
     * Book two values to be checked for numerical compatibility
     */
    void book(double a, double b);

    /**
     * Write to file.
     */
    void dump(const std::string& folder, const std::string& prefix,
	      const cPDVector& proc) const;

    /**
     * Write to persistent ostream
     */
    void persistentOutput(PersistentOStream&) const;

    /**
     * Read from persistent istream
     */
    void persistentInput(PersistentIStream&);

  };

  /**
   * Perform the check of epsilon pole cancellation.
   */
  void logPoles() const;

  /**
   * Return the virtual corrections
   */
  const vector<Ptr<MatchboxInsertionOperator>::ptr>& virtuals() const {
    return theVirtuals;
  }

  /**
   * Return the virtual corrections
   */
  vector<Ptr<MatchboxInsertionOperator>::ptr>& virtuals() {
    return theVirtuals;
  }

  /**
   * Instruct this matrix element to include one-loop corrections
   */
  void doOneLoop() { theOneLoop = true; }
      
  /**
   * Instruct this matrix element not to include one-loop corrections
   */
      
  void noOneLoop() { theOneLoop = false; }

  /**
   * Return true, if this matrix element includes one-loop corrections
   */
  bool oneLoop() const { return theOneLoop; }

  /**
   * Instruct this matrix element to include one-loop corrections but
   * no Born contributions
   */
  void doOneLoopNoBorn() { theOneLoop = true; theOneLoopNoBorn = true; }
      
  void noOneLoopNoBorn() { theOneLoop = false; theOneLoopNoBorn = false; }

  /**
   * Return true, if this matrix element includes one-loop corrections
   * but no Born contributions
   */
  bool oneLoopNoBorn() const { return theOneLoopNoBorn || onlyOneLoop(); }

  /**
   * Instruct this matrix element to include one-loop corrections but
   * no actual loop contributions
   */
  void doOneLoopNoLoops() { theOneLoop = true; theOneLoopNoLoops = true; }

  /**
   * Return true, if this matrix element includes one-loop corrections
   * but no actual loop contributions
   */
  bool oneLoopNoLoops() const { return theOneLoopNoLoops; }

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
	     const vector<Ptr<MatchboxMEBase>::ptr>&,bool slim=false) const;

  
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
   * Return the colour correlated matrix element squared in the
   * large-N approximation with respect to the given two partons as
   * appearing in mePartonData(), suitably scaled by sHat() to give a
   * dimension-less number.
   */
  virtual double largeNColourCorrelatedME2(pair<int,int> ij,
					   Ptr<ColourBasis>::tptr largeNBasis) const;

  /**
   * Return the colour and spin correlated matrix element squared for
   * the gluon indexed by the first argument using the given
   * correlation tensor.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;

  /**
   * Return the spin correlated matrix element squared for
   * the vector boson indexed by the first argument using the given
   * correlation tensor.
   */
  virtual double spinCorrelatedME2(pair<int,int> emitterSpectator,
				   const SpinCorrelationTensor& c) const;

  //@}

  /** @name Caching and diagnostic information */
  //@{

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches();

  /**
   * Return true, if verbose
   */
  bool verbose() const;

  /**
   * Return true, if verbose
   */
  bool initVerbose() const;

  /**
   * Dump the setup to an ostream
   */
  void print(ostream&) const;

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

  /**
   * Return the theMerger.
   */
  const MergerBasePtr merger() const;
    
  /**
   * Return the theMerger.
   */
  MergerBasePtr merger() ;
      
  /**
   * Set the theMerger.
   */
  void merger(MergerBasePtr v);
  
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
  void cloneDependencies(const std::string& prefix = "",bool slim = false );

  /**
   * Prepare an xcomb
   */
  void prepareXComb(MatchboxXCombData&) const;

  /**
   * For the given event generation setup return a xcomb object
   * appropriate to this matrix element.
   */
  virtual StdXCombPtr makeXComb(Energy newMaxEnergy, const cPDPair & inc,
				tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
				tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
				const PBPair & newPartonBins, tCutsPtr newCuts,
				const DiagramVector & newDiagrams, bool mir,
				const PartonPairVec& allPBins,
				tStdXCombPtr newHead = tStdXCombPtr(),
				tMEPtr newME = tMEPtr());

  /**
   * For the given event generation setup return a dependent xcomb object
   * appropriate to this matrix element.
   */
  virtual StdXCombPtr makeXComb(tStdXCombPtr newHead,
				const PBPair & newPartonBins,
				const DiagramVector & newDiagrams,
				tMEPtr newME = tMEPtr());

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The phase space generator to be used.
   */
  Ptr<MatchboxPhasespace>::ptr thePhasespace;

  /**
   * The amplitude to be used
   */
  Ptr<MatchboxAmplitude>::ptr theAmplitude;

  /**
   * The scale choice object
   */
  Ptr<MatchboxScaleChoice>::ptr theScaleChoice;

  /**
   * The virtual corrections.
   */
  vector<Ptr<MatchboxInsertionOperator>::ptr> theVirtuals;

  /**
   * A vector of reweight objects the sum of which
   * should be applied to reweight this matrix element
   */
  vector<Ptr<MatchboxReweightBase>::ptr> theReweights;

private:

  /**
   * The subprocess to be considered.
   */
  Process theSubprocess;

  /**
   * True, if this matrix element includes one-loop corrections
   */
  bool theOneLoop;

  /**
   * True, if this matrix element includes one-loop corrections
   * but no Born contributions
   */
  bool theOneLoopNoBorn;

  /**
   * True, if this matrix element includes one-loop corrections
   * but no actual loop contributions (e.g. finite collinear terms)
   */
  bool theOneLoopNoLoops;

  /**
   * The process index, if this is an OLP handled matrix element
   */
  vector<int> theOLPProcess;

  /**
   * Histograms of epsilon^2 pole cancellation
   */
  mutable map<cPDVector,AccuracyHistogram> epsilonSquarePoleHistograms;

  /**
   * Histograms of epsilon pole cancellation
   */
  mutable map<cPDVector,AccuracyHistogram> epsilonPoleHistograms;

  /**
   * True, if this is a real emission matrix element which does
   * not require colour correlators.
   */
  bool theNoCorrelations;

  /**
   * Flag which pdfs should be included.
   */
  mutable pair<bool,bool> theHavePDFs;

  /**
   * True, if already checked for which PDFs to include.
   */
  mutable bool checkedPDFs;
  
  /**
   * The merging helper to be used. 
   * Only the head ME has a pointer to this helper.
   */

  MergerBasePtr theMerger;


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxMEBase & operator=(const MatchboxMEBase &) = delete;

};

inline PersistentOStream& operator<<(PersistentOStream& os,
				     const MatchboxMEBase::AccuracyHistogram& h) {
  h.persistentOutput(os);
  return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is,
				     MatchboxMEBase::AccuracyHistogram& h) {
  h.persistentInput(is);
  return is;
}

}

#endif /* HERWIG_MatchboxMEBase_H */
