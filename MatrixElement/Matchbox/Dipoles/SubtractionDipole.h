// -*- C++ -*-
//
// SubtractionDipole.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SubtractionDipole_H
#define HERWIG_SubtractionDipole_H
//
// This is the declaration of the SubtractionDipole class.
//

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.fh"
#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.fh"

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximation.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief SubtractionDipole represents a dipole subtraction
 * term in the formalism of Catani and Seymour.
 *
 */
class SubtractionDipole: 
    public MEBase,  public LastMatchboxXCombInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SubtractionDipole();

  /**
   * The destructor.
   */
  virtual ~SubtractionDipole();
  //@}

public:

  /**
   * Return the factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr factory() const;

  /**
   * Set the factory which produced this matrix element
   */
  void factory(Ptr<MatchboxFactory>::tptr f);

  /** @name Subprocess and diagram information. */
  //@{

  /**
   * A helpre struct to communicate diagram merging and remapping
   * information
   */
  struct MergeInfo {

    /**
     * The merged emitter
     */
    int emitter;

    /**
     * The Born diagram
     */
    Ptr<Tree2toNDiagram>::ptr diagram;

    /**
     * The merging map
     */
    map<int,int> mergeLegs;

  };

  /**
   * Return true, if this dipole can possibly handle the indicated
   * emitter.
   */
  virtual bool canHandleEmitter(const cPDVector& partons, int emitter) const = 0;

  /**
   * Return true, if this dipole can possibly handle the indicated
   * splitting.
   */
  virtual bool canHandleSplitting(const cPDVector& partons, int emitter, int emission) const = 0;

  /**
   * Return true, if this dipole can possibly handle the indicated
   * spectator.
   */
  virtual bool canHandleSpectator(const cPDVector& partons, int spectator) const = 0;

  /**
   * Return true, if this dipole applies to the selected
   * configuration.
   */
  virtual bool canHandle(const cPDVector& partons,
			 int emitter, int emission, int spectator) const = 0;

  /**
   * Return true, if this dipole is symmetric with respect to emitter
   * and emission.
   */
  virtual bool isSymmetric() const { return false; }

  /**
   * If this is a dependent matrix element in a ME group, return true,
   * if it applies to the process set in lastXComb()
   */
  virtual bool apply() const { return theApply; }

  /**
   * Clear the bookkeeping
   */
  void clearBookkeeping();

  /**
   * Setup bookkeeping maps.
   */
  void setupBookkeeping(const map<Ptr<DiagramBase>::ptr,MergeInfo>& mergeInfo,bool slim);

  /**
   * Get bookkeeping information for the given
   * real emission diagram
   */
  void subtractionBookkeeping();

  /**
   * Determine bookkeeping information for
   * the underlying Born process supplied through
   * the lastHeadXComb() object.
   */
  void splittingBookkeeping();

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

  /**
   * Create a dependent xcomb object for the underlying
   * Born process, given a XComb driving the real emission
   */
  StdXCombPtr makeBornXComb(tStdXCombPtr realXC);

  /**
   * Create dependent xcomb objects for the real emission process,
   * given a XComb driving the underlying Born
   */
  vector<StdXCombPtr> makeRealXCombs(tStdXCombPtr bornXC);

  /**
   * Return true, if bookkeeping did not find a non-trivial setup.
   */
  bool empty() const { return theSplittingMap.empty()&&theMergingMap.empty(); }

  /**
   * Return the emitter as referred to by the real emission
   * matrix element.
   */
  int realEmitter() const { return theRealEmitter; }

  /**
   * Set the emitter as referred to by the real emission
   * matrix element.
   */
  void realEmitter(int id) { theRealEmitter = id; }

  /**
   * Return the emission as referred to by the real emission
   * matrix element.
   */
  int realEmission() const { return theRealEmission; }

  /**
   * Set the emission as referred to by the real emission
   * matrix element.
   */
  void realEmission(int id) { theRealEmission = id; }

  /**
   * Return the spectator as referred to by the real emission
   * matrix element.
   */
  int realSpectator() const { return theRealSpectator; }

  /**
   * Set the spectator as referred to by the real emission
   * matrix element.
   */
  void realSpectator(int id) { theRealSpectator = id; }

  /**
   * Return the emitter as referred to by the underlying
   * Born process.
   */
  int bornEmitter() const { return theBornEmitter; }

  /**
   * Set the emitter as referred to by the underlying
   * Born process.
   */
  void bornEmitter(int id) { theBornEmitter = id; }

  /**
   * Return the spectator as referred to by the underlying
   * Born process.
   */
  int bornSpectator() const { return theBornSpectator; }

  /**
   * Set the spectator as referred to by the underlying
   * Born process.
   */
  void bornSpectator(int id) { theBornSpectator = id; }

  /**
   * Define the real emission key type
   */
  typedef pair<pair<cPDVector,int>,pair<int,int> > RealEmissionKey;

  /**
   * Create a real emission key
   */
  static RealEmissionKey realEmissionKey(const cPDVector& proc, 
					 int em, int emm, int sp) {
    return make_pair(make_pair(proc,emm),make_pair(em,sp));
  }

  /**
   * Return the diagram of a real emission key
   */
  static const cPDVector& process(const RealEmissionKey& key) {
    return key.first.first;
  }

  /**
   * Return the emission id of a real emission key
   */
  static int emission(const RealEmissionKey& key) {
    return key.first.second;
  }

  /**
   * Return the emitter id of a real emission key
   */
  static int emitter(const RealEmissionKey& key) {
    return key.second.first;
  }

  /**
   * Return the spectator id of a real emission key
   */
  static int spectator(const RealEmissionKey& key) {
    return key.second.second;
  }

  /**
   * Define the underlying Born key type
   */
  typedef pair<cPDVector,pair<int,int> > UnderlyingBornKey;

  /**
   * Create a underlying Born key
   */
  static UnderlyingBornKey underlyingBornKey(const cPDVector& proc, 
					     int em, int sp) {
    return make_pair(proc,make_pair(em,sp));
  }

  /**
   * Return the diagram of a underlying Born key
   */
  static const cPDVector& process(const UnderlyingBornKey& key) {
    return key.first;
  }

  /**
   * Return the emitter id of a underlying Born key
   */
  static int emitter(const UnderlyingBornKey& key) {
    return key.second.first;
  }

  /**
   * Return the spectator id of a underlying Born key
   */
  static int spectator(const UnderlyingBornKey& key) {
    return key.second.second;
  }

  /**
   * Define real emission key and index dictionary
   * for partons not involved in the given dipole.
   */
  typedef pair<RealEmissionKey,map<int,int> > RealEmissionInfo;

  /**
   * Define underlying Born key and index dictionary
   * for partons not involved in the given dipole.
   */
  typedef pair<UnderlyingBornKey,map<int,int> > UnderlyingBornInfo;

  /**
   * Return the merging map
   */
  const map<RealEmissionKey,UnderlyingBornInfo>& mergingMap() const { return theMergingMap; }

  /**
   * Return the splitting map
   */
  const multimap<UnderlyingBornKey,RealEmissionInfo>& splittingMap() const { return theSplittingMap; }

  /**
   * Return the underlying Born diagrams to be considered
   * for the given real emission process.
   */
  const DiagramVector& underlyingBornDiagrams(const cPDVector& real) const;

  /**
   * Find the underlying Born diagram for the given real emission diagram
   */
  tcDiagPtr underlyingBornDiagram(tcDiagPtr realDiag) const;

  /**
   * Return the real emission diagrams to be considered
   * for the given Born process.
   */
  const DiagramVector& realEmissionDiagrams(const cPDVector& born) const;

  /**
   * Find the real emission diagram for the given underlying Born diagram
   */
  tcDiagPtr realEmissionDiagram(tcDiagPtr bornDiag) const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Return true, if this matrix element does not want to
   * make use of mirroring processes; in this case all
   * possible partonic subprocesses with a fixed assignment
   * of incoming particles need to be provided through the diagrams
   * added with the add(...) method.
   */
  virtual bool noMirror () const { return true; }

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * Select a ColpurLines geometry. The default version returns a
   * colour geometry selected among the ones returned from
   * colourGeometries(tcDiagPtr).
   */
  virtual const ColourLines &
  selectColourGeometry(tcDiagPtr diag) const;

  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix element
   * is given.
   */
  virtual unsigned int orderInAlphaS() const { return realEmissionME()->orderInAlphaS(); }

  /**
   * Return the order in \f$\alpha_{EM}\f$ in which this matrix
   * element is given. Returns 0.
   */
  virtual unsigned int orderInAlphaEW() const { return underlyingBornME()->orderInAlphaEW(); }

  //@}

  /** @name Phasespace generation */
  //@{

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr xc);

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * Generate internal degrees of freedom given nDim() uniform random
   * numbers in the interval ]0,1[. To help the phase space generator,
   * the 'dSigHatDR' should be a smooth function of these numbers,
   * although this is not strictly necessary. The return value should
   * be true of the generation succeeded. If so the generated momenta
   * should be stored in the meMomenta() vector.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * The number of internal degreed of freedom used in the matrix
   * element. This default version returns 0;
   */
  virtual int nDim() const;

  /**
   * Return true, if this matrix element expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS () const { return realEmissionME()->wantCMS(); }

  /**
   * Clear the information previously provided by a call to
   * setKinematics(...).
   */
  virtual void clearKinematics();

  /**
   * If this is a dependent matrix element in a ME group, return true,
   * if cuts should be ignored.
   */
  virtual bool ignoreCuts() const { return theIgnoreCuts; }

  /**
   * Indicate that cuts should be ignored
   */
  void doIgnoreCuts(bool is = true) { theIgnoreCuts = is; }

  //@}

  /** @name Tilde kinematics */
  //@{

  /**
   * Return the TildeKinematics object used
   */
  Ptr<TildeKinematics>::tcptr tildeKinematics() const { return theTildeKinematics; }

  /**
   * Set the TildeKinematics object used
   */
  void tildeKinematics(Ptr<TildeKinematics>::tptr);

  /**
   * Generate the tilde kinematics from real emission
   * kinematics accessible through the XComb's
   * head object and store it in meMomenta(). This default
   * implemenation uses the tildeKinematics() object.
   */
  virtual bool generateTildeKinematics();

  /**
   * Return the InvertedTildeKinematics object used
   */
  Ptr<InvertedTildeKinematics>::tcptr invertedTildeKinematics() const { return theInvertedTildeKinematics; }

  /**
   * Set the InvertedTildeKinematics object used
   */
  void invertedTildeKinematics(Ptr<InvertedTildeKinematics>::tptr);

  /**
   * Return the number of additional random numbers
   * needed to generate real emission kinematics off
   * the tilde kinematics previously supplied through
   * the XComb object. This default implementation
   * returns invertedTildeKinematics()->nDimRadiation()
   */
  virtual int nDimRadiation() const;

  /**
   * Generate the real emission kinematics
   * off the Born kinematics accessible through the XComb's
   * head object and store it in meMomenta(); store
   * the single particle phasespace in units of lastHeadXComb()->lastSHat()
   * in jacobian(). This default
   * implemenation uses the invertedTildeKinematics() object
   */
  virtual bool generateRadiationKinematics(const double *);

  /**
   * Set a pt cut when splitting
   */
  void ptCut(Energy cut);

  /**
   * Return the relevant dipole scale
   */
  Energy lastDipoleScale() const {
    return splitting() ? theLastSplittingScale : theLastSubtractionScale;
  }

  /**
   * Return the relevant pt
   */
  Energy lastPt() const {
    return splitting() ? theLastSplittingPt : theLastSubtractionPt;
  }

  /**
   * Return the relevant momentum fractions
   */
  double lastZ() const {
    return splitting() ? theLastSplittingZ : theLastSubtractionZ;
  }

  /**
   * Return true, if this dipole acts in splitting mode.
   */
  bool splitting() const { return theSplitting; }

  /**
   * Switch on splitting mode for this dipole.
   */  
  void doSplitting() { theSplitting = true; }

  /**
   * Switch off splitting mode for this dipole.
   */  
  void doSubtraction() { theSplitting = false; }

  /**
   * Return the subtraction parameters.
   */
  const vector<double>& subtractionParameters() const { return theSubtractionParameters; }

  /**
   * Access the subtraction parameters.
   */
  vector<double>& subtractionParameters() { return theSubtractionParameters; }

  /**
   * Return the shower hard scale encountered
   */
  Energy showerHardScale() const { return theShowerHardScale; }

  /**
   * Set the shower hard scale encountered
   */
  void showerHardScale(Energy s) { theShowerHardScale = s; }

  /**
   * Return the shower evolution scale encountered
   */
  Energy showerScale() const { return theShowerScale; }

  /**
   * Set the shower evolution scale encountered
   */
  void showerScale(Energy s) { theShowerScale = s; }

  /**
   * Return the shower splitting variables encountered
   */
  const vector<double>& showerParameters() const { return theShowerParameters; }

  /**
   * Access the shower splitting variables encountered
   */
  vector<double>& showerParameters() { return theShowerParameters; }

  /**
   * Return true, if this configuration is in the shower phase space
   */
  bool isInShowerPhasespace() const { return theIsInShowerPhasespace; }

  /**
   * Indicate whether this configuration is in the shower phase space
   */
  void isInShowerPhasespace(bool yes) { theIsInShowerPhasespace = yes; }

  /**
   * Return true, if this configuration is above the shower infrared cutoff
   */
  bool isAboveCutoff() const { return theIsAboveCutoff; }

  /**
   * Indicate whether this configuration is above the shower infrared cutoff
   */
  void isAboveCutoff(bool yes) { theIsAboveCutoff = yes; }

  //@}

  /** @name Scale choices, couplings and PDFs */
  //@{

  /**
   * Return true, if scales should be calculated from real emission kinematics
   */
  bool realEmissionScales() const { return theRealEmissionScales; }

  /**
   * Switch on or off that scales should be calculated from real emission kinematics
   */
  void doRealEmissionScales(bool on = true) { theRealEmissionScales = on; }

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const { 
    return realEmissionScales() ? 
      realEmissionME()->scale() :
      underlyingBornME()->scale();
  }

  /**
   * Return the value of \f$\alpha_S\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaS(scale()).
   */
  virtual double alphaS() const { 
    return realEmissionScales() ? 
      realEmissionME()->alphaS() :
      underlyingBornME()->alphaS();
  }

  /**
   * Return the value of \f$\alpha_EM\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaEM(scale()).
   */
  virtual double alphaEM() const { 
    return realEmissionScales() ? 
      realEmissionME()->alphaEM() :
      underlyingBornME()->alphaEM();
  }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the first incoming parton itself.
   */
  virtual bool havePDFWeight1() const { return realEmissionME()->havePDFWeight1(); }

  /**
   * Return true, if this matrix element provides the PDF
   * weight for the second incoming parton itself.
   */
  virtual bool havePDFWeight2() const { return realEmissionME()->havePDFWeight2(); }

  /**
   *  How to sample the z-distribution.
   *  FlatZ = 1
   *  OneOverZ = 2
   *  OneOverOneMinusZ = 3
   *  OneOverZOneMinusZ = 4
   */

  virtual int samplingZ() const {return 4;}
  //@}

  /** @name Matrix elements and evaluation */
  //@{

  /**
   * Return the real emission matrix element
   */
  Ptr<MatchboxMEBase>::tcptr realEmissionME() const { 
    return theRealEmissionME;
  }

  /**
   * Return the real emission matrix element
   */
  Ptr<MatchboxMEBase>::tptr realEmissionME() { 
    return theRealEmissionME;
  }

  /**
   * Set the real emission matrix element
   */
  void realEmissionME(Ptr<MatchboxMEBase>::tptr me) { theRealEmissionME = me; }

  /**
   * Return the underlying Born matrix element
   */
  Ptr<MatchboxMEBase>::tcptr underlyingBornME() const { 
    return theUnderlyingBornME;
  }

  /**
   * Return the underlying Born matrix element
   */
  Ptr<MatchboxMEBase>::tptr underlyingBornME() { 
    return theUnderlyingBornME;
  }

  /**
   * Set the underlying Born matrix element
   */
  void underlyingBornME(Ptr<MatchboxMEBase>::tptr me) { theUnderlyingBornME = me; }

  /**
   * Set the dipoles which have been found along with this dipole
   */
  void partnerDipoles(const vector<Ptr<SubtractionDipole>::tptr>& p) {
    thePartners = p;
  }

  /**
   * Return the dipoles which have been found along with this dipole
   */
  const vector<Ptr<SubtractionDipole>::tptr>& partnerDipoles() const {
    return thePartners;
  }

  /**
   * Return the matrix element averaged over spin correlations.
   */
  virtual double me2Avg(double ccme2) const = 0;

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR(Energy2 factorizationScale) const;

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const { return dSigHatDR(ZERO); }

      
        /// calculate the general prefactor for merging.
  CrossSection prefactor(Energy2 factorizationScale)const;
      
  /**
   *  Calculate the parton shower approximation for this dipole.
   **/
      
  CrossSection ps(Energy2 factorizationScale,Ptr<ColourBasis>::tptr largeNBasis) const;

  /**
   *  Calculate the dipole with clusterfsafe flag.
   **/

  CrossSection dip(Energy2 factorizationScale) const;

  

  /**
   *  Calculate the dipole dSigDR and the parton shower approximation for this dipole.
   **/

  pair<CrossSection,CrossSection> dipandPs(Energy2 factorizationScale,Ptr<ColourBasis>::tptr largeNBasis) const;

  //@}

  /** @name Methods relevant to matching */
  //@{

  /**
   * Set the shower approximation.
   */
  void showerApproximation(Ptr<ShowerApproximation>::tptr app) {
    theShowerApproximation = app;
  }

  /**
   * Return the shower approximation.
   */
  Ptr<ShowerApproximation>::tptr showerApproximation() const { return theShowerApproximation; }

  /**
   * Indicate that the shower real emission contribution should be subtracted.
   */
  void doRealShowerSubtraction() { theRealShowerSubtraction = true; }

  /**
   * Return true, if the shower real emission contribution should be subtracted.
   */
  bool realShowerSubtraction() const { return theRealShowerSubtraction; }

  /**
   * Indicate that the shower virtual contribution should be subtracted.
   */
  void doVirtualShowerSubtraction() { theVirtualShowerSubtraction = true; }

  /**
   * Return true, if the shower virtual contribution should be subtracted.
   */
  bool virtualShowerSubtraction() const { return theVirtualShowerSubtraction; }

  /**
   * Indicate that the loopsim matched virtual contribution should be subtracted.
   */
  void doLoopSimSubtraction() { theLoopSimSubtraction = true; }

  /**
   * Return true, if the loopsim matched virtual contribution should be subtracted.
   */
  bool loopSimSubtraction() const { return theLoopSimSubtraction; }

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
   * Indicate that the subtraction is being tested.
   */
  void doTestSubtraction() { theSubtractionTest = true; }

  /**
   * Return true, if the subtraction is being tested.
   */
  bool testSubtraction() const { return theSubtractionTest; }

  /**
   * Return true, if verbose
   */
  bool verbose() const { return realEmissionME()->verbose() || underlyingBornME()->verbose(); }

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
   * generateTildeKinematics
   */
  void logGenerateTildeKinematics() const;

  /**
   * Write out diagnostic information for
   * generateRadiationKinematics
   */
  void logGenerateRadiationKinematics(const double * r) const;

  /**
   * Write out diagnostic information for
   * me2 evaluation
   */
  void logME2() const;

  /**
   * Write out diagnostic information
   * for dsigdr evaluation
   */
  void logDSigHatDR(double effectiveJac) const;

  //@}

  /** @name Reweight objects */
  //@{

  /**
   * Insert a reweight object
   */
  void addReweight(Ptr<MatchboxReweightBase>::ptr rw) { theReweights.push_back(rw); }

  /**
   * Return the reweight objects
   */
  const vector<Ptr<MatchboxReweightBase>::ptr>& reweights() const { return theReweights; }

  /**
   * Access the reweight objects
   */
  vector<Ptr<MatchboxReweightBase>::ptr>& reweights() { return theReweights; }

  //@}

  /** @name Methods used to setup SubtractionDipole objects */
  //@{

  /**
   * Clone this dipole.
   */
  Ptr<SubtractionDipole>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<SubtractionDipole>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
      
  void cloneDependencies(const std::string& prefix = "", bool slim=false);

  //@}

  /** @name Methods required to setup the event record */
  //@{

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub);

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub, const ColourLines* cl);

  /**
   * Comlete a SubProcess object using the internal degrees of freedom
   * generated in the last generateKinematics() (and possible other
   * degrees of freedom which was intergated over in dSigHatDR(). This
   * default version does nothing. Will be made purely virtual in the
   * future.
   */
  virtual void generateSubCollision(SubProcess & sub);
      
  /**
   * Alpha parameter as in Nagy
   * (http://arxiv.org/pdf/hep-ph/0307268v2.pdf) to restrict dipole
   * phase space
   */
   double alpha() const;
      
   /*
    * True if phase space point is above the alpha cut for this dipole.
    */
      
   bool aboveAlpha() const;

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

private:

  /**
   * The factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr theFactory;

  /**
   * Wether or not this dipole acts in splitting mode.
   */
  bool theSplitting;

  /**
   * True, if should apply to process in the xcomb.
   */
  bool theApply;

  /**
   * True, if the subtraction is being tested.
   */
  bool theSubtractionTest;

  /**
   * True if cuts should be ignored
   */
  bool theIgnoreCuts;

  /**
   * The real emission matrix element to be considered
   */
  Ptr<MatchboxMEBase>::ptr theRealEmissionME;

  /**
   * The underlying Born matrix element
   */
  Ptr<MatchboxMEBase>::ptr theUnderlyingBornME;

  /**
   * The dipoles which have been found along with this dipole
   */
  vector<Ptr<SubtractionDipole>::tptr> thePartners;

  /**
   * The TildeKinematics to be used.
   */
  Ptr<TildeKinematics>::ptr theTildeKinematics;

  /**
   * The InvertedTildeKinematics to be used.
   */
  Ptr<InvertedTildeKinematics>::ptr theInvertedTildeKinematics;

  /**
   * A vector of reweight objects the sum of which
   * should be applied to reweight this dipole
   */
  vector<Ptr<MatchboxReweightBase>::ptr> theReweights;

  /**
   * The emitter as referred to by the real emission
   * matrix element.
   */
  int theRealEmitter;

  /**
   * The emission as referred to by the real emission
   * matrix element.
   */
  int theRealEmission;

  /**
   * The spectator as referred to by the real emission
   * matrix element.
   */
  int theRealSpectator;

  /**
   * The subtraction parameters
   */
  vector<double> theSubtractionParameters;

  /**
   * Map real emission diagrams to underlying Born diagrams 
   * and tilde emitter/spectator.
   */
  map<RealEmissionKey,UnderlyingBornInfo> theMergingMap;

  /**
   * Map underlying Born diagrams and tilde emitter/spectator
   * to real emission diagram containing the splitting.
   */
  multimap<UnderlyingBornKey,RealEmissionInfo> theSplittingMap;

  /**
   * Map underlying Born diagrams to emitter/spectator pairs
   */
  map<cPDVector,pair<int,int> > theIndexMap;

  /**
   * Map real emission processes to Born diagrams
   */
  map<cPDVector,DiagramVector> theUnderlyingBornDiagrams;

  /**
   * Map Born processes to real emission diagrams
   */
  map<cPDVector,DiagramVector> theRealEmissionDiagrams;

  /**
   * Map underlying Born diagrams to real emission diagrams.
   */
  map<tcDiagPtr,tcDiagPtr> theBornToRealDiagrams;

  /**
   * Map real emission diagrams to underlying Born diagrams.
   */
  map<tcDiagPtr,tcDiagPtr> theRealToBornDiagrams;

  /**
   * The last real emission key encountered
   */
  RealEmissionKey lastRealEmissionKey;

  /**
   * The last underlying Born key encountered
   */
  UnderlyingBornKey lastUnderlyingBornKey;

  /**
   * The last real emission info encountered
   */
  multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator lastRealEmissionInfo;

  /**
   * The emitter as referred to by the underlying Born
   * matrix element.
   */
  int theBornEmitter;

  /**
   * The spectator as referred to by the underlying Born
   * matrix element.
   */
  int theBornSpectator;

  /**
   * The last scale as generated from the tilde mapping
   */
  Energy theLastSubtractionScale;

  /**
   * The last scale as generated from the splitting mapping
   */
  Energy theLastSplittingScale;

  /**
   * The last pt as generated from the tilde mapping
   */
  Energy theLastSubtractionPt;

  /**
   * The last pt as generated from the splitting mapping
   */
  Energy theLastSplittingPt;

  /**
   * The last z as generated from the tilde mapping
   */
  double theLastSubtractionZ;

  /**
   * The last z as generated from the splitting mapping
   */
  double theLastSplittingZ;

  /**
   * The shower approximation.
   */
  Ptr<ShowerApproximation>::ptr theShowerApproximation;

  /**
   * True, if the shower real emission contribution should be subtracted.
   */
  bool theRealShowerSubtraction;

  /**
   * True, if the shower virtual contribution should be subtracted.
   */
  bool theVirtualShowerSubtraction;

  /**
   * True, if the loopsim matched virtual contribution should be subtracted.
   */
  bool theLoopSimSubtraction;

  /**
   * True, if scales should be calculated from real emission kinematics
   */
  bool theRealEmissionScales;

  /**
   * Return the shower hard scale encountered
   */
  Energy theShowerHardScale;

  /**
   * The shower evolution scale encountered
   */
  Energy theShowerScale;

  /**
   * The shower splitting variables encountered
   */
  vector<double> theShowerParameters;

  /**
   * True, if this configuration is in the shower phase space
   */
  bool theIsInShowerPhasespace;

  /**
   * True, if this configuration is above the shower infrared cutoff
   */
  bool theIsAboveCutoff;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SubtractionDipole & operator=(const SubtractionDipole &);

};

}

#endif /* HERWIG_SubtractionDipole_H */
