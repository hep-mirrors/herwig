// -*- C++ -*-
//
// SubtractionDipole.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SubtractionDipole_H
#define HERWIG_SubtractionDipole_H
//
// This is the declaration of the SubtractionDipole class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

namespace Herwig {

using namespace ThePEG;

class TildeKinematics;
class InvertedTildeKinematics;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief SubtractionDipole represents a dipole subtraction
 * term in the formalism of Catani and Seymour.
 *
 */
class SubtractionDipole: public MEBase {

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

  /** @name Subprocess and diagram information. */
  //@{

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
   * Setup bookkeeping maps.
   */
  void setupBookkeeping();

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
  bool empty() const { return theSplittingMap.empty(); }

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
  const map<UnderlyingBornKey,RealEmissionInfo>& splittingMap() const { return theSplittingMap; }

  /**
   * Return the underlying Born diagrams to be considered
   * for the given real emission process.
   */
  const DiagramVector& underlyingBornDiagrams(const cPDVector& real) const;

  /**
   * Return the real emission diagrams to be considered
   * for the given Born process.
   */
  const DiagramVector& realEmissionDiagrams(const cPDVector& born) const;

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
   * Select a diagram. Default version uses diagrams(const
   * DiagramVector &) to select a diagram according to the
   * weights. This is the only method used that should be outside of
   * MEBase.
   */
  virtual DiagramIndex diagram(const DiagramVector & dv) const;

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

  //@}

  /** @name Scale choices, couplings and PDFs */
  //@{

  /**
   * Return the scale associated with the phase space point provided
   * by the last call to setKinematics().
   */
  virtual Energy2 scale() const { return realEmissionME()->scale(); }

  /**
   * Return the value of \f$\alpha_S\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaS(scale()).
   */
  virtual double alphaS() const { return realEmissionME()->lastAlphaS(); }

  /**
   * Return the value of \f$\alpha_EM\f$ associated with the phase
   * space point provided by the last call to setKinematics(). This
   * versions returns SM().alphaEM(scale()).
   */
  virtual double alphaEM() const { return realEmissionME()->lastAlphaEM(); }

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
   * Return the matrix element averaged over spin correlations.
   */
  virtual double me2Avg(double ccme2) const = 0;

  /**
   * Return true, if the cross section should actually return the spin
   * averaged splitting function times the Born matrix element squared.
   */
  bool showerKernel() const { return theShowerKernel; }

  /**
   * Indicate that the cross section should actually return the spin
   * averaged splitting function times the Born matrix element squared.
   */
  void doShowerKernel(bool is = true) { theShowerKernel = is; }

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
  void cloneDependencies(const std::string& prefix = "");

  //@}

  /** @name Methods required to setup the event record */
  //@{

  /**
   * construct the spin information for the interaction
   */
  virtual void constructVertex(tSubProPtr sub);

  /**
   * Comlete a SubProcess object using the internal degrees of freedom
   * generated in the last generateKinematics() (and possible other
   * degrees of freedom which was intergated over in dSigHatDR(). This
   * default version does nothing. Will be made purely virtual in the
   * future.
   */
  virtual void generateSubCollision(SubProcess & sub);

  //@}

protected:

  /**
   * Handle integer powers appearing downstream.
   */
  double pow(double x, unsigned int p) const {
    return std::pow(x,(double)p);
  }

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
   * True, if the cross section should actually return the spin
   * averaged splitting function times the Born matrix element squared.
   */
  bool theShowerKernel;

  /**
   * The real emission matrix element to be considered
   */
  Ptr<MatchboxMEBase>::ptr theRealEmissionME;

  /**
   * The underlying Born matrix element
   */
  Ptr<MatchboxMEBase>::ptr theUnderlyingBornME;

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
  map<UnderlyingBornKey,RealEmissionInfo> theSplittingMap;

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
   * The last real emission key encountered
   */
  RealEmissionKey lastRealEmissionKey;

  /**
   * The last underlying Born key encountered
   */
  UnderlyingBornKey lastUnderlyingBornKey;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SubtractionDipole & operator=(const SubtractionDipole &);

};

}

#endif /* HERWIG_SubtractionDipole_H */
