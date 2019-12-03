// -*- C++ -*-
//
// MatchboxFactory.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxFactory_H
#define HERWIG_MatchboxFactory_H
//
// This is the declaration of the MatchboxFactory class.
//

#include "ThePEG/Handlers/SubProcessHandler.h"

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig/MatrixElement/Matchbox/Utility/Tree2toNGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxFactory automatically sets up a NLO
 * QCD calculation carried out in dipole subtraction.
 *
 * @see \ref MatchboxFactoryInterfaces "The interfaces"
 * defined for MatchboxFactory.
 */
class MatchboxFactory: public SubProcessHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxFactory();

  /**
   * The destructor.
   */
  virtual ~MatchboxFactory();
  //@}

public:

  /**
   * Pointer to the current factory object
   */
  static const Ptr<MatchboxFactory>::tptr currentFactory() {
    assert(theCurrentFactory);
    return theCurrentFactory;
  }

private:

  /**
   * Pointer to the current factory object
   */
  static Ptr<MatchboxFactory>::tptr theCurrentFactory;

public:

  /**
   * Flag to indicate that at least one MatchboxFactory object is in action
   */
  static bool isMatchboxRun() {
    return theIsMatchboxRun();
  }

  /** @name Process and diagram information */
  //@{

  /**
   * Return the diagram generator.
   */
  Ptr<Tree2toNGenerator>::tptr diagramGenerator() const { return theDiagramGenerator; }

  /**
   * Set the diagram generator.
   */
  void diagramGenerator(Ptr<Tree2toNGenerator>::ptr dg) { theDiagramGenerator = dg; }

  /**
   * Return the process data.
   */
  Ptr<ProcessData>::tptr processData() const { return theProcessData; }

  /**
   * Set the process data.
   */
  void processData(Ptr<ProcessData>::ptr pd) { theProcessData = pd; }

  /**
   * Return the number of light flavours, this matrix
   * element is calculated for.
   */
  unsigned int nLight() const { return theNLight; }

  /**
   * Set the number of light flavours, this matrix
   * element is calculated for.
   */
  void nLight(unsigned int n) { theNLight = n; }

  /**
   * Return the vector that contains the PDG ids of 
   * the light flavours, which are contained in the
   * jet particle group.
   */
  vector<long> nLightJetVec() const { return theNLightJetVec; }

  /**
   * Set the elements of the vector that contains the PDG
   * ids of the light flavours, which are contained in the
   * jet particle group.
   */
  void nLightJetVec(long n) { theNLightJetVec.push_back(n); }

  /**
   * Return the vector that contains the PDG ids of 
   * the heavy flavours, which are contained in the
   * jet particle group.
   */
  vector<long> nHeavyJetVec() const { return theNHeavyJetVec; }

  /**
   * Set the elements of the vector that contains the PDG
   * ids of the heavy flavours, which are contained in the
   * jet particle group.
   */
  void nHeavyJetVec(long n) { theNHeavyJetVec.push_back(n); }

  /**
   * Return the vector that contains the PDG ids of 
   * the light flavours, which are contained in the
   * proton particle group.
   */
  vector<long> nLightProtonVec() const { return theNLightProtonVec; }

  /**
   * Set the elements of the vector that contains the PDG
   * ids of the light flavours, which are contained in the
   * proton particle group.
   */
  void nLightProtonVec(long n) { theNLightProtonVec.push_back(n); }

  /**
   * Return the order in \f$\alpha_S\f$.
   */
  unsigned int orderInAlphaS() const { return theOrderInAlphaS; }

  /**
   * Set the order in \f$\alpha_S\f$.
   */
  void orderInAlphaS(unsigned int o) { theOrderInAlphaS = o; }

  /**
   * Return the order in \f$\alpha_{EM}\f$.
   */
  unsigned int orderInAlphaEW() const { return theOrderInAlphaEW; }

  /**
   * Set the order in \f$\alpha_{EM}\f$.
   */
  void orderInAlphaEW(unsigned int o) { theOrderInAlphaEW = o; }
 
  /**
   * The multiplicity of legs with virtual contributions.
   */
  size_t highestVirt() const {return theHighestVirtualSize;}

  /**
   * Set the highest 
   **/
  void setHighestVirt(size_t n){theHighestVirtualSize=n;}
 
  /**
   * Access the processes vector.
   */
   const vector<vector<string> > getProcesses() const {return processes;}

  /**
   * Return true, if all processes up to a maximum order are considered
   */
  bool allProcesses() const { return theAllProcesses; }

  /**
   * Switch on/off inclusino off all processes up to a maximum order
   */
  void setAllProcesses(bool on = true) { theAllProcesses = on; }

  /**
   * Return true, if Born contributions should be included.
   */
  bool bornContributions() const { return theBornContributions; }

  /**
   * Switch on or off Born contributions
   */
  void setBornContributions(bool on = true) { theBornContributions = on; }

  /**
   * Return true, if virtual contributions should be included.
   */
  bool virtualContributions() const { return theVirtualContributions; }

  /**
   * Switch on or off virtual contributions
   */
  void setVirtualContributions(bool on = true) { theVirtualContributions = on; }

  /**
   * Produce matrix element corrections, but no NLO
   */
  bool meCorrectionsOnly() const { return theMECorrectionsOnly; }

  /**
   * Switch to produce matrix element corrections, but no NLO
   */
  void setMECorrectionsOnly(bool on = true) { theMECorrectionsOnly = on; }

  /**
   * Produce matrix element corrections, with LoopSim NLO
   */
  bool loopSimCorrections() const { return theLoopSimCorrections; }

  /**
   * Switch to produce matrix element corrections, with LoopSim NLO
   */
  void setLoopSimCorrections(bool on = true) { theLoopSimCorrections = on; }

  /**
   * Return true, if subtracted real emission contributions should be included.
   */
  bool realContributions() const { return theRealContributions; }

  /**
   * Switch on or off subtracted real emission contributions
   */
  void setRealContributions(bool on = true) { theRealContributions = on; }

  /**
   * Return true, if virtual contributions should be treated as independent subprocesses
   */
  bool independentVirtuals() const { return theIndependentVirtuals; }

  /**
   * Switch on/off virtual contributions should be treated as independent subprocesses
   */
  void setIndependentVirtuals(bool on = true) { theIndependentVirtuals = on; }

  /**
   * Return true, if PK operator contributions should be treated as independent subprocesses
   */
  bool independentPKs() const { return theIndependentPKs; }

  /**
   * Switch on/off PK operator contributions should be treated as independent subprocesses
   */
  void setIndependentPKs(bool on = true) { theIndependentPKs = on; }

  /**
   * Return true, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  virtual bool subProcessGroups() const { return !showerApproximation(); }

  /**
   * Return true, if subtraction scales should be caluclated from real emission kinematics
   */
  bool realEmissionScales() const { return theRealEmissionScales; }

  /**
   * Switch on/off that subtraction scales should be caluclated from real emission kinematics
   */
  void setRealEmissionScales(bool on = true) { theRealEmissionScales = on; }

  /**
   * Set the shower approximation.
   */
  void showerApproximation(Ptr<ShowerApproximation>::tptr app) { theShowerApproximation = app; }

  /**
   * Return the shower approximation.
   */
  Ptr<ShowerApproximation>::tptr showerApproximation() const { return theShowerApproximation; }

  //@}

  /** @name Phasespace generation and scale choice */
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
   * Set the scale choice object
   */
  void scaleChoice(Ptr<MatchboxScaleChoice>::ptr sc) { theScaleChoice = sc; }

  /**
   * Return the scale choice object
   */
  Ptr<MatchboxScaleChoice>::tptr scaleChoice() const { return theScaleChoice; }

  /**
   * Get the factorization scale factor
   */
  double factorizationScaleFactor() const { return theFactorizationScaleFactor; }

  /**
   * Set the factorization scale factor
   */
  void factorizationScaleFactor(double f) { theFactorizationScaleFactor = f; }

  /**
   * Get the renormalization scale factor
   */
  double renormalizationScaleFactor() const { return theRenormalizationScaleFactor; }

  /**
   * Set the renormalization scale factor
   */
  void renormalizationScaleFactor(double f) { theRenormalizationScaleFactor = f; }

  /**
   * Return true, if fixed couplings are used.
   */
  bool fixedCouplings() const { return theFixedCouplings; }

  /**
   * Switch on fixed couplings.
   */
  void setFixedCouplings(bool on = true) { theFixedCouplings = on; }

  /**
   * Return true, if fixed couplings are used.
   */
  bool fixedQEDCouplings() const { return theFixedQEDCouplings; }

  /**
   * Switch on fixed couplings.
   */
  void setFixedQEDCouplings(bool on = true) { theFixedQEDCouplings = on; }

  /**
   * Return true, if veto scales should be set
   * for the real emission
   */
  bool vetoScales() const { return theVetoScales; }

  /**
   * Switch on setting veto scales
   */
  void doVetoScales() { theVetoScales = true; }

  /**
   * Switch off setting veto scales
   */
  void noVetoScales() { theVetoScales = true; }

  //@}

  /** @name Amplitudes and caching */
  //@{

  /**
   * Return the amplitudes to be considered
   */
  const vector<Ptr<MatchboxAmplitude>::ptr>& amplitudes() const { return theAmplitudes; }

  /**
   * Access the amplitudes to be considered
   */
  vector<Ptr<MatchboxAmplitude>::ptr>& amplitudes() { return theAmplitudes; }

  //@}

  /** @name Matrix element objects. */
  //@{

  /**
   * Return the Born matrix elements to be considered
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& bornMEs() const { return theBornMEs; }

  /**
   * Access the Born matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr>& bornMEs() { return theBornMEs; }

  /**
   * Return the loop induced matrix elements to be considered
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& loopInducedMEs() const { return theLoopInducedMEs; }

  /**
   * Access the loop induced matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr>& loopInducedMEs() { return theLoopInducedMEs; }

  /**
   * Return the processes to be ordered from an OLP
   */
  const map<Ptr<MatchboxAmplitude>::tptr,
	    map<pair<Process,int>,int> >&
  olpProcesses() const { return theOLPProcesses; }

  /**
   * Access the processes to be ordered from an OLP
   */
  map<Ptr<MatchboxAmplitude>::tptr,
      map<pair<Process,int>,int> >& 
  olpProcesses() { return theOLPProcesses; }

  /**
   * Order an OLP process and return its id
   */
  int orderOLPProcess(const Process& p,
		      Ptr<MatchboxAmplitude>::tptr amp,
		      int type);

  /**
   * Return the amplitudes which need external initialization
   */
  const set<Ptr<MatchboxAmplitude>::tptr>& externalAmplitudes() const {
    return theExternalAmplitudes;
  }

  /**
   * Access the amplitudes which need external initialization
   */
  set<Ptr<MatchboxAmplitude>::tptr>& externalAmplitudes() {
    return theExternalAmplitudes;
  }

  /**
   * Return the virtual corrections to be considered
   */
  const vector<Ptr<MatchboxInsertionOperator>::ptr>& virtuals() const { return theVirtuals; }

  /**
   * Access the virtual corrections to be considered
   */
  vector<Ptr<MatchboxInsertionOperator>::ptr>& virtuals() { return theVirtuals; }

  /**
   * Return the produced NLO matrix elements
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& bornVirtualMEs() const { return theBornVirtualMEs; }

  /**
   * Access the produced NLO matrix elements
   */
  vector<Ptr<MatchboxMEBase>::ptr>& bornVirtualMEs() { return theBornVirtualMEs; }

  /**
   * Return the real emission matrix elements to be considered
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& realEmissionMEs() const { return theRealEmissionMEs; }

  /**
   * Access the real emission matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr>& realEmissionMEs() { return theRealEmissionMEs; }

  /**
   * Return, which set of dipoles should be considered
   */
  int dipoleSet() const { return theDipoleSet; }

  /**
   * Return, which set of dipoles should be considered
   */
  void dipoleSet(int s) { theDipoleSet = s; }

  /**
   * Return the produced subtracted matrix elements
   */
  const vector<Ptr<SubtractedME>::ptr>& subtractedMEs() const { return theSubtractedMEs; }

  /**
   * Access the produced subtracted matrix elements
   */
  vector<Ptr<SubtractedME>::ptr>& subtractedMEs() { return theSubtractedMEs; }

  /**
   * Return the produced finite real emission matrix elements
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& finiteRealMEs() const { return theFiniteRealMEs; }

  /**
   * Access the produced finite real emission elements
   */
  vector<Ptr<MatchboxMEBase>::ptr>& finiteRealMEs() { return theFiniteRealMEs; }

  /**
   * Return the map of Born processes to splitting dipoles
   */
  const map<cPDVector,set<Ptr<SubtractionDipole>::ptr> >& splittingDipoles() const {
    return theSplittingDipoles;
  }

  /**
   * Identify a splitting channel
   */
  struct SplittingChannel {

    /**
     * The Born XComb
     */
    StdXCombPtr bornXComb;

    /**
     * The real XComb
     */
    StdXCombPtr realXComb;

    /**
     * The set of tilde XCombs to consider for the real xcomb
     */
    vector<StdXCombPtr> tildeXCombs;

    /**
     * The dipole in charge of the splitting
     */
    Ptr<SubtractionDipole>::ptr dipole;

    /**
     * Dump the setup
     */
    void print(ostream&) const;

  };

  /**
   * Generate all splitting channels for the Born process handled by
   * the given XComb
   */
  list<SplittingChannel> getSplittingChannels(tStdXCombPtr xc) const;

  /**
   * Return the reweight objects for matrix elements
   */
  const vector<ReweightPtr>& reweighters() const { return theReweighters; }

  /**
   * Access the reweight objects for matrix elements
   */
  vector<ReweightPtr>& reweighters() { return theReweighters; }

  /**
   * Return the preweight objects for matrix elements
   */
  const vector<ReweightPtr>& preweighters() const { return thePreweighters; }

  /**
   * Access the preweight objects for matrix elements
   */
  vector<ReweightPtr>& preweighters() { return thePreweighters; }

  //@}

  /** @name Setup the matrix elements */
  //@{

  /**
   * Return true if this object needs to be initialized before all
   * other objects (except those for which this function also returns
   * true).  This default version always returns false, but subclasses
   * may override it to return true.
   */
  virtual bool preInitialize() const { return true; }

  /**
   * Prepare a matrix element.
   */
  void prepareME(Ptr<MatchboxMEBase>::ptr);

  /**
   * Check consistency and switch to porduction mode.
   */
  virtual void productionMode();

  /**
   * Setup everything
   */
  virtual void setup();
  
   /**
   * The highest multiplicity of legs having virtual contributions.(needed for madgraph) 
   */

  size_t highestVirt(){return theHighestVirtualsize;}

  //@}

  /** @name Diagnostic information */
  //@{

  /**
   * Return true, if verbose
   */
  bool verbose() const { return theVerbose; }

  /**
   * Switch on diagnostic information.
   */
  void setVerbose(bool on = true) { theVerbose = on; }
  
  /**
   * Return true, if verbose while initializing
   */
  bool initVerbose() const { return theInitVerbose || verbose(); }

  /**
   * Switch on diagnostic information while initializing
   */
  void setInitVerbose(bool on = true) { theInitVerbose = on; }

  /**
   * Dump the setup
   */
  void print(ostream&) const;

  /**
   * Return the subtraction data prefix.
   */
  const string& subtractionData() const { return theSubtractionData; }

  /**
   * Set the subtraction data prefix.
   */
  void subtractionData(const string& s) { theSubtractionData = s; }

  /**
   * Return the subtraction plot type.
   */
  const int& subtractionPlotType() const { return theSubtractionPlotType; }

  /**
   * Set the subtraction plot type.
   */
  void subtractionPlotType(const int& t) { theSubtractionPlotType = t; }
  
  /**
   * Return whether subtraction data should be plotted for all phase space points individually
   */
  const bool& subtractionScatterPlot() const { return theSubtractionScatterPlot; }
  
  /**
   * Set whether subtraction data should be plotted for all phase space points individually
   */
  void subtractionScatterPlot(const bool& s) { theSubtractionScatterPlot = s; }
  
  /**
   * Return the pole data prefix.
   */
  const string& poleData() const { return thePoleData; }

  /**
   * Set the pole data prefix.
   */
  void poleData(const string& s) { thePoleData = s; }

  /**
   * Return true, if cancellationn of epsilon poles should be checked.
   */
  bool checkPoles() const { return poleData() != ""; }

  //@}

  /** @name Process generation */
  //@{

  /**
   * Return the particle groups.
   */
  const map<string,PDVector>& particleGroups() const { return theParticleGroups; }

  /**
   * Access the particle groups.
   */
  map<string,PDVector>& particleGroups() { return theParticleGroups; }

  /**
   * Return true, if the given particle is incoming
   */
  bool isIncoming(cPDPtr p) const {
    return theIncoming.find(p->id()) != theIncoming.end();
  }

  /**
   * Return true, if spin correlation information should be provided, if possible.
   */
  bool spinCorrelations() const { return theSpinCorrelations; }

  /**
   * Indicate that spin correlation information should be provided, if possible.
   */
  void setSpinCorrelations(bool yes) { theSpinCorrelations = yes; }

  //@}

  /** @name Truncated qtilde shower information */
  //@{

  /**
   * Return the subprocess of the real emission
   */
  tSubProPtr hardTreeSubprocess() { return theHardtreeSubprocess; }

  /**
   * Set the subprocess of the real emission for use in calculating the shower hardtree
   */
  void setHardTreeSubprocess(tSubProPtr hardTree) { theHardtreeSubprocess = hardTree; }

  /**
   * Return the born emitter 
   */
  int hardTreeEmitter() { return theHardtreeEmitter; }

  /**
   * Set the born emitter for use in calculating the shower hardtree
   */
  void setHardTreeEmitter(int emitter) { theHardtreeEmitter = emitter; }

  /**
   * Return the born spectator 
   */
  int hardTreeSpectator() { return theHardtreeSpectator; }

  /**
   * Set the born spectator for use in calculating the shower hardtree
   */
  void setHardTreeSpectator(int spectator) { theHardtreeSpectator = spectator; }

  //@}

  /** @name Data handling */
  //@{
  /**
   * Return (and possibly create) a directory to contain amplitude
   * information.
   */
  const string& buildStorage();

  /**
   * Return (and possibly create) a directory to contain integration grid
   * information.
   */
  const string& runStorage();
  
  /**
   *  alpha of http://arxiv.org/pdf/hep-ph/0307268v2.pdf to restrict 
   *  dipole phase space
   */
  double alphaParameter() const { return theAlphaParameter; }
  
  /**
   *  set the alpha parameter (needed for massive PK-Operator)
   */
  void setAlphaParameter(double a)const { theAlphaParameter = a; }
  
  //@}

public:

  /**
   * Print a summary of the parameters used
   */
  void summary(ostream&) const;

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
  //@}

private:

  /**
   * Flag to indicate that at least one MatchboxFactory object is in action
   */
  static bool& theIsMatchboxRun();

  /**
   * The diagram generator.
   */
  Ptr<Tree2toNGenerator>::ptr theDiagramGenerator;

  /**
   * The process data object to be used
   */
  Ptr<ProcessData>::ptr theProcessData;

  /**
   * The number of light flavours, this matrix
   * element is calculated for.
   */
  unsigned int theNLight;

  /**
   * Vector with the PDG ids of the light quark flavours,
   * which are contained in the jet particle group.
   */
  vector<long> theNLightJetVec;

  /**
   * Vector with the PDG ids of the heavy quark flavours,
   * which are contained in the jet particle group.
   */
  vector<long> theNHeavyJetVec;

  /**
   * Vector with the PDG ids of the light quark flavours,
   * which are contained in the proton particle group.
   */
  vector<long> theNLightProtonVec;

  /**
   * The order in \f$\alpha_S\f$.
   */
  unsigned int theOrderInAlphaS;

  /**
   * The order in \f$\alpha_{EM}\f$.
   */
  unsigned int theOrderInAlphaEW;

  /**
   * The maximum number of legs with virtual corrections.
   **/
  unsigned int theHighestVirtualSize;

  /**
   * Switch on or off Born contributions
   */
  bool theBornContributions;

  /**
   * Switch on or off virtual contributions
   */
  bool theVirtualContributions;

  /**
   * Switch on or off subtracted real emission contributions should be included.
   */
  bool theRealContributions;

  /**
   * True if virtual contributions should be treated as independent subprocesses
   */
  bool theIndependentVirtuals;

  /**
   * True if PK operator contributions should be treated as independent subprocesses
   */
  bool theIndependentPKs;

  /**
   * The phase space generator to be used.
   */
  Ptr<MatchboxPhasespace>::ptr thePhasespace;

  /**
   * The scale choice object
   */
  Ptr<MatchboxScaleChoice>::ptr theScaleChoice;

  /**
   * The factorization scale factor.
   */
  double theFactorizationScaleFactor;

  /**
   * The renormalization scale factor.
   */
  double theRenormalizationScaleFactor;

  /**
   * Use non-running couplings.
   */
  bool theFixedCouplings;

  /**
   * Use non-running couplings.
   */
  bool theFixedQEDCouplings;

  /**
   * True, if veto scales should be set
   * for the real emission
   */
  bool theVetoScales;

  /**
   * The amplitudes to be considered
   */
  vector<Ptr<MatchboxAmplitude>::ptr> theAmplitudes;

  /**
   * The Born matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr> theBornMEs;

  /**
   * The loop induced matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr> theLoopInducedMEs;

  /**
   * The virtual corrections to be considered
   */
  vector<Ptr<MatchboxInsertionOperator>::ptr> theVirtuals;

  /**
   * The real emission matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr> theRealEmissionMEs;

  /**
   * The produced NLO matrix elements
   */
  vector<Ptr<MatchboxMEBase>::ptr> theBornVirtualMEs;

  /**
   * The produced subtracted matrix elements
   */
  vector<Ptr<SubtractedME>::ptr> theSubtractedMEs;

  /**
   * The produced finite real emission matrix elements
   */
  vector<Ptr<MatchboxMEBase>::ptr> theFiniteRealMEs;

  /**
   * Which set of dipoles should be considered
   */
  int theDipoleSet;

  /**
   * Switch on or off verbosity
   */
  bool theVerbose;
  
  /**
   * True, if verbose while initializing
   */
  bool theInitVerbose;

  /**
   * Prefix for subtraction data
   */
  string theSubtractionData;

  /**
   * Set the type of plot that is to be generated for subtraction checking
   */
  int theSubtractionPlotType;
  
  /**
   * Set whether subtraction data should be plotted for all phase space points individually
   */
  bool theSubtractionScatterPlot;
  
  /**
   * Prefix for pole data.
   */
  string thePoleData;

  /**
   * Command to limit the real emission process to be considered.
   */
  string doSingleRealProcess(string);

  /**
   * The real emission process to be included; if empty, all possible
   * ones will be considered.
   */
  vector<vector<string> > realEmissionProcesses;

  /**
   * Particle groups.
   */
  map<string,PDVector> theParticleGroups;

  /**
   * Command to start a particle group.
   */
  string startParticleGroup(string);

  /**
   * The name of the particle group currently edited.
   */
  string particleGroupName;

  /**
   * The particle group currently edited.
   */
  PDVector particleGroup;

  /**
   * Command to end a particle group.
   */
  string endParticleGroup(string);

protected:
  
  /**
   * Parse a process description
   */
  virtual vector<string> parseProcess(string);

private:
  
  /**
   * Command to set the process.
   */
  string doProcess(string);

  /**
   * Command to set the process.
   */
  string doLoopInducedProcess(string);

  /**
   * The process to consider in terms of particle groups.
   */
  vector<vector<string> > processes;

  /**
   * The loop induced process to consider in terms of particle groups.
   */
  vector<vector<string> > loopInducedProcesses;

  /**
   * Generate subprocesses.
   */
  set<PDVector> makeSubProcesses(const vector<string>&) const;

 
public: 
  
  /**
   * Generate matrix element objects for the given process.
   */
  vector<Ptr<MatchboxMEBase>::ptr> makeMEs(const vector<string>&, 
					   unsigned int orderas,
					   bool virt);

  
private:
  /**
   * The shower approximation.
   */
  Ptr<ShowerApproximation>::ptr theShowerApproximation;

  /**
   * The map of Born processes to splitting dipoles
   */
  map<cPDVector,set<Ptr<SubtractionDipole>::ptr> > theSplittingDipoles;

  /**
   * True, if subtraction scales should be caluclated from real emission kinematics
   */
  bool theRealEmissionScales;

  /**
   * Consider all processes with order in couplings specifying the
   * maximum order.
   */
  bool theAllProcesses;

  /**
   * The processes to be ordered from an OLP
   */
  map<Ptr<MatchboxAmplitude>::tptr,map<pair<Process,int>,int> > theOLPProcesses;

  /**
   * Amplitudes which need external initialization
   */
  set<Ptr<MatchboxAmplitude>::tptr> theExternalAmplitudes;

  /**
   * Amplitudes to be selected on clashing responsibilities.
   */
  vector<Ptr<MatchboxAmplitude>::ptr> theSelectedAmplitudes;

  /**
   * Amplitudes to be deselected on clashing responsibilities.
   */
  vector<Ptr<MatchboxAmplitude>::ptr> theDeselectedAmplitudes;

  /**
   * Reweight objects for matrix elements
   */
  vector<ReweightPtr> theReweighters;

  /**
   * Preweight objects for matrix elements
   */
  vector<ReweightPtr> thePreweighters;

  /**
   * Produce matrix element corrections, but no NLO
   */
  bool theMECorrectionsOnly;

  /**
   * The highest multiplicity of legs having virtual contributions.(needed for madgraph) 
   */
  int theHighestVirtualsize;

  /**
   * Produce matrix element corrections, with LoopSim NLO
   */
  bool theLoopSimCorrections;

  /**
   * True, if the setup has already been run.
   */
  bool ranSetup;

  /**
   * PDG ids of incoming particles
   */
  set<long> theIncoming;

  /**
   * True, if first incoming partons originate from perturbative PDF
   */
  bool theFirstPerturbativePDF;

  /**
   * True, if second incoming partons originate from perturbative PDF
   */
  bool theSecondPerturbativePDF;

  /**
   * True, if this Factory is in production mode.
   */
  bool inProductionMode;

  /**
   * The real emission subprocess used when calculating the hardtree
   * in the truncated qtilde shower
   */
  tSubProPtr theHardtreeSubprocess;

  /**
   * The born emitter used when calculating the hardtree in
   * the truncated shower
   */
  int theHardtreeEmitter;

  /**
   * The born spectator used when calculating the hardtree in
   * the truncated shower
   */
  int theHardtreeSpectator;

  /**
   * True, if spin correlation information should be provided, if possible.
   */
  bool theSpinCorrelations;

  /**
   * The alpha parameter to be used for the dipole subtraction
   * JB: The parameter is muatble, since we need to be able to change it 
   * while calculating the difference of IPK with and without alpha.
   */  
  mutable double theAlphaParameter;

  /**
   * Wether or not charge conservation should be enforced for the processes
   * constructed.
   */
  bool theEnforceChargeConservation;

  /**
   * Wether or not colour conservation should be enforced for the processes
   * constructed.
   */
  bool theEnforceColourConservation;

  /**
   * Wether or not lepton number conservation should be enforced for the processes
   * constructed.
   */
  bool theEnforceLeptonNumberConservation;

  /**
   * Wether or not quark number conservation should be enforced for the processes
   * constructed.
   */
  bool theEnforceQuarkNumberConservation;

  /**
   * Assume flavour diagonal lepton interactions
   */
  bool theLeptonFlavourDiagonal;

  /**
   * Assume flavour diagonal quark interactions
   */
  bool theQuarkFlavourDiagonal;

  /**
   * Command for production mode
   */
  string doProductionMode(string) {
    productionMode(); return "";
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxFactory & operator=(const MatchboxFactory &) = delete;

};

}

#endif /* HERWIG_MatchboxFactory_H */
