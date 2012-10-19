// -*- C++ -*-
//
// MatchboxFactory.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxFactory_H
#define HERWIG_MatchboxFactory_H
//
// This is the declaration of the MatchboxFactory class.
//

#include "ThePEG/Handlers/SubProcessHandler.h"

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxAmplitude.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/Tree2toNGenerator.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/MatchboxScaleChoice.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/MatchboxMECache.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxNLOME.h"
#include "Herwig++/MatrixElement/Matchbox/Base/SubtractedME.h"

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
   * Return true, if subtracted real emission contributions should be included.
   */
  bool realContributions() const { return theRealContributions; }

  /**
   * Switch on or off subtracted real emission contributions
   */
  void setRealContributions(bool on = true) { theRealContributions = on; }

  /**
   * Return true, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  bool subProcessGroups() const { return theSubProcessGroups; }

  /**
   * Switch on or off producing subprocess groups.
   */
  void setSubProcessGroups(bool on = true) { theSubProcessGroups = on; }

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

  /**
   * Set the ME cache object
   */
  void cache(Ptr<MatchboxMECache>::ptr c) { theCache = c; }

  /**
   * Get the ME cache object
   */
  Ptr<MatchboxMECache>::tptr cache() const { return theCache; }

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
  const vector<Ptr<MatchboxNLOME>::ptr>& bornVirtualMEs() const { return theBornVirtualMEs; }

  /**
   * Access the produced NLO matrix elements
   */
  vector<Ptr<MatchboxNLOME>::ptr>& bornVirtualMEs() { return theBornVirtualMEs; }

  /**
   * Return the real emission matrix elements to be considered
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& realEmissionMEs() const { return theRealEmissionMEs; }

  /**
   * Access the real emission matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr>& realEmissionMEs() { return theRealEmissionMEs; }

  /**
   * Return the produced subtracted matrix elements
   */
  const vector<Ptr<SubtractedME>::ptr>& subtractedMEs() const { return theSubtractedMEs; }

  /**
   * Access the produced subtracted matrix elements
   */
  vector<Ptr<SubtractedME>::ptr>& subtractedMEs() { return theSubtractedMEs; }

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
  void prepareME(Ptr<MatchboxMEBase>::ptr) const;

  /**
   * Setup everything
   */
  void setup();

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
   * Return true, if cancellationn of epsilon poles should be checked.
   */
  bool checkPoles() const { return theCheckPoles; }

  /**
   * Switch on checking of epsilon pole cancellation.
   */
  void doCheckPoles() { theCheckPoles = true; }

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

private:

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
   * The order in \f$\alpha_S\f$.
   */
  unsigned int theOrderInAlphaS;

  /**
   * The order in \f$\alpha_{EM}\f$.
   */
  unsigned int theOrderInAlphaEW;

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
   * True, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  bool theSubProcessGroups;

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
   * The ME cache object
   */
  Ptr<MatchboxMECache>::ptr theCache;

  /**
   * The Born matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr> theBornMEs;

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
  vector<Ptr<MatchboxNLOME>::ptr> theBornVirtualMEs;

  /**
   * The produced subtracted matrix elements
   */
  vector<Ptr<SubtractedME>::ptr> theSubtractedMEs;

  /**
   * Switch on or off verbosity
   */
  bool theVerbose;

  /**
   * Prefix for subtraction data
   */
  string theSubtractionData;

  /**
   * Command to limit the real emission process to be considered.
   */
  string doSingleRealProcess(string);

  /**
   * The real emission process to be included; if empty, all possible
   * ones will be considered.
   */
  vector<string> realEmissionProcess;

  /**
   * True, if cancellationn of epsilon poles should be checked.
   */
  bool theCheckPoles;

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

  /**
   * Command to set the process.
   */
  string doProcess(string);

  /**
   * The process to consider in terms of particle groups.
   */
  vector<string> process;

  /**
   * Generate subprocesses.
   */
  set<PDVector> makeSubProcesses(const vector<string>&) const;

  /**
   * Generate matrix element objects for the given process.
   */
  vector<Ptr<MatchboxMEBase>::ptr> makeMEs(const vector<string>&, 
					   unsigned int orderas) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxFactory & operator=(const MatchboxFactory &);

};

}

#endif /* HERWIG_MatchboxFactory_H */
