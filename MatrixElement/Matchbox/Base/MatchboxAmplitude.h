// -*- C++ -*-
//
// MatchboxAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxAmplitude_H
#define HERWIG_MatchboxAmplitude_H
//
// This is the declaration of the MatchboxAmplitude class.
//

#include "ThePEG/MatrixElement/Amplitude.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/ColourBasis.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxAmplitude is the base class for amplitude
 * implementations inside Matchbox.
 *
 * @see \ref MatchboxAmplitudeInterfaces "The interfaces"
 * defined for MatchboxAmplitude.
 */
class MatchboxAmplitude: public Amplitude, public LastXCombInfo<StandardXComb> {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxAmplitude();

  /**
   * The destructor.
   */
  virtual ~MatchboxAmplitude();
  //@}

public:

  typedef map<vector<int>,CVector> AmplitudeMap;
  typedef map<vector<int>,CVector>::iterator AmplitudeIterator;
  typedef map<vector<int>,CVector>::const_iterator AmplitudeConstIterator;

  /**
   * Return the amplitude. Needs to be implemented from
   * ThePEG::Amplitude but is actually ill-defined, as colours of the
   * external particles are not specified. To this extent, this
   * implementation just asserts.
   */
  virtual Complex value(const tcPDVector & particles,
			const vector<Lorentz5Momentum> & momenta, 
			const vector<int> & helicities);

  /** @name Subprocess information */
  //@{

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const { return false; }

  /**
   * Provide the additional random numbers
   */
  void additionalKinematics(const double *);

  /**
   * Return the number of random numbers required to evaluate this
   * amplitude at a fixed phase space point.
   */
  virtual int nDimAdditional() const { return 0; }

  /**
   * Return a ME instance appropriate for this amplitude and the given
   * subprocesses
   */
  virtual Ptr<MatchboxMEBase>::ptr makeME(const vector<PDVector>&) const;

  /**
   * Return the process data.
   */
  Ptr<ProcessData>::tptr processData() const { return theProcessData; }

  /**
   * Set the process data.
   */
  void processData(Ptr<ProcessData>::ptr pd) { theProcessData = pd; }

  /**
   * Return the amplitude parton data.
   */
  const cPDVector& lastAmplitudePartonData() const { return theLastAmplitudePartonData->second; }

  /**
   * Access the amplitude parton data.
   */
  cPDVector& lastAmplitudePartonData() { return theLastAmplitudePartonData->second; }

  /**
   * Access the amplitude parton data.
   */
  map<tStdXCombPtr,cPDVector>& amplitudePartonData() { return processData()->amplitudePartonData(); }

  /**
   * Return the number of light flavours
   */
  unsigned int nLight() const { return theNLight; }

  /**
   * Set the number of light flavours
   */
  void nLight(unsigned int n) { theNLight = n; }

  /**
   * Set the (tree-level) order in \f$g_S\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGs(unsigned int) {}

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const = 0;

  /**
   * Set the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGem(unsigned int) {}

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const = 0;

  /**
   * Return the Herwig++ StandardModel object
   */
  Ptr<StandardModel>::tcptr standardModel() { 
    if ( !theStandardModel )
      theStandardModel = 
	dynamic_ptr_cast<Ptr<StandardModel>::tcptr>(HandlerBase::standardModel());
    return theStandardModel;
  }

  //@}

  /** @name Colour basis. */
  //@{

  /**
   * Return the colour basis.
   */
  Ptr<ColourBasis>::tptr colourBasis() const { return theColourBasis; }

  /**
   * Set the colour basis dimensionality.
   */
  void colourBasisDim(size_t dim) { theColourBasisDim = dim; }

  /**
   * Get the colour basis dimensionality.
   */
  size_t colourBasisDim() const { return theColourBasisDim; }

  /**
   * Return true, if this amplitude will not require colour correlations.
   */
  virtual bool noCorrelations() const { return !haveOneLoop(); }  

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { 
    return colourBasis() ? colourBasis()->haveColourFlows() : false;
  }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *> colourGeometries(tcDiagPtr diag) const {
    return 
      haveColourFlows() ? 
      theColourBasis->colourGeometries(diag,lastLargeNAmplitudes()) :
      Selector<const ColourLines *>();
  }

  /**
   * Return the colour crossing information as filled by the last call to
   * fillCrossingMap(...), mapping amplitude ids to colour basis ids.
   */
  const map<size_t,size_t>& lastAmplitudeToColourMap() const { return theLastAmplitudeToColourMap->second; }

  /**
   * Access the colour crossing information.
   */
  map<size_t,size_t>& lastAmplitudeToColourMap() { return theLastAmplitudeToColourMap->second; }  

  /**
   * Access the colour crossing information.
   */
  map<tStdXCombPtr,map<size_t,size_t> >& amplitudeToColourMap() { return processData()->amplitudeToColourMap(); }

  /**
   * Return the colour crossing information as filled by the last call to
   * fillCrossingMap(...), mapping amplitude ids to colour basis ids.
   */
  const map<size_t,size_t>& lastColourToAmplitudeMap() const { return theLastColourToAmplitudeMap->second; }

  /**
   * Access the colour crossing information.
   */
  map<size_t,size_t>& lastColourToAmplitudeMap() { return theLastColourToAmplitudeMap->second; }  

  /**
   * Access the colour crossing information.
   */
  map<tStdXCombPtr,map<size_t,size_t> >& colourToAmplitudeMap() { return processData()->colourToAmplitudeMap(); }

  /**
   * Return an ordering identifier for the current subprocess and
   * colour absis tensor index.
   */
  const string& colourOrdering(size_t id) const;

  //@}

  /** @name Phasespace point, crossing and helicities */
  //@{

  /**
   * Set the xcomb object.
   */
  virtual void setXComb(tStdXCombPtr xc) {
    theLastXComb = xc;
    fillCrossingMap();
  }

  /**
   * Return the momentum as crossed appropriate for this amplitude.
   */
  Lorentz5Momentum amplitudeMomentum(int) const;

  /**
   * Perform a normal ordering of external legs and fill the
   * crossing information as. This default implementation sorts
   * lexicographically in (abs(colour)/spin/abs(charge)), putting pairs
   * of particles/anti-particles where possible.
   */
  virtual void fillCrossingMap(size_t shift = 0);

  /**
   * Return the crossing sign.
   */
  double lastCrossingSign() const { return theLastCrossingSign; }

  /**
   * Set the crossing sign.
   */
  void lastCrossingSign(double s) { theLastCrossingSign = s; }

  /**
   * Return the crossing information as filled by the last call to
   * fillCrossingMap(...), mapping amplitude ids to process ids.
   */
  const vector<int>& lastCrossingMap() const { return theLastCrossingMap->second; }

  /**
   * Access the crossing information.
   */
  vector<int>& lastCrossingMap() { return theLastCrossingMap->second; }  

  /**
   * Access the crossing information.
   */
  map<tStdXCombPtr,vector<int> >& crossingMap() { return processData()->crossingMap(); }

  /**
   * Access the crossing signs.
   */
  map<tStdXCombPtr,double>& crossingSigns() { return processData()->crossingSigns(); }

  /**
   * Generate the helicity combinations.
   */
  virtual set<vector<int> > generateHelicities() const;

  //@}

  /** @name Tree-level amplitudes */
  //@{

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return last evaluated helicity amplitudes.
   */
  const AmplitudeMap& lastAmplitudes() const { return theLastAmplitudes->second; }

  /**
   * Access the last evaluated helicity amplitudes.
   */
  AmplitudeMap& lastAmplitudes() { return theLastAmplitudes->second; }

  /**
   * Return last evaluated, leading colour helicity amplitudes.
   */
  const AmplitudeMap& lastLargeNAmplitudes() const { return theLastLargeNAmplitudes->second; }

  /**
   * Access the last evaluated, leading colour helicity amplitudes.
   */
  AmplitudeMap& lastLargeNAmplitudes() { return theLastLargeNAmplitudes->second; }

  /**
   * Return the matrix element squared.
   */
  virtual double me2() const {
    return 
      lastCrossingSign()*colourBasis()->me2(mePartonData(),lastAmplitudes());
  }

  /**
   * Return the colour correlated matrix element.
   */
  virtual double colourCorrelatedME2(pair<int,int> ij) const;

  /**
   * Return a positive helicity polarization vector for a gluon of
   * momentum p (with reference vector n) to be used when evaluating
   * spin correlations.
   */
  virtual LorentzVector<Complex> plusPolarization(const Lorentz5Momentum& p,
						  const Lorentz5Momentum& n) const;

  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;


  /**
   * Return true, if tree-level contributions will be evaluated at amplitude level.
   */
  virtual bool treeAmplitudes() const { return true; }

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t, const vector<int>&, Complex&) { return 0.; }

  //@}

  /** @name One-loop amplitudes */
  //@{

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return false; }

  /**
   * Return true, if this amplitude only provides
   * one-loop (QCD) corrections.
   */
  virtual bool onlyOneLoop() const { return false; }

  /**
   * Return true, if one-loop contributions will be evaluated at amplitude level.
   */
  virtual bool oneLoopAmplitudes() const { return true; }

  /**
   * Return true, if one loop corrections have been calculated in
   * dimensional reduction. Otherwise conventional dimensional
   * regularization is assumed. Note that renormalization is always
   * assumed to be MSbar.
   */
  virtual bool isDR() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of the integrated dipoles.
   */
  virtual bool isCS() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */
  virtual bool isBDK() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of everything expanded.
   */
  virtual bool isExpanded() const { return false; }

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return 0.*GeV2; }

  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const { return 0.; }

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const { return 0.; }

  /**
   * Calculate the one-loop amplitudes for the phasespace point
   * stored in lastXComb, if provided.
   */
  virtual void prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return last evaluated one-loop helicity amplitudes.
   */
  const AmplitudeMap& lastOneLoopAmplitudes() const { return theLastOneLoopAmplitudes->second; }

  /**
   * Access the last evaluated one-loop helicity amplitudes.
   */
  AmplitudeMap& lastOneLoopAmplitudes() { return theLastOneLoopAmplitudes->second; }

  /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const {
    return 
      lastCrossingSign()*colourBasis()->interference(mePartonData(),
						     lastOneLoopAmplitudes(),lastAmplitudes());
  }

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluateOneLoop(size_t, const vector<int>&) { return 0.; }

  //@}

  /** @name Caching and helpers to setup amplitude objects. */
  //@{

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {
    calculateTrees = true;
    calculateLoops = true;
  }

  /**
   * Clone this amplitude.
   */
  Ptr<MatchboxAmplitude>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
  virtual void cloneDependencies(const std::string& prefix = "");

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /**
   * The additional random numbers requested by
   * this virtual correction.
   */
  vector<double> additionalRandomNumbers;

private:

  /**
   * Recursively generate helicities
   */
  void doGenerateHelicities(set<vector<int> >& res,
			    vector<int>& current,
			    size_t pos) const;

  /**
   * The Herwig++ StandardModel object
   */
  Ptr<StandardModel>::tcptr theStandardModel;

  /**
   * The process data object to be used
   */
  Ptr<ProcessData>::ptr theProcessData;

  /**
   * The number of light flavours to be used.
   */
  unsigned int theNLight;

  /**
   * The colour basis implementation to be used.
   */
  Ptr<ColourBasis>::ptr theColourBasis;

  /**
   * The dimensionality of the colour basis for the processes covered
   * by the colour basis.
   */
  size_t theColourBasisDim;

  /**
   * References to the amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  map<tStdXCombPtr,map<vector<int>,CVector> > theLastAmplitudeMap;

  /**
   * References to the leading N amplitude values which have been
   * contributing to the last call of prepareAmplitudes.
   */
  map<tStdXCombPtr,map<vector<int>,CVector> > theLastLargeNAmplitudeMap;

  /**
   * References to the one-loop amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  map<tStdXCombPtr,map<vector<int>,CVector> > theLastOneLoopAmplitudeMap;

  /**
   * References to the amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  map<tStdXCombPtr,map<vector<int>,CVector> >::iterator theLastAmplitudes;

  /**
   * References to the leading N amplitude values which have been
   * contributing to the last call of prepareAmplitudes.
   */
  map<tStdXCombPtr,map<vector<int>,CVector> >::iterator theLastLargeNAmplitudes;

  /**
   * References to the one-loop amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  map<tStdXCombPtr,map<vector<int>,CVector> >::iterator theLastOneLoopAmplitudes;

  /**
   * The crossing information as filled by the last call to
   * fillCrossingMap()
   */
  map<tStdXCombPtr,vector<int> >::iterator theLastCrossingMap;

  /**
   * The colour crossing information as filled by the last call to
   * fillCrossingMap()
   */
  map<tStdXCombPtr,map<size_t,size_t> >::iterator theLastAmplitudeToColourMap;

  /**
   * The colour crossing information as filled by the last call to
   * fillCrossingMap()
   */
  map<tStdXCombPtr,map<size_t,size_t> >::iterator theLastColourToAmplitudeMap;

  /**
   * The amplitude parton data.
   */
  map<tStdXCombPtr,cPDVector>::iterator theLastAmplitudePartonData;

  /**
   * The crossing sign.
   */
  double theLastCrossingSign;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxAmplitude & operator=(const MatchboxAmplitude &);

protected:

  /**
   * True, if tree amplitudes need to be recalculated.
   */
  bool calculateTrees;

  /**
   * True, if loop amplitudes need to be recalculated.
   */
  bool calculateLoops;

};

}

#endif /* HERWIG_MatchboxAmplitude_H */
