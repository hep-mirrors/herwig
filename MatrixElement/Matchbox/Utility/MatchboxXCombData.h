// -*- C++ -*-
//
// MatchboxXCombData.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxXCombData_H
#define Herwig_MatchboxXCombData_H
//
// This is the declaration of the MatchboxXCombData class.
//

// work around a Boost 1.64 bug where ublas headers would fail otherwise
#include <boost/version.hpp>
#if (BOOST_VERSION / 100 >= 1064)
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "ThePEG/MatrixElement/MEBase.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.fh"
#include "Herwig/Models/StandardModel/StandardModel.h"

#include "ThePEG/Persistency/PersistentOStream.fh"
#include "ThePEG/Persistency/PersistentIStream.fh"

namespace Herwig {

  using namespace ThePEG;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Define complex vector from boost::uBLAS
   */
  typedef boost::numeric::ublas::vector<Complex> CVector;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Define how amplitudes are stored
   */
  typedef map<vector<int>,CVector> AmplitudeMap;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Define amplitude iterators
   */
  typedef map<vector<int>,CVector>::iterator AmplitudeIterator;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Define amplitude const iterators
   */
  typedef map<vector<int>,CVector>::const_iterator AmplitudeConstIterator;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Matchbox extensions to StandardXComb
   */
  class MatchboxXCombData {

  public:

    /** @name Standard constructors and destructors. */
    //@{
    /**
     * Standard constructor.
     */
    explicit MatchboxXCombData(tMEPtr newME);

    /**
     * Default constructor.
     */
    MatchboxXCombData();

    /**
     * Destructor.
     */
    virtual ~MatchboxXCombData();

    //@}

  public:

    /**
     * Reset all cache flags
     */
    void flushCaches();

  public:

    /**
     * Get the factory
     */
    Ptr<MatchboxFactory>::tcptr factory() const;

    /**
     * Get the matrix element; may return null
     */
    Ptr<MatchboxMEBase>::tptr matchboxME() const;

    /**
     * Get the dipole; may return null
     */
    Ptr<SubtractionDipole>::tptr subtractionDipole() const;

    /**
     * The crossing information as filled by the last call to
     * fillCrossingMap()
     */
    const vector<int>& crossingMap() const { return theCrossingMap; }

    /**
     * The crossing information as filled by the last call to
     * fillCrossingMap()
     */
    vector<int>& crossingMap() { return theCrossingMap; }

    /**
     * The colour crossing information as filled by the last call to
     * fillCrossingMap()
     */
    const map<size_t,size_t>& amplitudeToColourMap() const { return theAmplitudeToColourMap; }

    /**
     * The colour crossing information as filled by the last call to
     * fillCrossingMap()
     */
    map<size_t,size_t>& amplitudeToColourMap() { return theAmplitudeToColourMap; }

    /**
     * The colour crossing information as filled by the last call to
     * fillCrossingMap()
     */
    const map<size_t,size_t>& colourToAmplitudeMap() const { return theColourToAmplitudeMap; }

    /**
     * The colour crossing information as filled by the last call to
     * fillCrossingMap()
     */
    map<size_t,size_t>& colourToAmplitudeMap() { return theColourToAmplitudeMap; }

    /**
     * The crossing sign as filled by the last call to
     * fillCrossingMap()
     */
    double crossingSign() const { return theCrossingSign; }

    /**
     * The crossing sign as filled by the last call to
     * fillCrossingMap()
     */
    void crossingSign(double c) { theCrossingSign = c; }

    /**
     * The last renormalization scale
     */
    Energy2 lastRenormalizationScale() const { return theLastRenormalizationScale; }

    /**
     * The last renormalization scale
     */
    void lastRenormalizationScale(Energy2 lrs) { theLastRenormalizationScale = lrs; }

    /**
     * The amplitude parton data.
     */
    const cPDVector& amplitudePartonData() const { return theAmplitudePartonData; }

    /**
     * The amplitude parton data.
     */
    cPDVector& amplitudePartonData() { return theAmplitudePartonData; }

    /**
     * The crossed momenta
     */
    const vector<Lorentz5Momentum>& amplitudeMomenta() const { return theAmplitudeMomenta; }

    /**
     * The crossed momenta
     */
    vector<Lorentz5Momentum>& amplitudeMomenta() { return theAmplitudeMomenta; }

    /**
     * True, if the the tree level amplitudes need to be calculated
     */
    bool calculateTreeAmplitudes() const { return theCalculateTreeAmplitudes; }

    /**
     * The amplitude values which have been contributing
     * to the last call of prepareAmplitudes.
     */
    const map<vector<int>,CVector>& lastAmplitudes() const { return theLastAmplitudes; }

    /**
     * True, if the the tree level amplitudes need to be calculated
     */
    void haveTreeAmplitudes(bool f = true) { theCalculateTreeAmplitudes = !f; }

    /**
     * The amplitude values which have been contributing
     * to the last call of prepareAmplitudes.
     */
    map<vector<int>,CVector>& lastAmplitudes() { return theLastAmplitudes; }

    /**
     * The leading N amplitude values which have been
     * contributing to the last call of prepareAmplitudes.
     */
    const map<vector<int>,CVector>& lastLargeNAmplitudes() const { return theLastLargeNAmplitudes; }

    /**
     * The leading N amplitude values which have been
     * contributing to the last call of prepareAmplitudes.
     */
    map<vector<int>,CVector>& lastLargeNAmplitudes() { return theLastLargeNAmplitudes; }

    /**
     * True, if the the one-loop amplitudes need to be calculated
     */
    bool calculateOneLoopAmplitudes() const { return theCalculateOneLoopAmplitudes; }

    /**
     * The one-loop amplitude values which have been contributing
     * to the last call of prepareAmplitudes.
     */
    const map<vector<int>,CVector>& lastOneLoopAmplitudes() const { return theLastOneLoopAmplitudes; }

    /**
     * True, if the the one-loop amplitudes need to be calculated
     */
    void haveOneLoopAmplitudes(bool f = true) { theCalculateOneLoopAmplitudes = !f; }

    /**
     * The one-loop amplitude values which have been contributing
     * to the last call of prepareAmplitudes.
     */
    map<vector<int>,CVector>& lastOneLoopAmplitudes() { return theLastOneLoopAmplitudes; }

    /**
     * True, if the tree-level matrix element squared needs to be
     * calculated.
     */
    bool calculateTreeME2() const { return theCalculateTreeME2; }

    /**
     * The last tree-level matrix element squared
     */
    double lastTreeME2() const { return theLastTreeME2; }

    /**
     * The last tree-level matrix element squared
     */
    void lastTreeME2(double v) { 
      theLastTreeME2 = v; theCalculateTreeME2 = false;
    }

    /**
     * True, if the tree-level matrix element squared needs to be
     * calculated.
     */
    bool calculateLargeNME2() const { return theCalculateLargeNME2; }

    /**
     * The last tree-level matrix element squared
     */
    double lastLargeNME2() const { return theLastLargeNME2; }

    /**
     * The last tree-level matrix element squared
     */
    void lastLargeNME2(double v) { 
      theLastLargeNME2 = v; theCalculateLargeNME2 = false;
    }

    /**
     * True, if the one-loop/tree-level interference.
     * be calculated.
     */
    bool calculateOneLoopInterference() const { return theCalculateOneLoopInterference; }

    /**
     * The last one-loop/tree-level interference.
     */
    double lastOneLoopInterference() const { return theLastOneLoopInterference; }

    /**
     * The last one-loop/tree-level interference.
     */
    void lastOneLoopInterference(double v) { 
      theLastOneLoopInterference = v; theCalculateOneLoopInterference = false;
    }

    /**
     * True, if the one-loop/tree-level interference.
     * be calculated.
     */
    bool calculateOneLoopPoles() const { return theCalculateOneLoopPoles; }

    /**
     * The last one-loop/tree-level interference.
     */
    pair<double,double> lastOneLoopPoles() const { return theLastOneLoopPoles; }

    /**
     * The last one-loop/tree-level interference.
     */
    void lastOneLoopPoles(pair<double,double> v) { 
      theLastOneLoopPoles = v; theCalculateOneLoopPoles = false;
    }

    /**
     * True, if the indexed colour correlated matrix element needs to be
     * calculated.
     */
    bool calculateColourCorrelator(pair<int,int> ij) const {
      if ( ij.first > ij.second )
	swap(ij.first,ij.second);
      map<pair<int,int>,bool>::const_iterator f =
	theCalculateColourCorrelators.find(ij);
      if ( f == theCalculateColourCorrelators.end() )
	return true;
      return f->second;
    }

    /**
     * The colour correlated matrix element.
     */
    double lastColourCorrelator(pair<int,int> ij) const {
      if ( ij.first > ij.second )
	swap(ij.first,ij.second);
      map<pair<int,int>,double>::const_iterator v =
	theColourCorrelators.find(ij);
      if ( v == theColourCorrelators.end() )
	return 0.;
      return v->second;
    }

    /**
     * The colour correlated matrix element.
     */
    void lastColourCorrelator(pair<int,int> ij, double v) {
      if ( ij.first > ij.second )
	swap(ij.first,ij.second);
      theColourCorrelators[ij] = v;
      theCalculateColourCorrelators[ij] = false;
    }

    /**
     * True, if the indexed large-N colour correlated matrix element needs to be
     * calculated.
     */
    bool calculateLargeNColourCorrelator(pair<int,int> ij) const {
      if ( ij.first > ij.second )
	swap(ij.first,ij.second);
      map<pair<int,int>,bool>::const_iterator f =
	theCalculateLargeNColourCorrelators.find(ij);
      if ( f == theCalculateLargeNColourCorrelators.end() )
	return true;
      return f->second;
    }

    /**
     * The large-N colour correlated matrix element.
     */
    double lastLargeNColourCorrelator(pair<int,int> ij) const {
      if ( ij.first > ij.second )
	swap(ij.first,ij.second);
      map<pair<int,int>,double>::const_iterator v =
	theLargeNColourCorrelators.find(ij);
      if ( v == theLargeNColourCorrelators.end() )
	return 0.;
      return v->second;
    }

    /**
     * The large-N colour correlated matrix element.
     */
    void lastLargeNColourCorrelator(pair<int,int> ij, double v) {
      if ( ij.first > ij.second )
	swap(ij.first,ij.second);
      theLargeNColourCorrelators[ij] = v;
      theCalculateLargeNColourCorrelators[ij] = false;
    }

    /**
     * True, if the indexed colour/spin correlated matrix element needs to be
     * calculated.
     */
    bool calculateColourSpinCorrelator(const pair<int,int>& ij) const {
      map<pair<int,int>,bool>::const_iterator f =
	theCalculateColourSpinCorrelators.find(ij);
      if ( f == theCalculateColourSpinCorrelators.end() )
	return true;
      return f->second;
    }

    /**
     * The colour/spin correlated matrix element.
     */
    Complex lastColourSpinCorrelator(const pair<int,int>& ij) const {
      map<pair<int,int>,Complex>::const_iterator v =
	theColourSpinCorrelators.find(ij);
      if ( v == theColourSpinCorrelators.end() )
	return 0.;
      return v->second;
    }

    /**
     * The colour/spin correlated matrix element.
     */
    void lastColourSpinCorrelator(const pair<int,int>& ij, Complex v) {
      theColourSpinCorrelators[ij] = v;
      theCalculateColourSpinCorrelators[ij] = false;
    }

    /**
     * True, if the indexed spin correlated matrix element needs to be
     * calculated.
     */
    bool calculateSpinCorrelator(const pair<int,int>& ij) const {
      map<pair<int,int>,bool>::const_iterator f =
	theCalculateSpinCorrelators.find(ij);
      if ( f == theCalculateSpinCorrelators.end() )
	return true;
      return f->second;
    }

    /**
     * The spin correlated matrix element.
     */
    Complex lastSpinCorrelator(const pair<int,int>& ij) const {
      map<pair<int,int>,Complex>::const_iterator v =
	theSpinCorrelators.find(ij);
      if ( v == theSpinCorrelators.end() )
	return 0.;
      return v->second;
    }

    /**
     * The spin correlated matrix element.
     */
    void lastSpinCorrelator(const pair<int,int>& ij, Complex v) {
      theSpinCorrelators[ij] = v;
      theCalculateSpinCorrelators[ij] = false;
    }

    /**
     * Return the number of light flavours to be considered for this process.
     */
    unsigned int nLight() const { return theNLight; }

    /**
     * Set the number of light flavours to be considered for this process.
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
     * Get the dimensionality of the colour basis for this process.
     */
    size_t colourBasisDim() const { return theColourBasisDim; }

    /**
     * Set the dimensionality of the colour basis for this process.
     */
    void colourBasisDim(size_t d) { theColourBasisDim = d; }

    /**
     * Return the number of degrees of freedom required by the phase space generator
     */
    int nDimPhasespace() const { return theNDimPhasespace; }

    /**
     * Set the number of degrees of freedom required by the phase space generator
     */
    void nDimPhasespace(int d) { theNDimPhasespace = d; }

    /**
     * Return the number of degrees of freedom required by the amplitude
     */
    int nDimAmplitude() const { return theNDimAmplitude; }

    /**
     * Set the number of degrees of freedom required by the amplitude
     */
    void nDimAmplitude(int d) { theNDimAmplitude = d; }

    /**
     * Return the number of degrees of freedom required by the insertion operators
     */
    int nDimInsertions() const { return theNDimInsertions; }

    /**
     * Set the number of degrees of freedom required by the insertion operators
     */
    void nDimInsertions(int d) { theNDimInsertions = d; }

    /**
     * Get the additional random numbers required by the amplitude
     */
    const vector<double>& amplitudeRandomNumbers() const { return theAmplitudeRandomNumbers; }

    /**
     * Access the additional random numbers required by the amplitude
     */
    vector<double>& amplitudeRandomNumbers() { return theAmplitudeRandomNumbers; }

    /**
     * Get the additional random numbers required by the insertion operator
     */
    const vector<double>& insertionRandomNumbers() const { return theInsertionRandomNumbers; }

    /**
     * Access the additional random numbers required by the insertion operator
     */
    vector<double>& insertionRandomNumbers() { return theInsertionRandomNumbers; }

    /**
     * Return the diagram weights indexed by diagram id.
     */
    const map<int,double>& diagramWeights() const { return theDiagramWeights; }

    /**
     * Access the diagram weights indexed by diagram id.
     */
    map<int,double>& diagramWeights() { return theDiagramWeights; }

    /**
     * Return the singular limits
     */
    const set<pair<size_t,size_t> >& singularLimits() const { return theSingularLimits; }

    /**
     * Access the singular limits
     */
    set<pair<size_t,size_t> >& singularLimits() { return theSingularLimits; }

    /**
     * Return the last matched singular limit.
     */
    const set<pair<size_t,size_t> >::const_iterator& lastSingularLimit() const { return theLastSingularLimit; }

    /**
     * Access the last matched singular limit.
     */
    set<pair<size_t,size_t> >::const_iterator& lastSingularLimit() { return theLastSingularLimit; }

    /**
     * Set the Herwig StandardModel object
     */
    void hwStandardModel(Ptr<StandardModel>::tcptr sm) { theStandardModel = sm; }

    /**
     * Get the Herwig StandardModel object
     */
    Ptr<StandardModel>::tcptr hwStandardModel() const { return theStandardModel; }

    /**
     * Return the symmetry factor
     */
    double symmetryFactor() const { return theSymmetryFactor; }

    /**
     * Set the symmetry factor
     */
    void symmetryFactor(double f) { theSymmetryFactor = f; }

    /**
     * Return the OLP process ids
     */
    const vector<int>& olpId() const { return theOLPId; }

    /**
     * Set the OLP process ids
     */
    void olpId(int pType, int id) {
      if ( theOLPId.empty() )
	theOLPId.resize(4,0);
      theOLPId[pType] = id;
    }

    /**
     * Set the OLP process ids
     */
    void olpId(const vector<int>& id) { 
      theOLPId = id;
    }

    /**
     * Return the olp momentum vector
     */
    double* olpMomenta() { return theOLPMomenta; }

    /**
     * Fill the olp momentum vector
     */
    void fillOLPMomenta(const vector<Lorentz5Momentum>& mm,
			const cPDVector& mePartonData = cPDVector(),
			const map<long,Energy>& reshuffleMap = map<long,Energy>());

    /**
     * Helper struct to define the reshuffling equation
     */
    struct ReshuffleEquation {

      ReshuffleEquation(Energy xq,
			cPDVector::const_iterator xdBegin,
			cPDVector::const_iterator xdEnd,
			vector<Lorentz5Momentum>::const_iterator xmBegin,
			const map<long,Energy>* xreshuffleMap)
	: q(xq), dBegin(xdBegin), dEnd(xdEnd), mBegin(xmBegin),
	  reshuffleMap(xreshuffleMap) {}

      typedef double ArgType;
      typedef double ValType;

      static double aUnit() { return 1.; }
      static double vUnit() { return 1.; }

      double operator() (double xi) const;

      Energy q;
      cPDVector::const_iterator dBegin;
      cPDVector::const_iterator dEnd;
      vector<Lorentz5Momentum>::const_iterator mBegin;
      const map<long,Energy>* reshuffleMap;

    };

    /**
     * Perform a reshuffling to the mass shells contained in the map
     */
    void reshuffle(vector<Lorentz5Momentum>& momenta,
		   const cPDVector& mePartonData,
		   const map<long,Energy>& reshuffleMap) const;

    /**
     * Return a generic process id to communicate with external codes
     */
    int externalId() const { return theExternalId; }

    /**
     * Set a generic process id to communicate with external codes
     */
    void externalId(int id) { theExternalId = id; }

    /**
     * True if the process has been initialized
     */
    bool initialized() const { return theInitialized; }

    /**
     * True if the process has been initialized
     */
    void isInitialized(bool is = true) { theInitialized = is; }

    /**
     * Return the external momentum vector
     */
    const vector<double*>& externalMomenta() const { return theExternalMomenta; }

    /**
     * Fill the external momentum vector
     */
    void fillExternalMomenta(const vector<Lorentz5Momentum>&);
    
    /**
     * Caching for the external madgraph colour structures
     */
    const map<vector<int>,vector < complex<double> > >& heljamp() const { return theHelJamp; }
    
    /**
     * Caching for the external madgraph colour structures (large N)
     */
    const map<vector<int>,vector < complex<double> > >& helLNjamp() const { return theLNHelJamp; }
    
    /**
     *  pushback the madgraph colour structures
     */
    void pushheljamp(const vector<int>& hel, const complex<double>& jamp) { theHelJamp[hel].push_back(jamp); }
    
    /**
     * clear the madgraph colour structures
     */
    void clearheljamp() { theHelJamp.clear(); }
        
    /**
     *  pushback the madgraph colour structures (large N)
     */
    void pushhelLNjamp(const vector<int>& hel, const complex<double>& jamp) { theLNHelJamp[hel].push_back(jamp); }
    
    /**
     * clear the madgraph colour structures (large N)
     */
    void clearhelLNjamp() { theLNHelJamp.clear(); }

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

    /**
     * Put a CVector to the persistent ostream
     */
    static void putCVector(PersistentOStream&, const CVector&);

    /**
     * Get a CVector from the persistent istream
     */
    static void getCVector(PersistentIStream&, CVector&);

    /**
     * Put an amplitude map to the persistent ostream
     */
    static void putAmplitudeMap(PersistentOStream&, const map<vector<int>,CVector>&);

    /**
     * Get an amplitude map from the persistent istream
     */
    static void getAmplitudeMap(PersistentIStream&, map<vector<int>,CVector>&);
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
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    MatchboxXCombData & operator=(const MatchboxXCombData &) = delete;

  private:

    /**
     * The matrix element
     */
    Ptr<MatchboxMEBase>::tptr theMatchboxME;

    /**
     * The dipole
     */
    Ptr<SubtractionDipole>::tptr theSubtractionDipole;

    /**
     * The crossing information as filled by the last call to
     * fillCrossingMap()
     */
    vector<int> theCrossingMap;

    /**
     * The colour crossing information as filled by the last call to
     * fillCrossingMap()
     */
    map<size_t,size_t> theAmplitudeToColourMap;

    /**
     * The colour crossing information as filled by the last call to
     * fillCrossingMap()
     */
    map<size_t,size_t> theColourToAmplitudeMap;

    /**
     * The crossing sign as filled by the last call to
     * fillCrossingMap()
     */
    double theCrossingSign;

    /**
     * The last renormalization scale
     */
    Energy2 theLastRenormalizationScale;

    /**
     * The amplitude parton data.
     */
    cPDVector theAmplitudePartonData;

    /**
     * The ccrossed momenta
     */
    vector<Lorentz5Momentum> theAmplitudeMomenta;

    /**
     * True, if the the tree level amplitudes need to be calculated
     */
    bool theCalculateTreeAmplitudes;

    /**
     * The amplitude values which have been contributing
     * to the last call of prepareAmplitudes.
     */
    map<vector<int>,CVector> theLastAmplitudes;

    /**
     * The leading N amplitude values which have been
     * contributing to the last call of prepareAmplitudes.
     */
    map<vector<int>,CVector> theLastLargeNAmplitudes;

    /**
     * True, if the the one-loop amplitudes need to be calculated
     */
    bool theCalculateOneLoopAmplitudes;

    /**
     * The one-loop amplitude values which have been contributing
     * to the last call of prepareAmplitudes.
     */
    map<vector<int>,CVector> theLastOneLoopAmplitudes;

    /**
     * True, if the tree-level matrix element squared needs to be
     * calculated.
     */
    bool theCalculateTreeME2;

    /**
     * The last tree-level matrix element squared
     */
    double theLastTreeME2;

    /**
     * True, if the tree-level matrix element squared needs to be
     * calculated.
     */
    bool theCalculateLargeNME2;

    /**
     * The last tree-level matrix element squared
     */
    double theLastLargeNME2;

    /**
     * True, if the one-loop/tree-level interference.
     * be calculated.
     */
    bool theCalculateOneLoopInterference;

    /**
     * The last one-loop/tree-level interference.
     */
    double theLastOneLoopInterference;

    /**
     * True, if the one-loop/tree-level interference.
     * be calculated.
     */
    bool theCalculateOneLoopPoles;

    /**
     * The last one-loop/tree-level interference.
     */
    pair<double,double> theLastOneLoopPoles;

    /**
     * True, if the indexed colour correlated matrix element needs to be
     * calculated.
     */
    map<pair<int,int>,bool> theCalculateColourCorrelators;

    /**
     * The colour correlated matrix element.
     */
    map<pair<int,int>,double> theColourCorrelators;

    /**
     * True, if the indexed large-N colour correlated matrix element needs to be
     * calculated.
     */
    map<pair<int,int>,bool> theCalculateLargeNColourCorrelators;

    /**
     * The large-N colour correlated matrix element.
     */
    map<pair<int,int>,double> theLargeNColourCorrelators;

    /**
     * True, if the indexed colour/spin correlated matrix element needs to be
     * calculated.
     */
    map<pair<int,int>,bool> theCalculateColourSpinCorrelators;

    /**
     * The colour/spin correlated matrix element.
     */
    map<pair<int,int>,Complex> theColourSpinCorrelators;

    /**
     * True, if the indexed spin correlated matrix element needs to be
     * calculated.
     */
    map<pair<int,int>,bool> theCalculateSpinCorrelators;

    /**
     * The spin correlated matrix element.
     */
    map<pair<int,int>,Complex> theSpinCorrelators;

    /**
     * The number of light flavours to be considered for this process.
     */
    static unsigned int theNLight;

    /**
     * Vector with the PDG ids of the light quark flavours,
     * which are contained in the jet particle group.
     */
    static vector<long> theNLightJetVec;

    /**
     * Vector with the PDG ids of the heavy quark flavours,
     * which are contained in the jet particle group.
     */
    static vector<long> theNHeavyJetVec;

    /**
     * Vector with the PDG ids of the light quark flavours,
     * which are contained in the proton particle group.
     */
    static vector<long> theNLightProtonVec;

    /**
     * The dimensionality of the colour basis for this process.
     */
    size_t theColourBasisDim;

    /**
     * The number of degrees of freedom required by the phase space generator
     */
    int theNDimPhasespace;

    /**
     * The number of degrees of freedom required by the amplitude
     */
    int theNDimAmplitude;

    /**
     * The number of degrees of freedom required by the insertion operators
     */
    int theNDimInsertions;

    /**
     * Additional random numbers required by the amplitude
     */
    vector<double> theAmplitudeRandomNumbers;

    /**
     * Additional random numbers required by the insertion operator
     */
    vector<double> theInsertionRandomNumbers;

    /**
     * The diagram weights indexed by diagram id.
     */
    map<int,double> theDiagramWeights;

    /**
     * If not empty, the entries here serve to limit phasespace
     * generation to the corresponding collinear limits, or soft limits
     * if both pair entries are equal.
     */
    set<pair<size_t,size_t> > theSingularLimits;

    /**
     * The last matched singular limit.
     */
    set<pair<size_t,size_t> >::const_iterator theLastSingularLimit;

    /**
     * The Herwig StandardModel object
     */
    Ptr<StandardModel>::tcptr theStandardModel;

    /**
     * The symmetry factor
     */
    double theSymmetryFactor;

    /**
     * The OLP process id
     */
    vector<int> theOLPId;

    /**
     * Return the olp momentum vector
     */
    double* theOLPMomenta;

    /**
     * True, if olp momenta have been filled
     */
    bool filledOLPMomenta;

    /**
     * A generic process id to communicate with external codes
     */
    int theExternalId;

    /**
     * True if the process has been initialized
     */
    bool theInitialized;

    /**
     * The external momenta
     */
    vector<double*> theExternalMomenta;

    /**
     * True, if external momenta have been filled
     */
    bool filledExternalMomenta;
    
    /**
     * caching of different colour structures (MadGraph-Interface)
     */
    map<vector<int>,vector < complex<double> > > theHelJamp;    
	
    /**
     * caching of different colour structures (MadGraph-Interface)
     */
    map<vector<int>,vector < complex<double> > > theLNHelJamp;

  };

}

#endif /* Herwig_MatchboxXCombData_H */
