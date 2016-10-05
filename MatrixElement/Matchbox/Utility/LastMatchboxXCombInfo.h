// -*- C++ -*-
//
// MatchboxXComb.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_LastMatchboxXCombInfo_H
#define Herwig_LastMatchboxXCombInfo_H
//
// This is the declaration of the MatchboxXComb class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/MatchboxXCombData.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Provide easy access to MatchboxXComb XComb extensions
 */
class LastMatchboxXCombInfo {

public:

  /**
   * Default constructor
   */
  LastMatchboxXCombInfo()
    : theLastMatchboxXComb(0), theLastHeadMatchboxXComb(0) {}

  /**
   * Return a pointer to the last selected XComb.
   */
  MatchboxXCombData* lastMatchboxXComb() const { return theLastMatchboxXComb; }

  /**
   * If the last selected XComb object belongs to a
   * group of XComb's return a pointer to the head 
   * XComb object for this group.
   */
  MatchboxXCombData* lastHeadMatchboxXComb() const { return theLastHeadMatchboxXComb; }

public:

  /**
   * The crossing information as filled by the last call to
   * fillCrossingMap()
   */
  const vector<int>& crossingMap() const { return lastMatchboxXComb()->crossingMap(); }

  /**
   * The colour crossing information as filled by the last call to
   * fillCrossingMap()
   */
  const map<size_t,size_t>& amplitudeToColourMap() const { return lastMatchboxXComb()->amplitudeToColourMap(); }

  /**
   * The colour crossing information as filled by the last call to
   * fillCrossingMap()
   */
  const map<size_t,size_t>& colourToAmplitudeMap() const { return lastMatchboxXComb()->colourToAmplitudeMap(); }

  /**
   * The crossing sign as filled by the last call to
   * fillCrossingMap()
   */
  double crossingSign() const { return lastMatchboxXComb()->crossingSign(); }

  /**
   * The last renormalization scale
   */
  Energy2 lastRenormalizationScale() const { return lastMatchboxXComb()->lastRenormalizationScale(); }

  /**
   * The amplitude parton data.
   */
  const cPDVector& amplitudePartonData() const { return lastMatchboxXComb()->amplitudePartonData(); }

  /**
   * The crossed momenta
   */
  const vector<Lorentz5Momentum>& amplitudeMomenta() const { return lastMatchboxXComb()->amplitudeMomenta(); }

  /**
   * True, if the the tree level amplitudes need to be calculated
   */
  bool calculateTreeAmplitudes() const { return lastMatchboxXComb()->calculateTreeAmplitudes(); }

  /**
   * The amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  const map<vector<int>,CVector>& lastAmplitudes() const { return lastMatchboxXComb()->lastAmplitudes(); }

  /**
   * The leading N amplitude values which have been
   * contributing to the last call of prepareAmplitudes.
   */
  const map<vector<int>,CVector>& lastLargeNAmplitudes() const { return lastMatchboxXComb()->lastLargeNAmplitudes(); }

  /**
   * True, if the the one-loop amplitudes need to be calculated
   */
  bool calculateOneLoopAmplitudes() const { return lastMatchboxXComb()->calculateOneLoopAmplitudes(); }

  /**
   * The one-loop amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  const map<vector<int>,CVector>& lastOneLoopAmplitudes() const { return lastMatchboxXComb()->lastOneLoopAmplitudes(); }

  /**
   * True, if the tree-level matrix element squared needs to be
   * calculated.
   */
  bool calculateTreeME2() const { return lastMatchboxXComb()->calculateTreeME2(); }

  /**
   * The last tree-level matrix element squared
   */
  double lastTreeME2() const { return lastMatchboxXComb()->lastTreeME2(); }

  /**
   * True, if the tree-level matrix element squared needs to be
   * calculated.
   */
  bool calculateLargeNME2() const { return lastMatchboxXComb()->calculateLargeNME2(); }

  /**
   * The last tree-level matrix element squared
   */
  double lastLargeNME2() const { return lastMatchboxXComb()->lastLargeNME2(); }

  /**
   * True, if the one-loop/tree-level interference.
   * be calculated.
   */
  bool calculateOneLoopInterference() const { return lastMatchboxXComb()->calculateOneLoopInterference(); }

  /**
   * The last one-loop/tree-level interference.
   */
  double lastOneLoopInterference() const { return lastMatchboxXComb()->lastOneLoopInterference(); }

  /**
   * True, if the one-loop/tree-level interference.
   * be calculated.
   */
  bool calculateOneLoopPoles() const { return lastMatchboxXComb()->calculateOneLoopPoles(); }

  /**
   * The last one-loop/tree-level interference.
   */
  pair<double,double> lastOneLoopPoles() const { return lastMatchboxXComb()->lastOneLoopPoles(); }


  /**
   * True, if the indexed colour correlated matrix element needs to be
   * calculated.
   */
  bool calculateColourCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->calculateColourCorrelator(ij); }

  /**
   * The colour correlated matrix element.
   */
  double lastColourCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->lastColourCorrelator(ij); }

  /**
   * True, if the indexed large-N colour correlated matrix element needs to be
   * calculated.
   */
  bool calculateLargeNColourCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->calculateLargeNColourCorrelator(ij); }

  /**
   * The large-N colour correlated matrix element.
   */
  double lastLargeNColourCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->lastLargeNColourCorrelator(ij); }

  /**
   * True, if the indexed colour/spin correlated matrix element needs to be
   * calculated.
   */
  bool calculateColourSpinCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->calculateColourSpinCorrelator(ij); }

  /**
   * The colour/spin correlated matrix element.
   */
  Complex lastColourSpinCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->lastColourSpinCorrelator(ij); }

  /**
   * True, if the indexed spin correlated matrix element needs to be
   * calculated.
   */
  bool calculateSpinCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->calculateSpinCorrelator(ij); }

  /**
   * The spin correlated matrix element.
   */
  Complex lastSpinCorrelator(const pair<int,int>& ij) const { return lastMatchboxXComb()->lastSpinCorrelator(ij); }

  /**
   * Return the number of light flavours to be considered for this process.
   */
  unsigned int nLight() const { return lastMatchboxXComb()->nLight(); }

  /**
   * Return the vector that contains the PDG ids of 
   * the light flavours, which are contained in the
   * jet particle group.
   */
  vector<int> nLightJetVec() const { return lastMatchboxXComb()->nLightJetVec(); }

  /**
   * Return the vector that contains the PDG ids of 
   * the heavy flavours, which are contained in the
   * jet particle group.
   */
  vector<int> nHeavyJetVec() const { return lastMatchboxXComb()->nHeavyJetVec(); }

  /**
   * Return the vector that contains the PDG ids of 
   * the light flavours, which are contained in the
   * proton particle group.
   */
  vector<int> nLightProtonVec() const { return lastMatchboxXComb()->nLightProtonVec(); }

  /**
   * Get the dimensionality of the colour basis for this process.
   */
  size_t colourBasisDim() const { return lastMatchboxXComb()->colourBasisDim(); }

  /**
   * Return the number of degrees of freedom required by the phase space generator
   */
  int nDimPhasespace() const { return lastMatchboxXComb()->nDimPhasespace(); }

  /**
   * Return the number of degrees of freedom required by the amplitude
   */
  int nDimAmplitude() const { return lastMatchboxXComb()->nDimAmplitude(); }

  /**
   * Return the number of degrees of freedom required by the insertion operators
   */
  int nDimInsertions() const { return lastMatchboxXComb()->nDimInsertions(); }

  /**
   * Get the additional random numbers required by the amplitude
   */
  const vector<double>& amplitudeRandomNumbers() const { return lastMatchboxXComb()->amplitudeRandomNumbers(); }

  /**
   * Get the additional random numbers required by the insertion operator
   */
  const vector<double>& insertionRandomNumbers() const { return lastMatchboxXComb()->insertionRandomNumbers(); }

  /**
   * Return the diagram weights indexed by diagram id.
   */
  const map<int,double>& diagramWeights() const { return lastMatchboxXComb()->diagramWeights(); }

  /**
   * Return the singular limits
   */
  const set<pair<size_t,size_t> >& singularLimits() const { return lastMatchboxXComb()->singularLimits(); }

  /**
   * Return the last matched singular limit.
   */
  const set<pair<size_t,size_t> >::const_iterator& lastSingularLimit() const { return lastMatchboxXComb()->lastSingularLimit(); }

  /**
   * Get the Herwig StandardModel object
   */
  Ptr<StandardModel>::tcptr hwStandardModel() const { return lastMatchboxXComb()->hwStandardModel(); }

  /**
   * Return the symmetry factor
   */
  double symmetryFactor() const { return lastMatchboxXComb()->symmetryFactor(); }
   
  /**
   * Return the OLP process id
   */
  const vector<int>& olpId() const { return lastMatchboxXComb()->olpId(); }

  /**
   * Return the olp momentum vector
   */
  double* olpMomenta() const { return lastMatchboxXComb()->olpMomenta(); }

  /**
   * Fill the olp momentum vector
   */
  void fillOLPMomenta(const vector<Lorentz5Momentum>& mm,
		      const cPDVector& mePartonData,
		      const map<long,Energy>& reshuffleMap) const { 
    lastMatchboxXComb()->fillOLPMomenta(mm,mePartonData,reshuffleMap);
  }

protected:

  /**
   * The crossing information as filled by the last call to
   * fillCrossingMap()
   */
  vector<int>& crossingMap() { return lastMatchboxXComb()->crossingMap(); }

  /**
   * The colour crossing information as filled by the last call to
   * fillCrossingMap()
   */
  map<size_t,size_t>& amplitudeToColourMap() { return lastMatchboxXComb()->amplitudeToColourMap(); }

  /**
   * The colour crossing information as filled by the last call to
   * fillCrossingMap()
   */
  map<size_t,size_t>& colourToAmplitudeMap() { return lastMatchboxXComb()->colourToAmplitudeMap(); }

  /**
   * The crossing sign as filled by the last call to
   * fillCrossingMap()
   */
  void crossingSign(double c) { lastMatchboxXComb()->crossingSign(c); }

  /**
   * The last renormalization scale
   */
  void lastRenormalizationScale(Energy2 lrs) { lastMatchboxXComb()->lastRenormalizationScale(lrs); }

  /**
   * The amplitude parton data.
   */
  cPDVector& amplitudePartonData() { return lastMatchboxXComb()->amplitudePartonData(); }

  /**
   * The crossed momenta
   */
  vector<Lorentz5Momentum>& amplitudeMomenta() { return lastMatchboxXComb()->amplitudeMomenta(); }

  /**
   * True, if the the tree level amplitudes need to be calculated
   */
  void haveTreeAmplitudes(bool f = true) { lastMatchboxXComb()->haveTreeAmplitudes(f); }

  /**
   * The amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  map<vector<int>,CVector>& lastAmplitudes() { return lastMatchboxXComb()->lastAmplitudes(); }

  /**
   * The leading N amplitude values which have been
   * contributing to the last call of prepareAmplitudes.
   */
  map<vector<int>,CVector>& lastLargeNAmplitudes() { return lastMatchboxXComb()->lastLargeNAmplitudes(); }

  /**
   * True, if the the one-loop amplitudes need to be calculated
   */
  void haveOneLoopAmplitudes(bool f = true) { lastMatchboxXComb()->haveOneLoopAmplitudes(f); }

  /**
   * The one-loop amplitude values which have been contributing
   * to the last call of prepareAmplitudes.
   */
  map<vector<int>,CVector>& lastOneLoopAmplitudes() { return lastMatchboxXComb()->lastOneLoopAmplitudes(); }

  /**
   * The last tree-level matrix element squared
   */
  void lastTreeME2(double v) const { lastMatchboxXComb()->lastTreeME2(v); }

  /**
   * The last tree-level matrix element squared
   */
  void lastLargeNME2(double v) const { lastMatchboxXComb()->lastLargeNME2(v); }

  /**
   * The last one-loop/tree-level interference.
   */
  void lastOneLoopInterference(double v) const { lastMatchboxXComb()->lastOneLoopInterference(v); }

  /**
   * The last one-loop/tree-level interference.
   */
  void lastOneLoopPoles(pair<double,double> v) const { lastMatchboxXComb()->lastOneLoopPoles(v); }

  /**
   * The colour correlated matrix element.
   */
  void lastColourCorrelator(const pair<int,int>& ij, double v) const { lastMatchboxXComb()->lastColourCorrelator(ij,v); }

  /**
   * The large-N colour correlated matrix element.
   */
  void lastLargeNColourCorrelator(const pair<int,int>& ij, double v) const { lastMatchboxXComb()->lastLargeNColourCorrelator(ij,v); }

  /**
   * The colour/spin correlated matrix element.
   */
  void lastColourSpinCorrelator(const pair<int,int>& ij, Complex v) const { lastMatchboxXComb()->lastColourSpinCorrelator(ij,v); }

  /**
   * The spin correlated matrix element.
   */
  void lastSpinCorrelator(const pair<int,int>& ij, Complex v) const { lastMatchboxXComb()->lastSpinCorrelator(ij,v); }

  /**
   * Set the number of light flavours to be considered for this process.
   */
  void nLight(unsigned int n) { lastMatchboxXComb()->nLight(n); }

  /**
   * Set the elements of the vector that contains the PDG
   * ids of the light flavours, which are contained in the
   * jet particle group.
   */
  void nLightJetVec(int n) { lastMatchboxXComb()->nLightJetVec(n); }

  /**
   * Set the elements of the vector that contains the PDG
   * ids of the heavy flavours, which are contained in the
   * jet particle group.
   */
  void nHeavyJetVec(int n) { lastMatchboxXComb()->nHeavyJetVec(n); }

  /**
   * Set the elements of the vector that contains the PDG
   * ids of the light flavours, which are contained in the
   * proton particle group.
   */
  void nLightProtonVec(int n) { lastMatchboxXComb()->nLightProtonVec(n); }

  /**
   * Set the dimensionality of the colour basis for this process.
   */
  void colourBasisDim(size_t d) { lastMatchboxXComb()->colourBasisDim(d); }

  /**
   * Set the number of degrees of freedom required by the phase space generator
   */
  void nDimPhasespace(int d) { lastMatchboxXComb()->nDimPhasespace(d); }

  /**
   * Set the number of degrees of freedom required by the amplitude
   */
  void nDimAmplitude(int d) { lastMatchboxXComb()->nDimAmplitude(d); }

  /**
   * Set the number of degrees of freedom required by the insertion operators
   */
  void nDimInsertions(int d) { lastMatchboxXComb()->nDimInsertions(d); }

  /**
   * Access the additional random numbers required by the amplitude
   */
  vector<double>& amplitudeRandomNumbers() { return lastMatchboxXComb()->amplitudeRandomNumbers(); }

  /**
   * Access the additional random numbers required by the insertion operator
   */
  vector<double>& insertionRandomNumbers() { return lastMatchboxXComb()->insertionRandomNumbers(); }

  /**
   * Access the diagram weights indexed by diagram id.
   */
  map<int,double>& diagramWeights() { return lastMatchboxXComb()->diagramWeights(); }

  /**
   * Access the singular limits
   */
  set<pair<size_t,size_t> >& singularLimits() { return lastMatchboxXComb()->singularLimits(); }

  /**
   * Access the last matched singular limit.
   */
  set<pair<size_t,size_t> >::const_iterator& lastSingularLimit() { return lastMatchboxXComb()->lastSingularLimit(); }

  /**
   * Set the Herwig StandardModel object
   */
  void hwStandardModel(Ptr<StandardModel>::tcptr sm) { lastMatchboxXComb()->hwStandardModel(sm); }

  /**
   * Set the symmetry factor
   */
  void symmetryFactor(double f) const { lastMatchboxXComb()->symmetryFactor(f); }

  /**
   * Set the OLP process id
   */
  void olpId(int pType, int id) { lastMatchboxXComb()->olpId(pType,id); }

protected:

  /**
   * Set the XComb pointer cast to MatchboxXComb
   */
  void lastMatchboxXComb(tStdXCombPtr xc) {
    theLastMatchboxXComb = xc ? &dynamic_cast<MatchboxXCombData&>(*xc) : 0;
    theLastHeadMatchboxXComb = 
      xc && xc->head() ? &dynamic_cast<MatchboxXCombData&>(*xc->head()) : 0;
  }

  /**
   * The XComb pointer cast to MatchboxXComb
   */
  MatchboxXCombData* theLastMatchboxXComb;

  /**
   * The head XComb pointer cast to MatchboxXComb
   */
  MatchboxXCombData* theLastHeadMatchboxXComb;

};

}

#endif // Herwig_LastMatchboxXCombInfo_H

