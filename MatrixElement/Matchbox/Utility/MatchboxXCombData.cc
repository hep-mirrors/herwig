// -*- C++ -*-
//
// MatchboxXCombData.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxXCombData class.
//

#include "MatchboxXCombData.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/Utilities/GSLBisection.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxXCombData::MatchboxXCombData() 
  : theCrossingSign(1.0), theCalculateTreeAmplitudes(true), 
    theCalculateOneLoopAmplitudes(true), theCalculateTreeME2(true), 
    theLastTreeME2(0.0), theCalculateLargeNME2(true), 
    theLastLargeNME2(0.0), theCalculateOneLoopInterference(true), 
    theLastOneLoopInterference(0.0), theCalculateOneLoopPoles(true), 
    theLastOneLoopPoles(0.0,0.0), 
    theColourBasisDim(0), theNDimPhasespace(0), 
    theNDimAmplitude(0), theNDimInsertions(0), 
    theSymmetryFactor(0.0), theOLPMomenta(0),
    filledOLPMomenta(false), theExternalId(0),
    theInitialized(false), filledExternalMomenta(false) {
  flushCaches();
}

unsigned int MatchboxXCombData::theNLight(0);
vector<long> MatchboxXCombData::theNLightJetVec=vector<long> ();
vector<long> MatchboxXCombData::theNHeavyJetVec=vector<long>() ;
vector<long> MatchboxXCombData::theNLightProtonVec=vector<long>() ;


MatchboxXCombData::~MatchboxXCombData() {
  if ( theOLPMomenta ) {
    delete[] theOLPMomenta;
    theOLPMomenta = 0;
  }
  for ( vector<double*>::iterator k = theExternalMomenta.begin();
	k != theExternalMomenta.end(); ++k ) {
    delete[] *k;
    *k = 0;
  }
  theExternalMomenta.clear();
}

MatchboxXCombData::MatchboxXCombData(tMEPtr newME) 
  : theCrossingSign(1.0), theCalculateTreeAmplitudes(true), 
    theCalculateOneLoopAmplitudes(true), theCalculateTreeME2(true), 
    theLastTreeME2(0.0), theCalculateLargeNME2(true), 
    theLastLargeNME2(0.0), theCalculateOneLoopInterference(true), 
    theLastOneLoopInterference(0.0), theCalculateOneLoopPoles(true), 
    theLastOneLoopPoles(0.0,0.0),
    theColourBasisDim(0), theNDimPhasespace(0), 
    theNDimAmplitude(0), theNDimInsertions(0), 
    theSymmetryFactor(0.0), theOLPMomenta(0),
    filledOLPMomenta(false), theExternalId(0),
    theInitialized(false), filledExternalMomenta(false) {
  flushCaches();
  theMatchboxME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(newME);
  theSubtractionDipole = dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(newME);
  Ptr<SubtractedME>::tptr subme = 
    dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(newME);
  if ( !theMatchboxME && !theSubtractionDipole && subme )
    theMatchboxME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(subme->head());
  assert(theMatchboxME || theSubtractionDipole);
  if ( theMatchboxME )
    theFactory = theMatchboxME->factory();
  else if ( theSubtractionDipole )
    theFactory = theSubtractionDipole->realEmissionME()->factory();
  assert(theFactory);
}

double MatchboxXCombData::ReshuffleEquation::operator() (double xi) const {

  double res = -q/GeV;

  cPDVector::const_iterator dit = dBegin;
  vector<Lorentz5Momentum>::const_iterator mit = mBegin;

  for ( ; dit != dEnd; ++dit, ++mit ) {
    map<long,Energy>::const_iterator rm =
      reshuffleMap->find((**dit).id());
    Energy2 tmass2 = 
      rm == reshuffleMap->end() ?
      mit->mass2() : sqr(rm->second);
    res += sqrt(sqr(xi)*(mit->vect().mag2()) + tmass2)/GeV;
  }

  return res;

}

void MatchboxXCombData::reshuffle(vector<Lorentz5Momentum>& momenta,
				  const cPDVector& mePartonData,
				  const map<long,Energy>& reshuffleMap) const {

  if ( momenta.size() == 3 ) // nothing to do; don't throw an exception
    return;

  bool needDoSomething = false;
  for ( cPDVector::const_iterator d = mePartonData.begin() + 2;
	d != mePartonData.end(); ++d )
    needDoSomething |= reshuffleMap.find((**d).id()) != reshuffleMap.end();

  if ( !needDoSomething )
    return;

  Lorentz5Momentum Q(ZERO,ZERO,ZERO,ZERO);
  for ( vector<Lorentz5Momentum>::const_iterator p = momenta.begin()+2;
	p != momenta.end(); ++p )
    Q += *p;

  Boost beta = Q.findBoostToCM();
  bool needBoost = (beta.mag2() > Constants::epsilon);

  if ( needBoost ) {
    for ( vector<Lorentz5Momentum>::iterator p = momenta.begin()+2;
	  p != momenta.end(); ++p )
      p->boost(beta);
  }

  ReshuffleEquation solve(Q.m(),
			  mePartonData.begin() + 2,
			  mePartonData.end(),
			  momenta.begin() + 2, &reshuffleMap);

  double xi = 1.;

  GSLBisection solver(1e-10,1e-8,10000);

  try {
    xi = solver.value(solve,0.0,1.1);
  } catch (GSLBisection::GSLerror) {
    throw Veto();
  } catch (GSLBisection::IntervalError) {
    throw Veto();
  }

  cPDVector::const_iterator d = mePartonData.begin() + 2;
  vector<Lorentz5Momentum>::iterator p = momenta.begin() + 2;
  for ( ; d != mePartonData.end(); ++d, ++p ) {
    p->setVect(xi*p->vect());
    map<long,Energy>::const_iterator mit = 
      reshuffleMap.find((**d).id());
    Energy2 newMass2 = mit == reshuffleMap.end() ? p->mass2() : sqr(mit->second);
    p->setE(sqrt(p->vect().mag2() + newMass2));
    p->rescaleMass();
  }

  if ( needBoost ) {
    for ( vector<Lorentz5Momentum>::iterator p = momenta.begin()+2;
	  p != momenta.end(); ++p )
      p->boost(-beta);
  }

}

void MatchboxXCombData::fillOLPMomenta(const vector<Lorentz5Momentum>& memomenta,
				       const cPDVector& mePartonData,
				       const map<long,Energy>& reshuffleMap) {
  if ( filledOLPMomenta )
    return;
  if ( !reshuffleMap.empty() && memomenta.size() > 3 ) {
    vector<Lorentz5Momentum> reshuffled = memomenta;
    reshuffle(reshuffled,mePartonData,reshuffleMap);
    fillOLPMomenta(reshuffled);
    return;
  }
  if ( !theOLPMomenta ) {
    theOLPMomenta = new double[5*memomenta.size()];
  }
  for ( size_t p = 0; p < memomenta.size(); ++p ) {
    theOLPMomenta[5*p] = memomenta[p].t()/GeV;
    theOLPMomenta[5*p+1] = memomenta[p].x()/GeV;
    theOLPMomenta[5*p+2] = memomenta[p].y()/GeV;
    theOLPMomenta[5*p+3] = memomenta[p].z()/GeV;
    theOLPMomenta[5*p+4] = memomenta[p].mass()/GeV;
  }
  filledOLPMomenta = true;
}

void MatchboxXCombData::fillExternalMomenta(const vector<Lorentz5Momentum>& memomenta) {
  if ( filledExternalMomenta )
    return;
  if ( theExternalMomenta.empty() ) {
    theExternalMomenta.resize(memomenta.size());
    for ( size_t k = 0; k < memomenta.size(); ++k )
      theExternalMomenta[k] = new double[4];
  }
  for ( size_t p = 0; p < memomenta.size(); ++p ) {
    theExternalMomenta[p][0] = memomenta[p].t()/GeV;
    theExternalMomenta[p][1] = memomenta[p].x()/GeV;
    theExternalMomenta[p][2] = memomenta[p].y()/GeV;
    theExternalMomenta[p][3] = memomenta[p].z()/GeV;
  }
  filledExternalMomenta = true;
}

Ptr<MatchboxFactory>::tcptr MatchboxXCombData::factory() const {
  return theFactory;
}

Ptr<MatchboxMEBase>::tptr MatchboxXCombData::matchboxME() const {
  return theMatchboxME;
}

Ptr<SubtractionDipole>::tptr MatchboxXCombData::subtractionDipole() const {
  return theSubtractionDipole;
}

void MatchboxXCombData::flushCaches() {
  theCalculateTreeAmplitudes = true;
  theCalculateOneLoopAmplitudes = true;
  theCalculateTreeME2 = true;
  theCalculateLargeNME2 = true;
  theCalculateOneLoopInterference = true;
  theCalculateOneLoopPoles = true;
  for ( auto & f : theCalculateColourCorrelators )
    f.second = true;
  for ( auto & f : theCalculateLargeNColourCorrelators )
    f.second = true;
  for ( auto & f : theCalculateColourSpinCorrelators )
    f.second = true;
  for ( auto & f : theCalculateSpinCorrelators )
    f.second = true;
  filledOLPMomenta = false;
  filledExternalMomenta = false;
    //theLastAmplitudes.clear();
    //theLastLargeNAmplitudes.clear();
    //theLastOneLoopAmplitudes.clear();
  theColourCorrelators.clear();
  theLargeNColourCorrelators.clear();
  theColourSpinCorrelators.clear();
  theSpinCorrelators.clear();
  theAmplitudeRandomNumbers.clear();
  theInsertionRandomNumbers.clear();
  theDiagramWeights.clear();
  theHelJamp.clear();
  theLNHelJamp.clear();
}

void MatchboxXCombData::putCVector(PersistentOStream& os, const CVector& v) {
  size_t n = v.size();
  os << n;
  for ( size_t k = 0; k < n; k++ )
    os << v(k);
}

void MatchboxXCombData::getCVector(PersistentIStream& is, CVector& v) {
  size_t n; is >> n;
  v.resize(n);
  Complex value;
  for ( size_t k = 0; k < n; k++ ) {
    is >> value; v(k) = value;
  }
}

void MatchboxXCombData::putAmplitudeMap(PersistentOStream& os, const map<vector<int>,CVector>& amps) {
  os << amps.size();
  for ( map<vector<int>,CVector>::const_iterator a = amps.begin();
	a != amps.end(); ++a ) {
    os << a->first;
    putCVector(os,a->second);
  }
}

void MatchboxXCombData::getAmplitudeMap(PersistentIStream& is, map<vector<int>,CVector>& amps) {
  size_t n; is >> n;
  for ( size_t k = 0; k < n; ++k ) {
    vector<int> hel; is >> hel;
    getCVector(is,amps[hel]);
  }
}

void MatchboxXCombData::persistentOutput(PersistentOStream & os) const {
  os << theFactory << theMatchboxME << theSubtractionDipole << theCrossingMap 
     << theAmplitudeToColourMap << theColourToAmplitudeMap 
     << theCrossingSign << theAmplitudePartonData << ounit(theAmplitudeMomenta,GeV) 
     << ounit(theLastRenormalizationScale,GeV2)
     << theCalculateTreeAmplitudes /* << theLastAmplitudes << theLastLargeNAmplitudes */
     << theCalculateOneLoopAmplitudes /* << theLastOneLoopAmplitudes */
     << theCalculateTreeME2 << theLastTreeME2 
     << theCalculateLargeNME2 << theLastLargeNME2 
     << theCalculateOneLoopInterference 
     << theLastOneLoopInterference << theCalculateOneLoopPoles
     << theLastOneLoopPoles << theCalculateColourCorrelators 
     << theColourCorrelators << theCalculateLargeNColourCorrelators 
     << theLargeNColourCorrelators << theCalculateColourSpinCorrelators 
     << theColourSpinCorrelators << theCalculateSpinCorrelators 
     << theSpinCorrelators << theNLight 
     << theNLightJetVec << theNHeavyJetVec << theNLightProtonVec 
     << theColourBasisDim 
     << theNDimPhasespace << theNDimAmplitude << theNDimInsertions 
     << theAmplitudeRandomNumbers << theInsertionRandomNumbers 
     << theDiagramWeights << theSingularLimits// << theLastSingularLimit 
     << theStandardModel << theSymmetryFactor
     << theOLPId << theExternalId;
  putAmplitudeMap(os,theLastAmplitudes);
  putAmplitudeMap(os,theLastLargeNAmplitudes);
  putAmplitudeMap(os,theLastOneLoopAmplitudes);
}

void MatchboxXCombData::persistentInput(PersistentIStream & is, int) {
  is >> theFactory >> theMatchboxME >> theSubtractionDipole >> theCrossingMap 
     >> theAmplitudeToColourMap >> theColourToAmplitudeMap 
     >> theCrossingSign >> theAmplitudePartonData >> iunit(theAmplitudeMomenta,GeV)
     >> iunit(theLastRenormalizationScale,GeV2)
     >> theCalculateTreeAmplitudes /* >> theLastAmplitudes >> theLastLargeNAmplitudes */
     >> theCalculateOneLoopAmplitudes /* >> theLastOneLoopAmplitudes */
     >> theCalculateTreeME2 >> theLastTreeME2 
     >> theCalculateLargeNME2 >> theLastLargeNME2 
     >> theCalculateOneLoopInterference 
     >> theLastOneLoopInterference >> theCalculateOneLoopPoles
     >> theLastOneLoopPoles >> theCalculateColourCorrelators 
     >> theColourCorrelators >> theCalculateLargeNColourCorrelators 
     >> theLargeNColourCorrelators >> theCalculateColourSpinCorrelators 
     >> theColourSpinCorrelators >> theCalculateSpinCorrelators 
     >> theSpinCorrelators >> theNLight 
     >> theNLightJetVec >> theNHeavyJetVec >> theNLightProtonVec 
     >> theColourBasisDim 
     >> theNDimPhasespace >> theNDimAmplitude >> theNDimInsertions 
     >> theAmplitudeRandomNumbers >> theInsertionRandomNumbers 
     >> theDiagramWeights >> theSingularLimits// >> theLastSingularLimit 
     >> theStandardModel >> theSymmetryFactor
     >> theOLPId >> theExternalId;
  getAmplitudeMap(is,theLastAmplitudes);
  getAmplitudeMap(is,theLastLargeNAmplitudes);
  getAmplitudeMap(is,theLastOneLoopAmplitudes);
}

