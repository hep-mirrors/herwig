// -*- C++ -*-
//
// MatchboxXCombData.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxXCombData class.
//

#include "MatchboxXCombData.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxXCombData::MatchboxXCombData() 
  : theCrossingSign(1.0), theCalculateTreeAmplitudes(true), 
    theCalculateOneLoopAmplitudes(true), theCalculateTreeME2(true), 
    theLastTreeME2(0.0), theCalculateOneLoopInterference(true), 
    theLastOneLoopInterference(0.0), theCalculateOneLoopPoles(true), 
    theLastOneLoopPoles(0.0,0.0), 
    theNLight(0), 
    theColourBasisDim(0), theNDimPhasespace(0), 
    theNDimAmplitude(0), theNDimInsertions(0), 
    theSymmetryFactor(0.0), theOLPMomenta(0),
    filledOLPMomenta(false), theExternalId(0),
    theInitialized(false), filledExternalMomenta(false) {
  flushCaches();
}

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
    theLastTreeME2(0.0), theCalculateOneLoopInterference(true), 
    theLastOneLoopInterference(0.0), theCalculateOneLoopPoles(true), 
    theLastOneLoopPoles(0.0,0.0), theNLight(0), 
    theColourBasisDim(0), theNDimPhasespace(0), 
    theNDimAmplitude(0), theNDimInsertions(0), 
    theSymmetryFactor(0.0), theOLPMomenta(0),
    filledOLPMomenta(false), theExternalId(0),
    theInitialized(false), filledExternalMomenta(false) {
  flushCaches();
  theMatchboxME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(newME);
  theSubtractionDipole = dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(newME);
  if ( theMatchboxME )
    theFactory = theMatchboxME->factory();
  else if ( theSubtractionDipole )
    theFactory = theSubtractionDipole->realEmissionME()->factory();
}

void MatchboxXCombData::fillOLPMomenta(const vector<Lorentz5Momentum>& memomenta) {
  if ( filledOLPMomenta )
    return;
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
  theCalculateOneLoopInterference = true;
  theCalculateOneLoopPoles = true;
  for ( map<pair<int,int>,bool>::iterator f = theCalculateColourCorrelators.begin();
	f != theCalculateColourCorrelators.end(); ++f )
    f->second = true;
  for ( map<pair<int,int>,bool>::iterator f = theCalculateLargeNColourCorrelators.begin();
	f != theCalculateLargeNColourCorrelators.end(); ++f )
    f->second = true;
  for ( map<pair<int,int>,bool>::iterator f = theCalculateColourSpinCorrelators.begin();
	f != theCalculateColourSpinCorrelators.end(); ++f )
    f->second = true;
  filledOLPMomenta = false;
  filledExternalMomenta = false;
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
     << theCalculateTreeAmplitudes /* << theLastAmplitudes << theLastLargeNAmplitudes */
     << theCalculateOneLoopAmplitudes /* << theLastOneLoopAmplitudes */
     << theCalculateTreeME2 << theLastTreeME2 << theCalculateOneLoopInterference 
     << theLastOneLoopInterference << theCalculateOneLoopPoles
     << theLastOneLoopPoles << theCalculateColourCorrelators 
     << theColourCorrelators << theCalculateLargeNColourCorrelators 
     << theLargeNColourCorrelators << theCalculateColourSpinCorrelators 
     << theColourSpinCorrelators << theNLight << theColourBasisDim 
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
     >> theCalculateTreeAmplitudes /* >> theLastAmplitudes >> theLastLargeNAmplitudes */
     >> theCalculateOneLoopAmplitudes /* >> theLastOneLoopAmplitudes */
     >> theCalculateTreeME2 >> theLastTreeME2 >> theCalculateOneLoopInterference 
     >> theLastOneLoopInterference >> theCalculateOneLoopPoles
     >> theLastOneLoopPoles >> theCalculateColourCorrelators 
     >> theColourCorrelators >> theCalculateLargeNColourCorrelators 
     >> theLargeNColourCorrelators >> theCalculateColourSpinCorrelators 
     >> theColourSpinCorrelators >> theNLight >> theColourBasisDim 
     >> theNDimPhasespace >> theNDimAmplitude >> theNDimInsertions 
     >> theAmplitudeRandomNumbers >> theInsertionRandomNumbers 
     >> theDiagramWeights >> theSingularLimits// >> theLastSingularLimit 
     >> theStandardModel >> theSymmetryFactor
     >> theOLPId >> theExternalId;
  getAmplitudeMap(is,theLastAmplitudes);
  getAmplitudeMap(is,theLastLargeNAmplitudes);
  getAmplitudeMap(is,theLastOneLoopAmplitudes);
}

