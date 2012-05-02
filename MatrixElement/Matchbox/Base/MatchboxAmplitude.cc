// -*- C++ -*-
//
// MatchboxAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitude class.
//

#include "MatchboxAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "MatchboxMEBase.h"

using namespace Herwig;

MatchboxAmplitude::MatchboxAmplitude() 
  : Amplitude(), theColourBasisDim(0),
    calculateTrees(true), calculateLoops(true) {}

MatchboxAmplitude::~MatchboxAmplitude() {}

void MatchboxAmplitude::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theColourBasis << theColourBasisDim
     << theCrossingMap << theColourMap << theCrossingSigns 
     << theAmplitudePartonData << theNLight;
}

void MatchboxAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theColourBasis >> theColourBasisDim
     >> theCrossingMap >> theColourMap >> theCrossingSigns
     >> theAmplitudePartonData >> theNLight;
}

void MatchboxAmplitude::cloneDependencies(const std::string&) {
}

Ptr<MatchboxMEBase>::ptr MatchboxAmplitude::makeME(const vector<PDVector>&) const {
  return new_ptr(MatchboxMEBase());
}

void MatchboxAmplitude::dumpInfo(const string& prefix) const {
  generator()->log() << prefix << fullName()
		     << " [" << this << "]\n";
  generator()->log() << prefix << "  | XComb " << lastXCombPtr()
		     << " for ";
  if ( lastXCombPtr() ) {
    for ( cPDVector::const_iterator p = lastXComb().mePartonData().begin();
	  p != lastXComb().mePartonData().end(); ++p ) {
      generator()->log() << (**p).PDGName() << " ";
    }
  }
  generator()->log() << "\n";
}

struct orderPartonData {

  bool operator()(const pair<tcPDPtr,int>& a, 
		  const pair<tcPDPtr,int>& b) const {

    if ( a.first == b.first )
      return a.second < b.second;

    int acolour = a.first->iColour();
    int bcolour = b.first->iColour();

    if ( abs(acolour) != abs(bcolour) )
      return abs(acolour) < abs(bcolour);

    if ( a.first->iSpin() != b.first->iSpin() )
      return a.first->iSpin() < b.first->iSpin();

    int acharge = a.first->iCharge();
    int bcharge = b.first->iCharge();

    if ( abs(acharge) != abs(bcharge) )
      return abs(acharge) < abs(bcharge);

    if ( abs(a.first->id()) != abs(b.first->id()) )
      return abs(a.first->id()) < abs(b.first->id());

    return a.first->id() > b.first->id();

  }

};

void MatchboxAmplitude::fillCrossingMap(size_t shift) {

  theLastCrossingMap = crossingMap().find(lastXCombPtr());

  if ( theLastCrossingMap != crossingMap().end() ) {
    assert(amplitudePartonData().find(lastXCombPtr()) != amplitudePartonData().end());
    assert(colourMap().find(lastXCombPtr()) != colourMap().end());
    theLastAmplitudePartonData = amplitudePartonData().find(lastXCombPtr());
    theLastColourMap = colourMap().find(lastXCombPtr());
    assert(crossingSigns().find(lastXCombPtr()) != crossingSigns().end());
    theLastCrossingSign = crossingSigns().find(lastXCombPtr())->second;
    return;
  } else {
    crossingMap()[lastXCombPtr()] = vector<int>(mePartonData().size());
    amplitudePartonData()[lastXCombPtr()] = cPDVector(mePartonData().size());
    colourMap()[lastXCombPtr()] = map<size_t,size_t>();
    theLastCrossingMap = crossingMap().find(lastXCombPtr());
    theLastAmplitudePartonData = amplitudePartonData().find(lastXCombPtr());
    theLastColourMap = colourMap().find(lastXCombPtr());
  }

  double csign = 1.;
  set<pair<tcPDPtr,int>,orderPartonData > processLegs;
  for ( unsigned int l = 0; l < mePartonData().size(); ++l ) {
    if ( l > 1 )
      processLegs.insert(make_pair(mePartonData()[l],l));
    else {
      if ( mePartonData()[l]->CC() ) {
	processLegs.insert(make_pair(mePartonData()[l]->CC(),l));
	if ( mePartonData()[l]->iSpin() == PDT::Spin1Half )
	  csign *= -1.;
      } else {
	processLegs.insert(make_pair(mePartonData()[l],l));
      }
    }
  }

  lastCrossingSign(csign);
  crossingSigns()[lastXCombPtr()] = csign;

  set<pair<tcPDPtr,int> > amplitudeLegs;

  int ampCount = 0;

  // process legs are already sorted, we only need to arrange for
  // adjacent particles and anti-particles
  while ( !processLegs.empty() ) {
    set<pair<tcPDPtr,int>,orderPartonData >::iterator next
      = processLegs.begin();
    while ( next->first->id() < 0 ) {
      if ( ++next == processLegs.end() )
	break;
    }
    assert(next != processLegs.end());
    lastCrossingMap()[ampCount] = next->second - shift;
    amplitudeLegs.insert(make_pair(next->first,ampCount));
    tcPDPtr check = next->first;
    processLegs.erase(next);
    ++ampCount;
    if ( check->CC() ) {
      set<pair<tcPDPtr,int>,orderPartonData>::iterator checkcc
	= processLegs.end();
      for ( set<pair<tcPDPtr,int>,orderPartonData>::iterator c = processLegs.begin();
	    c != processLegs.end(); ++c ) {
	if ( c->first == check->CC() ) {
	  checkcc = c; break;
	}
      }
      if ( checkcc == processLegs.end() )
	for ( set<pair<tcPDPtr,int>,orderPartonData>::iterator c = processLegs.begin();
	      c != processLegs.end(); ++c ) {
	  assert(SU2Helper::SU2CC(check)->CC());
	  if ( c->first == SU2Helper::SU2CC(check)->CC() ) {
	    checkcc = c; break;
	  }
	}
      assert(checkcc != processLegs.end());
      lastCrossingMap()[ampCount] = checkcc->second - shift;
      amplitudeLegs.insert(make_pair(checkcc->first,ampCount));
      processLegs.erase(checkcc);
      ++ampCount;
    }
  }

  for ( set<pair<tcPDPtr,int> >::const_iterator l = amplitudeLegs.begin();
	l != amplitudeLegs.end(); ++l )
    lastAmplitudePartonData()[l->second] = l->first;

  if ( colourBasis() ) {
    assert(colourBasis()->indexMap().find(mePartonData()) !=
	   colourBasis()->indexMap().end());
    const map<size_t,size_t> colourCross = 
      colourBasis()->indexMap().find(mePartonData())->second;
    for ( size_t k = 0; k < lastCrossingMap().size(); ++k ) {
      if ( colourCross.find(lastCrossingMap()[k]) !=
	   colourCross.end() ) {
	lastColourMap()[k] = colourCross.find(lastCrossingMap()[k])->second;
      }
    }
  }

}

Lorentz5Momentum MatchboxAmplitude::amplitudeMomentum(int i) const {
  int iCrossed = lastCrossingMap()[i];
  Lorentz5Momentum res = meMomenta()[iCrossed];
  return iCrossed > 1 ? res : -res;
}

set<vector<int> > MatchboxAmplitude::generateHelicities() const {
  set<vector<int> > res;
  vector<int> current(lastAmplitudePartonData().size());
  doGenerateHelicities(res,current,0);
  return res;
}

void MatchboxAmplitude::doGenerateHelicities(set<vector<int> >& res,
					     vector<int>& current,
					     size_t pos) const {

  if ( pos == lastAmplitudePartonData().size() ) {
    res.insert(current);
    return;
  }

  if ( lastAmplitudePartonData()[pos]->iSpin() == PDT::Spin0 ||
       ( lastAmplitudePartonData()[pos]->iSpin() == PDT::Spin1 &&
	 lastAmplitudePartonData()[pos]->mass() != ZERO ) ) {
    current[pos] = 0;
    doGenerateHelicities(res,current,pos+1);
  } else if ( lastAmplitudePartonData()[pos]->iSpin() == PDT::Spin1Half ||
	      lastAmplitudePartonData()[pos]->iSpin() == PDT::Spin1 ) {
    current[pos] = 1;
    doGenerateHelicities(res,current,pos+1);
    current[pos] = -1;
    doGenerateHelicities(res,current,pos+1);
  }

}

void MatchboxAmplitude::prepareAmplitudes() {

  if ( !calculateTrees )
    return;

  if ( lastAmplitudes().empty() ) {
    set<vector<int> > helicities = generateHelicities();
    for ( set<vector<int> >::const_iterator h = helicities.begin();
	  h != helicities.end(); ++h )
      lastAmplitudes().insert(make_pair(*h,CVector(colourBasisDim())));
    lastLargeNAmplitudes() = lastAmplitudes();
  }

  AmplitudeIterator amp = lastAmplitudes().begin();
  AmplitudeIterator lamp = lastLargeNAmplitudes().begin();
  for ( ;amp != lastAmplitudes().end(); ++amp, ++lamp ) {
    for ( size_t k = 0; k < colourBasisDim(); ++k )
      amp->second(k) = evaluate(k,amp->first,lamp->second(k));
  }

  calculateTrees = false;

}

void MatchboxAmplitude::prepareOneLoopAmplitudes() {

  if ( !calculateLoops )
    return;

  if ( lastOneLoopAmplitudes().empty() ) {
    set<vector<int> > helicities = generateHelicities();
    for ( set<vector<int> >::const_iterator h = helicities.begin();
	  h != helicities.end(); ++h )
      lastOneLoopAmplitudes().insert(make_pair(*h,CVector(colourBasisDim())));
  }

  for ( AmplitudeIterator amp = lastOneLoopAmplitudes().begin();
	amp != lastOneLoopAmplitudes().end(); ++amp ) {
    for ( size_t k = 0; k < colourBasisDim(); ++k )
      amp->second(k) = evaluateOneLoop(k,amp->first);
  }

  calculateLoops = false;

}

Complex MatchboxAmplitude::value(const tcPDVector&,
				 const vector<Lorentz5Momentum>&, 
				 const vector<int>&) {
  assert(false && "ThePEG::Amplitude interface is not sufficient at the moment.");
  throw Exception() << "ThePEG::Amplitude interface is not sufficient at the moment."
		    << Exception::abortnow;
  return 0.;
}

double MatchboxAmplitude::colourCorrelatedME2(pair<int,int> ij) const {
  double Nc = generator()->standardModel()->Nc();
  double cfac = mePartonData()[ij.first]->id() == ParticleID::g ? Nc : (sqr(Nc)-1.)/(2.*Nc);
  return colourBasis()->colourCorrelatedME2(ij,mePartonData(),lastAmplitudes())/cfac;
}

// compare int vectors modulo certain element
// which needs to differe between the two
bool equalsModulo(unsigned int i, const vector<int>& a, const vector<int>& b) {
  assert(a.size()==b.size());
  if ( a[i] == b[i] )
    return false;
  for ( unsigned int k = 0; k < a.size(); ++k ) {
    if ( k == i )
      continue;
    if ( a[k] != b[k] )
      return false;
  }
  return true;
}

double MatchboxAmplitude::spinColourCorrelatedME2(pair<int,int> ij,
						  const SpinCorrelationTensor& c) const {

  using namespace SpinorHelicity;

  Lorentz5Momentum p = meMomenta()[ij.first];
  Lorentz5Momentum n = meMomenta()[ij.second];

  LorentzVector<complex<Energy> > num =
    PlusSpinorCurrent(PlusConjugateSpinor(n),MinusSpinor(p)).eval();

  complex<Energy> den =
    sqrt(2.)*PlusSpinorProduct(PlusConjugateSpinor(n),PlusSpinor(p)).eval();

  LorentzVector<Complex> polarization(num.x()/den,num.y()/den,num.z()/den,num.t()/den);

  Complex pFactor = (polarization*c.momentum())/sqrt(abs(c.scale()));

  double avg =
    colourCorrelatedME2(ij)*(-c.diagonal()+ (c.scale() > ZERO ? 1. : -1.)*norm(pFactor));

  Complex corr = 0.;

  int iCrossed = -1;
  for ( unsigned int k = 0; k < lastCrossingMap().size(); ++k )
    if ( lastCrossingMap()[k] == ij.first ) {
      iCrossed = k;
      break;
    }      
  assert(iCrossed >= 0);

  set<const CVector*> done;
  for ( AmplitudeConstIterator a = theLastAmplitudes.begin();
	a != theLastAmplitudes.end(); ++a ) {
    if ( done.find(&(a->second)) != done.end() )
      continue;
    AmplitudeConstIterator b = theLastAmplitudes.begin();
    while ( !equalsModulo(iCrossed,a->first,b->first) )
      if ( ++b == theLastAmplitudes.end() )
	break;
    if ( b == theLastAmplitudes.end() || done.find(&(b->second)) != done.end() )
      continue;
    done.insert(&(a->second)); done.insert(&(b->second));
    if ( a->first[iCrossed] == 1 )
      swap(a,b);
    corr += colourBasis()->interference(mePartonData(),a->second,b->second);
  }

  double Nc = generator()->standardModel()->Nc();
  double cfac = mePartonData()[ij.first]->id() == ParticleID::g ? Nc : (sqr(Nc)-1.)/(2.*Nc);

  return avg + 2.*(c.scale() > ZERO ? 1. : -1.)*real(corr*sqr(pFactor))/cfac;

}

void MatchboxAmplitude::Init() {

  static ClassDocumentation<MatchboxAmplitude> documentation
    ("MatchboxAmplitude is the base class for amplitude "
     "implementations inside Matchbox.");

  static Reference<MatchboxAmplitude,ColourBasis> interfaceColourBasis
    ("ColourBasis",
     "Set the colour basis implementation.",
     &MatchboxAmplitude::theColourBasis, false, false, true, true, false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxAmplitude,Amplitude>
describeMatchboxAmplitude("Herwig::MatchboxAmplitude", "HwMatchbox.so");
