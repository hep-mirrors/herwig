// -*- C++ -*-
//
// SubtractionDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SubtractionDipole class.
//

#include "SubtractionDipole.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

SubtractionDipole::SubtractionDipole() 
  : MEBase(), theSplitting(false), theApply(true), theSubtractionTest(false),
    theIgnoreCuts(false),
    theRealEmitter(-1), theRealEmission(-1), theRealSpectator(-1), 
    lastRealEmissionKey(realEmissionKey(cPDVector(),-1,-1,-1)),
    lastUnderlyingBornKey(underlyingBornKey(cPDVector(),-1,-1)),
    theBornEmitter(-1), theBornSpectator(-1),
    theLastSubtractionScale(ZERO), theLastSplittingScale(ZERO),
    theLastSubtractionPt(ZERO), theLastSplittingPt(ZERO),
    theLastSubtractionZ(0.0), theLastSplittingZ(0.0),
    theRealShowerSubtraction(false), theVirtualShowerSubtraction(false),
    theLoopSimSubtraction(false), theRealEmissionScales(false),
    theShowerHardScale(ZERO), theShowerScale(ZERO), 
    theIsInShowerPhasespace(false), theIsAboveCutoff(false) {}

SubtractionDipole::~SubtractionDipole() {}

double SubtractionDipole::alpha() const{
  return factory()->alphaParameter();
}

void SubtractionDipole::clearBookkeeping() {
  theRealEmitter = -1;
  theRealEmission = -1;
  theRealSpectator = -1;
  theBornEmitter = -1;
  theBornSpectator = -1;
  theMergingMap.clear();
  theSplittingMap.clear();
  theIndexMap.clear();
  theUnderlyingBornDiagrams.clear();
  theRealEmissionDiagrams.clear();
  theBornToRealDiagrams.clear();
  theRealToBornDiagrams.clear();
}

void SubtractionDipole::setupBookkeeping(const map<Ptr<DiagramBase>::ptr,SubtractionDipole::MergeInfo>& mergeInfo,bool slim) {

  theMergingMap.clear();
  theSplittingMap.clear();

  theUnderlyingBornDiagrams.clear();
  theRealEmissionDiagrams.clear();

  theBornToRealDiagrams.clear();
  theRealToBornDiagrams.clear();

  int xemitter = -1;
  int xspectator = -1;
  map<int,int> mergeLegs;
  map<int,int> remapLegs;
  map<int,int> realBornMap;
  map<int,int> bornRealMap;

  set<Ptr<DiagramBase>::cptr> usedDiagrams;

  for ( map<Ptr<DiagramBase>::ptr,MergeInfo>::const_iterator mit = mergeInfo.begin();
	mit != mergeInfo.end(); ++mit ) {

    DiagramVector::const_iterator bd = 
      theUnderlyingBornME->diagrams().end();

    // work out the most similar underlying Born diagram
    map<int,int> xRemapLegs;
    int nomapScore = 0;
    for ( DiagramVector::const_iterator b = 
	    theUnderlyingBornME->diagrams().begin();
	  b != theUnderlyingBornME->diagrams().end(); ++b ) {
      map<int,int> theRemapLegs;
      if ( mit->second.diagram->isSame(*b,theRemapLegs) &&
	   usedDiagrams.find(*b) == usedDiagrams.end() ) {
	int theNomapScore = 0;
	for ( map<int,int>::const_iterator m = theRemapLegs.begin();
	      m != theRemapLegs.end(); ++m )
	  if ( m->first == m->second )
	    theNomapScore += 1;
	if ( theNomapScore >= nomapScore ) {
	  nomapScore = theNomapScore;
	  xRemapLegs = theRemapLegs;
	  bd = b;
	}
      }
    }

    // no underlying Born
    if ( bd == theUnderlyingBornME->diagrams().end() )
      continue;

    // as we deal with one splitting only we now mark this diagram as used
    // since we fixed the overall remapping of the process from the first
    // occurence, see below. TODO: This confuses this code even more, and
    // clearly calls for a cleanup. This is just grown historically and got
    // messed up with experiencing different processes and setups.
    usedDiagrams.insert(*bd);

    if ( xemitter == -1 ) {

      xemitter = mit->second.emitter;
      mergeLegs = mit->second.mergeLegs;
      remapLegs = xRemapLegs;

      assert(remapLegs.find(xemitter) != remapLegs.end());
      xemitter = remapLegs[xemitter];

      // work out the leg remapping real -> born
      for ( map<int,int>::const_iterator k = mergeLegs.begin();
	    k != mergeLegs.end(); ++k ) {
	assert(remapLegs.find(k->second) != remapLegs.end());
	realBornMap[k->first] = remapLegs[k->second];
      }

      // work out the leg remapping born -> real
      for ( map<int,int>::const_iterator k = realBornMap.begin();
	    k != realBornMap.end(); ++k ) {
	bornRealMap[k->second] = k->first;
      }

      // work out the spectator
      assert(mergeLegs.find(realSpectator()) != mergeLegs.end());
      assert(remapLegs.find(mergeLegs[realSpectator()]) != remapLegs.end());
      xspectator = realBornMap[realSpectator()];

    }

    RealEmissionKey realKey = realEmissionKey((*mit->first).partons(),realEmitter(),realEmission(),realSpectator());
    UnderlyingBornKey bornKey = underlyingBornKey((**bd).partons(),xemitter,xspectator);
    if ( theMergingMap.find(realKey) == theMergingMap.end() )
      theMergingMap.insert(make_pair(realKey,make_pair(bornKey,realBornMap)));
    RealEmissionInfo realInfo = make_pair(realKey,bornRealMap);
    bool gotit = false;
    typedef multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator spIterator;
    pair<spIterator,spIterator> range = theSplittingMap.equal_range(bornKey);
    for ( ; range.first != range.second; ++range.first )
      if ( range.first->second == realInfo ) {
	gotit = true;
	break;
      }
    if ( !gotit )
      theSplittingMap.insert(make_pair(bornKey,realInfo));
    theUnderlyingBornDiagrams[process(realKey)].push_back(*bd);
    theRealEmissionDiagrams[process(bornKey)].push_back(mit->first);
    theBornToRealDiagrams[*bd] = mit->first;
    theRealToBornDiagrams[mit->first] = *bd;

  }
  
  
  if (slim) {
    theIndexMap.clear();
    theSplittingMap.clear();
    theBornToRealDiagrams.clear();
    theRealEmissionDiagrams.clear();
  }

  if ( theSplittingMap.empty() )
    return;

  theIndexMap.clear();

  for ( multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator s =
	  theSplittingMap.begin(); s != theSplittingMap.end(); ++s ) {
    theIndexMap[process(s->first)] = make_pair(emitter(s->first),spectator(s->first));
  }

}

void SubtractionDipole::subtractionBookkeeping() {
  /*
  if ( theMergingMap.empty() )
    setupBookkeeping();
  */
  assert(!theMergingMap.empty());
  lastRealEmissionKey = 
    realEmissionKey(lastHeadXComb().mePartonData(),realEmitter(),realEmission(),realSpectator());
  map<RealEmissionKey,UnderlyingBornInfo>::const_iterator k =
    theMergingMap.find(lastRealEmissionKey);
  if ( k == theMergingMap.end() ) {
    theApply = false;
    return;
  }
  theApply = true;
  lastUnderlyingBornKey = k->second.first;
  bornEmitter(emitter(lastUnderlyingBornKey));
  bornSpectator(spectator(lastUnderlyingBornKey));
}

void SubtractionDipole::splittingBookkeeping() {
  /*
  if ( theMergingMap.empty() )
    setupBookkeeping();
  */
  assert(!theMergingMap.empty());
  map<cPDVector,pair<int,int> >::const_iterator esit =
    theIndexMap.find(lastHeadXComb().mePartonData());
  if ( esit == theIndexMap.end() ) {
    theApply = false;
    return;
  }
  theApply = true;
  pair<int,int> es = esit->second;
  bornEmitter(es.first);
  bornSpectator(es.second);
  lastUnderlyingBornKey = underlyingBornKey(lastHeadXComb().mePartonData(),bornEmitter(),bornSpectator());
  typedef multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator spit;
  pair<spit,spit> kr = theSplittingMap.equal_range(lastUnderlyingBornKey);
  assert(kr.first != kr.second);
  lastRealEmissionInfo = kr.first;
  for ( ; lastRealEmissionInfo != kr.second; ++lastRealEmissionInfo )
    if ( process(lastRealEmissionInfo->second.first) == lastXComb().mePartonData() )
      break;
  assert(lastRealEmissionInfo != kr.second);
  lastRealEmissionKey = lastRealEmissionInfo->second.first;
  realEmitter(emitter(lastRealEmissionKey));
  realEmission(emission(lastRealEmissionKey));
  realSpectator(spectator(lastRealEmissionKey));
}

StdXCombPtr SubtractionDipole::makeXComb(Energy newMaxEnergy, const cPDPair & inc,
					 tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
					 tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
					 const PBPair & newPartonBins, tCutsPtr newCuts,
					 const DiagramVector & newDiagrams, bool mir,
					 const PartonPairVec& allBins,
					 tStdXCombPtr newHead,
					 tMEPtr newME) {

  if ( !newME )
    newME = this;

  if ( !splitting() ) {
    return
      underlyingBornME()->makeXComb(newMaxEnergy, inc,
				    newEventHandler, newSubProcessHandler,
				    newExtractor, newCKKW,
				    newPartonBins, newCuts,
				    newDiagrams, mir, allBins,
				    newHead, newME);
  }

  return
    realEmissionME()->makeXComb(newMaxEnergy, inc,
				newEventHandler, newSubProcessHandler,
				newExtractor, newCKKW,
				newPartonBins, newCuts,
				newDiagrams, mir, allBins,
				newHead, newME);

}

StdXCombPtr SubtractionDipole::makeXComb(tStdXCombPtr newHead,
					 const PBPair & newPartonBins,
					 const DiagramVector & newDiagrams,
					 tMEPtr newME) {

  if ( !newME )
    newME = this;

  if ( !splitting() ) {
    return
      underlyingBornME()->makeXComb(newHead, newPartonBins,
				    newDiagrams, newME);
  }

  return
    realEmissionME()->makeXComb(newHead, newPartonBins,
				newDiagrams, newME);

}

StdXCombPtr SubtractionDipole::makeBornXComb(tStdXCombPtr realXC) {

  const cPDVector& proc = const_cast<const StandardXComb&>(*realXC).mePartonData();

  lastRealEmissionKey = 
    realEmissionKey(proc,realEmitter(),realEmission(),realSpectator());
  map<RealEmissionKey,UnderlyingBornInfo>::const_iterator k =
    theMergingMap.find(lastRealEmissionKey);

  if ( k == theMergingMap.end() )
    return StdXCombPtr();

  PartonPairVec pbs = realXC->pExtractor()->getPartons(realXC->maxEnergy(), 
						       realXC->particles(),
						       *(realXC->cuts()));

  DiagramVector bornDiags = underlyingBornDiagrams(proc);
  assert(!bornDiags.empty());

  PartonPairVec::iterator ppit = pbs.begin();
  for ( ; ppit != pbs.end(); ++ppit ) {
    if ( ppit->first->parton() == bornDiags.front()->partons()[0] &&
	 ppit->second->parton() == bornDiags.front()->partons()[1] )
      break;
  }

  assert(ppit != pbs.end());

  return
    underlyingBornME()->makeXComb(realXC,*ppit,bornDiags,this);

}

vector<StdXCombPtr> SubtractionDipole::makeRealXCombs(tStdXCombPtr bornXC) {

  const cPDVector& proc = const_cast<const StandardXComb&>(*bornXC).mePartonData();

  map<cPDVector,pair<int,int> >::const_iterator esit = theIndexMap.find(proc);
  if ( esit == theIndexMap.end() ) 
    return vector<StdXCombPtr>();
  pair<int,int> es = esit->second;
  bornEmitter(es.first);
  bornSpectator(es.second);
  lastUnderlyingBornKey = underlyingBornKey(proc,bornEmitter(),bornSpectator());

  if ( theSplittingMap.find(lastUnderlyingBornKey) == theSplittingMap.end() )
    return vector<StdXCombPtr>();

  PartonPairVec pbs = bornXC->pExtractor()->getPartons(bornXC->maxEnergy(), 
						       bornXC->particles(),
						       *(bornXC->cuts()));

  DiagramVector realDiags = realEmissionDiagrams(proc);
  assert(!realDiags.empty());

  vector<StdXCombPtr> res;

  map<cPDVector,DiagramVector> realProcs;

  for ( MEBase::DiagramVector::const_iterator d = realDiags.begin();
	d != realDiags.end(); ++d ) {
    realProcs[(**d).partons()].push_back(*d);
  }

  for ( map<cPDVector,DiagramVector>::const_iterator pr =
	  realProcs.begin(); pr != realProcs.end(); ++pr ) {

    PartonPairVec::iterator ppit = pbs.begin();
    for ( ; ppit != pbs.end(); ++ppit ) {
      if ( ppit->first->parton() == pr->second.front()->partons()[0] &&
	   ppit->second->parton() == pr->second.front()->partons()[1] )
	break;
    }

    assert(ppit != pbs.end());

    StdXCombPtr rxc =
      realEmissionME()->makeXComb(bornXC,*ppit,pr->second,this);

    res.push_back(rxc);

  }

  return res;

}

const MEBase::DiagramVector& SubtractionDipole::underlyingBornDiagrams(const cPDVector& real) const {
  static DiagramVector empty;
  map<cPDVector,DiagramVector>::const_iterator k = theUnderlyingBornDiagrams.find(real);
  if (k == theUnderlyingBornDiagrams.end() )
    return empty;
  return k->second;
}

tcDiagPtr SubtractionDipole::underlyingBornDiagram(tcDiagPtr realDiag) const {
  map<tcDiagPtr,tcDiagPtr>::const_iterator it = theRealToBornDiagrams.find(realDiag);
  assert(it != theRealToBornDiagrams.end());
  return it->second;
}

const MEBase::DiagramVector& SubtractionDipole::realEmissionDiagrams(const cPDVector& born) const {
  static DiagramVector empty;
  map<cPDVector,DiagramVector>::const_iterator k = theRealEmissionDiagrams.find(born);
  if ( k == theRealEmissionDiagrams.end() )
    return empty;
  return k->second;
}

tcDiagPtr SubtractionDipole::realEmissionDiagram(tcDiagPtr bornDiag) const {
  map<tcDiagPtr,tcDiagPtr>::const_iterator it = theBornToRealDiagrams.find(bornDiag);
  assert(it != theBornToRealDiagrams.end());
  return it->second;
}

void SubtractionDipole::getDiagrams() const {
  if ( splitting() ) {
    realEmissionME()->diagrams();
    useDiagrams(realEmissionME());
  } else {
    underlyingBornME()->diagrams();
    useDiagrams(underlyingBornME());
  }
}

Selector<MEBase::DiagramIndex> SubtractionDipole::diagrams(const DiagramVector & dv) const {
  Ptr<MatchboxMEBase>::tcptr me = 
    splitting() ?
    realEmissionME() :
    underlyingBornME();
  if ( me->phasespace() ) {
    me->phasespace()->setXComb(lastXCombPtr());
    me->phasespace()->fillDiagramWeights();
  }
  return 
    me->diagrams(dv);
}

Selector<const ColourLines *>
SubtractionDipole::colourGeometries(tcDiagPtr diag) const {
  return 
    splitting() ?
    realEmissionME()->colourGeometries(diag) :
    underlyingBornME()->colourGeometries(diag);
}

const ColourLines &
SubtractionDipole::selectColourGeometry(tcDiagPtr diag) const {
  return 
    splitting() ?
    realEmissionME()->selectColourGeometry(diag) :
    underlyingBornME()->selectColourGeometry(diag);
}

void SubtractionDipole::flushCaches() {
  theUnderlyingBornME->flushCaches();
  theRealEmissionME->flushCaches();
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator r =
	  reweights().begin(); r != reweights().end(); ++r ) {
    (**r).flushCaches();
  }
}

void SubtractionDipole::setXComb(tStdXCombPtr xc) {
  if ( !xc ) {
    theApply = false;
    return;
  } else {
    theApply = true;
  }
  lastMatchboxXComb(xc);
  MEBase::setXComb(xc); 
  if ( splitting() ) {
    realEmissionME()->setXComb(xc);
    underlyingBornME()->setXComb(xc->head());
    splittingBookkeeping();
  } else {
    realEmissionME()->setXComb(xc->head());
    underlyingBornME()->setXComb(xc);
    subtractionBookkeeping();
  }
  if ( !apply() )
    return;
}

void SubtractionDipole::setKinematics() {
  MEBase::setKinematics(); 
  if ( splitting() )
    realEmissionME()->setKinematics();
  else
    underlyingBornME()->setKinematics();
}

bool SubtractionDipole::generateKinematics(const double * r) {
  if ( lastXCombPtr()->kinematicsGenerated() )
    return true;
  if ( splitting() ) {
    if ( !generateRadiationKinematics(r) )
      return false;
    if( ! realEmissionME()->lastXCombPtr()->setIncomingPartons())
      return false;
    realEmissionME()->setScale();
    double jac = jacobian();
    jac *= pow(underlyingBornME()->lastXComb().lastSHat() / realEmissionME()->lastXComb().lastSHat(),
	       realEmissionME()->lastXComb().mePartonData().size()-4.);
    jacobian(jac);
    assert(lastXCombPtr() == realEmissionME()->lastXCombPtr());
    lastXCombPtr()->didGenerateKinematics();
    return true;
  }
  if ( !generateTildeKinematics() ){ return false;}
  if( ! underlyingBornME()->lastXCombPtr()->setIncomingPartons() )
    return false;
  underlyingBornME()->setScale();
  assert(lastXCombPtr() == underlyingBornME()->lastXCombPtr());
  if( ! underlyingBornME()->lastXCombPtr()->setIncomingPartons() )
      return false;
  // need to have the scale and x's available for checking shower phase space
  if ( showerApproximation() &&
       lastXCombPtr()->willPassCuts() )
    showerApproximation()->getShowerVariables();
  lastXCombPtr()->didGenerateKinematics();
  return true;
}

int SubtractionDipole::nDim() const {
  if ( !splitting() )
    return underlyingBornME()->nDim();
  return underlyingBornME()->nDim() + nDimRadiation();
}

void SubtractionDipole::clearKinematics() {
  MEBase::clearKinematics(); 
  if ( splitting() )
    realEmissionME()->clearKinematics();
  else
    underlyingBornME()->clearKinematics();
}

void SubtractionDipole::tildeKinematics(Ptr<TildeKinematics>::tptr tk) { 
  theTildeKinematics = tk;
}

bool SubtractionDipole::generateTildeKinematics() {

  assert(!splitting());

  Ptr<TildeKinematics>::tptr kinematics = theTildeKinematics;
  if ( showerApproximation() ) {
    showerApproximation()->setBornXComb(lastXCombPtr());
    showerApproximation()->setRealXComb(realEmissionME()->lastXCombPtr());
    showerApproximation()->setDipole(this);
    showerApproximation()->checkCutoff();
    if ( showerApproximation()->showerTildeKinematics() &&
	 isAboveCutoff() &&
	 realShowerSubtraction() )
      kinematics = showerApproximation()->showerTildeKinematics();
  }

  if ( !kinematics ) {
    jacobian(0.0);
    return false;
  }

  kinematics->prepare(lastHeadXCombPtr(),lastXCombPtr());

  if ( !kinematics->doMap() ) {
    jacobian(0.0);
    return false;
  }

  theLastSubtractionScale = kinematics->lastScale();
  theLastSubtractionPt = kinematics->lastPt();
  theLastSubtractionZ = kinematics->lastZ();

  meMomenta().resize(lastHeadXComb().meMomenta().size() - 1);

  assert(mergingMap().find(lastRealEmissionKey) != mergingMap().end());
  map<int,int>& trans = theMergingMap[lastRealEmissionKey].second;

  int n = lastHeadXComb().meMomenta().size();
  for ( int k = 0; k < n; ++k ) {
    if ( k == realEmitter() || k == realEmission() || k == realSpectator() )
      continue;
    meMomenta()[trans[k]] = lastHeadXComb().meMomenta()[k];
    if ( kinematics->doesTransform() && k > 1 )
      meMomenta()[trans[k]] = kinematics->transform(meMomenta()[trans[k]]);
  }

  meMomenta()[bornEmitter()] = 
    const_cast<const TildeKinematics&>(*kinematics).bornEmitterMomentum();
  meMomenta()[bornSpectator()] = 
    const_cast<const TildeKinematics&>(*kinematics).bornSpectatorMomentum();

  cPDVector::const_iterator pd = mePartonData().begin();
  vector<Lorentz5Momentum>::iterator p = meMomenta().begin();
  for ( ; pd != mePartonData().end(); ++pd, ++p ) {
    p->setMass((**pd).hardProcessMass());
    p->rescaleRho();
  }

  jacobian(realEmissionME()->lastXComb().jacobian());

  logGenerateTildeKinematics();

  return true;

}

void SubtractionDipole::invertedTildeKinematics(Ptr<InvertedTildeKinematics>::tptr itk) { 
  theInvertedTildeKinematics = itk;
}

int SubtractionDipole::nDimRadiation() const {
  return invertedTildeKinematics() ? 
    invertedTildeKinematics()->nDimRadiation() :
    0;
}

bool SubtractionDipole::generateRadiationKinematics(const double * r) {

  assert(splitting());

  Ptr<InvertedTildeKinematics>::tptr kinematics = theInvertedTildeKinematics;
  if ( showerApproximation() ) {
    showerApproximation()->setBornXComb(lastHeadXCombPtr());
    showerApproximation()->setRealXComb(lastXCombPtr());
    showerApproximation()->setDipole(this);
    if ( showerApproximation()->showerInvertedTildeKinematics() ) {
      kinematics = showerApproximation()->showerInvertedTildeKinematics();
    }
  }

  if ( !kinematics ) {
    jacobian(0.0);
    return false;
  }

  kinematics->prepare(lastXCombPtr(),lastHeadXCombPtr());

  if ( !kinematics->doMap(r) ) {
    jacobian(0.0);
    return false;
  }

  theLastSplittingScale = kinematics->lastScale();
  theLastSplittingPt = kinematics->lastPt();
  theLastSplittingZ = kinematics->lastZ();

  meMomenta().resize(lastHeadXComb().meMomenta().size() + 1);

  assert(splittingMap().find(lastUnderlyingBornKey) != splittingMap().end());
  map<int,int>& trans = const_cast<map<int,int>&>(lastRealEmissionInfo->second.second);

  int n = lastHeadXComb().meMomenta().size();
  for ( int k = 0; k < n; ++k ) {
    if ( k == bornEmitter() || k == bornSpectator() )
      continue;
    meMomenta()[trans[k]] = lastHeadXComb().meMomenta()[k];
    if ( kinematics->doesTransform() && k > 1 )
      meMomenta()[trans[k]] = kinematics->transform(meMomenta()[trans[k]]);
  }

  meMomenta()[realEmitter()] = 
    const_cast<const InvertedTildeKinematics&>(*kinematics).realEmitterMomentum();
  meMomenta()[realEmission()] = 
    const_cast<const InvertedTildeKinematics&>(*kinematics).realEmissionMomentum();
  meMomenta()[realSpectator()] = 
    const_cast<const InvertedTildeKinematics&>(*kinematics).realSpectatorMomentum();

  cPDVector::const_iterator pd = mePartonData().begin();
  vector<Lorentz5Momentum>::iterator p = meMomenta().begin();
  for ( ; pd != mePartonData().end(); ++pd, ++p ) {
    p->setMass((**pd).hardProcessMass());
    p->rescaleRho();
  }

  jacobian(underlyingBornME()->lastXComb().jacobian() *
	   kinematics->jacobian());

  logGenerateRadiationKinematics(r);

  return true;

}

void SubtractionDipole::ptCut(Energy cut) {
  theInvertedTildeKinematics->ptCut(cut);
}

CrossSection SubtractionDipole::dSigHatDR(Energy2 factorizationScale) const {

  double pdfweight = 1.;

  double jac = jacobian();

  if ( splitting() && jac == 0.0 ) {
    lastMECrossSection(ZERO);
    return ZERO;
  }

  if ( factorizationScale == ZERO ) {
    factorizationScale = underlyingBornME()->lastScale();
  }

  if ( havePDFWeight1() ) {
    pdfweight *= realEmissionME()->pdf1(factorizationScale);
  }
 
  if ( havePDFWeight2() ) {
    pdfweight *= realEmissionME()->pdf2(factorizationScale);
  }

  lastMEPDFWeight(pdfweight);

  bool needTheDipole = true;
  CrossSection shower = ZERO;

  double lastThetaMu = 1.0;

  double showerFactor = 1.;

  if ( showerApproximation() ) {
    assert(!splitting());
    showerApproximation()->setBornXComb(lastXCombPtr());
    showerApproximation()->setRealXComb(realEmissionME()->lastXCombPtr());
    showerApproximation()->setDipole(const_cast<SubtractionDipole*>(this));
    if ( !isAboveCutoff() ) {
      showerApproximation()->wasBelowCutoff();
      lastThetaMu = 0.0;
    } else {
      lastThetaMu = 1.0;
    }
    if ( lastThetaMu > 0.0 && isInShowerPhasespace() ) {
      if ( realShowerSubtraction() )
	shower = showerApproximation()->dSigHatDR()*lastThetaMu;
      if ( virtualShowerSubtraction() || loopSimSubtraction() )
	shower = -showerApproximation()->dSigHatDR()*lastThetaMu;
      if ( virtualShowerSubtraction() &&
	   isAboveCutoff() &&
	   showerApproximation()->showerTildeKinematics() ) {
	// map shower to dipole kinematics; we are always above the
	// cutoff in this case
	showerFactor *= 
	  showerApproximation()->showerTildeKinematics()->jacobianRatio();
      }
      shower *= showerFactor;
    }
    if ( realShowerSubtraction() && lastThetaMu == 1.0 )
      needTheDipole = false;
    if ( virtualShowerSubtraction() && lastThetaMu == 0.0 )
      needTheDipole = false;
    if ( factory()->loopSimCorrections() ||
	 factory()->meCorrectionsOnly() )
      needTheDipole = false;
  }

  double xme2 = 0.0;

  if ( needTheDipole )
    xme2 = me2();

  if ( factory()->loopSimCorrections() ||
       factory()->meCorrectionsOnly() ) {

    assert(showerApproximation());
    xme2 = realEmissionME()->me2() * showerApproximation()->channelWeight();

    double rws =
      pow(underlyingBornME()->lastXComb().lastAlphaS()/
	  realEmissionME()->lastXComb().lastAlphaS(),
	  realEmissionME()->orderInAlphaS());

    xme2 *= rws;

    double rwe =
      pow(underlyingBornME()->lastXComb().lastAlphaEM()/
	  realEmissionME()->lastXComb().lastAlphaEM(),
	  underlyingBornME()->orderInAlphaEW());

    xme2 *= rwe;

  }

  if ( realShowerSubtraction() )
    xme2 *= 1. - lastThetaMu;
  if ( virtualShowerSubtraction() || loopSimSubtraction() )
    xme2 *= lastThetaMu;

  double coupl = lastMECouplings();
  coupl *= underlyingBornME()->lastXComb().lastAlphaS();
  lastMECouplings(coupl);

  CrossSection res = 
    sqr(hbarc) * jac * pdfweight * xme2 /
    (2. * realEmissionME()->lastXComb().lastSHat());

  if ( !showerApproximation() && xme2 != 0.0 ) {
    double weight = 0.0;
    bool applied = false;
    for ( vector<Ptr<MatchboxReweightBase>::ptr>::const_iterator rw =
	    theReweights.begin(); rw != theReweights.end(); ++rw ) {
      (**rw).setXComb(theRealEmissionME->lastXCombPtr());
      if ( !(**rw).apply() )
	continue;
      weight += (**rw).evaluate();
      applied = true;
    }
    if ( applied )
      res *= weight;
  }

  lastMECrossSection(-res-shower);

  logDSigHatDR(jac);

  return lastMECrossSection();
}


bool SubtractionDipole::aboveAlpha() const{return theTildeKinematics->aboveAlpha();}



CrossSection SubtractionDipole::prefactor(Energy2 factorizationScale)const{
  
  const double jac = jacobian();
  assert( factorizationScale != ZERO );
  assert (! splitting());
  double pdfweight = 1.;
  if ( havePDFWeight1() ) pdfweight *= realEmissionME()->pdf1(factorizationScale);
  if ( havePDFWeight2() ) pdfweight *= realEmissionME()->pdf2(factorizationScale);

  
  return sqr(hbarc) * jac * pdfweight /  (2. * realEmissionME()->lastXComb().lastSHat());
  
}




CrossSection SubtractionDipole::ps(Energy2 factorizationScale,Ptr<ColourBasis>::tptr largeNBasis) const {

  double ccme2 =underlyingBornME()->me2()*
                underlyingBornME()->
                largeNColourCorrelatedME2(
                  make_pair(bornEmitter(),bornSpectator()),largeNBasis)/
                underlyingBornME()->largeNME2(largeNBasis);
  
  return prefactor(factorizationScale) * me2Avg(ccme2);
}




pair<CrossSection,CrossSection> SubtractionDipole::dipandPs(Energy2 factorizationScale,Ptr<ColourBasis>::tptr largeNBasis) const {
  
  CrossSection  factor= prefactor(factorizationScale);

  double ccme2 =underlyingBornME()->me2()*
                underlyingBornME()->
                largeNColourCorrelatedME2(
                              make_pair(bornEmitter(),bornSpectator()),largeNBasis)/
                underlyingBornME()->largeNME2(largeNBasis);

  double ps = me2Avg(ccme2);
  double dip = me2();
  
  return make_pair(factor*dip,factor*ps);
}

CrossSection SubtractionDipole::dip(Energy2 factorizationScale) const {
  CrossSection  factor= prefactor(factorizationScale);
  double dip = me2();
  return factor*dip;
}


void SubtractionDipole::print(ostream& os) const {

  os << "--- SubtractionDipole setup ----------------------------------------------------\n";

  os << " subtraction '" << name() << "'\n for real emission '"
     << theRealEmissionME->name() << "'\n using underlying Born '"
     << theUnderlyingBornME->name() << "'\n";

  os << " tilde kinematics are '"
     << (theTildeKinematics ? theTildeKinematics->name() : "") 
     << " '\n inverted tilde kinematics are '"
     << (theInvertedTildeKinematics ? theInvertedTildeKinematics->name() : "") << "'\n";

  os << " the following subtraction mappings have been found:\n";

  for ( map<RealEmissionKey,UnderlyingBornInfo>::const_iterator m =
	  theMergingMap.begin(); m != theMergingMap.end(); ++m ) {
    os << " " << process(m->second.first)[0]->PDGName() << " "
       << process(m->second.first)[1]->PDGName() << " -> ";
    for ( cPDVector::const_iterator p = process(m->second.first).begin() + 2;
	  p != process(m->second.first).end(); ++p ) {
      os << (**p).PDGName() << " ";
    }
    os << "[" << emitter(m->second.first) << "," << spectator(m->second.first) << "] <=> ";
    os << process(m->first)[0]->PDGName() << " "
       << process(m->first)[1]->PDGName() << " -> ";
    for ( cPDVector::const_iterator p = process(m->first).begin() + 2;
	  p != process(m->first).end(); ++p ) {
      os << (**p).PDGName() << " ";
    }
    os << "[(" << emitter(m->first) << "," << emission(m->first) << ")," << spectator(m->first) << "]\n"
       << " non-dipole momenta ( ";
    for ( map<int,int>::const_iterator k = m->second.second.begin();
	  k != m->second.second.end(); ++k ) {
      if ( k->first == spectator(m->first) )
	continue;
      os << k->second << " ";
    }
    os << ") <=> ( ";
    for ( map<int,int>::const_iterator k = m->second.second.begin();
	  k != m->second.second.end(); ++k ) {
      if ( k->first == spectator(m->first) )
	continue;
      os << k->first << " ";
    }
    os << ")\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void SubtractionDipole::printLastEvent(ostream& os) const {

  os << "--- SubtractionDipole last event information -----------------------------------\n";

  os << " for dipole '" << name() << "' applying [" 
     << bornEmitter() << "," << bornSpectator() << "] <=> [("
     << realEmitter() << "," << realEmission() << ")," << realSpectator() << "]\n"
     << " evaluated the cross section/nb " << (lastMECrossSection()/nanobarn) << "\n"
     << " with subtraction parameters x[0] = " << subtractionParameters()[0]
     << " x[1] = " << subtractionParameters()[1] << "\n";

  os << " the last real emission event was:\n";
  realEmissionME()->printLastEvent(os);

  os << " the last underlying Born event was:\n";
  underlyingBornME()->printLastEvent(os);


  os << "--- end SubtractionDipole last event information -------------------------------\n";

  os << flush;

}

void SubtractionDipole::logME2() const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  tcStdXCombPtr bornxc = splitting() ? lastHeadXCombPtr() : lastXCombPtr();
  tcStdXCombPtr realxc = splitting() ? lastXCombPtr() : lastHeadXCombPtr();

  generator()->log() << "'" << name() << "' evaluated me2 using\n"
		     << "Born XComb " << bornxc << " real XComb " << realxc << "\n";

  generator()->log() << "subtraction parameters: ";
  copy(subtractionParameters().begin(),subtractionParameters().end(),
       ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n";

  generator()->log() << "Born phase space point (in GeV):\n";

  vector<Lorentz5Momentum>::const_iterator pit = bornxc->meMomenta().begin();
  cPDVector::const_iterator dit = bornxc->mePartonData().begin();

  for ( ; pit != bornxc->meMomenta().end() ; ++pit, ++dit )
    generator()->log() << (**dit).PDGName() << " : "
		       << (*pit/GeV) << "\n";

  generator()->log() << "with x1 = " << bornxc->lastX1() << " x2 = " << bornxc->lastX2() << "\n"
		     << "sHat/GeV2 = " << (bornxc->lastSHat()/GeV2) << "\n";

  generator()->log() << "Real emission phase space point (in GeV):\n";

  pit = realxc->meMomenta().begin();
  dit = realxc->mePartonData().begin();

  for ( ; pit != realxc->meMomenta().end() ; ++pit, ++dit )
    generator()->log() << (**dit).PDGName() << " : "
		       << (*pit/GeV) << "\n";

  generator()->log() << "with x1 = " << realxc->lastX1() << " x2 = " << realxc->lastX2() << "\n"
		     << "sHat/GeV2 = " << (realxc->lastSHat()/GeV2) << "\n";

}

void SubtractionDipole::logDSigHatDR(double effectiveJac) const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  tcStdXCombPtr bornxc = splitting() ? lastHeadXCombPtr() : lastXCombPtr();
  tcStdXCombPtr realxc = splitting() ? lastXCombPtr() : lastHeadXCombPtr();

  generator()->log() << "'" << name() << "' evaluated cross section using\n"
		     << "Born XComb " << bornxc << " real XComb " << realxc << "\n"
		     << "Jacobian = " << jacobian()
		     << " effective Jacobian = " << effectiveJac << "\n"
		     << "Born sHat/GeV2 = " << (bornxc->lastSHat()/GeV2)
		     << " real sHat/GeV2 = " << (realxc->lastSHat()/GeV2)
		     << " dsig/nb = "
		     << (lastMECrossSection()/nanobarn) << "\n" << flush;

}

void SubtractionDipole::logGenerateTildeKinematics() const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  generator()->log() << "'" << name() << "' generating tilde kinematics.\n"
		     << "configuration: [" << bornEmitter() << ","
		     << bornSpectator() << "] => "
		     << "[(" << realEmitter() << "," << realEmission() << "),"
		     << realSpectator() << "]\n"
		     << "with real xcomb " << lastHeadXCombPtr() << " born xcomb "
		     << lastXCombPtr() << "\n"
		     << "from real emission phase space point:\n";
  Lorentz5Momentum rSum;
  vector<Lorentz5Momentum>::const_iterator pr = lastHeadXComb().meMomenta().begin();
  cPDVector::const_iterator dr = lastHeadXComb().mePartonData().begin();
  size_t count = 0;
  for ( ; pr != lastHeadXComb().meMomenta().end(); ++pr,++dr ) {
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";
    if ( count < 2 ) {
      rSum -= *pr;
    } else {
      rSum += *pr;
    }
    ++count;

  }
  generator()->log() << "sum : " << (rSum/GeV) << "\n";

  generator()->log() << "subtraction parameters: ";
  copy(subtractionParameters().begin(),subtractionParameters().end(),
       ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n"
		     << "with scale/GeV = " << (theLastSubtractionScale/GeV)
		     << "and pt/GeV = " << (theLastSubtractionPt/GeV) << "\n";

  generator()->log() << "generated tilde kinematics:\n";
  pr = lastXComb().meMomenta().begin();
  dr = lastXComb().mePartonData().begin();
  count = 0;
  Lorentz5Momentum bSum;
  for ( ; pr != lastXComb().meMomenta().end(); ++pr,++dr ) {
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";
    if ( count < 2 ) {
      bSum -= *pr;
    } else {
      bSum += *pr;
    }
    ++count;
  }
  generator()->log() << "sum : " << (bSum/GeV) << "\n";

  generator()->log() << "Jacobian = " << jacobian() << "\n" << flush;

}


void SubtractionDipole::logGenerateRadiationKinematics(const double * r) const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  generator()->log() << "'" << name() << "' generating radiation kinematics.\n"
		     << "configuration: [" << bornEmitter() << ","
		     << bornSpectator() << "] => "
		     << "[(" << realEmitter() << "," << realEmission() << "),"
		     << realSpectator() << "]\n"
		     << "with born xcomb " << lastHeadXCombPtr() << " real xcomb "
		     << lastXCombPtr() << "\n"
		     << "from random numbers:\n";
  copy(r,r+nDimRadiation(),ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n";
  generator()->log() << "and born phase space point:\n";
  vector<Lorentz5Momentum>::const_iterator pr = lastHeadXComb().meMomenta().begin();
  cPDVector::const_iterator dr = lastHeadXComb().mePartonData().begin();
  for ( ; pr != lastHeadXComb().meMomenta().end(); ++pr,++dr )
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";

  generator()->log() << "subtraction parameters: ";
  copy(subtractionParameters().begin(),subtractionParameters().end(),
       ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n" << flush;

  generator()->log() << "scales: scale/GeV = " << (theLastSplittingScale/GeV)
		     << " pt/GeV = " << (theLastSplittingPt/GeV) << "\n" << flush;

  generator()->log() << "generated real emission kinematics:\n";
  pr = lastXComb().meMomenta().begin();
  dr = lastXComb().mePartonData().begin();
  for ( ; pr != lastXComb().meMomenta().end(); ++pr,++dr )
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";

  generator()->log() << "Jacobian = "
		     << jacobian() << " = "
		     << underlyingBornME()->lastXComb().jacobian()
		     << "|Born * "
		     << invertedTildeKinematics()->jacobian()
		     << "|Radiation\n" << flush;

}

void SubtractionDipole::doinit() {
  MEBase::doinit();
  if ( underlyingBornME() ) {
    theUnderlyingBornME->init();
  }
  if ( realEmissionME() ) {
    theRealEmissionME->init();
  }
  if ( tildeKinematics() ) {
    theTildeKinematics->init();
  }
  if ( invertedTildeKinematics() ) {
    theInvertedTildeKinematics->init();
  }
  if ( showerApproximation() ) {
    theShowerApproximation->init();
  }
  for ( vector<Ptr<SubtractionDipole>::tptr>::iterator p = thePartners.begin();
	p != thePartners.end(); ++p ) {
    (**p).init();
  }
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    (**rw).init();
  }
}

void SubtractionDipole::doinitrun() {
  MEBase::doinitrun();
  if ( underlyingBornME() ) {
    theUnderlyingBornME->initrun();
  }
  if ( realEmissionME() ) {
    theRealEmissionME->initrun();
  }
  if ( tildeKinematics() ) {
    theTildeKinematics->initrun();
  }
  if ( invertedTildeKinematics() ) {
    theInvertedTildeKinematics->initrun();
  }
  if ( showerApproximation() ) {
    theShowerApproximation->initrun();
  }
  for ( vector<Ptr<SubtractionDipole>::tptr>::iterator p = thePartners.begin();
	p != thePartners.end(); ++p ) {
    (**p).initrun();
  }
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    (**rw).initrun();
  }
}

void SubtractionDipole::cloneDependencies(const std::string& prefix,bool slim) {

  if ( underlyingBornME() ) {
    Ptr<MatchboxMEBase>::ptr myUnderlyingBornME = underlyingBornME()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myUnderlyingBornME->name();
    if ( ! (generator()->preinitRegister(myUnderlyingBornME,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Matrix element " << pname.str() << " already existing." << Exception::runerror;
    myUnderlyingBornME->cloneDependencies(pname.str(),slim);
    underlyingBornME(myUnderlyingBornME);
  }

  if ( realEmissionME()&& !slim ) {
    Ptr<MatchboxMEBase>::ptr myRealEmissionME = realEmissionME()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myRealEmissionME->name();
    if ( ! (generator()->preinitRegister(myRealEmissionME,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Matrix element " << pname.str() << " already existing." << Exception::runerror;
    myRealEmissionME->cloneDependencies(pname.str());
    realEmissionME(myRealEmissionME);
  }

  if ( tildeKinematics() ) {
    Ptr<TildeKinematics>::ptr myTildeKinematics = tildeKinematics()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myTildeKinematics->name();
    if ( ! (generator()->preinitRegister(myTildeKinematics,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Tilde kinematics " << pname.str() << " already existing." << Exception::runerror;
    myTildeKinematics->dipole(this);
    tildeKinematics(myTildeKinematics);
  }

  if ( invertedTildeKinematics()&& !slim ) {
    Ptr<InvertedTildeKinematics>::ptr myInvertedTildeKinematics = invertedTildeKinematics()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myInvertedTildeKinematics->name();
    if ( ! (generator()->preinitRegister(myInvertedTildeKinematics,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Inverted tilde kinematics " << pname.str() << " already existing." << Exception::runerror;
    myInvertedTildeKinematics->dipole(this);
    invertedTildeKinematics(myInvertedTildeKinematics);
  }

  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    Ptr<MatchboxReweightBase>::ptr myReweight = (**rw).cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << (**rw).name();
    if ( ! (generator()->preinitRegister(myReweight,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Reweight " << pname.str() << " already existing." << Exception::runerror;
    myReweight->cloneDependencies(pname.str());
    *rw = myReweight;
  }

}

void SubtractionDipole::constructVertex(tSubProPtr sub) {
  if ( splitting() )
    realEmissionME()->constructVertex(sub);
  else 
    underlyingBornME()->constructVertex(sub);
}

void SubtractionDipole::constructVertex(tSubProPtr sub, const ColourLines* cl) {
  if ( splitting() )
    realEmissionME()->constructVertex(sub,cl);
  else 
    underlyingBornME()->constructVertex(sub,cl);
}

void SubtractionDipole::generateSubCollision(SubProcess & sub) {
  if ( splitting() )
    realEmissionME()->generateSubCollision(sub);
  else
    underlyingBornME()->generateSubCollision(sub);
}

void SubtractionDipole::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theSplitting << theApply << theSubtractionTest 
     << theIgnoreCuts << theRealEmissionME << theUnderlyingBornME 
     << thePartners << theTildeKinematics << theInvertedTildeKinematics 
     << theReweights << theRealEmitter << theRealEmission << theRealSpectator 
     << theSubtractionParameters << theMergingMap << theSplittingMap 
     << theIndexMap << theUnderlyingBornDiagrams << theRealEmissionDiagrams 
     << theBornToRealDiagrams << theRealToBornDiagrams
     << lastRealEmissionKey << lastUnderlyingBornKey 
     << theBornEmitter << theBornSpectator << ounit(theLastSubtractionScale,GeV) 
     << ounit(theLastSplittingScale,GeV) << ounit(theLastSubtractionPt,GeV) 
     << ounit(theLastSplittingPt,GeV) << theLastSubtractionZ
     << theLastSplittingZ << theShowerApproximation 
     << theRealShowerSubtraction << theVirtualShowerSubtraction 
     << theLoopSimSubtraction << theRealEmissionScales << theFactory
     << ounit(theShowerHardScale,GeV) << ounit(theShowerScale,GeV) 
     << theShowerParameters << theIsInShowerPhasespace << theIsAboveCutoff;
}

void SubtractionDipole::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theSplitting >> theApply >> theSubtractionTest 
     >> theIgnoreCuts >> theRealEmissionME >> theUnderlyingBornME 
     >> thePartners >> theTildeKinematics >> theInvertedTildeKinematics 
     >> theReweights >> theRealEmitter >> theRealEmission >> theRealSpectator 
     >> theSubtractionParameters >> theMergingMap >> theSplittingMap 
     >> theIndexMap >> theUnderlyingBornDiagrams >> theRealEmissionDiagrams 
     >> theBornToRealDiagrams >> theRealToBornDiagrams
     >> lastRealEmissionKey >> lastUnderlyingBornKey 
     >> theBornEmitter >> theBornSpectator >> iunit(theLastSubtractionScale,GeV) 
     >> iunit(theLastSplittingScale,GeV) >> iunit(theLastSubtractionPt,GeV) 
     >> iunit(theLastSplittingPt,GeV) >> theLastSubtractionZ
     >> theLastSplittingZ >> theShowerApproximation 
     >> theRealShowerSubtraction >> theVirtualShowerSubtraction 
     >> theLoopSimSubtraction >> theRealEmissionScales >> theFactory
     >> iunit(theShowerHardScale,GeV) >> iunit(theShowerScale,GeV) 
     >> theShowerParameters >> theIsInShowerPhasespace >> theIsAboveCutoff;
  lastMatchboxXComb(theLastXComb);
  typedef multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator spit;
  pair<spit,spit> kr = theSplittingMap.equal_range(lastUnderlyingBornKey);
  lastRealEmissionInfo = kr.first;
  for ( ; lastRealEmissionInfo != kr.second; ++lastRealEmissionInfo )
    if ( process(lastRealEmissionInfo->second.first) == lastXComb().mePartonData() )
      break;
}

Ptr<MatchboxFactory>::tptr SubtractionDipole::factory() const {
  return theFactory;
}

void SubtractionDipole::factory(Ptr<MatchboxFactory>::tptr f) {
  theFactory = f;
}


void SubtractionDipole::Init() {

  static ClassDocumentation<SubtractionDipole> documentation
    ("SubtractionDipole represents a dipole subtraction "
     "term in the formalism of Catani and Seymour.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<SubtractionDipole,MEBase>
describeSubtractionDipole("Herwig::SubtractionDipole", "Herwig.so");
