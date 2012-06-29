// -*- C++ -*-
//
// PowhegSplittingKernel.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegSplittingKernel class.
//

#include "PowhegSplittingKernel.h"
#include "PowhegSplittingGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxNLOME.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Handlers/StdXCombGroup.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

PowhegSplittingKernel::PowhegSplittingKernel()
  : ME2byDipoles(), theBornScreening(true), 
    thePresamplingPoints(10000), theMaxTry(100000),
    theScreeningScale(ZERO),
    theBornRandom(0,0), theRadiationRandom(0,0),
    thePresampling(false) {}

PowhegSplittingKernel::~PowhegSplittingKernel() {}

IBPtr PowhegSplittingKernel::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegSplittingKernel::fullclone() const {
  return new_ptr(*this);
}

void PowhegSplittingKernel::splittingGenerator(Ptr<PowhegSplittingGenerator>::tptr gen) { 
  theGenerator = gen; 
}

void PowhegSplittingKernel::setXComb(tStdXCombPtr real) {

  if ( real )
    ME2byDipoles::setXComb(real);
  projectionDipole()->setXComb(real);
  if ( !real )
    return;

  if ( theBornRandom.first == theBornRandom.second ) {

    theBornRandom.first = 
      projectionDipole()->lastHeadXComb().pExtractor()->
      nDims(projectionDipole()->lastHeadXComb().partonBins()).first;
    theBornRandom.second =
      theBornRandom.first + projectionDipole()->underlyingBornME()->nDim();

    theRadiationRandom.first = theBornRandom.second +
      ( projectionDipole()->lastHeadXComb().matrixElement()->nDim() - // includes radiation & the coll remainders
	projectionDipole()->underlyingBornME()->nDim() -
	projectionDipole()->nDimRadiation() );

    theRadiationRandom.second =
      theRadiationRandom.first + projectionDipole()->nDimRadiation();

  }

  if ( thePresamplingXCombs.find(real->head()) == thePresamplingXCombs.end() ) {

    tStdXCombPtr bornxc = real->head();

    thePresamplingXCombs[bornxc] =
      new_ptr(StandardXComb(bornxc->maxEnergy(),bornxc->particles(),
			    bornxc->eventHandlerPtr(),
			    const_ptr_cast<tSubHdlPtr>(bornxc->subProcessHandler()),
			    bornxc->pExtractor(),bornxc->CKKWHandler(),
			    bornxc->partonBins(),bornxc->cuts(),
			    projectionDipole()->underlyingBornME(),
			    bornxc->diagrams(),bornxc->mirror()));

  }

}

tSubProPtr PowhegSplittingKernel::construct(Energy pt) {

  // only cluster real emission

  if ( lastCuts().jetFinder() )
    lastCuts().jetFinder()->minOutgoing(lastXComb().mePartonData().size()-2);

  tSubProPtr subpro = lastXCombPtr()->construct();

  if ( projectionDipole()->realEmitter() == 0 ||
       projectionDipole()->realSpectator() == 0 ) {
    if ( subpro->incoming().first->vetoScale() < 0.0*GeV2 ||
	 subpro->incoming().first->vetoScale() > sqr(pt) )
      subpro->incoming().first->vetoScale(sqr(pt));
  }

  if ( projectionDipole()->realEmitter() == 1 ||
       projectionDipole()->realSpectator() == 1 ) {
    if ( subpro->incoming().second->vetoScale() < 0.0*GeV2 ||
	 subpro->incoming().second->vetoScale() > sqr(pt) )
      subpro->incoming().second->vetoScale(sqr(pt));
  }

  if ( projectionDipole()->realEmitter() > 1 ) {
    if ( subpro->outgoing()[projectionDipole()->realEmitter()-2]->vetoScale() < 0.0*GeV2 ||
	 subpro->outgoing()[projectionDipole()->realEmitter()-2]->vetoScale() > sqr(pt) )
      subpro->outgoing()[projectionDipole()->realEmitter()-2]->vetoScale(sqr(pt));
  }

  if ( projectionDipole()->realSpectator() > 1 ) {
    if ( subpro->outgoing()[projectionDipole()->realSpectator()-2]->vetoScale() < 0.0*GeV2 ||
	 subpro->outgoing()[projectionDipole()->realSpectator()-2]->vetoScale() > sqr(pt) )
      subpro->outgoing()[projectionDipole()->realSpectator()-2]->vetoScale(sqr(pt));
  }

  if ( subpro->outgoing()[projectionDipole()->realEmission()-2]->vetoScale() < 0.0*GeV2 ||
       subpro->outgoing()[projectionDipole()->realEmission()-2]->vetoScale() > sqr(pt) )
    subpro->outgoing()[projectionDipole()->realEmission()-2]->vetoScale(sqr(pt));  

  return subpro;

}

double PowhegSplittingKernel::evaluate() const {

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' evaluating\n";

  if ( !projectionDipole()->underlyingBornME()->
       lastXCombPtr()->willPassCuts() )
    return 0.;

  double fscaleFactor =
    projectionDipole()->realEmissionME()->factorizationScaleFactor();

  Energy2 fscale = 
    fscaleFactor * (sqr(projectionDipole()->lastPt()) + sqr(theScreeningScale));

  double dummy;
  double ratio = ME2byDipoles::evaluate(dummy);

  Energy2 bornSHat =
    projectionDipole()->underlyingBornME()->lastXComb().lastSHat();
  double bornJacobian =
    projectionDipole()->underlyingBornME()->lastXComb().jacobian();

  if ( projectionDipole()->verbose() )
    generator()->log() << "Born sHat/GeV2 = " << (bornSHat/GeV2)
		       << " Born Jacobian = " << bornJacobian << "\n" << flush;

  CrossSection born =
    sqr(hbarc) * scaledBorn(fscale) * bornJacobian / (2.*bornSHat);

  if ( born == ZERO ) {
    if ( projectionDipole()->verbose() )
      generator()->log() << "'" << name() << "' done evaluating\n";
    return 0.0;
  }

  if ( bornScreening() ) {
    born += sqr(hbarc) * scaledBornScreen() * bornJacobian / (2.*bornSHat);
  }

  ratio *= ( projectionDipole()->dSigHatDR(fscale) / born );

  double rscaleFactor =
    projectionDipole()->realEmissionME()->renormalizationScaleFactor();

  double runAlpha =
    SM().alphaS(rscaleFactor * (sqr(projectionDipole()->lastPt())
				+ sqr(theScreeningScale)));

  if ( projectionDipole()->verbose() )
    generator()->log() << "real emission alpha_s = " 
		       << projectionDipole()->realEmissionME()->lastXComb().lastAlphaS()
		       << " pt running alpha_s = " << runAlpha 
		       << " from pt/GeV = " << (projectionDipole()->lastPt()/GeV)
		       << "\n" << flush;

  ratio *= runAlpha / projectionDipole()->realEmissionME()->lastXComb().lastAlphaS();

  if ( projectionDipole()->verbose() )
    generator()->log() << "'" << name() << "' done evaluating\n";

  return ratio;

}

int PowhegSplittingKernel::nDim() const {
  return
    projectionDipole()->lastHeadXComb().pExtractor()->
    nDims(projectionDipole()->lastHeadXComb().partonBins()).first +
    projectionDipole()->lastHeadXComb().pExtractor()->
    nDims(projectionDipole()->lastHeadXComb().partonBins()).second +
    projectionDipole()->underlyingBornME()->nDim() +
    projectionDipole()->nDimRadiation();
}

const vector<double>& PowhegSplittingKernel::parameterPoint() {
  assert(!presampling());
  theLastParameterPoint.resize(nDim());
  copy(projectionDipole()->lastHeadXComb().lastRandomNumbers().begin(),
       projectionDipole()->lastHeadXComb().lastRandomNumbers().begin() + theBornRandom.second,
       theLastParameterPoint.begin());
  copy(projectionDipole()->lastHeadXComb().lastRandomNumbers().begin() + theRadiationRandom.first,
       projectionDipole()->lastHeadXComb().lastRandomNumbers().end(),
       theLastParameterPoint.begin() + theBornRandom.second);
  theLastParameterPoint[evolutionVariable()] = 1.;
  return theLastParameterPoint;
}

int PowhegSplittingKernel::evolutionVariable() const {
  return
    theBornRandom.second +
    projectionDipole()->invertedTildeKinematics()->evolutionVariable();
}

const vector<bool>& PowhegSplittingKernel::sampleFlags() {
  if ( !theFlags.empty() )
    return theFlags;
  theFlags.resize(nDim(),false);
  for ( int k = theBornRandom.second; 
	k < theBornRandom.second + projectionDipole()->nDimRadiation(); ++k )
    theFlags[k] = true;
  return theFlags;
}

const pair<vector<double>,vector<double> >& PowhegSplittingKernel::support() {
  if ( !theSupport.first.empty() )
    return theSupport;
  vector<double> l(nDim(),0.0);
  vector<double> u(nDim(),1.0);
  theSupport.first = l;
  theSupport.second = u;
  return theSupport;
}

void PowhegSplittingKernel::startPresampling() {
  thePresampling = true;
  theXCombBackup = lastHeadXCombPtr();
  lastXCombPtr()->head(thePresamplingXCombs[theXCombBackup]);
  ME2byDipoles::setXComb(lastXCombPtr());
  projectionDipole()->setXComb(lastXCombPtr());
  thePresamplingXCombs[theXCombBackup]->prepare(theXCombBackup->lastParticles());
}

void PowhegSplittingKernel::stopPresampling() {
  thePresampling = false;
  lastXCombPtr()->head(theXCombBackup);
  projectionDipole()->setXComb(lastXCombPtr());
  ME2byDipoles::setXComb(lastXCombPtr());
  theGenerator->setDiscardNext();
}

double PowhegSplittingKernel::evaluate(const vector<double>& p) {

  try {

    if ( projectionDipole()->verbose() )
      generator()->log() << "'" << name() << "' preparing\n";

    if ( presampling() ) {

      if ( projectionDipole()->verbose() )
	generator()->log() << "presampling\n";

      int bdim = thePresamplingXCombs[theXCombBackup]->nDim();

      thePresamplingPoint.resize(bdim);

      copy(p.begin(),p.begin()+theBornRandom.second,
	   thePresamplingPoint.begin());
      copy(p.begin()+theBornRandom.second+projectionDipole()->nDimRadiation(),p.end(),
	   thePresamplingPoint.begin()+theBornRandom.second);

      if ( thePresamplingXCombs[theXCombBackup]->
	   dSigDR(make_pair(0.0,0.0), bdim,
		  &thePresamplingPoint[0]) == ZERO ) {
	if ( projectionDipole()->verbose() )
	  generator()->log() << "Born outside phase space\n";
	return 0.0;
      }

    }

    const double * rad = &p[theBornRandom.second];

    tStdDependentXCombPtr depXComb =
      dynamic_ptr_cast<tStdDependentXCombPtr>(lastXCombPtr());  
    assert(depXComb);
    depXComb->setProcess();
    if ( !projectionDipole()->generateKinematics(rad) ) {
      return 0.0;
    }
    depXComb->remakeIncoming();
    depXComb->setIncomingPartons();

    return evaluate();

  } catch (...) {
    if ( presampling() )
      stopPresampling();
    throw;
  }

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PowhegSplittingKernel::persistentOutput(PersistentOStream & os) const {
  os << theBornScreening << thePresamplingPoints << theMaxTry << ounit(theScreeningScale,GeV)
     << theBornRandom << theRadiationRandom << theLastParameterPoint
     << thePresamplingPoint << thePresampling;
}

void PowhegSplittingKernel::persistentInput(PersistentIStream & is, int) {
  is >> theBornScreening >> thePresamplingPoints >> theMaxTry >> iunit(theScreeningScale,GeV)
     >> theBornRandom >> theRadiationRandom >> theLastParameterPoint
     >> thePresamplingPoint >> thePresampling;
}

void PowhegSplittingKernel::Init() {

  static ClassDocumentation<PowhegSplittingKernel> documentation
    ("PowhegSplittingKernel");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegSplittingKernel,ME2byDipoles>
describeHerwigPowhegSplittingKernel("Herwig::PowhegSplittingKernel", "HwMatchbox.so");
