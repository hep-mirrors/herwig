// -*- C++ -*-
//
// ShowerApproximationGenerator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerApproximationGenerator class.
//

#include "ShowerApproximationGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ShowerApproximationGenerator::ShowerApproximationGenerator() 
  : thePresamplingPoints(10000), theMaxTry(100000) {}

ShowerApproximationGenerator::~ShowerApproximationGenerator() {}

double ShowerApproximationGenerator::generateFraction(tcPDPtr pd, double r, double xmin) const {
  if ( pd->coloured() || pd->id() == ParticleID::gamma ) {
    return pow(xmin,r);
  }
  double x0 = 1.e-5;
  return 1. + x0 - x0*pow((1.+x0)/x0,r);
}

double ShowerApproximationGenerator::invertFraction(tcPDPtr pd, double x, double xmin) const {
  if ( pd->coloured() || pd->id() == ParticleID::gamma ) {
    return log(x)/log(xmin);
  }
  double x0 = 1.e-5;
  return log((1.-x+x0)/x0)/log((1.+x0)/x0);
}

bool ShowerApproximationGenerator::prepare() {

  tSubProPtr oldSub = lastIncomingXComb->subProcess();

  tcStdXCombPtr cIncomingXC = lastIncomingXComb;

  theLastMomenta.resize(cIncomingXC->mePartonData().size());
  theLastRandomNumbers.resize(thePhasespace->nDim(theLastMomenta.size()-2) + 2);

  double x1 =
    oldSub->incoming().first->momentum().plus() /
    lastIncomingXComb->lastParticles().first->momentum().plus();

  theLastRandomNumbers[0] = invertFraction(oldSub->incoming().first->dataPtr(),x1,
					   lastIncomingXComb->cuts()->x1Min());

  double x2 =
    oldSub->incoming().second->momentum().minus() /
    lastIncomingXComb->lastParticles().second->momentum().minus();

  theLastRandomNumbers[1] = invertFraction(oldSub->incoming().second->dataPtr(),x2,
					   lastIncomingXComb->cuts()->x2Min());

  Boost toCMS = 
    (oldSub->incoming().first->momentum() +
     oldSub->incoming().second->momentum()).findBoostToCM();

  theLastMomenta[0] = oldSub->incoming().first->momentum();
  theLastMomenta[0].boost(toCMS);
  theLastMomenta[1] = oldSub->incoming().second->momentum();
  theLastMomenta[1].boost(toCMS);

  ParticleVector::const_iterator out = oldSub->outgoing().begin();
  vector<Lorentz5Momentum>::iterator p = theLastMomenta.begin() + 2;
  for ( ; out != oldSub->outgoing().end(); ++out, ++p ) {
    *p = (**out).momentum();
    p->boost(toCMS);
  }

  theLastPartons.first = 
    oldSub->incoming().first->data().produceParticle(oldSub->incoming().first->momentum());
  theLastPartons.second = 
    oldSub->incoming().second->data().produceParticle(oldSub->incoming().second->momentum());

  thePhasespace->invertKinematics(theLastMomenta,&theLastRandomNumbers[2]);

  theLastBornXComb->clean();
  theLastBornXComb->fill(lastIncomingXComb->lastParticles(),theLastPartons,
			 theLastMomenta,theLastRandomNumbers);

  if ( !theLastBornXComb->cuts()->initSubProcess(theLastBornXComb->lastSHat(), 
						 theLastBornXComb->lastY(), 
						 theLastBornXComb->mirror()) )
    return false;

  theLastBornME->setXComb(theLastBornXComb);

  if ( !theLastBornME->generateKinematics(&theLastRandomNumbers[2]) )
    return false;

  CrossSection bornXS = theLastBornME->dSigHatDR();

  if ( bornXS == ZERO )
    return false;

  return true;

}

bool ShowerApproximationGenerator::generate(const vector<double>& r) {

  theLastBornXComb->clean();

  double x = generateFraction(theLastPartons.first->dataPtr(),r[0],
			      lastIncomingXComb->cuts()->x1Min());
  Energy Q = lastIncomingXComb->lastParticles().first->momentum().plus();
  Energy mass = theLastPartons.first->dataPtr()->mass();
  double xi = (4.*sqr(x*Q) - sqr(mass))/(4.*sqr(Q)*x);
  Lorentz5Momentum p1(ZERO,ZERO,xi*Q);
  p1.setMass(mass); p1.rescaleEnergy();
  theLastPartons.first->set5Momentum(p1);

  x = generateFraction(theLastPartons.second->dataPtr(),r[1],
		       lastIncomingXComb->cuts()->x2Min());
  Q = lastIncomingXComb->lastParticles().second->momentum().minus();
  mass = theLastPartons.second->dataPtr()->mass();
  xi = (4.*sqr(x*Q) - sqr(mass))/(4.*sqr(Q)*x);
  Lorentz5Momentum p2(ZERO,ZERO,-xi*Q);
  p2.setMass(mass); p2.rescaleEnergy();
  theLastPartons.second->set5Momentum(p2);

  theLastPresamplingMomenta.resize(theLastMomenta.size());

  Boost toCMS = 
    (theLastPartons.first->momentum() +
     theLastPartons.second->momentum()).findBoostToCM();

  theLastPresamplingMomenta[0] = theLastPartons.first->momentum();
  theLastPresamplingMomenta[0].boost(toCMS);
  theLastPresamplingMomenta[1] = theLastPartons.second->momentum();
  theLastPresamplingMomenta[1].boost(toCMS);

  theLastBornXComb->fill(lastIncomingXComb->lastParticles(),theLastPartons,
			 theLastPresamplingMomenta,r);

  if ( !theLastBornXComb->cuts()->initSubProcess(theLastBornXComb->lastSHat(), 
						 theLastBornXComb->lastY(), 
						 theLastBornXComb->mirror()) )
    return false;

  theLastBornME->setXComb(theLastBornXComb);

  if ( !theLastBornME->generateKinematics(&r[2]) )
    return false;

  CrossSection bornXS = theLastBornME->dSigHatDR();

  if ( bornXS == ZERO )
    return false;

  return true;

}

void ShowerApproximationGenerator::restore() {

  tSubProPtr oldSub = lastIncomingXComb->subProcess();

  theLastPartons.first->set5Momentum(oldSub->incoming().first->momentum());
  theLastPartons.second->set5Momentum(oldSub->incoming().second->momentum());

  theLastBornXComb->clean();
  theLastBornXComb->fill(lastIncomingXComb->lastParticles(),theLastPartons,
			 theLastMomenta,theLastRandomNumbers);

  theLastBornME->setXComb(theLastBornXComb);
  theLastBornME->generateKinematics(&theLastRandomNumbers[2]);
  theLastBornME->dSigHatDR();

}

void ShowerApproximationGenerator::
handle(EventHandler & eh, const tPVector &,
       const Hint &) {

  lastIncomingXComb = dynamic_ptr_cast<tStdXCombPtr>(eh.lastXCombPtr());
  if ( !lastIncomingXComb )
    throw Exception() << "expecting a standard event handler"
		      << Exception::abortnow;

  if ( lastIncomingXComb->lastProjector() )
    lastIncomingXComb = lastIncomingXComb->lastProjector();

  const StandardXComb& xc = *lastIncomingXComb;

  map<cPDVector,set<Ptr<ShowerApproximationKernel>::ptr> >::const_iterator
    kernelit = theKernelMap.find(xc.mePartonData());

  if ( kernelit == theKernelMap.end() ) {
    list<MatchboxFactory::SplittingChannel> channels =
      theFactory->getSplittingChannels(lastIncomingXComb);
    set<Ptr<ShowerApproximationKernel>::ptr> newKernels;
    for ( list<MatchboxFactory::SplittingChannel>::const_iterator c =
	    channels.begin(); c != channels.end(); ++c ) {
      Ptr<ShowerApproximationKernel>::ptr kernel =
	new_ptr(ShowerApproximationKernel());
      kernel->setBornXComb(c->bornXComb);
      kernel->setRealXComb(c->realXComb);
      kernel->setTildeXCombs(c->tildeXCombs);
      kernel->setDipole(c->dipole);
      kernel->showerApproximation(theShowerApproximation);
      kernel->presamplingPoints(thePresamplingPoints);
      kernel->maxtry(theMaxTry);
      kernel->showerApproximationGenerator(this);
      if ( kernel->dipole()->bornEmitter() > 1 &&
	   kernel->dipole()->bornSpectator() > 1 ) {
	kernel->ptCut(ffPtCut());
      } else if ( ( kernel->dipole()->bornEmitter() > 1 &&
		    kernel->dipole()->bornSpectator() < 2 ) ||
		  ( kernel->dipole()->bornEmitter() < 2 &&
		    kernel->dipole()->bornSpectator() > 1 ) ) {
	kernel->ptCut(fiPtCut());
      } else {
	assert(kernel->dipole()->bornEmitter() < 2 &&
	       kernel->dipole()->bornSpectator() < 2);
	kernel->ptCut(iiPtCut());
      }
      newKernels.insert(kernel);
    }
    theKernelMap[xc.mePartonData()] = newKernels;
    kernelit = theKernelMap.find(xc.mePartonData());
  }

  if ( kernelit->second.empty() )
    return;

  const set<Ptr<ShowerApproximationKernel>::ptr>& kernels = kernelit->second;

  theLastBornME = (**kernels.begin()).dipole()->underlyingBornME();
  theLastBornME->phasespace(thePhasespace);
  theLastBornXComb = (**kernels.begin()).bornXComb();

  if ( !prepare() )
    return;

  Energy winnerPt = ZERO;
  Ptr<ShowerApproximationKernel>::ptr winnerKernel;

  for ( set<Ptr<ShowerApproximationKernel>::ptr>::const_iterator k =
	  kernels.begin(); k != kernels.end(); ++k ) {
    if ( (**k).generate() != 0. && (*k)->dipole()->lastPt() > winnerPt){
      winnerKernel = *k;
      winnerPt = winnerKernel->dipole()->lastPt();
    }
  }

  if ( !winnerKernel || winnerPt == ZERO )
    return;

  SubProPtr oldSub = lastIncomingXComb->subProcess();
  SubProPtr newSub;

  try {
    newSub = winnerKernel->realXComb()->construct();
  } catch(Veto&) {
    return;
  }

  tParticleSet firstS = oldSub->incoming().first->siblings();
  assert(firstS.empty() || firstS.size() == 1);
  if ( !firstS.empty() ) {
    eh.currentStep()->removeParticle(*firstS.begin());
  }

  tParticleSet secondS = oldSub->incoming().second->siblings();
  assert(secondS.empty() || secondS.size() == 1);
  if ( !secondS.empty() ) {
    eh.currentStep()->removeParticle(*secondS.begin());
  }

  // prevent the colliding particles from disappearing
  // in the initial state and appearing
  // in the final state when we've cut off all their
  // (physical) children; only applies to the case
  // where we have a parton extractor not build from
  // noPDF, so check wether the incoming particle
  // doesnt equal the incoming parton -- this needs fixing in ThePEG
  PPtr dummy = new_ptr(Particle(getParticleData(ParticleID::gamma)));
  bool usedDummy = false;
  if ( eh.currentStep()->incoming().first != oldSub->incoming().first ) {
    eh.currentStep()->addDecayProduct(eh.currentStep()->incoming().first,dummy);
    usedDummy = true;
  }
  if ( eh.currentStep()->incoming().second != oldSub->incoming().second ) {
    eh.currentStep()->addDecayProduct(eh.currentStep()->incoming().second,dummy);
    usedDummy = true;
  }

  eh.currentStep()->removeSubProcess(oldSub);
  eh.currentStep()->addSubProcess(newSub);

  // get rid of the dummy
  if ( usedDummy ) {
    eh.currentStep()->removeParticle(dummy);
  }

  eh.select(winnerKernel->realXComb());

  winnerKernel->realXComb()->recreatePartonBinInstances(winnerKernel->realXComb()->lastScale());
  winnerKernel->realXComb()->refillPartonBinInstances(&(xc.lastRandomNumbers()[0]));

  winnerKernel->realXComb()->pExtractor()->constructRemnants(winnerKernel->realXComb()->partonBinInstances(),
							     newSub, eh.currentStep());

}

IBPtr ShowerApproximationGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerApproximationGenerator::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ShowerApproximationGenerator::persistentOutput(PersistentOStream & os) const {
  os << theShowerApproximation << thePhasespace << theFactory 
     << theKernelMap << thePresamplingPoints << theMaxTry 
     << lastIncomingXComb << theLastBornME << ounit(theLastMomenta,GeV) 
     << ounit(theLastPresamplingMomenta,GeV) << theLastRandomNumbers
     << theLastBornXComb << theLastPartons;
}

void ShowerApproximationGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theShowerApproximation >> thePhasespace >> theFactory 
     >> theKernelMap >> thePresamplingPoints >> theMaxTry 
     >> lastIncomingXComb >> theLastBornME >> iunit(theLastMomenta,GeV) 
     >> iunit(theLastPresamplingMomenta,GeV) >> theLastRandomNumbers
     >> theLastBornXComb >> theLastPartons;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ShowerApproximationGenerator,StepHandler>
  describeHerwigShowerApproximationGenerator("Herwig::ShowerApproximationGenerator", "HwMatchbox.so");

void ShowerApproximationGenerator::Init() {

  static ClassDocumentation<ShowerApproximationGenerator> documentation
    ("ShowerApproximationGenerator generates emissions according to a "
     "shower approximation entering a NLO matching.");


  static Reference<ShowerApproximationGenerator,ShowerApproximation> interfaceShowerApproximation
    ("ShowerApproximation",
     "Set the shower approximation to sample.",
     &ShowerApproximationGenerator::theShowerApproximation, false, false, true, false, false);


  static Reference<ShowerApproximationGenerator,MatchboxPhasespace> interfacePhasespace
    ("Phasespace",
     "The phase space generator to use.",
     &ShowerApproximationGenerator::thePhasespace, false, false, true, false, false);


  static Reference<ShowerApproximationGenerator,MatchboxFactory> interfaceFactory
    ("Factory",
     "The factory object to use.",
     &ShowerApproximationGenerator::theFactory, false, false, true, false, false);

  static Parameter<ShowerApproximationGenerator,unsigned long> interfacePresamplingPoints
    ("PresamplingPoints",
     "Set the number of presampling points.",
     &ShowerApproximationGenerator::thePresamplingPoints, 10000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximationGenerator,unsigned long> interfaceMaxTry
    ("MaxTry",
     "Set the number of maximum attempts.",
     &ShowerApproximationGenerator::theMaxTry, 100000, 1, 0,
     false, false, Interface::lowerlim);

}

