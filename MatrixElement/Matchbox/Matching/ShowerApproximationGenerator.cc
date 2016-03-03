// -*- C++ -*-
//
// ShowerApproximationGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerApproximationGenerator class.
//

#include <config.h>
#include "ShowerApproximationGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
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
  : thePresamplingPoints(2000), theMaxTry(100000), theFreezeGrid(500000),
    theDoCompensate(false) {}

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

bool ShowerApproximationGenerator::prepare(bool didproject) {

  tSubProPtr oldSub = lastIncomingXComb->subProcess();

  tcStdXCombPtr cIncomingXC = lastIncomingXComb;

  bool hasFractions =
    thePhasespace->haveX1X2() ||
    cIncomingXC->mePartonData().size() == 3;

  theLastMomenta.resize(cIncomingXC->mePartonData().size());

  if ( !hasFractions )
    theLastRandomNumbers.resize(thePhasespace->nDim(cIncomingXC->mePartonData()) + 2);
  else
    theLastRandomNumbers.resize(thePhasespace->nDim(cIncomingXC->mePartonData()));

  if ( !hasFractions ) {

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

  }

  theLastMomenta = cIncomingXC->meMomenta();

  theLastPartons.first = 
    oldSub->incoming().first->data().produceParticle(oldSub->incoming().first->momentum());
  theLastPartons.second = 
    oldSub->incoming().second->data().produceParticle(oldSub->incoming().second->momentum());

  thePhasespace->setXComb(lastIncomingXComb);

  // this is a brute force fix for private ticket #241 ; only done to get fixed
  // for the release but will need to be looked at in more detail later on by
  // cleaning up the XCombs for these cases
  if ( theLastMomenta.size() == 3 && didproject ) {
    // boost them where they belong so invertKinematics is doing something sensible
    Boost toLab = (lastIncomingXComb->lastPartons().first->momentum() + 
		   lastIncomingXComb->lastPartons().second->momentum()).boostVector();
    for ( int i = 0; i < 3; ++i )
      theLastMomenta[i].boost(toLab);
  }

  thePhasespace->invertKinematics(theLastMomenta,
				  !hasFractions ? &theLastRandomNumbers[2] : &theLastRandomNumbers[0]);

  theLastBornXComb->clean();
  theLastBornXComb->fill(lastIncomingXComb->lastParticles(),theLastPartons,
			 theLastMomenta,theLastRandomNumbers);

  if ( !theLastBornXComb->cuts()->initSubProcess(theLastBornXComb->lastSHat(), 
						 theLastBornXComb->lastY(), 
						 theLastBornXComb->mirror()) )
    return false;

  theLastBornME->setXComb(theLastBornXComb);

  if ( !theLastBornME->generateKinematics(!hasFractions ? &theLastRandomNumbers[2] : &theLastRandomNumbers[0]) )
    return false;

  CrossSection bornXS = theLastBornME->dSigHatDR();

  if ( bornXS == ZERO )
    return false;

  return true;

}

bool ShowerApproximationGenerator::generate(const vector<double>& r) {

  theLastBornXComb->clean();

  bool hasFractions =
    thePhasespace->haveX1X2() ||
    theLastBornXComb->mePartonData().size() == 3;

  if ( !hasFractions ) {

    double x = generateFraction(theLastPartons.first->dataPtr(),r[0],
				lastIncomingXComb->cuts()->x1Min());
    Energy Q = lastIncomingXComb->lastParticles().first->momentum().plus();
    Energy mass = theLastPartons.first->dataPtr()->mass();
    double xi = (sqr(x*Q) - sqr(mass))/(sqr(Q)*x);
    Lorentz5Momentum p1(ZERO,ZERO,xi*Q/2.);
    p1.setMass(mass); p1.rescaleEnergy();
    theLastPartons.first->set5Momentum(p1);

    x = generateFraction(theLastPartons.second->dataPtr(),r[1],
			 lastIncomingXComb->cuts()->x2Min());
    Q = lastIncomingXComb->lastParticles().second->momentum().minus();
    mass = theLastPartons.second->dataPtr()->mass();
    xi = (sqr(x*Q) - sqr(mass))/(sqr(Q)*x);
    Lorentz5Momentum p2(ZERO,ZERO,-xi*Q/2.);
    p2.setMass(mass); p2.rescaleEnergy();
    theLastPartons.second->set5Momentum(p2);

  } else {

    theLastBornME->setXComb(theLastBornXComb);
    theLastBornXComb->lastParticles(lastIncomingXComb->lastParticles());
    theLastBornXComb->lastP1P2(make_pair(0.0, 0.0));
    theLastBornXComb->lastS(lastIncomingXComb->lastS());

    if ( !theLastBornME->generateKinematics(&r[0]) )
      return false;

    theLastPartons.first->set5Momentum(theLastBornME->lastMEMomenta()[0]);
    theLastPartons.second->set5Momentum(theLastBornME->lastMEMomenta()[1]);

  }

  theLastPresamplingMomenta.resize(theLastMomenta.size());

  Boost toCMS = 
    (theLastPartons.first->momentum() +
     theLastPartons.second->momentum()).findBoostToCM();

  theLastPresamplingMomenta[0] = theLastPartons.first->momentum();
  if ( !hasFractions )
    theLastPresamplingMomenta[0].boost(toCMS);
  theLastPresamplingMomenta[1] = theLastPartons.second->momentum();
  if ( !hasFractions )
    theLastPresamplingMomenta[1].boost(toCMS);

  if ( hasFractions ) {
    for ( size_t k = 2; k < theLastBornME->lastMEMomenta().size(); ++k )
      theLastPresamplingMomenta[k] = theLastBornME->lastMEMomenta()[k];
  }

  theLastBornXComb->fill(lastIncomingXComb->lastParticles(),theLastPartons,
			 theLastPresamplingMomenta,r);

  if ( !theLastBornXComb->cuts()->initSubProcess(theLastBornXComb->lastSHat(), 
						 theLastBornXComb->lastY(), 
						 theLastBornXComb->mirror()) )
    return false;

  if ( !hasFractions ) {
    theLastBornME->setXComb(theLastBornXComb);
    if ( !theLastBornME->generateKinematics(&r[2]) )
      return false;
  }

  CrossSection bornXS = theLastBornME->dSigHatDR();

  if ( bornXS == ZERO )
    return false;

  return true;

}

void ShowerApproximationGenerator::restore() {

  theLastBornXComb->clean();

  bool hasFractions =
    thePhasespace->haveX1X2() ||
    theLastBornXComb->mePartonData().size() == 3;

  if ( !hasFractions ) {

    tSubProPtr oldSub = lastIncomingXComb->subProcess();
    
    theLastPartons.first->set5Momentum(oldSub->incoming().first->momentum());
    theLastPartons.second->set5Momentum(oldSub->incoming().second->momentum());

  } else {

    theLastBornME->setXComb(theLastBornXComb);
    theLastBornXComb->lastParticles(lastIncomingXComb->lastParticles());
    theLastBornXComb->lastP1P2(make_pair(0.0, 0.0));
    theLastBornXComb->lastS(lastIncomingXComb->lastS());

    theLastBornME->generateKinematics(&theLastRandomNumbers[0]);

    theLastPartons.first->set5Momentum(theLastBornME->lastMEMomenta()[0]);
    theLastPartons.second->set5Momentum(theLastBornME->lastMEMomenta()[1]);

  }

  theLastBornXComb->fill(lastIncomingXComb->lastParticles(),theLastPartons,
			 theLastMomenta,theLastRandomNumbers);

  if ( !hasFractions ) {
    theLastBornME->setXComb(theLastBornXComb);
    theLastBornME->generateKinematics(&theLastRandomNumbers[2]);
  }

  theLastBornME->dSigHatDR();

}

void ShowerApproximationGenerator::
handle(EventHandler & eh, const tPVector &,
       const Hint &) {
  theFactory->setHardTreeEmitter(-1);
  theFactory->setHardTreeSpectator(-1);
  theFactory->setHardTreeSubprocess(SubProPtr());


  lastIncomingXComb = dynamic_ptr_cast<tStdXCombPtr>(eh.lastXCombPtr());
  if ( !lastIncomingXComb )
    throw Exception() << "ShowerApproximationGenerator::handle(): Expecting a standard event handler."
		      << Exception::runerror;

  bool didproject = false;

  if ( lastIncomingXComb->lastProjector() ) {
    lastIncomingXComb = lastIncomingXComb->lastProjector();
    didproject = true;
  }

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
      kernel->freezeGrid(theFreezeGrid);
      kernel->showerApproximationGenerator(this);
      kernel->doCompensate(theDoCompensate);
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
  if ( theLastBornME->phasespace()->wantCMS() != thePhasespace->wantCMS() ) {
    throw Exception() << "Mismatch in centre-of-mass-system requirements of hard matrix element phasespace ("
                      << (theLastBornME->phasespace()->wantCMS()?"true":"false")
                      << ") and shower approximation phasespace ("
                      << (thePhasespace->wantCMS()?"true":"false") << ")"
                      << Exception::abortnow;
  }
  theLastBornME->phasespace(thePhasespace);
  theLastBornXComb = (**kernels.begin()).bornXComb();

  if ( !prepare(didproject) )
    return;

  Energy winnerPt = ZERO;
  Ptr<ShowerApproximationKernel>::ptr winnerKernel;

  try {
    for ( set<Ptr<ShowerApproximationKernel>::ptr>::const_iterator k =
	    kernels.begin(); k != kernels.end(); ++k ) {
      if ( (**k).generate() != 0. && (*k)->dipole()->lastPt() > winnerPt){
	winnerKernel = *k;
	winnerPt = winnerKernel->dipole()->lastPt();
      }
    }
  } catch(ShowerApproximationKernel::MaxTryException&) {
    throw Exception() << "Too many tries needed to generate the matrix element correction in '"
		      << name() << "'" << Exception::eventerror;
  }

  if ( !winnerKernel || winnerPt == ZERO )
    return;
  
  //Hardest emission should be this one.

  winnerKernel->realXComb()->lastShowerScale(sqr(winnerPt));
  winnerKernel->bornXComb()->lastShowerScale(sqr(winnerPt));
  lastIncomingXComb->lastShowerScale(sqr(winnerPt));


  SubProPtr oldSub = lastIncomingXComb->subProcess();
  SubProPtr newSub;

  try {
    tcDiagPtr bornDiag = lastIncomingXComb->lastDiagram();
    tcDiagPtr realDiag = 
      winnerKernel->dipole()->realEmissionDiagram(bornDiag);
    winnerKernel->realXComb()->externalDiagram(realDiag);
    newSub = winnerKernel->realXComb()->construct();
  } catch(Veto&) {
    return;
  }

  if ( !theShowerApproximation->needsTruncatedShower() ){
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
  else{
    theFactory->setHardTreeSubprocess(newSub);
    theFactory->setHardTreeEmitter(winnerKernel->dipole()->bornEmitter()); 
    theFactory->setHardTreeSpectator(winnerKernel->dipole()->bornSpectator()); 
  }

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
     << theKernelMap << thePresamplingPoints << theMaxTry << theFreezeGrid
     << lastIncomingXComb << theLastBornME << ounit(theLastMomenta,GeV) 
     << ounit(theLastPresamplingMomenta,GeV) << theLastRandomNumbers
     << theLastBornXComb << theLastPartons << theDoCompensate;
}

void ShowerApproximationGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theShowerApproximation >> thePhasespace >> theFactory 
     >> theKernelMap >> thePresamplingPoints >> theMaxTry >> theFreezeGrid
     >> lastIncomingXComb >> theLastBornME >> iunit(theLastMomenta,GeV) 
     >> iunit(theLastPresamplingMomenta,GeV) >> theLastRandomNumbers
     >> theLastBornXComb >> theLastPartons >> theDoCompensate;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ShowerApproximationGenerator,StepHandler>
  describeHerwigShowerApproximationGenerator("Herwig::ShowerApproximationGenerator", "Herwig.so");

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
     &ShowerApproximationGenerator::thePresamplingPoints, 2000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximationGenerator,unsigned long> interfaceMaxTry
    ("MaxTry",
     "Set the number of maximum attempts.",
     &ShowerApproximationGenerator::theMaxTry, 100000, 1, 0,
     false, false, Interface::lowerlim);

  static Parameter<ShowerApproximationGenerator,unsigned long> interfaceFreezeGrid
    ("FreezeGrid",
     "",
     &ShowerApproximationGenerator::theFreezeGrid, 500000, 1, 0,
     false, false, Interface::lowerlim);

  static Switch<ShowerApproximationGenerator,bool> interfaceDoCompensate
    ("DoCompensate",
     "",
     &ShowerApproximationGenerator::theDoCompensate, false, false, false);
  static SwitchOption interfaceDoCompensateYes
    (interfaceDoCompensate,
     "Yes",
     "",
     true);
  static SwitchOption interfaceDoCompensateNo
    (interfaceDoCompensate,
     "No",
     "",
     false);

}

