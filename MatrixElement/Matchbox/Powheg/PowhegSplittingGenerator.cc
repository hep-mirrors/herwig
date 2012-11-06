// -*- C++ -*-
//
// PowhegSplittingGenerator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegSplittingGenerator class.
//

#include "PowhegSplittingGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Handlers/StdXCombGroup.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDT/EnumParticles.h"

#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "PowhegInclusiveME.h"

#include "ThePEG/PDF/PDF.h"

using namespace Herwig;

PowhegSplittingGenerator::PowhegSplittingGenerator() 
  : StepHandler(), 
    theFFPtCut(1.0*GeV), theFFScreeningScale(ZERO),
    theFIPtCut(1.0*GeV), theFIScreeningScale(ZERO),
    theIIPtCut(1.0*GeV), theIIScreeningScale(ZERO),
    discardNoEmissions(false), theVerbose(false),
    theDiscardNext(false) {}

PowhegSplittingGenerator::~PowhegSplittingGenerator() {
  for ( GeneratorMap::iterator g = theGeneratorMap.begin();
	g != theGeneratorMap.end(); ++g )
    delete g->second.second;
  theGeneratorMap.clear();
}

void PowhegSplittingGenerator::
handle(EventHandler & eh, const tPVector &, const Hint &) {

  if ( theVerbose ) {
    generator()->log() << "PowhegSplittingGenerator generating real emission off the sub-process\n"
		       << (*(eh.lastXCombPtr()->subProcess())) << "\n"
		       << "with x1 = " << eh.lastXCombPtr()->lastX1() 
		       << " x2 = " << eh.lastXCombPtr()->lastX2() << "\n" << flush;
  }

  if ( !generate(eh) ) {
    if ( theVerbose ) {
      generator()->log() << "PowhegSplittingGenerator did not select radiation above the IR cutoff\n" << flush;
    }
    if ( discardNext() ) {
      setDiscardNext(false);
      if ( theVerbose ) {
	generator()->log() << "Splitting kernels have been presampled, will discard this event.\n" << flush;
      }
      throw Veto();
    }
    if ( discardNoEmissions )
      throw Veto();
    veto(eh);
    return;
  }

  if ( theVerbose ) {
    generator()->log() << "PowhegSplittingGenerator selected the kernel '"
		       << lastSplitting()->name() << "' to generate radiation\n" << flush;
  }

  if ( discardNext() ) {
    setDiscardNext(false);
    if ( theVerbose ) {
      generator()->log() << "Splitting kernels have been presampled, will discard this event.\n" << flush;
    }
    throw Veto();
  }

  SubProPtr oldSub =
    lastSplitting()->bornSubProcess();

  SubProPtr newSub;

  try {
    Energy pt = lastSplitting()->projectionDipole()->lastPt();
    newSub = lastSplitting()->construct(pt);
  } catch(Veto&) {
    if ( theVerbose ) {
      generator()->log() << "The generated real emission process did not pass the cuts.\n" << flush;
    }
    veto(eh);
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

  eh.select(lastSplitting()->lastXCombPtr());
  lastSplitting()->lastXCombPtr()->
    recreatePartonBinInstances(lastSplitting()->lastXCombPtr()->lastScale());
  eh.lastExtractor()->constructRemnants(lastSplitting()->lastXCombPtr()->partonBinInstances(),
					newSub, eh.currentStep());

  if ( theVerbose ) {
    generator()->log() << "PowhegSplittingGenerator generated the real emission sub-process\n"
		       << (*(eh.lastXCombPtr()->subProcess())) << "\n"
		       << "with x1 = " << eh.lastXCombPtr()->lastX1() 
		       << " x2 = " << eh.lastXCombPtr()->lastX2() << "\n" << flush;
  }

}

bool PowhegSplittingGenerator::generate(EventHandler & eh) {

  pair<GeneratorMap::iterator,GeneratorMap::iterator> generators
    = getGenerators(eh);

  Energy winnerPt = 0.*GeV;
  Energy pt;
  GeneratorMap::iterator winner = theGeneratorMap.end();

  for ( GeneratorMap::iterator gen = generators.first; gen != generators.second; ++gen ) {
    assert(gen->second.first->lastHeadXCombPtr() == eh.lastXCombPtr());
    assert(gen->second.first->projectionDipole()->lastHeadXCombPtr() == eh.lastXCombPtr());
    assert(gen->second.first->projectionDipole()->lastXCombPtr() ==
	   gen->second.first->lastXCombPtr());
    pt = generate(gen->second);
    if ( pt > winnerPt ) {
      winnerPt = pt;
      winner = gen;
    }
  }

  if ( winner == theGeneratorMap.end() ) {
    theLastSplitting = Ptr<PowhegSplittingKernel>::tptr();
    return false;
  }

  theLastSplitting = winner->second.first;

  return true;

}

void PowhegSplittingGenerator::veto(EventHandler & eh) const {

  tSubProPtr sub = eh.currentStep()->subProcesses().front();

  if ( sub->incoming().first->coloured() ) {
    sub->incoming().first->vetoScale(ZERO);
  }

  if ( sub->incoming().second->coloured() ) {
    sub->incoming().first->vetoScale(ZERO);
  }

  for ( ParticleVector::const_iterator p = sub->outgoing().begin();
	p != sub->outgoing().end(); ++p )
    if ( (**p).coloured() )
      (**p).vetoScale(ZERO);

}

IBPtr PowhegSplittingGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegSplittingGenerator::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PowhegSplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(theFFPtCut,GeV) << ounit(theFFScreeningScale,GeV) 
     << ounit(theFIPtCut,GeV) << ounit(theFIScreeningScale,GeV) 
     << ounit(theIIPtCut,GeV) << ounit(theIIScreeningScale,GeV) 
     << discardNoEmissions << theVerbose;
}

void PowhegSplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theFFPtCut,GeV) >> iunit(theFFScreeningScale,GeV) 
     >> iunit(theFIPtCut,GeV) >> iunit(theFIScreeningScale,GeV) 
     >> iunit(theIIPtCut,GeV) >> iunit(theIIScreeningScale,GeV) 
     >> discardNoEmissions >> theVerbose;
}

pair<PowhegSplittingGenerator::GeneratorMap::iterator,
     PowhegSplittingGenerator::GeneratorMap::iterator>
PowhegSplittingGenerator::getGenerators(EventHandler& eh) {

  tXCombPtr xc = eh.lastXCombPtr();

  const ThePEG::StandardXComb& xcref = *dynamic_ptr_cast<tStdXCombPtr>(xc);

  if ( theVerbose ) {
    generator()->log() << "getting splitting generators for xcomb "
		       << xc << " and process ";
    generator()->log() << xcref.mePartonData()[0]->PDGName() << " "
		       << xcref.mePartonData()[1]->PDGName() << " -> ";
    for ( ThePEG::cPDVector::const_iterator pid =
	    xcref.mePartonData().begin() + 2;
	  pid != xcref.mePartonData().end(); ++pid )
      generator()->log() << (**pid).PDGName() << " ";
    generator()->log() << "\n" << flush;
  }

  pair<GeneratorMap::iterator,GeneratorMap::iterator> res =
    theGeneratorMap.equal_range(xc);

  if ( res.first != res.second ) {
    if ( theVerbose )
      generator()->log() << "generators already known\n" << flush;
    return res;
  }

  tStdXCombGroupPtr xcGrp =
    dynamic_ptr_cast<tStdXCombGroupPtr>(xc);

  if ( !xcGrp ) {
    if ( theVerbose )
      generator()->log() << "xcomb is not an xcomb group\n" << flush;
    return make_pair(theGeneratorMap.end(),theGeneratorMap.end());
  }

  Ptr<PowhegInclusiveME>::tptr me =
    dynamic_ptr_cast<Ptr<PowhegInclusiveME>::tptr>((*xcGrp).matrixElement());

  if ( !me ) {
    if ( theVerbose )
      generator()->log() << "matrix element is not a powheg inclusive me\n" << flush;
    return make_pair(theGeneratorMap.end(),theGeneratorMap.end());
  }

  for ( vector<Ptr<PowhegSplittingKernel>::ptr>::iterator k =
	  me->splittingKernels().begin(); k != me->splittingKernels().end(); ++k ) {

    if ( !(**k).apply() )
      continue;

    if ( theVerbose )
      generator()->log() << "initializing generator for kernel '" << (**k).name() << "'\n" << flush;

    assert((**k).lastHeadXCombPtr() == eh.lastXCombPtr());
    assert((**k).projectionDipole()->lastHeadXCombPtr() == eh.lastXCombPtr());
    assert((**k).lastXCombPtr() == (**k).projectionDipole()->lastXCombPtr());

    (**k).splittingGenerator(this);
    if ( (**k).projectionDipole()->realEmitter() > 1 &&
	 (**k).projectionDipole()->realSpectator() > 1 ) {
      (**k).ptCut(theFFPtCut);
      (**k).screeningScale(theFFScreeningScale);
    } else if ( (**k).projectionDipole()->realEmitter() < 2 &&
		(**k).projectionDipole()->realSpectator() < 2 ) {
      (**k).ptCut(theIIPtCut);
      (**k).screeningScale(theIIScreeningScale);
    } else {
      (**k).ptCut(theFIPtCut);
      (**k).screeningScale(theFIScreeningScale);
    }
    ExponentialGeneratorPtr gen = new ExponentialGenerator();
    gen->sampling_parameters().maxtry = (**k).maxtry();
    gen->sampling_parameters().presampling_points = (**k).presamplingPoints();
    gen->function(*k);
    gen->initialize();

    theGeneratorMap.insert(make_pair(xc,make_pair(*k,gen)));

  }

  return getGenerators(eh);

}

Energy PowhegSplittingGenerator::generate(pair<Ptr<PowhegSplittingKernel>::ptr,ExponentialGeneratorPtr>& gen) {

  double res = 0.;
  gen.first->splittingGenerator(this);

  while (true) {
    try {
      res = gen.second->generate();
    } catch (exsample::exponential_regenerate&) {
      continue;
    } catch (exsample::hit_and_miss_maxtry&) {
      throw Veto();
    } catch (exsample::selection_maxtry&) {
      throw Veto();
    } 
    break;
  }

  if ( theVerbose ) {
    generator()->log() << "Generating splitting from '" << gen.first->projectionDipole()->name() << "'.\n" << flush;
    if ( res == 0. )
      generator()->log() << "Below infrared cutoff.\n" << flush;
    else
      generator()->log() << "pt/GeV = " << (gen.first->projectionDipole()->lastPt()/GeV) << ".\n" << flush;
  }

  if ( res == 0. )
    return 0.*GeV;

  return gen.first->projectionDipole()->lastPt();

}

void PowhegSplittingGenerator::Init() {

  static ClassDocumentation<PowhegSplittingGenerator> documentation
    ("PowhegSplittingGenerator");


  static Parameter<PowhegSplittingGenerator,Energy> interfaceFFPtCut
    ("FFPtCut",
     "Set the pt infrared cutoff",
     &PowhegSplittingGenerator::theFFPtCut, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<PowhegSplittingGenerator,Energy> interfaceFFScreeningScale
    ("FFScreeningScale",
     "Set the screening scale",
     &PowhegSplittingGenerator::theFFScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<PowhegSplittingGenerator,Energy> interfaceFIPtCut
    ("FIPtCut",
     "Set the pt infrared cutoff",
     &PowhegSplittingGenerator::theFIPtCut, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<PowhegSplittingGenerator,Energy> interfaceFIScreeningScale
    ("FIScreeningScale",
     "Set the screening scale",
     &PowhegSplittingGenerator::theFIScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<PowhegSplittingGenerator,Energy> interfaceIIPtCut
    ("IIPtCut",
     "Set the pt infrared cutoff",
     &PowhegSplittingGenerator::theIIPtCut, GeV, 1.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<PowhegSplittingGenerator,Energy> interfaceIIScreeningScale
    ("IIScreeningScale",
     "Set the screening scale",
     &PowhegSplittingGenerator::theIIScreeningScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<PowhegSplittingGenerator,bool> interfaceVerbose
    ("Verbose",
     "",
     &PowhegSplittingGenerator::theVerbose, false, false, false);
  static SwitchOption interfaceVerboseOn
    (interfaceVerbose,
     "On",
     "",
     true);
  static SwitchOption interfaceVerboseOff
    (interfaceVerbose,
     "Off",
     "",
     false);

  static Switch<PowhegSplittingGenerator,bool> interfaceDiscardNoEmissions
    ("DiscardNoEmissions",
     "Discard events without radiation.",
     &PowhegSplittingGenerator::discardNoEmissions, false, false, false);
  static SwitchOption interfaceDiscardNoEmissionsOn
    (interfaceDiscardNoEmissions,
     "On",
     "Discard events without radiation.",
     true);
  static SwitchOption interfaceDiscardNoEmissionsOff 
    (interfaceDiscardNoEmissions,
     "Off",
     "Do not discard events without radiation.",
     false);

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegSplittingGenerator,StepHandler>
describeHerwigPowhegSplittingGenerator("Herwig::PowhegSplittingGenerator", "HwMatchbox.so");
