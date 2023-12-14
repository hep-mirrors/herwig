// -*- C++ -*-
//
// PartonicDecayerBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonicDecayerBase class.
//

#include "PartonicDecayerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Hadronization/CluHadConfig.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"

using namespace Herwig;

PartonicDecayerBase::PartonicDecayerBase() : _exclusive(true), 
					     _partontries(100), _inter(false),
					     _shower(false)
{}

void PartonicDecayerBase::doinit() {
  HwDecayerBase::doinit();
  if(_shower) {
    if(!_partnerFinder) {
      throw InitException() << "PartonicDecayerBase::doinit() " << fullName()
			    << " must have the PartonFinder set if"
			    << " showers are being performecd in partonic decays"
			    << Exception::runerror;
    }
    if(!_splittingGenerator) {
      throw InitException() << "PartonicDecayerBase::doinit() " << fullName()
			    << " must have the SplittingGenerator set if"
			    << " showers are being performecd in partonic decays"
			    << Exception::runerror;
    }
    if(!_reconstructor) {
      throw InitException() << "PartonicDecayerBase::doinit() " << fullName()
			    << " must have the KinematicsReconstructor set if"
			    << " showers are being performecd in partonic decays"
			    << Exception::runerror;
    }

  }
}

ParticleVector PartonicDecayerBase::decay(const DecayMode & dm,
					  const Particle & p) const {
  // handling of the decay including the special features of the
  // DecayMode  
  // get the primary products
  tPDVector products=dm.orderedProducts();
  // add products for which the decay mode is all ready specified
  if(!dm.cascadeProducts().empty()) 
    throw Exception() << "PartonicDecayerBase::decay() cannot handle"
		      << " cascadeProducts() from the DecayMode "
		      << Exception::runerror;
  unsigned int ptry(0);
  bool hadronized(false);
  ParticleVector outpart;
  ParticleVector outhad;
  // copy of particle to act as parent so hadronization konws about it
  PPtr ptemp(new_ptr(Particle(p)));
  do {
    hadronized=false;
    // increment loop counter
    ++ptry;
    // perform the primary decay
    ParticleVector partons=decay(p,products);
    // must add partons are children so reshuffling vs leptons will work
    for(unsigned int ix=0;ix<partons.size();++ix) ptemp->addChild(partons[ix]);
    PVector currentlist = partons;
    // perform the shower if needed
    if(_shower) {
      currentlist = shower(p, partons);
      if(currentlist.empty()) continue;
    }
    // split the gluons
    _partonSplitter->split(currentlist);
    // form the clusters
    ClusterVector clusters = _clusterFinder->formClusters(currentlist);
    _clusterFinder->reduceToTwoComponents(clusters);
    tPVector finalHadrons = _clusterFissioner->fission(clusters,false);
    bool lightOK = _lightClusterDecayer->decay(clusters,finalHadrons);
    // abandon child here so always done
    for(unsigned int ix=0;ix<partons.size();++ix) ptemp->abandonChild(partons[ix]);
    // try again if can't reshuffle
    if(!lightOK) continue;
    // decay the remaining clusters
    _clusterDecayer->decay(clusters,finalHadrons);
    hadronized = !duplicateMode(p,finalHadrons);
    if(hadronized) {
      outhad  = ParticleVector(finalHadrons.begin(),finalHadrons.end());
      outpart = partons;
    }
  }
  while(!hadronized&&ptry<_partontries);
  if(ptry>=_partontries) return ParticleVector();
  // return decay products
  if(_inter) {
    return outpart;
  }
  else {
    for(unsigned int ix=0;ix<outhad.size();++ix) {
      tParticleVector parents=outhad[ix]->parents();
      for(unsigned int iy=0;iy<parents.size();++iy) {
	parents[iy]->abandonChild(outhad[ix]);
      }
    }
    for(unsigned int ix=0;ix<outpart.size();++ix) {
      if(!outpart[ix]->coloured()) outhad.push_back(outpart[ix]);
    }
    return outhad;
  }
}

void PartonicDecayerBase::persistentOutput(PersistentOStream & os) const {
  os << _partonSplitter << _clusterFinder << _clusterFissioner
    << _lightClusterDecayer << _clusterDecayer << _exclusive << _partontries
     << _inter << _shower << _splittingGenerator << _partnerFinder << _reconstructor;
}

void PartonicDecayerBase::persistentInput(PersistentIStream & is, int) {
  is >> _partonSplitter >> _clusterFinder >> _clusterFissioner
    >> _lightClusterDecayer >> _clusterDecayer >> _exclusive >> _partontries
     >> _inter >> _shower >> _splittingGenerator >> _partnerFinder >> _reconstructor;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<PartonicDecayerBase,HwDecayerBase>
describeHerwigPartonicDecayerBase("Herwig::PartonicDecayerBase", "HwShower.so HwPartonicDecay.so");

void PartonicDecayerBase::Init() {

  static ClassDocumentation<PartonicDecayerBase> documentation
    ("The PartonicDecayerBase class is the base class for partonic decays"
     " and handles the hadronization of the decays");

  static Reference<PartonicDecayerBase,PartonSplitter> 
    interfacePartonSplitter("PartonSplitter", 
		      "A reference to the PartonSplitter object", 
		      &Herwig::PartonicDecayerBase::_partonSplitter,
		      false, false, true, false);

  static Reference<PartonicDecayerBase,ClusterFinder> 
    interfaceClusterFinder("ClusterFinder", 
		      "A reference to the ClusterFinder object", 
		      &Herwig::PartonicDecayerBase::_clusterFinder,
		      false, false, true, false);

  static Reference<PartonicDecayerBase,ClusterFissioner> 
    interfaceClusterFissioner("ClusterFissioner", 
		      "A reference to the ClusterFissioner object", 
		      &Herwig::PartonicDecayerBase::_clusterFissioner,
		      false, false, true, false);

  static Reference<PartonicDecayerBase,LightClusterDecayer> 
    interfaceLightClusterDecayer("LightClusterDecayer", 
		    "A reference to the LightClusterDecayer object", 
		    &Herwig::PartonicDecayerBase::_lightClusterDecayer,
		    false, false, true, false);

  static Reference<PartonicDecayerBase,ClusterDecayer> 
    interfaceClusterDecayer("ClusterDecayer", 
		       "A reference to the ClusterDecayer object", 
		       &Herwig::PartonicDecayerBase::_clusterDecayer,
		       false, false, true, false);

  static Switch<PartonicDecayerBase,bool> interface_exclusive
    ("Exclusive",
     "Ensure that the hadrons produced in the partonic decays of bottom"
     " and charm baryons do not duplicate the inclusive modes.",
     &PartonicDecayerBase::_exclusive, true, false, false);
  static SwitchOption interface_exclusiveNoDuplication
    (interface_exclusive,
     "Yes",
     "Forbid duplication",
     true);
  static SwitchOption interface_exclusiveDuplication
    (interface_exclusive,
     "No",
     "Duplication allowed",
     false);
  
  static Switch<PartonicDecayerBase,bool> interfaceIntermediates
    ("Intermediates",
     "Whether or not to include the intermediate particles produced by the"
     " cluster alogorithm in the event record.",
     &PartonicDecayerBase::_inter, false, false, false);
  static SwitchOption interfaceIntermediatesIntermediates
    (interfaceIntermediates,
     "Yes",
     "Include the intermediates",
     true);
  static SwitchOption interfaceIntermediatesNoIntermediates
    (interfaceIntermediates,
     "No",
     "Don't include the intermediates.",
     false);

  static Parameter<PartonicDecayerBase,unsigned int> interfacePartonic_Tries
    ("Partonic_Tries",
     "Number of attempts to generator the hadronisation of the decay",
     &PartonicDecayerBase::_partontries, 100, 1, 1000,
     false, false, Interface::limited);
  
  static Switch<PartonicDecayerBase,bool> interfaceShower
    ("Shower",
     "Whether or not to perform the parton shower before the hadronization",
     &PartonicDecayerBase::_shower, false, false, false);
  static SwitchOption interfaceShowerYes
    (interfaceShower,
     "Yes",
     "Perform the shower",
     true);
  static SwitchOption interfaceShowerNo
    (interfaceShower,
     "No",
     "Don't perform the shower",
     false);

  static Reference<PartonicDecayerBase,PartnerFinder> interfacePartnerFinder
    ("PartnerFinder",
     "The parton finder",
     &PartonicDecayerBase::_partnerFinder, false, false, true, true, false);

  static Reference<PartonicDecayerBase,SplittingGenerator> interfaceSplittingGenerator
    ("SplittingGenerator",
     "The splitting generator",
     &PartonicDecayerBase::_splittingGenerator, false, false, true, true, false);
  
  static Reference<PartonicDecayerBase,KinematicsReconstructor> interfaceKinematicsReconstructor
    ("KinematicsReconstructor",
     "The kinematics reconstructor",
     &PartonicDecayerBase::_reconstructor, false, false, true, true, false);

}

bool PartonicDecayerBase::duplicateMode(const Particle & parent,
					const vector<tPPtr> & hadrons) const {
  // if not exclusive return
  if(!_exclusive) return false;
  // now find the hadrons and check them
  cParticleMSet hadronsb;
  bool found(false);
  for (unsigned ix = 0; ix < hadrons.size(); ++ix)
    hadronsb.insert(hadrons[ix]->dataPtr());
  // now check particle's decay modes 
  Selector<tDMPtr>::const_iterator modeptr 
    = parent.dataPtr()->decaySelector().begin();
  Selector<tDMPtr>::const_iterator end     
    = parent.dataPtr()->decaySelector().end();
  // check not a duplicate of a known mode
  for(;modeptr!=end;++modeptr) {
    tcDMPtr mode=(*modeptr).second;
    // check same number of products
    if(mode->products().size() != hadronsb.size()) continue;
    ParticleMSet::const_iterator dit;
    cParticleMSet::const_iterator pit;
    for(dit=mode->products().begin(), pit=hadronsb.begin();
	dit!=mode->products().end(); ++dit,++pit) {
      if((*dit)!=(*pit)) break;
    }
    if(dit != mode->products().end()) continue;
    found = true;
    break;
  }
  return found;
}

void PartonicDecayerBase::dataBaseOutput(ofstream & output,bool header) const {
  // header for MySQL
  if(header) output << "update decayers set parameters=\"";
  // parameters
  output << "newdef  " << name() << ":PartonSplitter " 
	 << _partonSplitter->name() << " \n";
  output << "newdef  " << name() << ":ClusterFinder " 
	 << _clusterFinder->name() << " \n";
  output << "newdef  " << name() << ":ClusterFissioner " 
	 << _clusterFissioner->name() << " \n";
  output << "newdef  " << name() << ":LightClusterDecayer " 
	 << _lightClusterDecayer->name() << " \n";
  output << "newdef  " << name() << ":ClusterDecayer " 
	 << _clusterDecayer->name() << " \n";
  output << "newdef  " << name() << ":Exclusive " <<  _exclusive<< " \n";
  output << "newdef  " << name() << ":Intermediates " << _inter << " \n";
  output << "newdef  " << name() << ":Shower " << _shower << " \n";
  if(_splittingGenerator)
    output << "newdef  " << name() << ":SplittingGenerator " << _splittingGenerator->fullName() << " \n";
  if(_partnerFinder)
    output << "newdef  " << name() << ":PartonFinder " << _partnerFinder->fullName() << " \n";
  if(_reconstructor)
    output << "newdef  " << name() << ":KinematicsReconstructor " << _reconstructor->fullName() << " \n";
  output << "newdef  " << name() << ":Shower " << _shower << " \n";
  output << "newdef  " << name() << ":Partonic_Tries " << _partontries << " \n";
  // footer for MySQL
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

PVector PartonicDecayerBase::shower(const Particle & parent,
				    ParticleVector & partons) const {
  try {
    // set up the shower
    // create the shower particles
    ShowerParticleVector particles;
    for(auto p : partons) {
      particles.push_back(new_ptr(ShowerParticle(*p,2,true)));
    }
    // set up the partners and evolution scales
    _partnerFinder->setInitialEvolutionScales(particles,true,ShowerInteraction::QCD,false);
    // do the shower
    PVector output;
    for(tShowerParticlePtr p : particles) {
      p->initializeFinalState();
      timeLikeShower(p,Branching(),true,output);
    }
    // check if in CMF frame
    Boost beta_cm = parent.momentum().findBoostToCM();
    bool gottaBoost(false);
    LorentzRotation trans = LorentzRotation();
    if(beta_cm.mag() > 1e-12) {
      gottaBoost = true;
      trans.boost(beta_cm);
    }
    bool radiated(false);
    vector<JetKinStruct> jetKinematics;
    for(unsigned int ix=0;ix<particles.size();++ix) {
      JetKinStruct tempJetKin;      
      tempJetKin.parent = particles[ix];
      if(gottaBoost) {
	_reconstructor->deepTransform(tempJetKin.parent,trans);
      }
      tempJetKin.p = particles[ix]->momentum();
      radiated |= _reconstructor->reconstructTimeLikeJet(particles[ix],particles[ix]);
      tempJetKin.q = particles[ix]->momentum();
      jetKinematics.push_back(tempJetKin);
    }
    // find the rescaling factor
    double k = 0.0;
    if(radiated) {
      k = _reconstructor->solveKfactor(parent.mass(), jetKinematics);
      // perform the rescaling and boosts
      for(const JetKinStruct & jet : jetKinematics) {
	LorentzRotation Trafo = _reconstructor->solveBoost(k, jet.q, jet.p);
	_reconstructor->deepTransform(jet.parent,Trafo);
      }
    }
    if(gottaBoost) {
      LorentzRotation rinv = trans.inverse();
      for(const JetKinStruct & jet : jetKinematics) {
	_reconstructor->deepTransform(jet.parent,rinv);
      }
    }
    // add shower particles as children
    for(unsigned int ix=0;ix<partons.size();++ix)
      partons[ix]->addChild(particles[ix]);
    // return final-state of shower
    return output;
  }
  catch (KinematicsReconstructionVeto & ) {
    return PVector();
  }
}

bool PartonicDecayerBase::timeLikeShower(tShowerParticlePtr particle,
					 Branching fb, bool first, PVector & output) const {
  // generate the emission
  if(!fb.kinematics) {
    fb=_splittingGenerator->chooseForwardBranching(*particle,1.,ShowerInteraction::QCD);
    if(fb.kinematics) fb.hard = false;
  }
  // no emission, return
  if(!fb.kinematics) {
    if(particle->spinInfo()) particle->spinInfo()->develop();
    output.push_back(particle);
    return false;
  }
  Branching fc[2] = {Branching(),Branching()};
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->showerKinematics(fb.kinematics);
  // create the children
  ShowerParticleVector children;
  children.reserve(2);
  for(unsigned int ix=0;ix<fb.ids.size()-1;++ix) {
    children.push_back(new_ptr(ShowerParticle(fb.ids[ix+1],true)));
    children[ix]->set5Momentum(Lorentz5Momentum(fb.ids[ix+1]->mass()));
  }
  // update the children
  particle->showerKinematics()->updateChildren(particle, children,1,fb.type);
  // select branchings for children
  for(unsigned int ix=0;ix<children.size()-1;++ix) {
    fc[ix] = _splittingGenerator->chooseForwardBranching(*children[ix],1.,ShowerInteraction::QCD);
  }
  // shower the children
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(fc[ix].kinematics) {
      timeLikeShower(children[ix],fc[ix],false,output);
    }
    else {
      output.push_back(children[ix]);
    }
    if(children[ix]->spinInfo()) children[ix]->spinInfo()->develop();
  }
  particle->showerKinematics()->updateParent(particle, children,1,fb.type);
  // branching has happened
  if(first&&!children.empty())
    particle->showerKinematics()->resetChildren(particle,children);
  if(particle->spinInfo()) particle->spinInfo()->develop();
  return true;
}
