// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonEvolver class.
//

#include "NasonEvolver.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/ShowerHandler.h"
#include "DefaultEmissionGenerator.h"

using namespace Herwig;

void NasonEvolver::persistentOutput(PersistentOStream & os) const {
  os << _hardgenerator;
}

void NasonEvolver::persistentInput(PersistentIStream & is, int) {
  is >> _hardgenerator;
}

ClassDescription<NasonEvolver> NasonEvolver::initNasonEvolver;
// Definition of the static class description member.

void NasonEvolver::Init() {

  static ClassDocumentation<NasonEvolver> documentation
    ("The NasonEvolver implements the Nason approach to MC@NLO");

  static RefVector<NasonEvolver,HardestEmissionGenerator> interfaceHardGenerator
    ("HardGenerator",
     "The objects responsible for generating the hardestr emission",
     &NasonEvolver::_hardgenerator, -1, false, false, true, false, false);

}

void NasonEvolver::showerDecay(ShowerTreePtr tree) {
  // set the tree
  currentTree(tree);
  // set up the shower
  setupShower(false);
  cerr << "testing in nason decay shower\n";
  exit(0);



}

vector<ShowerProgenitorPtr> NasonEvolver::setupShower(bool hard) {
  // set the colour partners
  setColourPartners(hard);
  // generate the hardest emission
  hardestEmission();
  // get the particles to be showered
  map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
  vector<ShowerProgenitorPtr> particlesToShower;
  // incoming particles
  for(cit=currentTree()->incomingLines().begin();
      cit!=currentTree()->incomingLines().end();++cit)
    particlesToShower.push_back(((*cit).first));
  assert((particlesToShower.size()==1&&!hard)||(particlesToShower.size()==2&&hard));
  // outgoing particles
  for(cit=currentTree()->outgoingLines().begin();
      cit!=currentTree()->outgoingLines().end();++cit)
    particlesToShower.push_back(((*cit).first));
  // remake the colour partners if needed
  if(currentTree()->hardMatrixElementCorrection()) {
    setColourPartners(hard);
    currentTree()->resetShowerProducts();
  }
  return particlesToShower;
}

void NasonEvolver::hardestEmission() {
  // see if there is an appropriate hard emission generator
  HardestEmissionGeneratorPtr currenthard=HardestEmissionGeneratorPtr();
  cerr << "testing before loop " << _hardgenerator.size() << "\n";
  for(unsigned int ix=0;ix<_hardgenerator.size();++ix) {
    if(!_hardgenerator[ix]->canHandle(currentTree())) continue;
    if(currenthard) {
      if(dynamic_ptr_cast<DefaultEmissionGeneratorPtr>(currenthard)) {
	currenthard=_hardgenerator[ix];
      }
      else if(!dynamic_ptr_cast<DefaultEmissionGeneratorPtr>(_hardgenerator[ix])) {
	ostringstream output;
	output << "There is more than one possible hard emission generator"
	       << " which could be used for ";
	map<ShowerProgenitorPtr,tShowerParticlePtr>::const_iterator cit;
	for(cit=currentTree()->incomingLines().begin();
	    cit!=currentTree()->incomingLines().end();++cit)
	  {output << cit->first->progenitor()->PDGName() << " ";}
	output << " -> ";
	for(cit=currentTree()->outgoingLines().begin();
	    cit!=currentTree()->outgoingLines().end();++cit)
	  {output << cit->first->progenitor()->PDGName() << " ";}
	output << "in NasonEvolver::hardestEmission()\n";
	throw Exception() << output << Exception::runerror;
      }
    }
    else {
      currenthard=_hardgenerator[ix];
    }
  }
  // if no suitable generator return
  if(!currenthard) return;
  // generate the hardest emission
  currenthard->generateHardest(currentTree());
}

