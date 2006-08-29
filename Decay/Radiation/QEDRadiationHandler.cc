// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QEDRadiationHandler class.
//

#include "QEDRadiationHandler.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "DecayRadiationGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Handlers/EventHandler.h"

using namespace Herwig;

QEDRadiationHandler::QEDRadiationHandler() 
{
  // only include electroweak gauge bosons
  _decayingParticles.push_back( 22);
  _decayingParticles.push_back( 23);
  _decayingParticles.push_back(-23);
  _decayingParticles.push_back( 24);
  // only include the charged leptons
  _decayProducts.push_back( 11);
  _decayProducts.push_back( 13);
  _decayProducts.push_back( 15);
  _decayProducts.push_back(-11);
  _decayProducts.push_back(-13);
  _decayProducts.push_back(-15);
}

void QEDRadiationHandler::
handle(EventHandler & eh, const tPVector & tagged,
       const Hint & hint) throw(Veto, Stop, Exception) {
  // find the potential decaying particles to be considered
  set<tPPtr> parents;
  for(unsigned int ix=0;ix<tagged.size();++ix) {
    long  id=tagged[ix]->id();
    if(find(_decayProducts.begin(),_decayProducts.end(),id)!=_decayProducts.end()) {
      PPtr par=tagged[ix]->parents()[0];
      id=par->id();
      if(tagged[ix]->parents()[0]->mass()>0.&&
	 find(_decayingParticles.begin(),_decayingParticles.end(),id)!=
	 _decayingParticles.end()) parents.insert(par);
    }
  }
  set<tPPtr>::const_iterator sit;
  StepPtr step=eh.currentStep();
  // loop over parents
  for(sit=parents.begin();sit!=parents.end();++sit) {
    // extract children
    ParticleVector children=(**sit).children();
    // store number of children
    unsigned int initsize  =children.size();
    // generate photons
    ParticleVector newchildren=_generator->generatePhotons(**sit,children);
    // if photons produced add as children and to step
    for(unsigned int ix=initsize;ix<newchildren.size();++ix) {
      (**sit).addChild(newchildren[ix]);
      step->addDecayProduct(newchildren[ix]);
    }
  }
}

void QEDRadiationHandler::persistentOutput(PersistentOStream & os) const {
  os << _generator << _decayingParticles << _decayProducts;
}

void QEDRadiationHandler::persistentInput(PersistentIStream & is, int) {
  is >> _generator >> _decayingParticles >> _decayProducts;
}

ClassDescription<QEDRadiationHandler> QEDRadiationHandler::initQEDRadiationHandler;
// Definition of the static class description member.

void QEDRadiationHandler::Init() {

  static ClassDocumentation<QEDRadiationHandler> documentation
    ("The QEDRadiationHandler class is designed to be used as a PostSubProcessHandler"
     "so that the same approach as for radiation in decays can be used for resonances"
     "produced as part of the hard process");

  static Reference<QEDRadiationHandler,DecayRadiationGenerator>
    interfaceRadiationGenerator
    ("RadiationGenerator",
     "Reference to the DecayRadiationGenerator",
     &QEDRadiationHandler::_generator, false, false, true, false, false);
 
  static ParVector<QEDRadiationHandler,long> interfaceDecayingParticles
    ("DecayingParticles",
     "List of PDF codes of the particles which should have radiation"
     " generated for them.",
     &QEDRadiationHandler::_decayingParticles, -1, long(24), 0, 0,
     false, false, Interface::nolimits);

  static ParVector<QEDRadiationHandler,long> interfaceDecayProducts
    ("DecayProducts",
     "List of PDG codes of the particles which should be present"
     " as decay products for the radiation to be generated.",
     &QEDRadiationHandler::_decayProducts, -1, long(11), 0, 0,
     false, false, Interface::nolimits);
  
}

