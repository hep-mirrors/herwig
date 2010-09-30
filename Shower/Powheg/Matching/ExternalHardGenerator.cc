// -*- C++ -*-//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExternalHardGenerator class.
//

#include "ExternalHardGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerTree.h"
#include "Herwig++/Shower/Base/Evolver.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/CKKWVeto.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

void ExternalHardGenerator::persistentOutput(PersistentOStream & os) const {
  os << _CKKWh << _hardGen;
}

void ExternalHardGenerator::persistentInput(PersistentIStream & is, int) {
  is >> _CKKWh >>  _hardGen;
}

ClassDescription<ExternalHardGenerator> ExternalHardGenerator::initExternalHardGenerator;
// Definition of the static class description member.

void ExternalHardGenerator::Init() {

  static ClassDocumentation<ExternalHardGenerator> documentation
    ("The ExternalHardGenerator class generates the hardest emission for"
     "vector boson decays to fermion-antifermion events in the POWHEG approach");

  static Reference<ExternalHardGenerator, PowhegHandler> interfaceCKKWHandler
    ("CKKWHandler",
     "The cascade handler for Powheg",
     &ExternalHardGenerator::_CKKWh, false, false,
     true, false);
  
  static Reference<ExternalHardGenerator,HardestEmissionGenerator> interfaceHardGenerator
    ( "HardGenerator",
      "The hardest emission generator used to fill the dead zone",
      &ExternalHardGenerator::_hardGen, false, false, true, true );
}

HardTreePtr ExternalHardGenerator::generateHardest(ShowerTreePtr tree) {

  ShowerProgenitorPtr 
    QProgenitor    = tree->outgoingLines().begin()->first,
    QbarProgenitor = tree->outgoingLines().rbegin()->first;
  if(QProgenitor->id()<0) swap(QProgenitor, QbarProgenitor);

  //get com energy from the progenitors
  Lorentz5Momentum tot_momentum = QProgenitor->progenitor()->momentum() + QbarProgenitor->progenitor()->momentum();
  Energy2 s = tot_momentum.m2();
  // Get the HardTree from the CKKW handler.
  HardTreePtr hardtree = _CKKWh->getHardTree();
  Energy veto_pt = sqrt( s );
  //evolver vetoes should be set to use shower pt in highest multiplicity case
  if( _CKKWh->highestMult() && hardtree ) veto_pt = hardtree->lowestPt( 1, s );

  for( map< ShowerProgenitorPtr, tShowerParticlePtr >::iterator it
	 = tree->outgoingLines().begin(); 
       it != tree->outgoingLines().end(); ++it ){
    if( ! it->second->coloured() ) continue;
    it->first->maximumpT( veto_pt );
  }     
  for( map< ShowerProgenitorPtr, ShowerParticlePtr >::iterator it
	 = tree->incomingLines().begin(); 
       it != tree->incomingLines().end(); ++it ){
    if( ! it->second->coloured() ) continue;
    it->first->maximumpT( veto_pt );
  }   

  if( ! hardtree ) return HardTreePtr();
  // check the trees match
  if( ! hardtree->connect(tree) ) return HardTreePtr();
  // Return the HardTree
  //fill the dead zone by generating an emission with the hard generator
  //updating the _theHardTree and vetoeing the event if the emission pt 
  //(from lowest pt) is greater than the merge scale
  if(  _CKKWh->lowestMult() && _hardGen ){ 
    //need access to the evolver to generate the hard emission
    tEvolverPtr theEvolver = _CKKWh->getEvolver();
    _hardGen->setEvolver( theEvolver );
    hardtree = _hardGen->generateHardest( tree );
    if( hardtree ){
      hardtree->findNodes();
      //reject the event if newhardtree->lowestPt is greater than pt_cut
      if( hardtree->lowestPt( 1, s ) > _CKKWh->theVeto()->getVeto() ) 
	throw Veto();
    }
    //maximum pt of shower progenitors will also have been 
    //updated to something below the veto scale
    //set the CKKWVeto to use the maximum pt scale as veto scale 
    //rather than the fixed veto scale
    //this means doing it as if it is a highest multiplicity emission 
    _CKKWh->theVeto()->setHighest( true );
  }
  return hardtree;
}

bool ExternalHardGenerator::canHandle(ShowerTreePtr) {
  return true;
}
