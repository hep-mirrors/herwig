// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartnerFinder class.
//

#include "PartnerFinder.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
// #include "Pythia7/Interface/Parameter.h" 
#include "Pythia7/Interface/Reference.h" 
#include "Pythia7/Repository/EventGenerator.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ShowerConstrainer.h"
#include "ShowerParticle.h"

using namespace Herwig;


PartnerFinder::~PartnerFinder() {}


void PartnerFinder::persistentOutput(PersistentOStream & os) const {
  // os << ;
}


void PartnerFinder::persistentInput(PersistentIStream & is, int) {
  // is >> ;
}


ClassDescription<PartnerFinder> PartnerFinder::initPartnerFinder;
// Definition of the static class description member.


void PartnerFinder::Init() {

  static ClassDocumentation<PartnerFinder> documentation
    ("This class is responsible for finding the partners for each interaction types ",
     "and within the evolution scale range specified by the ShowerConstrainer ",
     "then to determine the initial evolution scales for each pair of partners.");

}

//-----------------------------------------------------------------------------

bool PartnerFinder::setQCDInitialEvolutionScales( const tShoConstrPtr showerConstrainer,
						  const CollecShoParPtr particles,
						  const bool isDecayCase ) {
  bool isOK = true;

  // Loop over  particles  and consider only coloured particles which don't
  // have already their colour particles fixed and that don't have children
  // (the latter requirement is relaxed in the case isDecayCase is true). 
  for ( CollecShoParPtr::const_iterator cit = particles.begin(); 
	cit != particles.end(); ++cit ) {
    if ( (*cit)->data().coloured()  &&  ! (*cit)->partners()[ ShowerIndex::QCD ]  &&
	 ( (*cit)->children().size() == 0  ||  isDecayCase ) ) {

      //***LOOKHERE*** For simplicity we are making some strong simplifying
      //               assumptions, that MUST be removed later on
      //               (see Hadronization/ClusterFinder for a more correct
      //                treatment of a very similar problem) :
      //               --- no baryon violation in R-parity violating Susy;
      //               --- for gluons which have colour and anticolour partners
      //                   we randomly pick only one of them.
      bool hasAntiColourPartner = false;
      if ( (*cit)->colourLine() ) hasAntiColourPartner = true; 
      bool hasColourPartner = false;
      if ( (*cit)->antiColourLine() ) hasColourPartner = true;
      if ( hasAntiColourPartner && hasColourPartner ) {
	if ( rndbool() ) {
	  hasAntiColourPartner = false;
	} else {
	  hasColourPartner = false;
	}
      }
      CollecShoParPtr::const_iterator cjt = cit; 
      for ( ++cjt; cjt != particles.end(); ++cjt ) {
	if ( ( hasAntiColourPartner &&
	       (*cit)->colourLine() == (*cjt)->antiColourLine() ) ||
	     ( hasColourPartner &&
	       (*cit)->antiColourLine() == (*cjt)->colourLine() ) ) {
	  pair<Energy,Energy> pairScales = 
	    calculateInitialEvolutionScales( pair<tShoParPtr,tShoParPtr>( *cit, *cjt ) );
	  (*cit)->setEvolutionScale( ShowerIndex::QCD, pairScales.first );
	  (*cit)->setPartner( ShowerIndex::QCD, *cjt );
	  (*cjt)->setEvolutionScale( ShowerIndex::QCD, pairScales.second );
	  (*cjt)->setPartner( ShowerIndex::QCD, *cit );
	} 
      }
    } 
  }
  
  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "PartnerFinder::debuggingInfo "
		       << " ===> START DEBUGGING <=== " << endl;
    for ( CollecShoParPtr::const_iterator cit = particles.begin(); 
	  cit != particles.end(); ++cit ) {
      if ( (*cit)->data().coloured()  && 
	   ( (*cit)->children().size() == 0  ||  isDecayCase ) ) {
	generator()->log() << "\t particle = " << (*cit)->data().PDGName()
			   << "   partner = " << ( (*cit)->partners()[ ShowerIndex::QCD ] ?
						   (*cit)->partners()[ ShowerIndex::QCD ]
						   ->data().PDGName() : 0 ) << endl;
      }
    }
    generator()->log() << "PartnerFinder::debuggingInfo "
		       << " ===> END DEBUGGING <=== " << endl;
  }
  
  return isOK;
}


bool PartnerFinder::setQEDInitialEvolutionScales( const tShoConstrPtr showerConstrainer,
						  const CollecShoParPtr particles,
						  const bool isDecayCase ) {
  bool isOK = true;

  //***LOOKHERE*** To be implemented only if you want to have electromagnetic
  //               bremsstrahlung. You should use the data() method of
  //               the ShowerParticle objects to access to their properties
  //               like the electric charge, and then find the partners.
  //               For each find pair of partners, call
  //                    calculateInitialEvolutionScales
  //               and then fills:
  //                    particle1->evolutionScales()[ ShowerIndex::QED ] = scale; 
  //                    particle2->evolutionScales()[ ShowerIndex::QED ] = scale; 
  //                    particle1->partners()[ ShowerIndex::QED ] = particle2; 
  //                    particle2->partners()[ ShowerIndex::QED ] = particle1; 

  return isOK;
}


bool PartnerFinder::setEWKInitialEvolutionScales( const tShoConstrPtr showerConstrainer,
						  const CollecShoParPtr particles,
						  const bool isDecayCase ) {

  bool isOK = true;

  //***LOOKHERE*** To be implemented only if you want to have electroweak
  //               bremsstrahlung. You should use the data() method of
  //               the ShowerParticle objects to access to their properties
  //               like the electroweak charge, and then find the partners.
  //               For each find pair of partners, call
  //                    calculateInitialEvolutionScales
  //               and then fills:
  //                    particle1->evolutionScales()[ ShowerIndex::EWK ] = scale; 
  //                    particle2->evolutionScales()[ ShowerIndex::EWK ] = scale; 
  //                    particle1->partners()[ ShowerIndex::EWK ] = particle2; 
  //                    particle2->partners()[ ShowerIndex::EWK ] = particle1; 

  return isOK;
}


pair<Energy,Energy> PartnerFinder::
calculateInitialEvolutionScales( const pair<tShoParPtr,tShoParPtr> & particlePair ) {
  Energy firstQ = Energy();
  Energy secondQ = Energy();

  //***LOOKHERE***  Use the kinematical relationships between the 4-momenta
  //                of the pair of particles to calculate their initial
  //                evolution scales: firstQ and secondQ. 
  //                BELOW IT IS JUST A SIMPLE FAKE, JUST TO GET SOME VALUES

  double angle = particlePair.first->momentum().vect().angle( particlePair.second->momentum().vect() );
  firstQ  = particlePair.first->momentum().e() * angle;
  secondQ = particlePair.second->momentum().e() * angle;
  
  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "PartnerFinder::calculateInitialEvolutionScales " << endl
		       << "\t first  : " << particlePair.first->data().PDGName() 
                       << "\t p = " << particlePair.first->momentum() << endl
		       << "\t second : " << particlePair.second->data().PDGName() 
                       << "\t p = " << particlePair.second->momentum() << endl
                       << "\t angle = " << angle 
		       << "\t firstQ = " << firstQ / GeV 
		       << "\t secondQ = " << secondQ / GeV << "  (GeV) " << endl;
  }
  
  return pair<Energy,Energy>(firstQ, secondQ);
}



