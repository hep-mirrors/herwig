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
#include "Pythia7/Repository/UseRandom.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ShowerConstrainer.h"
#include "ShowerParticle.h"
#include "Pythia7/PDT/EnumParticles.h"

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
  // have already their colour partner fixed and that don't have children
  // (the latter requirement is relaxed in the case isDecayCase is true). 
  // Build a map which has as key one of these particles (i.e. a pointer 
  // to a ShowerParticle object) and as a corresponding value the vector
  // of all its possible *normal* candidate colour partners, defined as follows:
  //  --- have colour, and no children (this is not required in the case 
  //      isDecayCase is true);
  //  --- if both are initial (incoming) state particles, or if both are 
  //      final (outgoing) state particles, then both must have the
  //      same (non-null) colourLine() or the same (non-null) antiColourLine();
  //  --- if one is an initial (incoming) state particle and the other is
  //      a final (outgoing) state particle, then the (non-null) colourLine() 
  //      of one of them must match the (non-null) antiColourLine() of the other.
  // Notice that this definition exclude the special case of baryon-violating
  // processes (as in R-parity Susy), which will show up as particles
  // without candidate colour partners, and that we will be treated a part later
  // (this means that no modifications of the following loop is needed!)
  map<tShoParPtr, tCollecShoParPtr> mapParticleToCandidatePartners;
  for ( CollecShoParPtr::const_iterator cit = particles.begin(); 
	cit != particles.end(); ++cit ) {
    if ( (*cit)->data().coloured()  &&  ! (*cit)->partners()[ ShowerIndex::QCD ]  &&
	 ( (*cit)->children().size() == 0  ||  isDecayCase ) ) {
      tCollecShoParPtr partners;
      for ( CollecShoParPtr::const_iterator cjt = particles.begin(); 
	    cjt != particles.end(); ++cjt ) {
	if ( cit != cjt  &&  (*cjt)->data().coloured() &&
	     ( (*cit)->children().size() == 0  ||  isDecayCase ) ) { 
	  bool isPartner = false;
	  if ( ( (*cit)->isFinalState()  &&  ! (*cjt)->isFinalState() )  ||
	       ( ! (*cit)->isFinalState()  &&  (*cjt)->isFinalState() ) ) {
	    if ( ( (*cit)->colourLine()  &&  
		   (*cit)->colourLine() == (*cjt)->colourLine() )  ||
		 ( (*cit)->antiColourLine()  &&  
		   (*cit)->antiColourLine() == (*cjt)->antiColourLine() ) ) {
	      isPartner = true;
	    }
	  } else { 
	    if ( ( (*cit)->colourLine()  && 
		   (*cit)->colourLine() == (*cjt)->antiColourLine() )  ||
		 ( (*cit)->antiColourLine()  && 
		   (*cit)->antiColourLine() == (*cjt)->colourLine() ) ) {
	      isPartner = true;
	    }
	  } 
	  if ( isPartner ) {
	    partners.push_back( *cjt );
	  }
	}
      }
      mapParticleToCandidatePartners.insert( pair<tShoParPtr, tCollecShoParPtr>( *cit, partners ) );
    }
  }

  // Loop now over the map we have just filled (mapParticleToCandidatePartners)
  // and treat only those particles that have some candidate colour partners,
  // whereas those that don't have any are collected in the vector specialCases
  // and will be treated later.
  tCollecShoParPtr specialCases; 
  for ( map<tShoParPtr, tCollecShoParPtr>::const_iterator cit = 
	  mapParticleToCandidatePartners.begin(); 
	cit != mapParticleToCandidatePartners.end(); ++cit ) {
    if ( cit->second.size() ) {

      //***LOOKHERE*** In the case of more than one candidate colour partners,
      //               our treatment is based on two assumptions:
      //               1) the choice of which is the colour partner is done
      //                  *randomly* between the available candidates.
      //               2) the choice of which is the colour partner is done
      //                  *independently* from each particle: in other words,
      //                  if for a particle "i" its selected colour partner is  
      //                  the particle "j", then the colour partner of "j" 
      //                  does not have to be necessarily "i".
      //               To me (A.R.) both assumptions seem reasonable, but
      //               I am not 100% sure! Be careful...
      int position = UseRandom::irnd( cit->second.size() );

      pair<Energy,Energy> pairScales = 
	calculateInitialEvolutionScales( pair<tShoParPtr,tShoParPtr>( cit->first, 
								      cit->second[position] ) );
      cit->first->setEvolutionScale( ShowerIndex::QCD, pairScales.first );
      cit->first->setPartner( ShowerIndex::QCD, cit->second[position] );
    } else {
      specialCases.push_back( cit->first );
    }
  } 

  if ( specialCases.size() ) {

    //***LOOKHERE*** ADD TREATMENT OF SPECIAL CASES (BARYON-VIOLATING PROCESSES)
    //               AT THE MOMENT ONLY A WARNING IS ISSUED.
    generator()->logWarning( Exception("PartnerFinder::setQCDInitialEvolutionScales "
				       "***Special cases not yet treated*** ", 
				       Exception::warning) );
  }
  
  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "PartnerFinder::debuggingInfo "
		       << " ===> START DEBUGGING <=== " << endl;
    // To have a nice output, it is convenient to define two maps.
    // One numbers the ShowerParticles objects from 1; the other,
    // numbers the ShowerColourLine object from 501 (just to avoid
    // confusion between the two numbering). 
    map<tShoParPtr,int> mapShoPar;
    map<tShoColinePtr,int> mapShoColine;
    int countShoPar = 1, countShoColine = 501;
    for ( CollecShoParPtr::const_iterator cit = particles.begin(); 
	  cit != particles.end(); ++cit ) {
      mapShoPar.insert( pair<tShoParPtr,int>( *cit, countShoPar++ ) );
      tShoColinePtr shoColine = tShoColinePtr();
      for ( int i = 0; i < 2; i++ ) {
	if ( i == 0 ) {
	  shoColine = (*cit)->colourLine();
	} else {
	  shoColine = (*cit)->antiColourLine();
	}   
	if ( shoColine  &&  mapShoColine.find( shoColine ) == mapShoColine.end() ) {
	  mapShoColine.insert( pair<tShoColinePtr,int>( shoColine, countShoColine++ ) );
	}
      }
    }
    for ( CollecShoParPtr::const_iterator cit = particles.begin(); 
	  cit != particles.end(); ++cit ) {
      if ( (*cit)->data().coloured()  && 
	   ( (*cit)->children().size() == 0  ||  isDecayCase ) ) {
	int numShoParticle = 0;
	if ( mapShoPar.find( *cit ) != mapShoPar.end() ) {
	  numShoParticle = mapShoPar.find( *cit )->second;
	}
	int numShoColine = 0;
	if ( (*cit)->colourLine()  &&
	     mapShoColine.find( (*cit)->colourLine() ) != mapShoColine.end() ) {
	  numShoColine = mapShoColine.find( (*cit)->colourLine() )->second;
	}
	int numShoAntiColine = 0;
	if ( (*cit)->antiColourLine()  &&
	     mapShoColine.find( (*cit)->antiColourLine() ) != mapShoColine.end() ) {
	  numShoAntiColine = mapShoColine.find( (*cit)->antiColourLine() )->second;
	}
	int numShoPartner = 0;
	if ( (*cit)->partners()[ ShowerIndex::QCD ]  &&
              mapShoPar.find( (*cit)->partners()[ ShowerIndex::QCD ] ) != mapShoPar.end() ) {
	  numShoPartner = mapShoPar.find( (*cit)->partners()[ ShowerIndex::QCD ] )->second;   
	}
	generator()->log() << "\t" << "Particle number = " << numShoParticle  
	                   << "   PDGName= " << (*cit)->data().PDGName() 
			   << ( (*cit)->isFinalState() ? "   OUT" : "   IN" ) << endl
                           << "\t \t ColourLine number     = " << numShoColine << endl 
                           << "\t \t antiColourLine number = " << numShoAntiColine << endl 
                           << "\t \t Partner = " << numShoPartner << endl; 
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



