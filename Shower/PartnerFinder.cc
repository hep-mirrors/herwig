// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartnerFinder class.
//

#include "PartnerFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/Interface/Parameter.h" 
#include "ThePEG/Interface/Reference.h" 
#include "ThePEG/Repository/EventGenerator.h" 
#include "ThePEG/Repository/UseRandom.h" 
#include "Herwig++/Utilities/HwDebug.h"
#include "ShowerConstrainer.h"
#include "ShowerParticle.h"
#include "ThePEG/PDT/EnumParticles.h"

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
						  const ShowerParticleVector particles,
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
  map<tShowerParticlePtr, tShowerParticleVector> mapParticleToCandidatePartners;
  for ( ShowerParticleVector::const_iterator cit = particles.begin(); 
	cit != particles.end(); ++cit ) {
    if ( (*cit)->data().coloured()  &&  ! (*cit)->partners()[ ShowerIndex::QCD ]  &&
	 ( (*cit)->children().size() == 0  ||  isDecayCase ) ) {
      tShowerParticleVector partners;
      for ( ShowerParticleVector::const_iterator cjt = particles.begin(); 
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
      mapParticleToCandidatePartners.insert( pair<tShowerParticlePtr, tShowerParticleVector>( *cit, partners ) );
    }
  }

  // Loop now over the map we have just filled (mapParticleToCandidatePartners)
  // and treat only those particles that have some candidate colour partners,
  // whereas those that don't have any are collected in the vector specialCases
  // and will be treated later.
  tShowerParticleVector specialCases; 
  for ( map<tShowerParticlePtr, tShowerParticleVector>::const_iterator cit = 
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
	calculateInitialEvolutionScales( pair<tShowerParticlePtr,tShowerParticlePtr>( cit->first, cit->second[position] ), showerConstrainer);
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
    generator()->log() << "PartnerFinder::debuggingInfo full ______________________________________________" << endl;
    // To have a nice output, it is convenient to define two maps.
    // One numbers the ShowerParticles objects from 1; the other,
    // numbers the ShowerColourLine object from 501 (just to avoid
    // confusion between the two numbering). 
    map<tShowerParticlePtr,int> mapShoPar;
    map<tShoColinePtr,int> mapShoColine;
    int countShoPar = 1, countShoColine = 501;
    for ( ShowerParticleVector::const_iterator cit = particles.begin(); 
	  cit != particles.end(); ++cit ) {
      mapShoPar.insert( pair<tShowerParticlePtr,int>( *cit, countShoPar++ ) );
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
    for ( ShowerParticleVector::const_iterator cit = particles.begin(); 
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
//     generator()->log() << "PartnerFinder::debuggingInfo "
// 		       << " ===> END DEBUGGING <=== " << endl;
  }
  
  return isOK;
}


bool PartnerFinder::setQEDInitialEvolutionScales( const tShoConstrPtr showerConstrainer,
						  const ShowerParticleVector particles,
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

  // ***ACHTUNG!*** just copy the scales and partners from the QCD method 
  // in order to test the possibility of multiple interactions that compete in the 
  // shower.  The partners have to be determined according to some physical 
  // idea, i.e. they should be CHARGE partners and not colour partners. 

  for ( ShowerParticleVector::const_iterator cit = particles.begin(); 
	cit != particles.end(); ++cit ) {
    (*cit)->setEvolutionScale(ShowerIndex::QED, (*cit)->evolutionScales()[ ShowerIndex::QCD ]); 
    (*cit)->setPartner(ShowerIndex::QED, (*cit)->partners()[ ShowerIndex::QCD ]); 
  }

  return isOK;
}


bool PartnerFinder::setEWKInitialEvolutionScales( const tShoConstrPtr showerConstrainer,
						  const ShowerParticleVector particles,
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
calculateInitialEvolutionScales( const pair<tShowerParticlePtr,tShowerParticlePtr> & particlePair, const tShoConstrPtr showerConstrainer) {
  Energy firstQ = Energy();
  Energy secondQ = Energy();

  //***LOOKHERE***  Use the kinematical relationships between the 4-momenta
  //                of the pair of particles to calculate their initial
  //                evolution scales: firstQ and secondQ. 
  //                BELOW IT IS JUST A SIMPLE FAKE, JUST TO GET SOME VALUES

//   double angle = particlePair.first->momentum().vect().angle( particlePair.second->momentum().vect() );
//   firstQ  = particlePair.first->momentum().e() * angle;
//   secondQ = particlePair.second->momentum().e() * angle;

  Lorentz5Momentum p1, p2; 
  Lorentz5Momentum p, n; 
  p1 = particlePair.first->momentum(); 
  p2 = particlePair.second->momentum(); 
  p = p1; 
  p.boost((p1+p2).findBoostToCM());
  n = Lorentz5Momentum(0.0, -p.vect()); 
  firstQ = sqrt(2.*p*(p+n)); 
  p = p2; 
  p.boost((p1+p2).findBoostToCM());
  n = Lorentz5Momentum(0.0, - p.vect()); 
  secondQ = sqrt(2.*p*(p+n));   

  ////////////////////////////////////////////////////////////////////////////
  // only a hack for the moment! comment/uncomment for normal/asymmetric
  // if desired. 
  // get asymmetric distribution in x, xbar plane: 
  if (showerConstrainer->asyPS()) {
    Energy Q = sqrt(sqr(p1+p2)); 
    double r = p1.m()/Q; 
    //double v = sqrt(1.-sqr(r)); 
    if (particlePair.first->id() < 6 && particlePair.first->id() > 0) { 
      firstQ = 2.0*Q*sqrt(1.-2.*r); 
      secondQ = sqr(Q)/firstQ; 
      //secondQ = Q*(1.+v)*(firstQ/Q*(1.-v) +
      //v*(1.+v))/(4.*firstQ/Q+sqr(v)-1.);
    } else if (particlePair.first->id() > -6 
	       && particlePair.first->id() < 0) { 
      secondQ = 2.0*Q*sqrt(1.-2.*r); 
      firstQ = sqr(Q)/secondQ;    
      //firstQ = Q*(1.+v)*(secondQ/Q*(1.-v) +
      //v*(1.+v))/(4.*secondQ/Q+sqr(v)-1.);
    }
  }
  // swap the above values for even eventnumbers, 
  // (ie uncorrelated or randomly)
  if (showerConstrainer->rndPS()) {  
    Energy temp; 
    if (generator()->currentEventNumber()%2) {
      temp = firstQ; firstQ = secondQ; secondQ = temp;     
    }
  }
  ////////////////////////////////////////////////////////////////////////////

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "PartnerFinder::calculateInitialEvolutionScales full ____________________________" << endl
		       << "\t first:   " << particlePair.first->data().PDGName() 
                       << "\t p = " << particlePair.first->momentum() << endl
		       << "\t second:  " << particlePair.second->data().PDGName() 
                       << "\t p = " << particlePair.second->momentum() << endl
      //    << "\t angle = " << angle 
		       << "\t firstQ = " << firstQ / GeV 
		       << " GeV, secondQ = " << secondQ / GeV << " GeV " << endl;
  }
  
  return pair<Energy,Energy>(firstQ, secondQ);
}



