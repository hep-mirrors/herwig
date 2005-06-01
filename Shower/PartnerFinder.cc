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
#include "ShowerVariables.h"
#include "ShowerParticle.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;


PartnerFinder::~PartnerFinder() {}

void PartnerFinder::persistentOutput(PersistentOStream & os) const {}

void PartnerFinder::persistentInput(PersistentIStream & is, int) {}

ClassDescription<PartnerFinder> PartnerFinder::initPartnerFinder;
// Definition of the static class description member.

void PartnerFinder::Init() {

  static ClassDocumentation<PartnerFinder> documentation
    ("This class is responsible for finding the partners for each interaction types ",
     "and within the evolution scale range specified by the ShowerVariables ",
     "then to determine the initial evolution scales for each pair of partners.");
  static Parameter<PartnerFinder,int> approach
    ("Approximation",
     "This is a test variable to consider the different approaches of "
     "which colour dipoles of a hard process will shower.",
     &PartnerFinder::_approach, 0, 1, 0,false,false,false);

}

//-----------------------------------------------------------------------------

bool PartnerFinder::setQCDInitialEvolutionScales(const tShowerVarsPtr showerVariables,
						 const ShowerParticleVector particles,
						 const bool isDecayCase) {
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
  typedef map<tShowerParticlePtr, tShowerParticleVector> PartnerMap;
  PartnerMap candidatePartners;
  ShowerParticleVector::const_iterator cit, cjt;
  // Confusing...not correct either... changed by PJS 6/2/04
  /*for(cit = particles.begin(); cit != particles.end(); ++cit) {
    if((*cit)->data().coloured() && !(*cit)->partners()[ShowerIndex::QCD] &&
       ((*cit)->children().size() == 0 || isDecayCase )) {
      tShowerParticleVector partners;
      for(cjt = particles.begin(); cjt != particles.end(); cjt++) {
	if(cit != cjt && (*cjt)->data().coloured() &&
	   ((*cit)->children().size() == 0 || isDecayCase)) { 
	  bool isPartner = false;
	  if(((*cit)->isFinalState() && !(*cjt)->isFinalState()) ||
	       (!(*cit)->isFinalState() && (*cjt)->isFinalState())) {
	    if(((*cit)->colourLine() &&  
		(*cit)->colourLine() == (*cjt)->colourLine()) ||
	       ((*cit)->antiColourLine()  &&  
		(*cit)->antiColourLine() == (*cjt)->antiColourLine())) {
	      isPartner = true;
	    }
	  } else { 
	    if(((*cit)->colourLine()  && 
		(*cit)->colourLine() == (*cjt)->antiColourLine())  ||
	       ((*cit)->antiColourLine()  && 
		(*cit)->antiColourLine() == (*cjt)->colourLine())) {
	      isPartner = true;
	    }
	  } 
	  if(isPartner) partners.push_back(*cjt); 
	}
      }
      candidatePartners[*cit] = partners;
    }
  }*/
  for(cit = particles.begin(); cit != particles.end(); ++cit) {
    if(!(*cit)->data().coloured()) continue;
//     cout << "PF: cit = " << (*cit)->id() << ", " << *cit 
// 	 << " [" << (*cit)->colourLine() << ", " 
// 	 << (*cit)->antiColourLine() << "]" << endl; 
    // We now have a coloured particle
    tShowerParticleVector partners;
    for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
//       cout << "  PF: cjt = " << (*cjt)->id() << ", " << *cjt 
// 	   << " [" << (*cjt)->colourLine() << ", " 
// 	   << (*cjt)->antiColourLine() << "]" << endl; 
      if(!(*cjt)->data().coloured()) continue;
      bool isPartner = false;
#define FS(a) (*a)->isFinalState()
#define CL(a) (*a)->colourLine()
#define ACL(a) (*a)->antiColourLine()
      if((FS(cit) && !FS(cjt)) || (!FS(cit) && FS(cjt)) &&
         ((CL(cit) && CL(cit)==CL(cjt)) || (ACL(cit) && ACL(cit)==ACL(cjt))))
	 isPartner = true;
      else if((CL(cit) && CL(cit)==ACL(cjt)) || (ACL(cit) && ACL(cit)==CL(cjt)))
	isPartner = true;
      if(isPartner) partners.push_back(*cjt);
    }
//     cout << "Initial conditions: "
// 	 << (*cit) << " (" << (*cit)->id() 
// 	 << ") has " << partners.size() << " partners" << endl; 
    if (partners.size() == 0) {
//       cout << "No Partners among the showerParticles, " 
// 	   << "looking among the remnants..." << endl
// 	   << (*cit)->colourLine() << endl
// 	   << (*cit)->antiColourLine() << endl << flush;
      if ((*cit)->colourLine()) 
// 	cout << (*cit)->colourLine()->startParticle() << endl
// 	     << (*cit)->colourLine()->endParticle() << endl << flush;
      if ((*cit)->antiColourLine()) 
// 	cout << (*cit)->antiColourLine()->startParticle() << endl
// 	     << (*cit)->antiColourLine()->endParticle() << endl << flush;
      // set colourpartners 'by hand' from colourLine information
      // this should also give us access to the remnants
      if (CL(cit)) {
	if (CL(cit)->startParticle() == (*cit)) {
	  if (CL(cit)->endParticle()) {
	    partners.push_back(dynamic_ptr_cast<ShowerParticlePtr>(CL(cit)->endParticle()));
	  }
	} else if (CL(cit)->endParticle() == (*cit)) {
	  if (CL(cit)->startParticle()) {
	    partners.push_back(dynamic_ptr_cast<ShowerParticlePtr>(CL(cit)->startParticle()));
	  }
	}
      }
      if (ACL(cit)) {
	if (ACL(cit)->startParticle() == (*cit)) {
	  if (ACL(cit)->endParticle()) {
	    partners.push_back(dynamic_ptr_cast<ShowerParticlePtr>(ACL(cit)->endParticle()));
	  }
	} else if (ACL(cit)->endParticle() == (*cit)) {
	  if (ACL(cit)->startParticle()) {
	    partners.push_back(dynamic_ptr_cast<ShowerParticlePtr>(ACL(cit)->startParticle()));
	  }
	}
      }
      for(unsigned int i = 0; i < partners.size(); i++) {
	//	cout << "set partner by hand: " << partners[i] << endl;
      }
    }
    candidatePartners[*cit] = partners;
  }
#undef FS
#undef CL
#undef ACL


  // Loop now over the map we have just filled (mapParticleToCandidatePartners)
  // and treat only those particles that have some candidate colour partners,
  // whereas those that don't have any are collected in the vector specialCases
  // and will be treated later.
  tShowerParticleVector specialCases; 
  PartnerMap::const_iterator it;
  for(it = candidatePartners.begin(); it != candidatePartners.end(); ++it) {
    if(it->second.size()) {

      //***LOOKHERE*** In the case of more than one candidate colour partners,
      //               our treatment is based on two assumptions:
      //               1) the choice of which is the colour partner is done
      //                  *randomly* between the available candidates.
      //               2) the choice of which is the colour partner is done
      //                  *independently* from each particle: in other words,
      //                  if for a particle "i" its selected colour partner is  
      //                  the particle "j", then the colour partner of "j" 
      //                  does not have to be necessarily "i".
      int position = UseRandom::irnd(it->second.size());

      pair<Energy,Energy> pairScales = 
	calculateInitialEvolutionScales(ShowerPPair(it->first, 
						    it->second[position]),
			                showerVariables);
      switch(_approach) {
        case 0: // Totally random
         it->first->setEvolutionScale(ShowerIndex::QCD, pairScales.first);
         it->first->setPartner(ShowerIndex::QCD, it->second[position]);
	 break;
	case 1: // Partner is also set, if it has already been set, pick 50/50
	 if(!it->first->partners()[ShowerIndex::QCD] || UseRandom::rndbool()) {
	   it->first->setEvolutionScale(ShowerIndex::QCD, pairScales.first);
	   it->first->setPartner(ShowerIndex::QCD, it->second[position]);
	 }
	 if(!it->second[position]->partners()[ShowerIndex::QCD] ||
	    UseRandom::rndbool()) {
	   it->second[position]->setEvolutionScale(ShowerIndex::QCD, 
			                           pairScales.second);
	   it->second[position]->setPartner(ShowerIndex::QCD, it->first);
	 }
	 break;
      }
    } else specialCases.push_back(it->first);
  } 

  if(specialCases.size()) {

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
    generator()->log() << "PartnerFinder::debuggingInfo "
		       << "end ___________________________________" << endl;
  }
  
  return isOK;
}


bool PartnerFinder::setQEDInitialEvolutionScales(const tShowerVarsPtr showerVariables,
						 const ShowerParticleVector particles,
						 const bool isDecayCase) {
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
  ShowerParticleVector::const_iterator cit;
  for(cit = particles.begin(); cit != particles.end(); ++cit) {
    (*cit)->setEvolutionScale(ShowerIndex::QED, (*cit)->evolutionScales()[ShowerIndex::QCD]); 
    (*cit)->setPartner(ShowerIndex::QED, (*cit)->partners()[ShowerIndex::QCD]); 
  }

  return isOK;
}


bool PartnerFinder::setEWKInitialEvolutionScales(const tShowerVarsPtr showerVariables,
						 const ShowerParticleVector particles,
						 const bool isDecayCase) {

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
calculateInitialEvolutionScales(const ShowerPPair &particlePair, 
		                const tShowerVarsPtr showerVariables) {
  if(particlePair.first->isFinalState() && particlePair.second->isFinalState())
    return calculateFinalFinalScales(particlePair, showerVariables);
  else if(particlePair.first->isFinalState() && 
          !particlePair.second->isFinalState()) {
    ShowerPPair a(particlePair.second, particlePair.first);
    pair<Energy,Energy> rval = calculateInitialFinalScales(a,showerVariables);
    return pair<Energy,Energy>(rval.second,rval.first);
  }
  else if(!particlePair.first->isFinalState() && 
	  particlePair.second->isFinalState())
    return calculateInitialFinalScales(particlePair,showerVariables);
  else
    return calculateInitialInitialScales(particlePair,showerVariables);
}

pair<Energy,Energy> PartnerFinder::
calculateInitialFinalScales(const ShowerPPair &ppair, const tShowerVarsPtr s) {
  /********
   * In this case from JHEP 12(2003)045 we find the conditiongs
   * ktilda_b = (1+c) and ktilda_c = (1+2c)
   * We also find that c = m_c^2/Q^2. The process is a+b->c where
   * particle a is not colour connected (considered as a colour singlet).
   * Therefore we simply find that q_b = sqrt(Q^2+m_c^2) and 
   * q_c = sqrt(Q^2+2 m_c^2)
   * We also assume that the first particle in the pair is the initial
   * state particle and the second is the final state one (c)
   *********/
  Lorentz5Momentum p1, p2, p;
  p1 = ppair.first->momentum();
  p2 = ppair.second->momentum();
  p = p1+p2;
  p.boost(p.findBoostToCM());
  Energy2 mc2 = ppair.second->mass();
  Energy Q2 = p*p;
  return pair<Energy,Energy>(sqrt(Q2+mc2), sqrt(Q2+2*mc2));
}

pair<Energy,Energy> PartnerFinder::
calculateInitialInitialScales(const ShowerPPair &ppair, const tShowerVarsPtr s)
{
  /*******
   * This case is quite simple. From JHEP 12(2003)045 we find the condition
   * that ktilda_b = ktilda_c = 1. In this case we have the process
   * b+c->a so we need merely boost to the CM frame of the two incoming
   * particles and then qtilda is equal to the energy in that frame
   **********/
  Lorentz5Momentum p1, p2, p;
  p1 = ppair.first->momentum();
  p2 = ppair.second->momentum();
  p = p1+p2;
  p.boost((p1+p2).findBoostToCM());
  Energy Q = sqrt(p*p);
  return pair<Energy,Energy>(Q,Q);
}

pair<Energy,Energy> PartnerFinder::
calculateFinalFinalScales(const ShowerPPair &particlePair,
		          const tShowerVarsPtr showerVariables) {
  Energy firstQ = Energy();
  Energy secondQ = Energy();

  /********
   * Using JHEP 12(2003)045 we find that we need ktilda = 1/2(1+b-c+lambda)
   * ktilda = qtilda^2/Q^2 therefore qtilda = sqrt(ktilda*Q^2)
   * This is found from
   *   x -> a abar
   *   Px = 1/2 Q(1,0,0), 
   *   Pa = 1/2 Q(1+b-c,0,lambda), 
   *   Pabar = 1/2 Q(1-b+c,0,-lambda)
   *   so nabar = 1/2 Q(lambda,0,-lambda)
   * and we find (Pa+Pabar)(Pa+nabar) = 1/2 Q^2(1+b-c+lambda) which
   * is exactly the condition we want for qtilda^2.
   * We also find that this also applies for the ktilda for the abar
   * particle, where we flip b and c.
   *************/

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
  if (showerVariables->asyPS()) {
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
  if (showerVariables->rndPS()) {  
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



