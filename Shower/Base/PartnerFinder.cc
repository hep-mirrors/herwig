// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartnerFinder class.
//

#include "PartnerFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ShowerParticle.h"
#include "ThePEG/Utilities/Timer.h"
#include "ThePEG/Repository/UseRandom.h" 

using namespace Herwig;

// some useful functions to avoid using #define
namespace {

// return bool if final-state particle
inline bool FS(const tShowerParticlePtr a) {
  return a->isFinalState();
}

// return colour line pointer 
inline Ptr<ThePEG::ColourLine>::transient_pointer  
CL(const ShowerParticleVector::const_iterator & a) {
  return (*a)->colourLine();
}
  
inline Ptr<ThePEG::ColourLine>::transient_pointer
ACL(const ShowerParticleVector::const_iterator & a) {
  return (*a)->antiColourLine();
}

}

void PartnerFinder::persistentOutput(PersistentOStream & os) const {
  os << _approach;
}

void PartnerFinder::persistentInput(PersistentIStream & is, int) {
  is >> _approach;
}

AbstractClassDescription<PartnerFinder> PartnerFinder::initPartnerFinder;
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


bool PartnerFinder::setQCDInitialEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase) {
  Timer<1300> timer("PartnerFinder::setQCDInitialEvolutionScales");
  //  bool isOK = true;

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
  // (this means that no modifications of the following loop is needed!)
  ShowerParticleVector::const_iterator cit, cjt;
  for(cit = particles.begin(); cit != particles.end(); ++cit) {
    if(!(*cit)->data().coloured()) continue;
    // We now have a coloured particle
    tShowerParticleVector partners;
    for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
      if(!(*cjt)->data().coloured()||cit==cjt) continue;
      bool isPartner = false;
      if(FS(*cit) != FS(*cjt) &&
         ((CL(cit) && CL(cit)==CL(cjt)) || (ACL(cit) && ACL(cit)==ACL(cjt))))
	 isPartner = true;
      else if((CL(cit) && CL(cit)==ACL(cjt)) || (ACL(cit) && ACL(cit)==CL(cjt)))
	isPartner = true;
      if(isPartner) partners.push_back(*cjt);
    }
    if (partners.empty()) {
      // special for RPV
      tColinePtr col = CL(cit); 
      if(FS(*cit)&&col&&col->sourceNeighbours().first) {
	tColinePair cpair = col->sourceNeighbours();
	for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	  if((FS(*cjt) && ( CL(cjt) == cpair.first || CL(cjt)  == cpair.second))||
	     (             ACL(cjt) == cpair.first || ACL(cjt) == cpair.second )) {
	    partners.push_back(*cjt);
	  }
	}
      }
      else if(col&&col->sinkNeighbours().first) {
	throw Exception() << "PartnerFinder::setQCDInitialEvolutionScale() does not yet"
			  << "support a colour line connected to a sink"
			  << Exception::runerror;
      }
      col = ACL(cit);
      if(FS(*cit)&&col&&col->sinkNeighbours().first) {
	tColinePair cpair = col->sinkNeighbours();
	for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	  if((FS(*cjt) && ( CL(cjt) == cpair.first || CL(cjt)  == cpair.second))||
	     (             ACL(cjt) == cpair.first || ACL(cjt) == cpair.second )) {
	    partners.push_back(*cjt);
	  }
	}
      }
      else if(col&&col->sourceNeighbours().first) {
	throw Exception() << "PartnerFinder::setQCDInitialEvolutionScale() does not yet"
			  << "support an anti-colour line connected to a source"
			  << Exception::runerror;
      }
      if(partners.empty()) {
	throw Exception() << "`Failed to make colour connections in " 
			  << "PartnerFinder::setQCDInitialEvolutionScales"
			  << (**cit)
			  << Exception::eventerror;
      }
    }
    // In the case of more than one candidate colour partners,
    //               our treatment is based on two assumptions:
    //               1) the choice of which is the colour partner is done
    //                  *randomly* between the available candidates.
    //               2) the choice of which is the colour partner is done
    //                  *independently* from each particle: in other words,
    //                  if for a particle "i" its selected colour partner is  
    //                  the particle "j", then the colour partner of "j" 
    //                  does not have to be necessarily "i".
    int position = UseRandom::irnd(partners.size());
    pair<Energy,Energy> pairScales = 
      calculateInitialEvolutionScales(ShowerPPair(*cit,partners[position]),
                                      isDecayCase);
    switch(_approach) {
    case 0: // Totally random
      (*cit)->setEvolutionScale(ShowerIndex::QCD, pairScales.first);
      (*cit)->setPartner(ShowerIndex::QCD, partners[position]);
      break;
    case 1: // Partner is also set, if it has already been set, pick 50/50
      if(!(*cit)->partners()[ShowerIndex::QCD] || UseRandom::rndbool()) {
	(*cit)->setEvolutionScale(ShowerIndex::QCD, pairScales.first);
	(*cit)->setPartner(ShowerIndex::QCD, partners[position]);
      }
      if(!partners[position]->partners()[ShowerIndex::QCD] ||
	 UseRandom::rndbool()) {
	partners[position]->setEvolutionScale(ShowerIndex::QCD, 
						pairScales.second);
	partners[position]->setPartner(ShowerIndex::QCD, *cit);
      }
      break;
    default:
      throw Exception() << "Invalid approach for setting colour partner in"
			<< " PartnerFinder::setQCDInitialEvolutionScale()"
			<< Exception::abortnow;
    }
  } 
  return true;
}

bool PartnerFinder::setQEDInitialEvolutionScales(const ShowerParticleVector &,
						 const bool) {

  // ***LOOKHERE*** To be implemented only if you want to have electromagnetic
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
  throw Exception() << "PartnerFinder::setQEDInitialEvolutionScales "
		    << "implementation is not correct.\nMust match charge "
		    << "partners, not colour partners.\n"
		    << "Turn off QED in Shower.in" 
		    << Exception::runerror;


  return false;
}

bool PartnerFinder::setEWKInitialEvolutionScales(const ShowerParticleVector &,
						 const bool ) {
  // ***LOOKHERE*** To be implemented only if you want to have electroweak
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
  throw Exception() << "PartnerFinder::setEWKInitialEvolutionScales not "
		    << "implemented.\nTurn off EWK in Shower.in" 
		    << Exception::runerror;
  return false;
}

pair<Energy,Energy> PartnerFinder::
calculateInitialEvolutionScales(const ShowerPPair &particlePair, const bool isDecayCase) {
  bool FS1=FS(particlePair.first),FS2= FS(particlePair.second);

  if(FS1 && FS2)
    return calculateFinalFinalScales(particlePair);
  else if(FS1 && !FS2) {
    ShowerPPair a(particlePair.second, particlePair.first);
    pair<Energy,Energy> rval = calculateInitialFinalScales(a,isDecayCase);
    return pair<Energy,Energy>(rval.second,rval.first);
  }
  else if(!FS1 &&FS2)
    return calculateInitialFinalScales(particlePair,isDecayCase);
  else
    return calculateInitialInitialScales(particlePair);
}
