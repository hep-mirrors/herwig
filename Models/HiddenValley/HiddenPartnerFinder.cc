// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiddenPartnerFinder class.
//

#include "HiddenPartnerFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"
#include "ThePEG/Repository/UseRandom.h" 
#include "AdditionalGaugeParticleData.h"

using namespace Herwig;

IBPtr HiddenPartnerFinder::clone() const {
  return new_ptr(*this);
}

IBPtr HiddenPartnerFinder::fullclone() const {
  return new_ptr(*this);
}

void HiddenPartnerFinder::persistentOutput(PersistentOStream & os) const {
}

void HiddenPartnerFinder::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<HiddenPartnerFinder> HiddenPartnerFinder::initHiddenPartnerFinder;
// Definition of the static class description member.

void HiddenPartnerFinder::Init() {

  static ClassDocumentation<HiddenPartnerFinder> documentation
    ("There is no documentation for the HiddenPartnerFinder class");

}


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

inline bool hiddenColoured(tcPDPtr part) {
  tcAdditionalGaugeParticleDataPtr hp = 
    dynamic_ptr_cast<tcAdditionalGaugeParticleDataPtr>(part);
  if(!hp) return false;
  return hp->hiddenColoured();
}
}

bool HiddenPartnerFinder::setInitialEvolutionScales(const ShowerParticleVector &particles,
					      const bool isDecayCase,
					      const bool setPartners) {
  // set the partners and the scales
  if(setPartners) {
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
    ShowerParticleVector::const_iterator cit, cjt;
    for(cit = particles.begin(); cit != particles.end(); ++cit) {
      if(!hiddenColoured((*cit)->dataPtr())) continue;
      // We now have a coloured particle
      tShowerParticleVector partners;
      for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
	if(!hiddenColoured((*cjt)->dataPtr())||cit==cjt) continue;
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
	    if(( FS(*cjt) && ( CL(cjt) == cpair.first || CL(cjt)  == cpair.second))||
	       (!FS(*cjt) && (ACL(cjt) == cpair.first || ACL(cjt) == cpair.second ))) {
	      partners.push_back(*cjt);
	    }
	  }
	}
	else if(col&&col->sinkNeighbours().first) {
	  tColinePair cpair = col->sinkNeighbours();
	  for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	    if(( FS(*cjt) && (ACL(cjt) == cpair.first || ACL(cjt)  == cpair.second))||
	       (!FS(*cjt) && ( CL(cjt) == cpair.first ||  CL(cjt) == cpair.second))) {
	      partners.push_back(*cjt);
	    }
	  }
	}
	col = ACL(cit);
	if(FS(*cit)&&col&&col->sinkNeighbours().first) {
	  tColinePair cpair = col->sinkNeighbours();
	  for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	    if(( FS(*cjt) && (ACL(cjt) == cpair.first || ACL(cjt)  == cpair.second))||
	       (!FS(*cjt) && ( CL(cjt) == cpair.first ||  CL(cjt) == cpair.second ))) {
	      partners.push_back(*cjt);
	    }
	  }
	}
	else if(col&&col->sourceNeighbours().first) {
	  tColinePair cpair = col->sourceNeighbours();
	  for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	    if(( FS(*cjt) && ( CL(cjt) == cpair.first || CL(cjt) == cpair.second))||
	       (!FS(*cjt) && (ACL(cjt) == cpair.first ||ACL(cjt) == cpair.second))) {
	      partners.push_back(*cjt);
	    }
	  }
	}
	if(partners.empty()) {
	  throw Exception() << "`Failed to make colour connections in " 
			    << "HiddenPartnerFinder::setQCDInitialEvolutionScales"
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
      switch(approach()) {
      case 0: // Totally random
	(*cit)->setEvolutionScale(pairScales.first);
	(*cit)->setPartner(partners[position]);
	break;
      case 1: // Partner is also set, if it has already been set, pick 50/50
	if(!(*cit)->partner() || UseRandom::rndbool()) {
	  (*cit)->setEvolutionScale(pairScales.first);
	  (*cit)->setPartner(partners[position]);
	}
	if(!partners[position]->partner() || UseRandom::rndbool()) {
	  partners[position]->setEvolutionScale(pairScales.second);
	  partners[position]->setPartner(*cit);
	}
	break;
      default:
	exit(2);
	throw Exception() << "Invalid approach for setting colour partner in"
			  << " HiddenPartnerFinder::setQCDInitialEvolutionScale()"
			  << Exception::abortnow;
      }
    }
  }
  // partners all ready set only do the scales
  else {
    for(ShowerParticleVector::const_iterator cit = particles.begin();
	cit != particles.end(); ++cit) {
      if(!hiddenColoured((**cit).dataPtr())) continue;
      tShowerParticlePtr partner = (**cit).partner();
      pair<Energy,Energy> pairScales = 
	calculateInitialEvolutionScales(ShowerPPair(*cit,partner),
					isDecayCase);
      (*cit)->setEvolutionScale(pairScales.first);
    }
  }
  return true;
}
