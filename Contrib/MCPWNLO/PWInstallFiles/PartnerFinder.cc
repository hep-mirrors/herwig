// -*- C++ -*-
//
// PartnerFinder.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartnerFinder class.
//

#include "PartnerFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ShowerParticle.h"
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
  os << _approach << _PHpfinder;
}

void PartnerFinder::persistentInput(PersistentIStream & is, int) {
  is >> _approach >> _PHpfinder;
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

  static Switch<PartnerFinder, unsigned int> ifacePHpfinder
    ("PHPartnerFinder",
     "Choice of the Partner Finder",
     &PartnerFinder::_PHpfinder, 0, false, false);
   static SwitchOption PHpfinderoff
    (ifacePHpfinder,"No","PHpfinder off", 0);
   static SwitchOption PHpfinderon
    (ifacePHpfinder,"Yes","PHpfinder on", 1);
 

}

bool PartnerFinder::setInitialEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase,
						 const bool setPartners) {
  int particlenumber = 0; 

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


  /**************************************************************************/
      // int position = UseRandom::irnd(partners.size());
//       pair<Energy,Energy> pairScales = 
// 	calculateInitialEvolutionScales(ShowerPPair(*cit,partners[position]),
// 					isDecayCase);
      // pair<Energy,Energy> pairScales;
      
    pair<Energy,Energy> testscale,pairScales(make_pair(1e10*GeV,1e10*GeV));
    int position(0);
     int gluoncount = 0;
     int particlecount = 0;
     // count the number of gluons in shower tree.
     bool hardprocess = false;
     ShowerParticleVector::const_iterator ckt;
     for(ckt=particles.begin();ckt!=particles.end();++ckt) {
      if (fabs((**ckt).id()) > 5 && (**ckt).id() != 21) {hardprocess = true;} // dealing with hard process for W, Z or top pair production
            if ((**ckt).id() == 21){
      	++gluoncount;
       }
	    if(fabs((**ckt).id()) < 7 || (**ckt).id() == 21){
	      ++particlecount;}
}
       ++particlenumber; //update particle number for showered particles

       cout << "particleinfo" << particlenumber << "\t" << particlecount << endl; 

        if (PHpfinder() && hardprocess) {
	  
                  for(unsigned int ix=0;ix<partners.size();++ix) {
       testscale =
     calculateInitialEvolutionScales(ShowerPPair(*cit,partners[ix]),
                                                   isDecayCase);
                         if(testscale.first<pairScales.first) {
         pairScales=testscale;
         position=ix;
       
     }
     } 
       pair<Energy,Energy> testscale2(make_pair(0*GeV,0*GeV));

       //set scale for truncated gluon to larger of the 2 scales.

	   if (gluoncount > 1 && partners.size() == 2) {
		if(particlenumber == particlecount - 1) {
      	for(unsigned int ix=0;ix<partners.size();++ix) {
		   testscale2 =
	     calculateInitialEvolutionScales(ShowerPPair(*cit,partners[ix]),
     			     isDecayCase);
		   if(testscale2.first > pairScales.first) {
	    pairScales=testscale2;
	  position=ix;
       }
	}}
	   }}
           
  else {

    position = UseRandom::irnd(partners.size());
    pairScales = 
    calculateInitialEvolutionScales(ShowerPPair(*cit,partners[position]),
                                        isDecayCase);
     }
        /**************************************************************************/
      switch(_approach) {
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
			  << " PartnerFinder::setQCDInitialEvolutionScale()"
			  << Exception::abortnow;
      }
    }
  }
  // partners all ready set only do the scales
  else {
    for(ShowerParticleVector::const_iterator cit = particles.begin();
	cit != particles.end(); ++cit) {
      if(!(**cit).dataPtr()->coloured()) continue;
      tShowerParticlePtr partner = (**cit).partner();
      pair<Energy,Energy> pairScales = 
	calculateInitialEvolutionScales(ShowerPPair(*cit,partner),
					isDecayCase);
      (*cit)->setEvolutionScale(pairScales.first);
    }
  }
  return true;
}

pair<Energy,Energy> PartnerFinder::
calculateInitialEvolutionScales(const ShowerPPair &particlePair,
				const bool isDecayCase) {
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
