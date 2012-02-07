// -*- C++ -*-
//
// PartnerFinder.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
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
#include "ThePEG/Repository/UseRandom.h" 
#include "ThePEG/Interface/Switch.h"

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
  os << _approach << _partnerMethod;
}

void PartnerFinder::persistentInput(PersistentIStream & is, int) {
  is >> _approach >> _partnerMethod;
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

  static Switch<PartnerFinder,int> interfacePartnerMethod
    ("PartnerMethod",
     "Choice of partner finding method for gluon evolution.",
     &PartnerFinder::_partnerMethod, 0, false, false);
  static SwitchOption interfacePartnerMethodRandom
    (interfacePartnerMethod,
     "Random",
     "Choose partners of a gluon randomly.",
     0);
  static SwitchOption interfacePartnerMethodMaximum
    (interfacePartnerMethod,
     "Maximum",
     "Choose partner of gluon with largest angle.",
     1);

  static Switch<PartnerFinder,int> interfaceQEDPartner
    ("QEDPartner",
     "Control of which particles to use as the partner for QED radiation",
     &PartnerFinder::QEDPartner_, 0, false, false);
  static SwitchOption interfaceQEDPartnerAll
    (interfaceQEDPartner,
     "All",
     "Consider all possible choices which give a positive contribution"
     " in the soft limit.",
     0);
  static SwitchOption interfaceQEDPartnerIIandFF
    (interfaceQEDPartner,
     "IIandFF",
     "Only allow initial-initial or final-final combinations",
     1);
  static SwitchOption interfaceQEDPartnerIF
    (interfaceQEDPartner,
     "IF",
     "Only allow initial-final combinations",
     2);

}

bool PartnerFinder::setInitialEvolutionScales(const ShowerParticleVector &particles,
					      const bool isDecayCase,
					      ShowerInteraction::Type type,
					      const bool setPartners) {
  if(type==ShowerInteraction::QCD) {
    return setInitialQCDEvolutionScales(particles,isDecayCase,setPartners);
  }
  else if(type==ShowerInteraction::QED) {
    return setInitialQEDEvolutionScales(particles,isDecayCase,setPartners);
  }
  else {
    throw Exception() << "Must be either QCD or QED in "
		      << "PartnerFinder::setInitialEvolutionScales()\n"
		      << Exception::runerror;
  }
}

bool PartnerFinder::setInitialQCDEvolutionScales(const ShowerParticleVector &particles,
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
    // Define variables needed for angular and radiation line analysis
    double angle[2] = { 0.0, 0.0 }; 	
    int radiationLine[2] = { 0, 0 };    
    int ait = 0;    
    for(cit = particles.begin(); cit != particles.end(); ++cit) {
      if(!(*cit)->data().coloured()) continue;
      // We now have a coloured particle
      tShowerParticleVector partners;
      if (_partnerMethod == 1){
 	// Angular analysis need to be reset for each particle     
	angle[0] = 0;
	angle[1] = 0;	
	radiationLine[0] = 0;
	radiationLine[1] = 0;
	ait = 0;	
      }
      for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
	if(!(*cjt)->data().coloured()||cit==cjt) continue;
	bool isPartner = false;
	if(FS(*cit) != FS(*cjt) &&
	   ((CL(cit) && CL(cit)==CL(cjt)) || (ACL(cit) && ACL(cit)==ACL(cjt))))
	  isPartner = true;
	else if((CL(cit) && CL(cit)==ACL(cjt)) || (ACL(cit) && ACL(cit)==CL(cjt)))
	  isPartner = true;
	if(isPartner){
	  if (_partnerMethod == 1){
	    // Find the angle to the potential partner
	    angle[ait] = (*cit)->momentum().vect().angle((*cjt)->momentum().vect());	 
	    // Find out which colour line is connected to the partner
	    // Colour line has value 1, anti-colour line has value 2
	    if( CL(cit) == CL(cjt) ||  CL(cit) == ACL(cjt)) radiationLine[ait] = 1;
	    if(ACL(cit) == CL(cjt) || ACL(cit) == ACL(cjt)) radiationLine[ait] = 2;	
	    ait++;
	  }
	partners.push_back(*cjt);
	}
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
      //	       there are now two approaches to choosing the partner. The
      //	       first method is based on two assumptions:
      //               1) the choice of which is the colour partner is done
      //                  *randomly* between the available candidates.
      //               2) the choice of which is the colour partner is done
      //                  *independently* from each particle: in other words,
      //                  if for a particle "i" its selected colour partner is  
      //                  the particle "j", then the colour partner of "j" 
      //                  does not have to be necessarily "i".
      // 	      The second method always chooses the furthest partner
      //	      for hard gluons and gluinos.
      
      // First make a random choice of partner
      int position = UseRandom::irnd(partners.size());   
      // Set the scale from the random partner
      pair<Energy,Energy> pairScales = 
	calculateInitialEvolutionScales(ShowerPPair(*cit,partners[position]),
					isDecayCase);
      pair<Energy,Energy> pairScales2;
      if (_partnerMethod == 1){
        // Override choice of partner  
      	if ((*cit)->perturbative() == 1 && getParticleData((*cit)->id())->iColour()==PDT::Colour8){     
          // Parton is a hard octet
	  // Determine largest angle
      	  if (angle[1]>0 && angle[0]<angle[1]) position=1;
          else position=0;
	}      
        // Override scale for hard octet
        if (position==1){
	  pairScales = calculateInitialEvolutionScales(ShowerPPair(*cit,partners[1]),
					isDecayCase);
	  if (getParticleData((*cit)->id())->iColour()==PDT::Colour8){
	    // Set the second scale for hard octets			
	    pairScales2 = calculateInitialEvolutionScales(ShowerPPair(*cit,partners[0]),
						isDecayCase);  
	  }
	}
   	else {
	  pairScales = calculateInitialEvolutionScales(ShowerPPair(*cit,partners[0]),
					isDecayCase);
	  if (getParticleData((*cit)->id())->iColour()==PDT::Colour8){	
	    // Set the second scale for hard octets
	    pairScales2 = calculateInitialEvolutionScales(ShowerPPair(*cit,partners[1]),
						isDecayCase);	
	  }				
        }
        if ((*cit)->perturbative() == 1 && getParticleData((*cit)->id())->iColour()==PDT::Colour8){
        // Set radiation lines for hard octets
	  (*cit)->setRadiationLine(radiationLine[position]);
	  if( !(*cit)->progenitor() ){
	    // Set the hard partons to be the progenitors of the shower
	    (*cit)->setProgenitor(*cit);
	    // Set the second evolution scale of the progenitor
	    (*cit)->setEvolutionScale2(pairScales2.first);
	 }		
       }
       else if ((*cit)->perturbative() == 1 && (getParticleData((*cit)->id())->iColour()==PDT::Colour3 ||
       						getParticleData((*cit)->id())->iColour()==PDT::Colour3bar)){
       // Set radiation lines for hard triplets
	  if( !(*cit)->progenitor() ){
	    // Set the hard partons to be the progenitors of the shower
	    (*cit)->setProgenitor(*cit);
	    // Set the second evolution scale of the progenitor
	    (*cit)->setEvolutionScale2(pairScales2.first);
	    // Set the radiation line
	    (*cit)->setRadiationLine(0);
	 }		         
       }
     }
      switch(_approach) {
      case 0: // Totally random (unless chosen above)
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

bool PartnerFinder::setInitialQEDEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase,
						 const bool setPartners) {
  // set the partners and the scales
  if(setPartners) {
    ShowerParticleVector::const_iterator cit, cjt;
    for(cit = particles.begin(); cit != particles.end(); ++cit) {
      if(!(*cit)->data().charged()) continue;
      // We now have a charged particle
      vector<pair<double,tShowerParticlePtr> > partners;
      for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
	if(!(*cjt)->data().charged()||cit==cjt) continue;
	double charge = double((*cit)->data().iCharge()*(*cjt)->data().iCharge());
 	if( FS(*cit) != FS(*cjt) ) charge *=-1.;
	if( QEDPartner_ != 0 ) {
	  // only include II and FF is requested
	  if( QEDPartner_ == 1 && FS(*cit) != FS(*cjt) )
	    continue;
	  // ony include IF is requested
	  else if(QEDPartner_ == 2 && FS(*cit) == FS(*cjt) )
	    continue;
	}
	// only keep positive dipoles
	if(charge<0.) partners.push_back(make_pair(-charge,*cjt));
      }
      if(partners.empty()) {
	throw Exception() << "Failed to partner in " 
			  << "PartnerFinder::setQEDInitialEvolutionScales"
			  << (**cit) << Exception::eventerror;
      }
      double prob(0.);
      for(unsigned int ix=0;ix<partners.size();++ix) prob +=partners[ix].first;
      prob *= UseRandom::rnd();
      tShowerParticlePtr partner;
      for(unsigned int ix=0;ix<partners.size();++ix) {
	if(partners[ix].first>prob) {
	  partner = partners[ix].second;
	  break;
	}
	prob -= partners[ix].first;
      }
      if(!partner) 
	throw Exception() << "Failed to partner in " 
			  << "PartnerFinder::setQEDInitialEvolutionScales"
			  << (**cit) << Exception::eventerror;
      pair<Energy,Energy> pairScales = 
	calculateInitialEvolutionScales(ShowerPPair(*cit,partner),
					isDecayCase);
      (*cit)->setEvolutionScale(pairScales.first);
      (*cit)->setPartner(partner);
    }
  }
  // partners all ready set only do the scales
  else {
    for(ShowerParticleVector::const_iterator cit = particles.begin();
	cit != particles.end(); ++cit) {
      if(!(**cit).dataPtr()->charged()) continue;
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
