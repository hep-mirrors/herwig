// -*- C++ -*-
//
// PartnerFinder.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/Shower/Core/Base/ShowerParticle.h"
#include "ThePEG/Repository/UseRandom.h" 
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeAbstractClass<PartnerFinder,Interfaced>
describePartnerFinder ("Herwig::PartnerFinder","HwShower.so");

// some useful functions to avoid using #define
namespace {

  // return bool if final-state particle
  inline bool FS(const tShowerParticlePtr a) {
    return a->isFinalState();
  }

  // return colour line pointer 
  inline Ptr<ThePEG::ColourLine>::transient_pointer  
  CL(const tShowerParticlePtr a, unsigned int index=0) {
    return a->colourInfo()->colourLines().empty() ? ThePEG::tColinePtr() :
      const_ptr_cast<ThePEG::tColinePtr>(a->colourInfo()->colourLines()[index]);
  }

  // return colour line size
  inline unsigned int
  CLSIZE(const tShowerParticlePtr a) {
    return a->colourInfo()->colourLines().size();
  }
  
  inline Ptr<ThePEG::ColourLine>::transient_pointer
  ACL(const tShowerParticlePtr a, unsigned int index=0) {
    return a->colourInfo()->antiColourLines().empty() ?  ThePEG::tColinePtr() :
      const_ptr_cast<ThePEG::tColinePtr>(a->colourInfo()->antiColourLines()[index]);
  }

  inline unsigned int
  ACLSIZE(const tShowerParticlePtr a) {
    return a->colourInfo()->antiColourLines().size();
  }
}

void PartnerFinder::persistentOutput(PersistentOStream & os) const {
  os << partnerMethod_ << QEDPartner_ << scaleChoice_;
}

void PartnerFinder::persistentInput(PersistentIStream & is, int) {
  is >> partnerMethod_ >> QEDPartner_ >> scaleChoice_;
}

void PartnerFinder::Init() {

  static ClassDocumentation<PartnerFinder> documentation
    ("This class is responsible for finding the partners for each interaction types ",
     "and within the evolution scale range specified by the ShowerVariables ",
     "then to determine the initial evolution scales for each pair of partners.");

  static Switch<PartnerFinder,int> interfacePartnerMethod
    ("PartnerMethod",
     "Choice of partner finding method for gluon evolution.",
     &PartnerFinder::partnerMethod_, 0, false, false);
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

  static Switch<PartnerFinder,int> interfaceScaleChoice
    ("ScaleChoice",
     "The choice of the evolution scales",
     &PartnerFinder::scaleChoice_, 0, false, false);
  static SwitchOption interfaceScaleChoicePartner
    (interfaceScaleChoice,
     "Partner",
     "Scale of all interactions is that of the evolution partner",
     0);
  static SwitchOption interfaceScaleChoiceDifferent
    (interfaceScaleChoice,
     "Different",
     "Allow each interaction to have different scales",
     1);

}

void PartnerFinder::setInitialEvolutionScales(const ShowerParticleVector &particles,
					      const bool isDecayCase,
					      ShowerInteraction type,
					      const bool setPartners) {
  // clear the existing partners
  for(ShowerParticleVector::const_iterator cit = particles.begin();
      cit != particles.end(); ++cit) (*cit)->clearPartners();
  // set them
  if(type==ShowerInteraction::QCD) {
    setInitialQCDEvolutionScales(particles,isDecayCase,setPartners);
  }
  else if(type==ShowerInteraction::QED) {
    setInitialQEDEvolutionScales(particles,isDecayCase,setPartners);
  }
  else {
    setInitialQCDEvolutionScales(particles,isDecayCase,setPartners);
    setInitialQEDEvolutionScales(particles,isDecayCase,false);
  }
  // print out for debugging
  if(Debug::level>=10) {
    for(ShowerParticleVector::const_iterator cit = particles.begin();
	cit != particles.end(); ++cit) {
      generator()->log() << "Particle: " << **cit << "\n";
      if(!(**cit).partner()) continue;
      generator()->log() << "Primary partner: " << *(**cit).partner() << "\n";
      for(vector<ShowerParticle::EvolutionPartner>::const_iterator it= (**cit).partners().begin();
	  it!=(**cit).partners().end();++it) {
	generator()->log() << static_cast<long>(it->type) << " "
			   << it->weight << " " 
			   << it->scale/GeV << " " 
 			   << *(it->partner) 
			   << "\n";
      }
    }
    generator()->log() << flush;
  }
}

void PartnerFinder::setInitialQCDEvolutionScales(const ShowerParticleVector &particles,
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
    //  --- if both are initial (incoming) state particles, then the (non-null) colourLine() 
    //      of one of them must match the (non-null) antiColourLine() of the other.
    //  --- if one is an initial (incoming) state particle and the other is
    //      a final (outgoing) state particle,  then both must have the
    //      same (non-null) colourLine() or the same (non-null) antiColourLine();
    // Notice that this definition exclude the special case of baryon-violating
    // processes (as in R-parity Susy), which will show up as particles
    // without candidate colour partners, and that we will be treated a part later
    // (this means that no modifications of the following loop is needed!)
    ShowerParticleVector::const_iterator cit, cjt;
    for(cit = particles.begin(); cit != particles.end(); ++cit) {
      // Skip colourless particles
      if(!(*cit)->data().coloured()) continue;
      // find the partners
      vector< pair<ShowerPartnerType, tShowerParticlePtr> > partners = 
	findQCDPartners(*cit,particles);
      // must have a partner
      if(partners.empty()) {
	throw Exception() << "`Failed to make colour connections in " 
                          << "PartnerFinder::setQCDInitialEvolutionScales"
                          << (**cit)
                          << Exception::eventerror;
      }
      // Calculate the evolution scales for all possible pairs of of particles
      vector<pair<Energy,Energy> > scales;
      for(unsigned int ix=0;ix< partners.size();++ix) {
	scales.push_back(calculateInitialEvolutionScales(ShowerPPair(*cit,partners[ix].second),
							 isDecayCase));
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
      // store the choice
      int position(-1);
      // random choice
      if( partnerMethod_ == 0 ) {
	// random choice of partner
	position = UseRandom::irnd(partners.size());
      }
      // take the one with largest angle
      else if (partnerMethod_ == 1 ) {
	if ((*cit)->perturbative() == 1 && 
	    (*cit)->dataPtr()->iColour()==PDT::Colour8 ) {
	  assert(partners.size()==2);
	  // Determine largest angle
	  double maxAngle(0.);
	  for(unsigned int ix=0;ix<partners.size();++ix) {
	    double angle = (*cit)->momentum().vect().
	      angle(partners[ix].second->momentum().vect());
	    if(angle>maxAngle) {
	      maxAngle = angle;
	      position = ix;
	    }
	  }
	}
	else
	  position = UseRandom::irnd(partners.size());
      }
      else
	assert(false);
      // set the evolution partner
      (*cit)->partner(partners[position].second);
      for(unsigned int ix=0;ix<partners.size();++ix) {
	(**cit).addPartner(ShowerParticle::EvolutionPartner(partners[ix].second,
							    1.,partners[ix].first,
							    scales[ix].first));
      }
      // set scales for all interactions to that of the partner, default
      Energy scale = scales[position].first;
      for(unsigned int ix=0;ix<partners.size();++ix) {
	if(partners[ix].first==ShowerPartnerType::QCDColourLine) {
	  (**cit).scales().QCD_c = 
	    (**cit).scales().QCD_c_noAO = 
	    (scaleChoice_==0 ? scale : scales[ix].first);
	}
	else if(partners[ix].first==ShowerPartnerType::QCDAntiColourLine) {
	  (**cit).scales().QCD_ac = 
	    (**cit).scales().QCD_ac_noAO =
	    (scaleChoice_==0 ? scale : scales[ix].first);
	}
	else
	  assert(false);
      }
    }
  }
  // primary partner set, set the others and do the scale
  else {
    for(ShowerParticleVector::const_iterator cit = particles.begin();
        cit != particles.end(); ++cit) {
      // Skip colourless particles
      if(!(*cit)->data().coloured()) continue;
      // find the partners
      vector< pair<ShowerPartnerType, tShowerParticlePtr> > partners = 
	findQCDPartners(*cit,particles);
      // must have a partner
      if(partners.empty()) {
        throw Exception() << "`Failed to make colour connections in " 
                          << "PartnerFinder::setQCDInitialEvolutionScales"
                          << (**cit)
                          << Exception::eventerror;
      }
      // Calculate the evolution scales for all possible pairs of of particles
      vector<pair<Energy,Energy> > scales;
      int position(-1);
      for(unsigned int ix=0;ix< partners.size();++ix) {
	if(partners[ix].second) position = ix;
	scales.push_back(calculateInitialEvolutionScales(ShowerPPair(*cit,partners[ix].second),
							 isDecayCase));
      }
      assert(position>=0);
      for(unsigned int ix=0;ix<partners.size();++ix) {
	(**cit).addPartner(ShowerParticle::EvolutionPartner(partners[ix].second,
							    1.,partners[ix].first,
							    scales[ix].first));
      }
      // set scales for all interactions to that of the partner, default
      Energy scale = scales[position].first;
      for(unsigned int ix=0;ix<partners.size();++ix) {
	if(partners[ix].first==ShowerPartnerType::QCDColourLine) {
	  (**cit).scales().QCD_c = 
	    (**cit).scales().QCD_c_noAO = 
	    (scaleChoice_==0 ? scale : scales[ix].first);
	}
	else if(partners[ix].first==ShowerPartnerType::QCDAntiColourLine) {
	  (**cit).scales().QCD_ac = 
	    (**cit).scales().QCD_ac_noAO =
	    (scaleChoice_==0 ? scale : scales[ix].first);
	}
	else {
	  assert(false);
	}
      }
    }
  }
}

void PartnerFinder::setInitialQEDEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase,
						 const bool setPartners) {
  // loop over all the particles
  for(ShowerParticleVector::const_iterator cit = particles.begin();
      cit != particles.end(); ++cit) {
    // not charged or photon continue
    if(!(**cit).dataPtr()->charged()) continue;
    // find the potential partners
    vector<pair<double,tShowerParticlePtr> > partners = findQEDPartners(*cit,particles,isDecayCase);
    if(partners.empty()) {
      throw Exception() << "Failed to find partner in " 
			<< "PartnerFinder::setQEDInitialEvolutionScales"
			<< (**cit) << Exception::eventerror;
    }
    // calculate the probabilities
    double prob(0.);
    for(unsigned int ix=0;ix<partners.size();++ix) prob += partners[ix].first;
    // normalise
    for(unsigned int ix=0;ix<partners.size();++ix) partners[ix].first /= prob;
    // set the partner if required
    int position(-1);
    // use QCD partner if set
    if(!setPartners&&(*cit)->partner()) {
      for(unsigned int ix=0;ix<partners.size();++ix) {
	if((*cit)->partner()==partners[ix].second) {
	  position = ix;
	  break;
	}
      }
    }
    // set the partner
    if(setPartners||!(*cit)->partner()||position<0) {
      prob = UseRandom::rnd();
      for(unsigned int ix=0;ix<partners.size();++ix) {
 	if(partners[ix].first>prob) {
	  position = ix;
	  break;
	}
	prob -= partners[ix].first;
      }
      if(position>=0&&(setPartners||!(*cit)->partner())) {
	(*cit)->partner(partners[position].second);
      }
    }
    // must have a partner
    if(position<0) throw Exception() << "Failed to find partner in " 
				 << "PartnerFinder::setQEDInitialEvolutionScales"
				 << (**cit) << Exception::eventerror; 
    // Calculate the evolution scales for all possible pairs of of particles
    vector<pair<Energy,Energy> > scales;
    for(unsigned int ix=0;ix< partners.size();++ix) {
      scales.push_back(calculateInitialEvolutionScales(ShowerPPair(*cit,partners[ix].second),
						       isDecayCase));
    }
    // store all the possible partners
    for(unsigned int ix=0;ix<partners.size();++ix) {
      (**cit).addPartner(ShowerParticle::EvolutionPartner(partners[ix].second,
							  partners[ix].first,
							  ShowerPartnerType::QED,
							  scales[ix].first));
    }
    // set scales
    (**cit).scales().QED      = scales[position].first;
    (**cit).scales().QED_noAO = scales[position].first;
  }
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

vector< pair<ShowerPartnerType, tShowerParticlePtr> > 
PartnerFinder::findQCDPartners(tShowerParticlePtr particle,
			       const ShowerParticleVector &particles) {
  vector< pair<ShowerPartnerType, tShowerParticlePtr> > partners;
  ShowerParticleVector::const_iterator cjt;
  for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
    if(!(*cjt)->data().coloured() || particle==*cjt) continue;
    // one initial-state and one final-state particle
    if(FS(particle) != FS(*cjt)) {
      // loop over all the colours of both particles
      for(unsigned int ix=0; ix<CLSIZE(particle); ++ix) {
	for(unsigned int jx=0; jx<CLSIZE(*cjt); ++jx) {
	  if((CL(particle,ix) && CL(particle,ix)==CL(*cjt,jx))) {
	    partners.push_back(make_pair(ShowerPartnerType::    QCDColourLine,*cjt));
	  }
	}
      }
      //loop over all the anti-colours of both particles
      for(unsigned int ix=0; ix<ACLSIZE(particle); ++ix) {
	for(unsigned int jx=0; jx<ACLSIZE(*cjt); ++jx) {
	  if((ACL(particle,ix) && ACL(particle,ix)==ACL(*cjt,jx))) {
	    partners.push_back(make_pair(ShowerPartnerType::QCDAntiColourLine,*cjt));
	  }
	}
      }
    }
    // two initial-state or two final-state particles
    else {
      //loop over the colours of the first particle and the anti-colours of the other
      for(unsigned int ix=0; ix<CLSIZE(particle); ++ix){
	for(unsigned int jx=0; jx<ACLSIZE(*cjt); ++jx){
	  if(CL(particle,ix) && CL(particle,ix)==ACL(*cjt,jx)) {
	    partners.push_back(make_pair(ShowerPartnerType::    QCDColourLine,*cjt));
	  }
	}
      }
      //loop over the anti-colours of the first particle and the colours of the other
      for(unsigned int ix=0; ix<ACLSIZE(particle); ++ix){
	for(unsigned int jx=0; jx<CLSIZE(*cjt); jx++){
	  if(ACL(particle,ix) && ACL(particle,ix)==CL(*cjt,jx)) {
	    partners.push_back(make_pair(ShowerPartnerType::QCDAntiColourLine,*cjt));
	  }
	}
      }
    }
  }
  // if we haven't found any partners look for RPV
  if (partners.empty()) {
    // special for RPV
    tColinePtr col = CL(particle); 
    if(FS(particle)&&col&&col->sourceNeighbours().first) {
      tColinePair cpair = col->sourceNeighbours();
      for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	if(( FS(*cjt) && ( CL(*cjt) == cpair.first || CL(*cjt)  == cpair.second))||
	   (!FS(*cjt) && (ACL(*cjt) == cpair.first || ACL(*cjt) == cpair.second ))) {
	  partners.push_back(make_pair(ShowerPartnerType::    QCDColourLine,*cjt));
	}
      }
    }
    else if(col&&col->sinkNeighbours().first) {
      tColinePair cpair = col->sinkNeighbours();
      for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	if(( FS(*cjt) && (ACL(*cjt) == cpair.first || ACL(*cjt)  == cpair.second))||
	   (!FS(*cjt) && ( CL(*cjt) == cpair.first ||  CL(*cjt) == cpair.second))) {
	  partners.push_back(make_pair(ShowerPartnerType::    QCDColourLine,*cjt));    
	}
      }
    }
    col = ACL(particle);
    if(FS(particle)&&col&&col->sinkNeighbours().first) {
      tColinePair cpair = col->sinkNeighbours();
      for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	if(( FS(*cjt) && (ACL(*cjt) == cpair.first || ACL(*cjt)  == cpair.second))||
	   (!FS(*cjt) && ( CL(*cjt) == cpair.first ||  CL(*cjt) == cpair.second ))) {
	  partners.push_back(make_pair(ShowerPartnerType::QCDAntiColourLine,*cjt));  
	}
      }
    }
    else if(col&&col->sourceNeighbours().first) {
      tColinePair cpair = col->sourceNeighbours();
      for(cjt=particles.begin();cjt!=particles.end();++cjt) {
	if(( FS(*cjt) && ( CL(*cjt) == cpair.first || CL(*cjt) == cpair.second))||
	   (!FS(*cjt) && (ACL(*cjt) == cpair.first ||ACL(*cjt) == cpair.second))) {
	  partners.push_back(make_pair(ShowerPartnerType::QCDAntiColourLine,*cjt));     
	}
      }
    }
  }
  // return the partners
  return partners;
}

vector< pair<double, tShowerParticlePtr> > 
PartnerFinder::findQEDPartners(tShowerParticlePtr particle,
			       const ShowerParticleVector &particles,
			       const bool isDecayCase) {
  vector< pair<double, tShowerParticlePtr> > partners;
  ShowerParticleVector::const_iterator cjt;
  double pcharge = particle->id()==ParticleID::gamma ? 1 : double(particle->data().iCharge());
  vector< pair<double, tShowerParticlePtr> > photons;
  for(cjt = particles.begin(); cjt != particles.end(); ++cjt) {
    if(particle == *cjt) continue;
    if((**cjt).id()==ParticleID::gamma) photons.push_back(make_pair(1.,*cjt));
    if(!(*cjt)->data().charged() ) continue;
    double charge = pcharge*double((*cjt)->data().iCharge());
    if( FS(particle) != FS(*cjt) ) charge *=-1.;
    if( QEDPartner_ != 0 && !isDecayCase ) {
      // only include II and FF as requested
      if( QEDPartner_ == 1 && FS(particle) != FS(*cjt) )
	continue;
      // only include IF is requested
      else if(QEDPartner_ == 2 && FS(particle) == FS(*cjt) )
	continue;
    }
    if(particle->id()==ParticleID::gamma) charge = -abs(charge);
    // only keep positive dipoles
    if(charge<0.) partners.push_back(make_pair(-charge,*cjt));
  }
  if(particle->id()==ParticleID::gamma&& partners.empty()) {
    return photons;
  }
  return partners;
}
