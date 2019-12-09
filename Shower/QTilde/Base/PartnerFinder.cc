// -*- C++ -*-
//
// PartnerFinder.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ThePEG/Repository/UseRandom.h" 
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeClass<PartnerFinder,Interfaced>
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
  inline size_t
  CLSIZE(const tShowerParticlePtr a) {
    return a->colourInfo()->colourLines().size();
  }
  
  inline Ptr<ThePEG::ColourLine>::transient_pointer
  ACL(const tShowerParticlePtr a, unsigned int index=0) {
    return a->colourInfo()->antiColourLines().empty() ?  ThePEG::tColinePtr() :
      const_ptr_cast<ThePEG::tColinePtr>(a->colourInfo()->antiColourLines()[index]);
  }

  inline size_t
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
  // (this means that no modifications of the following loop is needed!

	for ( const auto & sp : particles ) {
		  // Skip colourless particles
      if(!sp->data().coloured()) continue;
			// find the partners
      auto partners = findQCDPartners(sp,particles);
      // must have a partner
      if(partners.empty()) {
				throw Exception() << "`Failed to make colour connections in " 
                          << "PartnerFinder::setQCDInitialEvolutionScales"
                          << *sp
                          << Exception::eventerror;
      }
      // Calculate the evolution scales for all possible pairs of of particles
			vector<pair<Energy,Energy> > scales;
			int position = -1;
			for(size_t ix=0; ix< partners.size(); ++ix) {
					scales.push_back(
					  calculateInitialEvolutionScales(
					      ShowerPPair(sp, partners[ix].second),
							  isDecayCase
						)
					);
					if (!setPartners && partners[ix].second) position = ix; 
			}
			assert(setPartners || position >= 0);

			// set partners if required
			if (setPartners) {
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
	      // random choice
	      if( partnerMethod_ == 0 ) {
					// random choice of partner
					position = UseRandom::irnd(partners.size());
	      }
	      // take the one with largest angle
	      else if (partnerMethod_ == 1 ) {
					if (sp->perturbative() == 1 && 
		    			sp->dataPtr()->iColour()==PDT::Colour8 ) {
		  			assert(partners.size()==2);
		  			// Determine largest angle
		  			double maxAngle(0.);
		  			for(unsigned int ix=0;ix<partners.size();++ix) {
		    			double angle = sp->momentum().vect().
		      										angle(partners[ix].second->momentum().vect());
		    			if(angle>maxAngle) {
		      			maxAngle = angle;
		      			position = ix;
		    			}
		  			}
					}
					else position = UseRandom::irnd(partners.size());
	      }
	      else assert(false);
		    // set the evolution partner
		    sp->partner(partners[position].second);
			}


      // primary partner set, set the others and do the scale
		  for(size_t ix=0; ix<partners.size(); ++ix) {
          sp->addPartner(
             ShowerParticle::EvolutionPartner(
              partners[ix].second,
					    1.,
					    partners[ix].first,
					    scales[ix].first
					  )
					);
      }

      // set scales for all interactions to that of the partner, default
      Energy scale = scales[position].first;
      for(unsigned int ix=0;ix<partners.size();++ix) {
				if(partners[ix].first==ShowerPartnerType::QCDColourLine) {
	  			sp->scales().QCD_c = 
	    			sp->scales().QCD_c_noAO = 
	    			(scaleChoice_==0 ? scale : scales[ix].first);
				}
				else if(partners[ix].first==ShowerPartnerType::QCDAntiColourLine) {
	  			sp->scales().QCD_ac = 
	    			sp->scales().QCD_ac_noAO =
	    			(scaleChoice_==0 ? scale : scales[ix].first);
				}
				else assert(false);
      }
  }
}




void PartnerFinder::setInitialQEDEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase,
						 const bool setPartners) {
  // loop over all the particles
  for(const auto & sp : particles) {
    // not charged or photon continue
    if(!sp->dataPtr()->charged()) continue;
    // find the potential partners
    vector<pair<double,tShowerParticlePtr> > partners = findQEDPartners(sp,particles,isDecayCase);
    if(partners.empty()) {
      throw Exception() << "Failed to find partner in " 
			<< "PartnerFinder::setQEDInitialEvolutionScales"
			<< *sp << Exception::eventerror;
    }
    // calculate the probabilities
    double prob(0.);
    for(unsigned int ix=0;ix<partners.size();++ix) prob += partners[ix].first;
    // normalise
    for(unsigned int ix=0;ix<partners.size();++ix) partners[ix].first /= prob;
    // set the partner if required
    int position(-1);
    // use QCD partner if set
    if(!setPartners&&sp->partner()) {
      for(unsigned int ix=0;ix<partners.size();++ix) {
	if(sp->partner()==partners[ix].second) {
	  position = ix;
	  break;
	}
      }
    }
    // set the partner
    if(setPartners||!sp->partner()||position<0) {
      prob = UseRandom::rnd();
      for(unsigned int ix=0;ix<partners.size();++ix) {
 	if(partners[ix].first>prob) {
	  position = ix;
	  break;
	}
	prob -= partners[ix].first;
      }
      if(position>=0&&(setPartners||!sp->partner())) {
	sp->partner(partners[position].second);
      }
    }
    // must have a partner
    if(position<0) throw Exception() << "Failed to find partner in " 
				 << "PartnerFinder::setQEDInitialEvolutionScales"
				 << *sp << Exception::eventerror; 
    // Calculate the evolution scales for all possible pairs of of particles
    vector<pair<Energy,Energy> > scales;
    for(unsigned int ix=0;ix< partners.size();++ix) {
      scales.push_back(calculateInitialEvolutionScales(ShowerPPair(sp,partners[ix].second),
						       isDecayCase));
    }
    // store all the possible partners
    for(unsigned int ix=0;ix<partners.size();++ix) {
      sp->addPartner(ShowerParticle::EvolutionPartner(partners[ix].second,
							  partners[ix].first,
							  ShowerPartnerType::QED,
							  scales[ix].first));
    }
    // set scales
    sp->scales().QED      = scales[position].first;
    sp->scales().QED_noAO = scales[position].first;
  }
}

pair<Energy,Energy> PartnerFinder::
calculateInitialEvolutionScales(const ShowerPPair &particlePair,
                                const bool isDecayCase) {
  bool FS1=FS(particlePair.first),FS2= FS(particlePair.second);

  if(FS1 && FS2)
    return calculateFinalFinalScales(particlePair.first->momentum(),
                                     particlePair.second->momentum());
  else if(FS1 && !FS2) {
    pair<Energy,Energy> rval = calculateInitialFinalScales(particlePair.second->momentum(),
                                                           particlePair.first->momentum(),
                                                           isDecayCase);
    return { rval.second, rval.first };
  }
  else if(!FS1 &&FS2)
    return calculateInitialFinalScales(particlePair.first->momentum(),particlePair.second->momentum(),isDecayCase);
  else
    return calculateInitialInitialScales(particlePair.first->momentum(),particlePair.second->momentum());
}

vector< pair<ShowerPartnerType, tShowerParticlePtr> > 
PartnerFinder::findQCDPartners(tShowerParticlePtr particle,
			       const ShowerParticleVector &particles) {
  vector< pair<ShowerPartnerType, tShowerParticlePtr> > partners;
  for(const auto & sp : particles) {
    if(!sp->data().coloured() || particle==sp) continue;
    // one initial-state and one final-state particle
    if(FS(particle) != FS(sp)) {
      // loop over all the colours of both particles
      for(size_t ix=0; ix<CLSIZE(particle); ++ix) {
	for(size_t jx=0; jx<CLSIZE(sp); ++jx) {
	  if((CL(particle,ix) && CL(particle,ix)==CL(sp,jx))) {
	    partners.push_back({ ShowerPartnerType::    QCDColourLine, sp });
	  }
	}
      }
      //loop over all the anti-colours of both particles
      for(size_t ix=0; ix<ACLSIZE(particle); ++ix) {
	for(size_t jx=0; jx<ACLSIZE(sp); ++jx) {
	  if((ACL(particle,ix) && ACL(particle,ix)==ACL(sp,jx))) {
	    partners.push_back({ ShowerPartnerType::QCDAntiColourLine, sp });
	  }
	}
      }
    }
    // two initial-state or two final-state particles
    else {
      //loop over the colours of the first particle and the anti-colours of the other
      for(size_t ix=0; ix<CLSIZE(particle); ++ix){
	for(size_t jx=0; jx<ACLSIZE(sp); ++jx){
	  if(CL(particle,ix) && CL(particle,ix)==ACL(sp,jx)) {
	    partners.push_back({ ShowerPartnerType::    QCDColourLine, sp });
	  }
	}
      }
      //loop over the anti-colours of the first particle and the colours of the other
      for(size_t ix=0; ix<ACLSIZE(particle); ++ix){
	for(size_t jx=0; jx<CLSIZE(sp); jx++){
	  if(ACL(particle,ix) && ACL(particle,ix)==CL(sp,jx)) {
	    partners.push_back({ ShowerPartnerType::QCDAntiColourLine, sp });
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
      for(const auto & sp : particles) {
	if(( FS(sp) && ( CL(sp) == cpair.first || CL(sp)  == cpair.second))||
	   (!FS(sp) && (ACL(sp) == cpair.first || ACL(sp) == cpair.second ))) {
	  partners.push_back({ ShowerPartnerType::    QCDColourLine, sp });
	}
      }
    }
    else if(col&&col->sinkNeighbours().first) {
      tColinePair cpair = col->sinkNeighbours();
      for(const auto & sp : particles) {
	if(( FS(sp) && (ACL(sp) == cpair.first || ACL(sp)  == cpair.second))||
	   (!FS(sp) && ( CL(sp) == cpair.first ||  CL(sp) == cpair.second))) {
	  partners.push_back({ ShowerPartnerType::    QCDColourLine, sp });
	}
      }
    }
    col = ACL(particle);
    if(FS(particle)&&col&&col->sinkNeighbours().first) {
      tColinePair cpair = col->sinkNeighbours();
      for(const auto & sp : particles) {
	if(( FS(sp) && (ACL(sp) == cpair.first || ACL(sp)  == cpair.second))||
	   (!FS(sp) && ( CL(sp) == cpair.first ||  CL(sp) == cpair.second ))) {
	  partners.push_back({ ShowerPartnerType::QCDAntiColourLine, sp });
	}
      }
    }
    else if(col&&col->sourceNeighbours().first) {
      tColinePair cpair = col->sourceNeighbours();
      for(const auto & sp : particles) {
	if(( FS(sp) && ( CL(sp) == cpair.first || CL(sp) == cpair.second))||
	   (!FS(sp) && (ACL(sp) == cpair.first ||ACL(sp) == cpair.second))) {
	  partners.push_back({ ShowerPartnerType::QCDAntiColourLine, sp });
	}
      }
    }
  }
  // return the partners
  return partners;
}

vector< pair<double, tShowerParticlePtr> > 
PartnerFinder::findQEDPartners(tShowerParticlePtr particle,
			       const ShowerParticleVector & particles,
			       const bool isDecayCase) {
  vector< pair<double, tShowerParticlePtr> > partners;
  const double pcharge = 
    particle->id()==ParticleID::gamma ? 1 : double(particle->data().iCharge());
  vector< pair<double, tShowerParticlePtr> > photons;
  for (const auto & sp : particles) {
    if (particle == sp) continue;
    if (sp->id()==ParticleID::gamma) photons.push_back(make_pair(1.,sp));
    if (!sp->data().charged() ) continue;
    double charge = pcharge*double((sp)->data().iCharge());
    if ( FS(particle) != FS(sp) ) charge *=-1.;
    if ( QEDPartner_ != 0 && !isDecayCase ) {
      // only include II and FF as requested
      if ( QEDPartner_ == 1 && FS(particle) != FS(sp) )
	continue;
      // only include IF is requested
      else if (QEDPartner_ == 2 && FS(particle) == FS(sp) )
	continue;
    }
    if (particle->id()==ParticleID::gamma) charge = -abs(charge);
    // only keep positive dipoles
    if (charge<0.) partners.push_back({ -charge, sp });
  }
  if (particle->id()==ParticleID::gamma && partners.empty()) {
    return photons;
  }
  return partners;
}

pair<Energy,Energy> 
PartnerFinder::calculateFinalFinalScales(
        const Lorentz5Momentum & p1,
        const Lorentz5Momentum & p2) 
{
  static const double eps=1e-7;
  // Using JHEP 12(2003)045 we find that we need ktilde = 1/2(1+b-c+lambda)
  // ktilde = qtilde^2/Q^2 therefore qtilde = sqrt(ktilde*Q^2)
  // find momenta in rest frame of system
  // calculate quantities for the scales
  Energy2 Q2 = (p1+p2).m2();
  double b = p1.mass2()/Q2;
  double c = p2.mass2()/Q2;
  if(b<0.) {
    if(b<-eps) {
      throw Exception() << "Negative Mass squared b = " << b
			<< "in PartnerFinder::calculateFinalFinalScales()"
			<< Exception::eventerror;
    }
    b = 0.;
  }
  if(c<0.) {
    if(c<-eps) {
      throw Exception() << "Negative Mass squared c = " << c
			<< "in PartnerFinder::calculateFinalFinalScales()"
			<< Exception::eventerror;
    }
    c = 0.;
  }
  // KMH & PR - 16 May 2008 - swapped lambda calculation from 
  // double lam=2.*p1.vect().mag()/Q; to sqrt(kallen(1,b,c)), 
  // which should be identical for p1 & p2 onshell in their COM
  // but in the inverse construction for the Nason method, this
  // was not the case, leading to misuse. 
  const double lam=sqrt((1.+sqrt(b)+sqrt(c))*(1.-sqrt(b)-sqrt(c))
                   *(sqrt(b)-1.-sqrt(c))*(sqrt(c)-1.-sqrt(b)));
  // symmetric case
  return { sqrt(0.5*Q2*(1.+b-c+lam)), sqrt(0.5*Q2*(1.-b+c+lam)) };
}


pair<Energy,Energy>
PartnerFinder::calculateInitialFinalScales(const Lorentz5Momentum& pb, const Lorentz5Momentum& pc,
			    const bool isDecayCase) {
  if(!isDecayCase) { 
    // In this case from JHEP 12(2003)045 we find the conditions
    // ktilde_b = (1+c) and ktilde_c = (1+2c)
    // We also find that c = m_c^2/Q^2. The process is a+b->c where
    // particle a is not colour connected (considered as a colour singlet).
    // Therefore we simply find that q_b = sqrt(Q^2+m_c^2) and 
    // q_c = sqrt(Q^2+2 m_c^2)
    // We also assume that the first particle in the pair is the initial
    // state particle and the second is the final state one c 
    const Energy2  mc2 = sqr(pc.mass());
    const Energy2  Q2  = -(pb-pc).m2();
    return { sqrt(Q2+mc2), sqrt(Q2+2*mc2) };
  }
  else {    
    // In this case from JHEP 12(2003)045 we find, for the decay
    // process b->c+a(neutral), the condition
    // (ktilde_b-1)*(ktilde_c-c)=(1/4)*sqr(1-a+c+lambda). 
    // We also assume that the first particle in the pair is the initial
    // state particle (b) and the second is the final state one (c).
    //  - We find maximal phase space coverage through emissions from 
    //    c if we set ktilde_c = 4.*(sqr(1.-sqrt(a))-c)
    //  - We find the most 'symmetric' way to populate the phase space
    //    occurs for (ktilde_b-1)=(ktilde_c-c)=(1/2)*(1-a+c+lambda) 
    //  - We find the most 'smooth' way to populate the phase space
    //    occurs for...
    Energy2 mb2(sqr(pb.mass()));
    double a=(pb-pc).m2()/mb2;
    double c=sqr(pc.mass())/mb2;
    double lambda   = 1. + a*a + c*c - 2.*a - 2.*c - 2.*a*c;
    lambda = sqrt(max(lambda,0.));
    const double PROD     = 0.25*sqr(1. - a + c + lambda);
    const double ktilde_c = 0.5*(1-a+c+lambda) + c ;
    const double ktilde_b = 1.+PROD/(ktilde_c-c)   ;
    return { sqrt(mb2*ktilde_b), sqrt(mb2*ktilde_c) };
  }
}

pair<Energy,Energy>
PartnerFinder::calculateInitialInitialScales(const Lorentz5Momentum& p1, const Lorentz5Momentum& p2) {
  // This case is quite simple. From JHEP 12(2003)045 we find the condition
  // that ktilde_b = ktilde_c = 1. In this case we have the process
  // b+c->a so we need merely boost to the CM frame of the two incoming
  // particles and then qtilde is equal to the energy in that frame
  const Energy Q = sqrt((p1+p2).m2());
  return {Q,Q};
}
