// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartnerFinder class.
//

#include "PartnerFinder.h"
#include "ShowerVariables.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PartnerFinder.tcc"
#endif

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower2/Kinematics/ShowerParticle.h"
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
  os << _approach << _showerVariables;
}

void PartnerFinder::persistentInput(PersistentIStream & is, int) {
  is >> _approach >> _showerVariables;
}

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
    if (partners.empty())
      throw Exception() << "Failed to make colour connections in " 
			<< "PartnerFinder::setQCDInitialEvolutionScales"
			<< (**cit)
			<< Exception::eventerror;
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

bool PartnerFinder::setQEDInitialEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase) {

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

bool PartnerFinder::setEWKInitialEvolutionScales(const ShowerParticleVector &particles,
						 const bool isDecayCase) {
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

pair<Energy,Energy> PartnerFinder::
calculateInitialFinalScales(const ShowerPPair &ppair, const bool isDecayCase) {

  Lorentz5Momentum pa, pb, pc;
  pb = ppair.first->momentum();
  pc = ppair.second->momentum();

  if(!isDecayCase) {
    /********
     * In this case from JHEP 12(2003)045 we find the conditions
     * ktilda_b = (1+c) and ktilda_c = (1+2c)
     * We also find that c = m_c^2/Q^2. The process is a+b->c where
     * particle a is not colour connected (considered as a colour singlet).
     * Therefore we simply find that q_b = sqrt(Q^2+m_c^2) and 
     * q_c = sqrt(Q^2+2 m_c^2)
     * We also assume that the first particle in the pair is the initial
     * state particle and the second is the final state one (c)
     *********/ 
    pa = pb+pc;
    pa.boost(pa.findBoostToCM());
    Energy2  mc2 = sqr(ppair.second->mass());
    Energy2  Q2  = pa*pa;
    return pair<Energy,Energy>(sqrt(Q2+mc2), sqrt(Q2+2*mc2));
  }
  else {
    /********
     * In this case from JHEP 12(2003)045 we find, for the decay
     * process b->c+a(neutral), the condition
     * (ktilda_b-1)*(ktilda_c-c)=(1/4)*sqr(1-a+c+lambda). 
     * We also assume that the first particle in the pair is the initial
     * state particle (b) and the second is the final state one (c).
     *  - We find maximal phase space coverage through emissions from 
     *    c if we set ktilde_c = 4.*(sqr(1.-sqrt(a))-c)
     *  - We find the most 'symmetric' way to populate the phase space
     *    occurs for (ktilda_b-1)=(ktilda_c-c)=(1/2)*(1-a+c+lambda) 
     *  - We find the most 'smooth' way to populate the phase space
     *    occurs for... 
     *********/ 
    unsigned int phase_space = _showerVariables->decay_shower_partition();
    pa              = pb-pc               ;
    Energy ma       = sqrt(pa*pa)         ;      
    Energy mb       = ppair.first->mass() ;      
    Energy mc       = ppair.second->mass();
    double a        = sqr(ma/mb)          ;
    double c        = sqr(mc/mb)          ;
    double lambda   = sqrt(1. + a*a + c*c - 2.*a - 2.*c - 2.*a*c);
    double PROD     = 0.25*sqr(1. - a + c + lambda);
    double ktilde_b, ktilde_c,cosi(0.);
    switch(phase_space) {
    case 0: // the 'symmetric' choice
	ktilde_c = 0.5*(1-a+c+lambda) + c ;
	ktilde_b = 1.+PROD/(ktilde_c-c)   ;
      break;
    case 1:  // the 'maximal' choice
	ktilde_c = 4.0*(sqr(1.-sqrt(a))-c);
	ktilde_b = 1.+PROD/(ktilde_c-c)   ;
      break;
    case 2:  // the 'smooth' choice
	cosi = (sqr(1-sqrt(c))-a)/lambda;
	ktilde_b = 2.0/(1.0-cosi);
	ktilde_c = (1.0-a+c+lambda)*(1.0+c-a-lambda*cosi)/(2.0*(1.0+cosi));
      break;
    default:
      throw Exception() << "Invalid option for decay shower's phase space"
			<< " PartnerFinder::calculateInitialFinalScales"
			<< Exception::abortnow;
    }
    return pair<Energy,Energy>(mb*sqrt(ktilde_b),mb*sqrt(ktilde_c));
  }
}

pair<Energy,Energy> PartnerFinder::
calculateInitialInitialScales(const ShowerPPair &ppair)
{
  /*******
   * This case is quite simple. From JHEP 12(2003)045 we find the condition
   * that ktilda_b = ktilda_c = 1. In this case we have the process
   * b+c->a so we need merely boost to the CM frame of the two incoming
   * particles and then qtilda is equal to the energy in that frame
   **********/
  Lorentz5Momentum p(ppair.first->momentum()+ppair.second->momentum());
  p.boost(p.findBoostToCM());
  Energy Q = sqrt(p.m2());
  return pair<Energy,Energy>(Q,Q);
}

pair<Energy,Energy> PartnerFinder::
calculateFinalFinalScales(const ShowerPPair &particlePair) {
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

  Lorentz5Momentum p1 = particlePair.first->momentum(); 
  Lorentz5Momentum p2 = particlePair.second->momentum(); 
  Lorentz5Momentum p = p1; 
  p.boost((p1+p2).findBoostToCM());
  Lorentz5Momentum n(0.0, -p.vect()); 
  Energy firstQ = sqrt(2.*p*(p+n)); 
  p = p2;
  p.boost((p1+p2).findBoostToCM());
  n = Lorentz5Momentum(0.0, - p.vect()); 
  Energy secondQ = sqrt(2.*p*(p+n));
  // get asymmetric distribution in x, xbar plane: 
  if (_showerVariables->asyPS()) {
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
    }
  }
  // random swap if needed
  if (_showerVariables->rndPS()&&UseRandom::rndbool()) 
    swap(firstQ,secondQ);
  return pair<Energy,Energy>(firstQ, secondQ);
}
