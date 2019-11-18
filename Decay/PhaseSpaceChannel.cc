// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhaseSpaceChannel class.
//

#include "PhaseSpaceChannel.h"
#include "PhaseSpaceMode.h"
#include "Herwig/Utilities/Kinematics.h"

/** 
 *  Constructor with incoming particles
 */
PhaseSpaceChannel::PhaseSpaceChannel(tPhaseSpaceModePtr inm, bool skip) : mode_(inm), weight_(1.),
									  initialized_(false),
									  skipFirst_(skip) {
  if(!inm->incoming().second)
     intermediates_.push_back(PhaseSpaceResonance(inm->incoming().first));
}

void PhaseSpaceChannel::init(tPhaseSpaceModePtr mode) {
  mode_=mode;
  if(initialized_) return;
  initialized_=true;
  // find the descendents
  for(PhaseSpaceResonance & res : intermediates_)
    findChildren(res,res.descendents);
  // ensure intermediates either have the width set, or
  // can't possibly be on-shell
  // first the maximum energy release
  Energy massmax = mode->eMax_;
  for(tcPDPtr part : mode->outgoing_)
    massmax -= mode->testOnShell_ ? part->mass() : part->massMin();
  for(PhaseSpaceResonance & res : intermediates_) {
    if(!res.particle || res.mWidth!=ZERO ||
       res.jacobian != PhaseSpaceResonance::BreitWigner) continue;
    Energy massmin(ZERO);
    for(const int & ext : res.descendents)
      massmin += mode->testOnShell_ ? 
	mode->outgoing_[ext-1]->mass() : mode->outgoing_[ext-1]->massMin();
    // check if can be on-shell
    Energy mass = sqrt(res.mass2);
    if(mass>=massmin&&mass<=massmax+massmin) {
      string modeout = mode->incoming_.first->PDGName() + " ";
      if(mode->incoming_.second)
	modeout += mode->incoming_.second->PDGName() + " ";
      for( tcPDPtr out : mode->outgoing_)
	modeout += out->PDGName() + " ";
      throw InitException() << "Width zero for " << res.particle->PDGName()
 			      << " in PhaseSpaceChannel::init() "
			    << modeout << Exception::runerror;
    }
  }
}

ostream & Herwig::operator<<(ostream & os, const PhaseSpaceChannel & channel) {
  // output of the external particles
  if(!channel.mode_->incoming().second) {
    os << "Channel for the decay of " << channel.mode_->incoming().first->PDGName() << " -> ";
  }
  else 
    os << "Channel for " << channel.mode_->incoming().first->PDGName()
       << ", " << channel.mode_->incoming().second->PDGName()
       << " -> ";
  for(tcPDPtr part : channel.mode_->outgoing())
    os << part->PDGName() << " ";
  os << endl;
  os << "Proceeds in following steps ";
  for(unsigned int ix=0;ix<channel.intermediates_.size();++ix) {
    os << channel.intermediates_[ix].particle->PDGName() << " -> ";
    if(channel.intermediates_[ix].children.first>0) {
      os << channel.mode_->outgoing()[channel.intermediates_[ix].children.first-1]->PDGName()  
  	 << "(" << channel.intermediates_[ix].children.first<< ") ";
    }
    else {
      os << channel.intermediates_[-channel.intermediates_[ix].children.first].particle->PDGName() 
     	 << "(" << channel.intermediates_[ix].children.first<< ") ";
    }
    if(channel.intermediates_[ix].children.second>0) {
      os << channel.mode_->outgoing()[channel.intermediates_[ix].children.second-1]->PDGName()  
  	 << "(" << channel.intermediates_[ix].children.second<< ") ";
    }
    else {
      os << channel.intermediates_[-channel.intermediates_[ix].children.second].particle->PDGName() 
     	 << "(" << channel.intermediates_[ix].children.second<< ") ";
    }
    os << endl;
  }
  return os;
}
void PhaseSpaceChannel::initrun(tPhaseSpaceModePtr mode) {
  mode_=mode;
  if(!mode_->testOnShell_) return;
  // ensure intermediates either have the width set, or
  // can't possibly be on-shell
  Energy massmax = mode_->incoming().first->massMax();
  for(tPDPtr out : mode_->outgoing()) 
    massmax -= out->massMin();
  for(unsigned int ix=1;ix<intermediates_.size();++ix) { 
    if(intermediates_[ix].mWidth==ZERO && intermediates_[ix].jacobian==PhaseSpaceResonance::BreitWigner) {
      Energy massmin(ZERO);
      for(const int & iloc : intermediates_[ix].descendents) 
	massmin += mode_->outgoing()[iloc-1]->massMin();
      // check if can be on-shell
      if(intermediates_[ix].mass2>=sqr(massmin)&&
	 intermediates_[ix].mass2<=sqr(massmax+massmin)) {
	string modeout = mode_->incoming().first->PDGName() + " -> ";
	for(tPDPtr out : mode_->outgoing()) 
	  modeout += out->PDGName() + " ";
	throw Exception() << "Width zero for " << intermediates_[ix].particle->PDGName()
			  << " in PhaseSpaceChannel::initrun() "
			  << modeout
			  << Exception::runerror;
      }
    }
  }
}

// generate the momenta of the external particles
vector<Lorentz5Momentum> 
PhaseSpaceChannel::generateMomenta(const Lorentz5Momentum & pin,
				   const vector<Energy> & massext) const {
  // storage of the momenta of the external particles
  vector<Lorentz5Momentum> pexternal(massext.size());
  // and the internal particles
  vector<Lorentz5Momentum> pinter(1,pin);
  pinter.resize(intermediates_.size());
  // masses of the intermediate particles
  vector<Energy> massint(intermediates_.size());
  massint[0]=pin.mass();
  // generate all the decays in the chain
  for(unsigned int ix=0;ix<intermediates_.size();++ix) { 
    int idau[2] = {abs(intermediates_[ix].children.first),
		   abs(intermediates_[ix].children.second)};
    // if both decay products off-shell
    if(intermediates_[ix].children.first<0&&intermediates_[ix].children.second<0) {
      double rnd  = mode_->rStack_.top();
      mode_->rStack_.pop();
      Energy lowerb[2];
       // lower limits on the masses of the two resonances
      for(unsigned int iy=0;iy<2;++iy) {
   	lowerb[iy]=ZERO;
   	bool massless=true;
   	for(const int & des : intermediates_[idau[iy]].descendents) {
   	  if(massext[des-1]!=ZERO) massless = false;
   	  lowerb[iy]+=massext[des-1];
  	}
   	if(massless) lowerb[iy] = mode_->epsilonPS(); 
      }
      double rnd2 = mode_->rStack_.top();
      mode_->rStack_.pop();
      // randomize the order
      if(rnd<0.5) {
	rnd*=2.;
  	// mass of the first resonance
  	Energy upper = massint[ix]-lowerb[1];
  	Energy lower = lowerb[0];
   	massint[idau[0]]=generateMass(intermediates_[idau[0]],lower,upper,rnd );
   	// mass of the second resonance
   	upper = massint[ix]-massint[idau[0]];
   	lower = lowerb[1];
   	massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper,rnd2);
      }
      else {
	rnd = 2.*rnd-1.;
  	// mass of the second resonance
  	Energy upper = massint[ix]-lowerb[0];
  	Energy lower = lowerb[1];
   	massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper,rnd2);
   	// mass of the first resonance
   	upper = massint[ix]-massint[idau[1]];
   	lower = lowerb[0];
   	massint[idau[0]]=generateMass(intermediates_[idau[0]],lower,upper,rnd );
      }
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massint[idau[1]],
   		   pinter[idau[0]],pinter[idau[1]]);
    }
    // only first off-shell
    else if(intermediates_[ix].children.first<0) {
      double rnd  = mode_->rStack_.top();
      mode_->rStack_.pop();
      // compute the limits of integration
      Energy upper = massint[ix]-massext[idau[1]-1];
      Energy lower = ZERO;
      bool massless=true;
      for(const int & des : intermediates_[idau[0]].descendents) {
   	if(massext[des-1]!=ZERO) massless = false;
   	lower+=massext[des-1];
      }
      if(massless) lower = mode_->epsilonPS();
      massint[idau[0]] = generateMass(intermediates_[idau[0]],lower,upper,rnd);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massext[idau[1]-1], 
		   pinter[idau[0]],pexternal[idau[1]-1]);
    }
    // only second off-shell
    else if(intermediates_[ix].children.second<0) {
      double rnd  = mode_->rStack_.top();
      mode_->rStack_.pop();
      // compute the limits of integration
      Energy upper = massint[ix]-massext[idau[0]-1];
      Energy lower = ZERO;
      bool massless=true;
      for(const int & des : intermediates_[idau[1]].descendents) {	
 	if(massext[des-1]!=ZERO) massless = false;
 	lower+=massext[des-1];
      }
      if(massless) lower = mode_->epsilonPS();
      massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper,rnd);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massext[idau[0]-1],massint[idau[1]], 
   		   pexternal[idau[0]-1],pinter[idau[1]]);
    }
    // both on-shell
    else {
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massext[idau[0]-1],massext[idau[1]-1], 
   		   pexternal[idau[0]-1],pexternal[idau[1]-1]);
    }
  }
  // return the external momenta
  return pexternal;
}
void PhaseSpaceChannel::twoBodyDecay(const Lorentz5Momentum & p,
				     const Energy m1, const Energy m2,
				     Lorentz5Momentum & p1,
				     Lorentz5Momentum & p2 ) const {
  static const double eps=1e-6;
  double ctheta = 2.*mode_->rStack_.top()-1.;
  mode_->rStack_.pop();
  double phi = Constants::twopi*mode_->rStack_.top();
  mode_->rStack_.pop();
  Kinematics::generateAngles(ctheta,phi);
  Axis unitDir1=Kinematics::unitDirection(ctheta,phi);
  Momentum3 pstarVector;
  Energy min=p.mass();
  if ( min >= m1 + m2  &&  m1 >= ZERO  &&  m2 >= ZERO  ) {
    pstarVector = unitDir1 * Kinematics::pstarTwoBodyDecay(min,m1,m2);
  }
  else if( m1 >= ZERO  &&  m2 >= ZERO && (m1+m2-min)/(min+m1+m2)<eps) {
    pstarVector = Momentum3();
  }
  else {
    throw PhaseSpaceError() << "Two body decay cannot proceed "
				 << "p = " << p / GeV 
				 << " p.m() = " << min / GeV
				 << " -> " << m1/GeV 
				 << ' ' << m2/GeV << Exception::eventerror;
  }
  p1 = Lorentz5Momentum(m1, pstarVector);
  p2 = Lorentz5Momentum(m2,-pstarVector);
  // boost from CM to LAB
  Boost bv = p.boostVector();
  double gammarest = p.e()/p.mass();
  p1.boost( bv , gammarest );
  p2.boost( bv , gammarest );
}
double PhaseSpaceChannel::atanhelper(const PhaseSpaceResonance & res,
				     Energy limit) const {
  return atan2( sqr(limit) - res.mass2, res.mWidth );
}

// return the weight for a given resonance
InvEnergy2 PhaseSpaceChannel::massWeight(const PhaseSpaceResonance & res,
					 Energy moff, Energy lower,
					 Energy upper) const {
  InvEnergy2 wgt = ZERO;
  if(lower>upper) {
    string modestring = mode_->incoming().first->PDGName() + " -> ";
    for(tPDPtr part :mode_->outgoing()) modestring += part->PDGName() + " ";
    throw PhaseSpaceError() << "PhaseSpaceChannel::massWeight not allowed in"
			    << modestring << " "
			    << res.particle->PDGName() << "   " 
			    << moff/GeV << " " << lower/GeV << " " << upper/GeV
			    << Exception::eventerror;
  }
  // use a Breit-Wigner 
  if ( res.jacobian == PhaseSpaceResonance::BreitWigner ) {
    double rhomin  = atanhelper(res,lower);
    double rhomax  = atanhelper(res,upper) - rhomin;
    if ( rhomax != 0.0 ) {
      Energy2 moff2=moff*moff-res.mass2;
      wgt = res.mWidth/rhomax/(moff2*moff2+res.mWidth*res.mWidth);
    }
    else {
      wgt = 1./((sqr(upper)-sqr(lower))*sqr(sqr(moff)-res.mass2)/
 		(sqr(lower)-res.mass2)/(sqr(upper)-res.mass2));
    }
  }
  // power law
  else if(res.jacobian == PhaseSpaceResonance::Power) {
    double rhomin = pow(sqr(lower/MeV),res.power+1.);
    double rhomax = pow(sqr(upper/MeV),res.power+1.)-rhomin;
    wgt = (res.power+1.)/rhomax*pow(sqr(moff/MeV),res.power)
      /MeV/MeV;
  }
  else if(res.jacobian == PhaseSpaceResonance::OnShell ) {
    wgt = 1./Constants::pi/res.mWidth;
  } 
  else {
    throw PhaseSpaceError() << "Unknown type of Jacobian in " 
			    << "PhaseSpaceChannel::massWeight"
			    << Exception::eventerror;
  }
  return wgt;
}

Energy PhaseSpaceChannel::generateMass(const PhaseSpaceResonance & res,
				       Energy lower,Energy upper,
				       const double & rnd) const {
  static const Energy eps=1e-9*MeV;
  if(lower<eps) lower=eps;
  Energy mass=ZERO;
  if(lower>upper) {
    string modestring = mode_->incoming().first->PDGName() + " -> ";
    for(tPDPtr part :mode_->outgoing()) modestring += part->PDGName() + " ";
    throw PhaseSpaceError() << "PhaseSpaceChannel::generateMass"
			    << " not allowed in"
			    << modestring << " "
			    << res.particle->PDGName()
			    << " " << lower/GeV << " " << upper/GeV
			    << Exception::eventerror;
  }
  if(abs(lower-upper)/(lower+upper)>2e-10) {
    lower +=1e-10*(lower+upper);
    upper -=1e-10*(lower+upper);
  }
  else 
    return 0.5*(lower+upper);
  // use a Breit-Wigner
  if(res.jacobian==PhaseSpaceResonance::BreitWigner) {
    if(res.mWidth!=ZERO) {
      Energy2 lower2 = sqr(lower);
      Energy2 upper2 = sqr(upper);
      
      double rhomin = atan2((lower2 - res.mass2),res.mWidth);
      double rhomax = atan2((upper2 - res.mass2),res.mWidth)-rhomin;
      double rho = rhomin+rhomax*rnd;
      Energy2 mass2 = max(lower2,min(upper2,res.mass2+res.mWidth*tan(rho)));
      if(mass2<ZERO) mass2 = ZERO;
      mass = sqrt(mass2);
    }
    else {
      mass = sqrt(res.mass2+
  		  (sqr(lower)-res.mass2)*(sqr(upper)-res.mass2)/
  		  (sqr(lower)-res.mass2-rnd*(sqr(lower)-sqr(upper))));
    }
  }
  // use a power-law
  else if(res.jacobian == PhaseSpaceResonance::Power) {
    double rhomin = pow(sqr(lower/MeV),res.power+1.);
    double rhomax = pow(sqr(upper/MeV),res.power+1.)-rhomin;
    double rho = rhomin+rhomax*rnd;
    mass = pow(rho,0.5/(res.power+1.))*MeV;
  }
  else if(res.jacobian == PhaseSpaceResonance::OnShell) {
    mass = sqrt(res.mass2);
  } 
  else {
    throw PhaseSpaceError() << "Unknown type of Jacobian in " 
			    << "PhaseSpaceChannel::generateMass" 
			    << Exception::eventerror;
  }
  if(mass<lower+1e-10*(lower+upper))      mass=lower+1e-10*(lower+upper);
  else if(mass>upper-1e-10*(lower+upper)) mass=upper-1e-10*(lower+upper);
  return mass;
}

// generate the weight for this channel given a phase space configuration
double PhaseSpaceChannel::generateWeight(const vector<Lorentz5Momentum> & output) const {
  using Constants::pi;
  // include the prefactor due to the weight of the channel
  double wgt = weight_;
  // work out the masses of the intermediate particles
  vector<Energy> intmass;
  for(const PhaseSpaceResonance & res: intermediates_) {
    Lorentz5Momentum pinter;
    for(const int & des : res.descendents) pinter += output[des-1];
    pinter.rescaleMass();
    intmass.push_back( pinter.mass() );
  }
  Energy2 scale(sqr(intmass[0]));
  // calculate the terms for each of the decays
  for(unsigned int ix=0;ix<intermediates_.size();++ix) {
    int idau[2] = {abs(intermediates_[ix].children.first),
		   abs(intermediates_[ix].children.second)};
    // if both decay products off-shell
    if(intermediates_[ix].children.first<0&&intermediates_[ix].children.second<0) {
      // lower limits on the masses of the two resonances
      Energy lowerb[2];
      for(unsigned int iy=0;iy<2;++iy) {
       	lowerb[iy]=ZERO;
   	for(const int & des : intermediates_[idau[iy]].descendents) {
   	  lowerb[iy]+=output[des-1].mass();
  	}
      }
      // undo effect of randomising
      // weight for the first order
      // contribution of first resonance
      Energy upper = intmass[ix]-lowerb[1];
      Energy lower = lowerb[0];
      InvEnergy2 wgta=massWeight(intermediates_[idau[0]],
				 intmass[idau[0]],lower,upper);
      // contribution of second resonance
      upper = intmass[ix]-intmass[idau[0]];
      lower = lowerb[1];
      InvEnergy4 wgta2 = wgta*massWeight(intermediates_[idau[1]],
					 intmass[idau[1]],lower,upper);
      // weight for the second order
      upper = intmass[ix]-lowerb[0];
      lower = lowerb[1];
      InvEnergy2 wgtb=massWeight(intermediates_[idau[1]],
				 intmass[idau[1]],lower,upper);
      upper = intmass[ix]-intmass[idau[1]];
      lower = lowerb[0];
      InvEnergy4 wgtb2=wgtb*massWeight(intermediates_[idau[0]],
				       intmass[idau[0]],lower,upper);
      // weight factor
      wgt *=0.5*sqr(scale)*(wgta2+wgtb2);
      // factor for the kinematics
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[0]],
						 intmass[idau[1]]);
      if(pcm>ZERO)
      	wgt *= intmass[ix]*8.*pi*pi/pcm;
      else
      	wgt = 0.;
    }
    // only first off-shell
    else if(intermediates_[ix].children.first<0) {
      // compute the limits of integration
      Energy upper = intmass[ix]-output[idau[1]-1].mass();
      Energy lower = ZERO;
      for(const int & des : intermediates_[idau[0]].descendents) {
	lower += output[des-1].mass();
      }
      wgt *=scale*massWeight(intermediates_[idau[0]],intmass[idau[0]],lower,upper);
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[0]],
   					  output[idau[1]-1].mass());
      if(pcm>ZERO)
   	wgt *= intmass[ix]*8.*pi*pi/pcm;
      else
   	wgt = 0.;
    }
    // only second off-shell
    else if(intermediates_[ix].children.second<0) {
      // compute the limits of integration
      Energy upper = intmass[ix]-output[idau[0]-1].mass(); 
      Energy lower = ZERO;
      for(const int & des : intermediates_[idau[1]].descendents) {
	lower += output[des-1].mass();
      }
      wgt *=scale*massWeight(intermediates_[idau[1]],intmass[idau[1]],lower,upper);
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[1]],
						 output[idau[0]-1].mass());
      if(pcm>ZERO)
      	wgt *=intmass[ix]*8.*pi*pi/pcm;
      else
      	wgt=0.;
    }
    // both on-shell
    else {
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],output[idau[1]-1].mass(),
						 output[idau[0]-1].mass());
      if(pcm>ZERO)
	wgt *=intmass[ix]*8.*pi*pi/pcm;
      else
	wgt = 0.;
    }
  }
  // finally the overall factor
  wgt /= pi;
  // return the answer
  return wgt;
}

// generate the final-state particles including the intermediate resonances
void PhaseSpaceChannel::generateIntermediates(bool cc, const Particle & in,
					      ParticleVector & out) {
  // create the particles
  // incoming particle
  ParticleVector external;
  external.push_back(const_ptr_cast<tPPtr>(&in));
  // outgoing
  for(unsigned int ix=0;ix<out.size();++ix)
    external.push_back(out[ix]);
  out.clear();
  // now create the intermediates
  ParticleVector resonance;
  resonance.push_back(external[0]);
  // Lorentz5Momentum pinter;
  for(unsigned ix=1;ix<intermediates_.size();++ix) {
    Lorentz5Momentum pinter;
    for(const int & des : intermediates_[ix].descendents)
      pinter += external[des]->momentum();
    pinter.rescaleMass();
    PPtr respart = (cc&&intermediates_[ix].particle->CC()) ? 
      intermediates_[ix].particle->CC()->produceParticle(pinter) : 
      intermediates_[ix].particle      ->produceParticle(pinter);
    resonance.push_back(respart);
  }
  // set up the mother daughter relations
  for(unsigned int ix=1;ix<intermediates_.size();++ix) {
    resonance[ix]->addChild( intermediates_[ix].children.first <0 ? 
   			     resonance[-intermediates_[ix].children.first ] :
			     external[intermediates_[ix].children.first   ]);
    
    resonance[ix]->addChild( intermediates_[ix].children.second<0 ? 
   			     resonance[-intermediates_[ix].children.second] :
			     external[intermediates_[ix].children.second  ]);
    if(resonance[ix]->dataPtr()->stable())
      resonance[ix]->setLifeLength(Lorentz5Distance());
  }
  // construct the output with the particles in the first step
  out.push_back( intermediates_[0].children.first >0 ?
		 external[intermediates_[0].children.first] :
		 resonance[-intermediates_[0].children.first ]);
  out.push_back( intermediates_[0].children.second>0 ?
		 external[intermediates_[0].children.second] :
		 resonance[-intermediates_[0].children.second]);
}
 
ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer PhaseSpaceChannel::createDiagram() const {
  assert(mode_->incoming().second);
  // create the diagram with incoming and s chnnel particle
  ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer diag =
    new_ptr((Tree2toNDiagram(2), mode_->incoming().first,
	     mode_->incoming().second,1,intermediates_[0].particle));
  map<int,int> children;
  map<unsigned int,int> res;
  int ires=3;
  res[0]=3;
  // add the intermediates
  for(unsigned int ix=0;ix<intermediates_.size();++ix) {
    if(intermediates_[ix].children.first>0)
      children[intermediates_[ix].children.first]= res[ix];
    else {
      ires+=1;
      diag = new_ptr((*diag,res[ix],intermediates_[-intermediates_[ix].children.first].particle));
      res[-intermediates_[ix].children.first]=ires;
    }
    if(intermediates_[ix].children.second>0)
      children[intermediates_[ix].children.second]= res[ix];
    else {
      ires+=1;
      diag = new_ptr((*diag,res[ix],intermediates_[-intermediates_[ix].children.second].particle));
      res[-intermediates_[ix].children.second]=ires;
    }
  }
  // add the children in the corret order
  for(map<int,int>::const_iterator it=children.begin();it!=children.end();++it) {
    diag = new_ptr((*diag,it->second,mode_->outgoing()[it->first-1]));
  }
  return diag;
}

bool PhaseSpaceChannel::checkKinematics() {
  // recalculate the masses and widths of the resonances
  for(auto inter : intermediates_) {
    inter.mass2 = sqr(inter.particle->mass());
    inter.mWidth = inter.particle->mass()*inter.particle->width();
  }
  Energy massmax = mode_->incoming().first->massMax();
  for(tPDPtr out : mode_->outgoing()) 
    massmax -= out->massMin();
  for(unsigned int ix=1;ix<intermediates_.size();++ix) {
    Energy massmin(ZERO);
    for(const int & iloc : intermediates_[ix].descendents) 
      massmin += mode_->outgoing()[iloc-1]->massMin();
    if(intermediates_[ix].mass2>=sqr(massmin)&&
       intermediates_[ix].mass2<=sqr(massmax+massmin)&&
       intermediates_[ix].mWidth==ZERO)
      return false;
  }
  return true;
}
