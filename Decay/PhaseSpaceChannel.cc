// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhaseSpaceChannel class.
//

#include "PhaseSpaceChannel.h"
#include "PhaseSpaceMode.h"

/** 
 *  Constructor with incoming particles
 */
PhaseSpaceChannel::PhaseSpaceChannel(tPhaseSpaceModePtr inm) : mode_(inm), weight_(1.) {
  if(!inm->incoming().second)
     intermediates_.push_back(PhaseSpaceResonance(inm->incoming().first));
}

void PhaseSpaceChannel::init(tPhaseSpaceModePtr mode) {
  mode_=mode;
  // find the descentents
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
	mode->outgoing_[ext]->mass() : mode->outgoing_[ext]->massMin();
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
  if(!channel.mode_->incoming().second)
    os << "Channel for the decay of " << channel.mode_->incoming().first->PDGName() << " -> ";
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
  assert(false);
  // _intmass.clear();
  // _intwidth.clear();
  // _intmass2.clear();
  // _intmwidth.clear();
  // // masses and widths of the intermediate particles
  // for(unsigned int ix=0;ix<_intpart.size();++ix) {
  //   _intmass.push_back(_intpart[ix]->mass());
  //   _intwidth.push_back(_intpart[ix]->width());
  //   _intmass2.push_back(_intpart[ix]->mass()*_intpart[ix]->mass());
  //   _intmwidth.push_back(_intpart[ix]->mass()*_intpart[ix]->width());
  // }
  // // ensure intermediates either have the width set, or
  // // can't possibly be on-shell
  // Energy massmax = _mode->externalParticles(0)->massMax();
  // for(unsigned int ix=1;ix<_mode->numberofParticles();++ix) 
  //   massmax -= _mode->externalParticles(ix)->massMin();
  // for(unsigned int ix=0;ix<_intpart.size();++ix) {
  //   if(_intwidth[ix]==0.*MeV && ix>0 && _jactype[ix]==0 ) {
  //     Energy massmin(0.*GeV);
  //     for(unsigned int iy=0;iy<_intext[ix].size();++iy)
  // 	massmin += _mode->externalParticles(_intext[ix][iy])->massMin();
  //     // check if can be on-shell
  //     if(_intmass[ix]>=massmin&&_intmass[ix]<=massmax+massmin) {
  // 	string modeout;
  // 	for(unsigned int iy=0;iy<_mode->numberofParticles();++iy) {
  // 	  modeout += _mode->externalParticles(iy)->PDGName() + " ";
  // 	}
  // 	throw Exception() << "Width zero for " << _intpart[ix]->PDGName()
  // 			  << " in PhaseSpaceChannel::doinitrun() "
  // 			  << modeout
  // 			  << Exception::runerror;
  //     }
  //   }
  // }
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
      Energy lowerb[2];
       // lower limits on the masses of the two resonances
      for(unsigned int iy=0;iy<2;++iy) {
   	lowerb[iy]=ZERO;
   	bool massless=true;
   	for(const int & des : intermediates_[idau[iy]].descendents) {
   	  if(massext[des-1]!=ZERO) massless = false;
   	  lowerb[iy]+=massext[des-1];
  	}
  // 	if(massless) lowerb[iy] = _mode->epsilonPS(); 
      }
      // randomize the order
      if(UseRandom::rnd()<0.5) {
  	// mass of the first resonance
  	Energy upper = massint[ix]-lowerb[1];
  	Energy lower = lowerb[0];
   	massint[idau[0]]=generateMass(intermediates_[idau[0]],lower,upper);
   	// mass of the second resonance
   	upper = massint[ix]-massint[idau[0]];
   	lower = lowerb[1];
   	massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper);
      }
      else {
  	// mass of the second resonance
  	Energy upper = massint[ix]-lowerb[0];
  	Energy lower = lowerb[1];
   	massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper);
   	// mass of the first resonance
   	upper = massint[ix]-massint[idau[1]];
   	lower = lowerb[0];
   	massint[idau[0]]=generateMass(intermediates_[idau[0]],lower,upper);
      }
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massint[idau[1]],
   		   pinter[idau[0]],pinter[idau[1]]);
    }
    // only first off-shell
    else if(intermediates_[ix].children.first<0) {
      // compute the limits of integration
      Energy upper = massint[ix]-massext[idau[1]-1];
      Energy lower = ZERO;
      bool massless=true;
      for(const int & des : intermediates_[idau[0]].descendents) {
   	if(massext[des-1]!=ZERO) massless = false;
   	lower+=massext[des-1];
      }
      //     if(massless) lower = _mode->epsilonPS();
      massint[idau[0]] = generateMass(intermediates_[idau[0]],lower,upper);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massext[idau[1]-1], 
		   pinter[idau[0]],pexternal[idau[1]-1]);
    }
    // only second off-shell
    else if(intermediates_[ix].children.second<0) {
      // compute the limits of integration
      Energy upper = massint[ix]-massext[idau[0]];
      Energy lower = ZERO;
      bool massless=true;
      for(const int & des : intermediates_[idau[1]].descendents) {	
 	if(massext[des-1]!=ZERO) massless = false;
 	lower+=massext[des-1];
      }
      //     if(massless) lower = _mode->epsilonPS();
      massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper);
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
  double ctheta,phi;
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
    cerr << "testing in mass weight " << lower/GeV << " " << upper/GeV << "\n";
    assert(false);
    throw PhaseSpaceError() << "PhaseSpaceChannel::massWeight not allowed " 
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
				       Energy lower,Energy upper) const {
  static const Energy eps=1e-9*MeV;
  if(lower<eps) lower=eps;
  Energy mass=ZERO;
  if(lower>upper) throw PhaseSpaceError() << "PhaseSpaceChannel::generateMass"
					  << " not allowed" 
					  << Exception::eventerror;
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
      double rho = rhomin+rhomax*UseRandom::rnd();
      Energy2 mass2 = max(lower2,min(upper2,res.mass2+res.mWidth*tan(rho)));
      if(mass2<ZERO) mass2 = ZERO;
      mass = sqrt(mass2);
    }
    else {
      mass = sqrt(res.mass2+
  		  (sqr(lower)-res.mass2)*(sqr(upper)-res.mass2)/
  		  (sqr(lower)-res.mass2-UseRandom::rnd()*(sqr(lower)-sqr(upper))));
    }
  }
  // use a power-law
  else if(res.jacobian == PhaseSpaceResonance::Power) {
    double rhomin = pow(sqr(lower/MeV),res.power+1.);
    double rhomax = pow(sqr(upper/MeV),res.power+1.)-rhomin;
    double rho = rhomin+rhomax*UseRandom::rnd();
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
      Energy upper = intmass[ix]-output[idau[0]].mass(); 
      Energy lower = ZERO;
      for(const int & des : intermediates_[idau[1]].descendents) {
	lower += output[des-1].mass();
      }
      wgt *=scale*massWeight(intermediates_[idau[1]],intmass[idau[1]],lower,upper);
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[1]],
						 output[idau[0]].mass());
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
