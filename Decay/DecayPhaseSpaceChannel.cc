// -*- C++ -*-
//
// DecayPhaseSpaceChannel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayPhaseSpaceChannel class.
//
// Author: Peter Richardson
// 

#include "DecayPhaseSpaceChannel.h"
#include "DecayPhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"

#include <ThePEG/Repository/CurrentGenerator.h>

using namespace Herwig;

  
DecayPhaseSpaceChannel::DecayPhaseSpaceChannel(tcDecayPhaseSpaceModePtr in) 
  : _mode(in) {}

void DecayPhaseSpaceChannel::persistentOutput(PersistentOStream & os) const {
  os << _intpart << _jactype << ounit(_intmass,GeV) << ounit(_intwidth,GeV)
     << ounit(_intmass2,GeV2) << ounit(_intmwidth,GeV2) << _intpower 
     << _intdau1 << _intdau2 << _intext << _mode;
}
  
void DecayPhaseSpaceChannel::persistentInput(PersistentIStream & is, int) {
  is >> _intpart >> _jactype >> iunit(_intmass,GeV) >> iunit(_intwidth,GeV) 
     >> iunit(_intmass2,GeV2) >> iunit(_intmwidth,GeV2) >> _intpower
     >> _intdau1 >> _intdau2 >> _intext >> _mode;
}
  
ClassDescription<DecayPhaseSpaceChannel> 
DecayPhaseSpaceChannel::initDecayPhaseSpaceChannel;
// Definition of the static class description member.

void DecayPhaseSpaceChannel::Init() {
    
  static RefVector<DecayPhaseSpaceChannel,ParticleData> interfaceIntermediateParticles
    ("IntermediateParticles",
     "The intermediate particles in the decay chain.",
     &DecayPhaseSpaceChannel::_intpart, 0, false, false, true, false);
  
  static ParVector<DecayPhaseSpaceChannel,int> interfacejactype
    ("Jacobian",
     "The type of Jacobian to use for the intermediate particle",
     &DecayPhaseSpaceChannel::_jactype,
     0, 0, 0, 0, 1, false, false, true);
  
  static ParVector<DecayPhaseSpaceChannel,double> interfaceIntermediatePower
    ("IntermediatePower",
     "The power to use in the Jacobian",
     &DecayPhaseSpaceChannel::_intpower,
     0, 0, 0, -10, 10, false, false, true);
  
  static ParVector<DecayPhaseSpaceChannel,int> interfaceIntermediateDau1
    ("IntermediateDaughter1",
     "First Daughter of the intermediate",
     &DecayPhaseSpaceChannel::_intdau1,
     0, 0, 0, -10, 10, false, false, true);
  
  static ParVector<DecayPhaseSpaceChannel,int> interfaceIntermediateDau2
    ("IntermediateDaughter2",
     "Second Daughter of the intermediate",
     &DecayPhaseSpaceChannel::_intdau2,
     0, 0, 0, -10, 10, false, false, true);
  
  static ClassDocumentation<DecayPhaseSpaceChannel> documentation
    ("The DecayPhaseSpaceChannel class defines a channel"
     " for the multichannel integration of the phase space for a decay.");
  
}

// generate the momenta of the external particles
vector<Lorentz5Momentum> 
DecayPhaseSpaceChannel::generateMomenta(const Lorentz5Momentum & pin,
					const vector<Energy> & massext) {
  // integers for loops
  unsigned int ix,iy,idau[2],iz;
  // storage of the momenta of the external particles
  vector<Lorentz5Momentum> pexternal;
  // and the internal particles
  vector<Lorentz5Momentum> pinter; 
  // copy the momentum of the incoming particle
  pexternal.push_back(pin); pinter.push_back(pin);
  pexternal.resize(_mode->numberofParticles());
  pinter.resize(_intpart.size());
  // masses of the intermediate particles
  vector<Energy> massint(_intpart.size());
  massint[0]=pin.mass();
  // generate all the decays in the chain
  Energy lower,upper,lowerb[2];
  for(ix=0;ix<_intpart.size();++ix) {
    idau[0] = abs(_intdau1[ix]);
    idau[1] = abs(_intdau2[ix]);
    // if both decay products off-shell
    if(_intdau1[ix]<0&&_intdau2[ix]<0) {
      // lower limits on the masses of the two resonances
      for(iy=0;iy<2;++iy) {
	lowerb[iy]=ZERO;
	for(iz=0;iz<_intext[idau[iy]].size();++iz) {
	  lowerb[iy]+=massext[_intext[idau[iy]][iz]];
	}
      }
      // randomize the order
      if(UseRandom::rnd()<0.5) {
	// mass of the first resonance
	upper = massint[ix]-lowerb[1];
	lower = lowerb[0];
	massint[idau[0]]=generateMass(idau[0],lower,upper);
	// mass of the second resonance
	upper = massint[ix]-massint[idau[0]];
	lower = lowerb[1];
	massint[idau[1]]=generateMass(idau[1],lower,upper);
      }
      else {
	// mass of the second resonance
	upper = massint[ix]-lowerb[0];
	lower = lowerb[1];
	massint[idau[1]]=generateMass(idau[1],lower,upper);
	// mass of the first resonance
	upper = massint[ix]-massint[idau[1]];
	lower = lowerb[0];
	massint[idau[0]]=generateMass(idau[0],lower,upper);
      }
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massint[idau[1]],
		   pinter[idau[0]],pinter[idau[1]]);
    }
    // only first off-shell
    else if(_intdau1[ix]<0) {
      // compute the limits of integration
      upper = massint[ix]-massext[idau[1]];
      lower = ZERO;
      for(iy=0;iy<_intext[idau[0]].size();++iy) {
	lower+=massext[_intext[idau[0]][iy]];
      }
      massint[idau[0]]=generateMass(idau[0],lower,upper);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massext[idau[1]], 
		   pinter[idau[0]],pexternal[idau[1]]);
    }
    // only second off-shell
    else if(_intdau2[ix]<0) {
      // compute the limits of integration
      upper = massint[ix]-massext[idau[0]];
      lower = ZERO;
      for(iy=0;iy<_intext[idau[1]].size();++iy) {
	lower+=massext[_intext[idau[1]][iy]];
      }
      massint[idau[1]]=generateMass(idau[1],lower,upper);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massext[idau[0]],massint[idau[1]], 
		   pexternal[idau[0]],pinter[idau[1]]);
    }
    // both on-shell
    else {
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massext[idau[0]],massext[idau[1]], 
		   pexternal[idau[0]],pexternal[idau[1]]);
    }
  }
  // return the external momenta
  return pexternal;
}

// generate the weight for this channel given a phase space configuration
double DecayPhaseSpaceChannel::generateWeight(const vector<Lorentz5Momentum> & output) {
  using Constants::pi;
  // integers for loops
  unsigned int ix,iy,idau[2],iz;
  // include the prefactor due to the weight of the channel
  double wgt=1.;
  // work out the masses of the intermediate particles
  static vector<Energy> intmass;
  intmass.clear();
  Lorentz5Momentum pinter;
  for(ix=0;ix<_intpart.size();++ix) {
    pinter=output[_intext[ix][0]];
    for(iz=1;iz<_intext[ix].size();++iz) pinter+=output[_intext[ix][iz]];
    pinter.rescaleMass();
    intmass.push_back( pinter.mass() );
  }
  Energy2 scale(sqr(intmass[0]));
  // calculate the terms for each of the decays
  Energy lower,upper,lowerb[2];
  for(ix=0;ix<_intpart.size();++ix) {
    idau[0] = abs(_intdau1[ix]);
    idau[1] = abs(_intdau2[ix]);
    // if both decay products off-shell
    Energy pcm;
    if(_intdau1[ix]<0&&_intdau2[ix]<0) {
      // lower limits on the masses of the two resonances
      for(iy=0;iy<2;++iy) {
	lowerb[iy]=ZERO;
	for(iz=0;iz<_intext[idau[iy]].size();++iz)
	  lowerb[iy]+=output[_intext[idau[iy]][iz]].mass();
      }
      // undo effect of randomising
      // weight for the first order
      // contribution of first resonance
      upper = intmass[ix]-lowerb[1];
      lower = lowerb[0];
      InvEnergy2 wgta=massWeight(idau[0],intmass[idau[0]],lower,upper);
      // contribution of second resonance
      upper = intmass[ix]-intmass[idau[0]];
      lower = lowerb[1];
      InvEnergy4 wgta2 = wgta*massWeight(idau[1],intmass[idau[1]],lower,upper);
      // weight for the second order
      upper = intmass[ix]-lowerb[0];
      lower = lowerb[1];
      InvEnergy2 wgtb=massWeight(idau[1],intmass[idau[1]],lower,upper);
      upper = intmass[ix]-intmass[idau[1]];
      lower = lowerb[0];
      InvEnergy4 wgtb2=wgtb*massWeight(idau[0],intmass[idau[0]],lower,upper);
      // weight factor
      wgt *=0.5*sqr(scale)*(wgta2+wgtb2);
      // factor for the kinematics
      pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[0]],
				   intmass[idau[1]]);
      wgt *= intmass[ix]*8.*pi*pi/pcm;
    }
    // only first off-shell
    else if(_intdau1[ix]<0) {
      // compute the limits of integration
      upper = intmass[ix]-output[idau[1]].mass();
      lower = ZERO;
      for(iy=0;iy<_intext[idau[0]].size();++iy)
	lower+=output[_intext[idau[0]][iy]].mass();
      wgt *=scale*massWeight(idau[0],intmass[idau[0]],lower,upper);
      pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[0]],
				   output[idau[1]].mass());
      wgt *= intmass[ix]*8.*pi*pi/pcm;
    }
    // only second off-shell
    else if(_intdau2[ix]<0) {
      // compute the limits of integration
      upper = intmass[ix]-output[idau[0]].mass(); 
      lower = ZERO;
      for(iy=0;iy<_intext[idau[1]].size();++iy)
	lower+=output[_intext[idau[1]][iy]].mass();
      wgt *=scale*massWeight(idau[1],intmass[idau[1]],lower,upper);
      pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[1]],
				   output[idau[0]].mass());
      wgt *=intmass[ix]*8.*pi*pi/pcm;
    }
    // both on-shell
    else {
      pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],output[idau[1]].mass(),
					  output[idau[0]].mass());
      if(pcm!=ZERO)
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

// output the information to a stream
ostream & Herwig::operator<<(ostream & os, const DecayPhaseSpaceChannel & channel) {
  // output of the external particles
  os << "Channel for the decay of " << channel._mode->externalParticles(0)->PDGName() 
     << " -> ";
  for(unsigned int ix=1;ix<channel._mode->numberofParticles();++ix)
    os << channel._mode->externalParticles(ix)->PDGName() << " ";
  os << endl;
  os << "Decay proceeds in following steps ";
  for(unsigned int ix=0;ix<channel._intpart.size();++ix) {
    os << channel._intpart[ix]->PDGName() << " -> ";
    if(channel._intdau1[ix]>0) {
      os << channel._mode->externalParticles(channel._intdau1[ix])->PDGName()  
	 << "(" << channel._intdau1[ix]<< ") ";
    }
    else {
      os << channel._intpart[-channel._intdau1[ix]]->PDGName() 
	 << "(" << channel._intdau1[ix]<< ") ";
    }
    if(channel._intdau2[ix]>0) {
      os << channel._mode->externalParticles(channel._intdau2[ix])->PDGName()  
	 << "(" <<channel._intdau2[ix] << ") ";
    }
    else{
      os << channel._intpart[-channel._intdau2[ix]]->PDGName() 
	 << "(" <<channel._intdau2[ix] << ") ";
    }
    os << endl;
  }
  return os;
}

// doinit method  
void DecayPhaseSpaceChannel::doinit() {
  Interfaced::doinit();
  // check if the mode pointer exists
  if(!_mode){throw InitException() << "DecayPhaseSpaceChannel::doinit() the " 
				   << "channel must have a pointer to a decay mode "
				   << Exception::abortnow;}
  // masses and widths of the intermediate particles
  for(unsigned int ix=0;ix<_intpart.size();++ix) {
    _intmass.push_back(_intpart[ix]->mass());
    _intwidth.push_back(_intpart[ix]->width());
    _intmass2.push_back(_intpart[ix]->mass()*_intpart[ix]->mass());
    _intmwidth.push_back(_intpart[ix]->mass()*_intpart[ix]->width());
  }
  // external particles for each intermediate particle
  vector<int> temp;
  _intext.resize(_intpart.size());
  // loop over the intermediate particles
  for(int ix=_intpart.size()-1;ix>=0;--ix) {
    temp.clear();
    // add the first daughter
    if(_intdau1[ix]>=0) {
      temp.push_back(_intdau1[ix]);
    }
    else {
      int iy = -_intdau1[ix];
      vector<int>::iterator istart=_intext[iy].begin();
      vector<int>::iterator iend=_intext[iy].end();
      for(;istart!=iend;++istart) temp.push_back(*istart);
    }
    // add the second daughter
    if(_intdau2[ix]>=0) {
      temp.push_back(_intdau2[ix]);
    }
    else {
      int iy = -_intdau2[ix];
      vector<int>::iterator istart=_intext[iy].begin();
      vector<int>::iterator iend=_intext[iy].end();
      for(;istart!=iend;++istart) temp.push_back(*istart);
    }
    _intext[ix]=temp;
  }
  // ensure intermediates either have the width set, or
  // can't possibly be on-shell
  Energy massmax;
  if(_mode->testOnShell()) {
    massmax = _mode->externalParticles(0)->mass();
    for(unsigned int ix=1;ix<_mode->numberofParticles();++ix) 
      massmax -= _mode->externalParticles(ix)->mass();
  }
  else { 
    massmax = _mode->externalParticles(0)->massMax();
    for(unsigned int ix=1;ix<_mode->numberofParticles();++ix) 
      massmax -= _mode->externalParticles(ix)->massMin();
  }
  for(unsigned int ix=0;ix<_intpart.size();++ix) {
    if(_intwidth[ix]==ZERO && ix>0 && _jactype[ix]==0 ) {
      Energy massmin(ZERO);
      for(unsigned int iy=0;iy<_intext[ix].size();++iy)
	massmin += _mode->testOnShell() ? 
	  _mode->externalParticles(_intext[ix][iy])->mass() :
	  _mode->externalParticles(_intext[ix][iy])->massMin();
      // check if can be on-shell
      if(_intmass[ix]>=massmin&&_intmass[ix]<=massmax+massmin) {
	string modeout;
	for(unsigned int iy=0;iy<_mode->numberofParticles();++iy) {
	  modeout += _mode->externalParticles(iy)->PDGName() + " ";
	}
	throw InitException() << "Width zero for " << _intpart[ix]->PDGName()
			      << " in DecayPhaseSpaceChannel::doinit() "
			      << modeout
			      << Exception::runerror;
      }
    }
  }
}

void DecayPhaseSpaceChannel::doinitrun() {
  Interfaced::doinitrun();
  if(!_mode->testOnShell()) return;
  _intmass.clear();
  _intwidth.clear();
  _intmass2.clear();
  _intmwidth.clear();
  // masses and widths of the intermediate particles
  for(unsigned int ix=0;ix<_intpart.size();++ix) {
    _intmass.push_back(_intpart[ix]->mass());
    _intwidth.push_back(_intpart[ix]->width());
    _intmass2.push_back(_intpart[ix]->mass()*_intpart[ix]->mass());
    _intmwidth.push_back(_intpart[ix]->mass()*_intpart[ix]->width());
  }
  // ensure intermediates either have the width set, or
  // can't possibly be on-shell
  Energy massmax = _mode->externalParticles(0)->massMax();
  for(unsigned int ix=1;ix<_mode->numberofParticles();++ix) 
    massmax -= _mode->externalParticles(ix)->massMin();
  for(unsigned int ix=0;ix<_intpart.size();++ix) {
    if(_intwidth[ix]==0.*MeV && ix>0 && _jactype[ix]==0 ) {
      Energy massmin(0.*GeV);
      for(unsigned int iy=0;iy<_intext[ix].size();++iy)
	massmin += _mode->externalParticles(_intext[ix][iy])->massMin();
      // check if can be on-shell
      if(_intmass[ix]>=massmin&&_intmass[ix]<=massmax+massmin) {
	string modeout;
	for(unsigned int iy=0;iy<_mode->numberofParticles();++iy) {
	  modeout += _mode->externalParticles(iy)->PDGName() + " ";
	}
	throw Exception() << "Width zero for " << _intpart[ix]->PDGName()
			  << " in DecayPhaseSpaceChannel::doinitrun() "
			  << modeout
			  << Exception::runerror;
      }
    }
  }
}

// generate the final-state particles including the intermediate resonances
void DecayPhaseSpaceChannel::generateIntermediates(bool cc, const Particle & in,
						   ParticleVector & out) {
  // integers for the loops
  unsigned int ix,iz;
  // create the particles
  // incoming particle
  ParticleVector external;
  external.push_back(const_ptr_cast<tPPtr>(&in));
  // outgoing
  for(ix=0;ix<out.size();++ix) external.push_back(out[ix]);
  out.clear();
  // now create the intermediates
  ParticleVector resonance;
  resonance.push_back(external[0]);
  PPtr respart;
  tcPDPtr parttemp;
  Lorentz5Momentum pinter;
  for(ix=1;ix<_intpart.size();++ix) {
    pinter=external[_intext[ix][0]]->momentum();
    for(iz=1;iz<_intext[ix].size();++iz) 
      pinter+=external[_intext[ix][iz]]->momentum();
    pinter.rescaleMass();
    respart = (cc&&_intpart[ix]->CC()) ? 
      _intpart[ix]->CC()->produceParticle(pinter) : 
      _intpart[ix]      ->produceParticle(pinter);
    resonance.push_back(respart);
  }
  // set up the mother daughter relations
  for(ix=1;ix<_intpart.size();++ix) {
    resonance[ix]->addChild( _intdau1[ix]<0 ? 
			     resonance[-_intdau1[ix]] : external[_intdau1[ix]]);
    resonance[ix]->addChild( _intdau2[ix]<0 ? 
			     resonance[-_intdau2[ix]] : external[_intdau2[ix]]);
    if(resonance[ix]->dataPtr()->stable())
      resonance[ix]->setLifeLength(Lorentz5Distance());
  }
  // construct the output with the particles in the first step
  out.push_back( _intdau1[0]>0 ? external[_intdau1[0]] : resonance[-_intdau1[0]]);
  out.push_back( _intdau2[0]>0 ? external[_intdau2[0]] : resonance[-_intdau2[0]]);
}
 
double DecayPhaseSpaceChannel::atanhelper_(int ires, Energy limit) {
  return atan2( limit*limit-_intmass2[ires], _intmwidth[ires] );
}

// return the weight for a given resonance
InvEnergy2 DecayPhaseSpaceChannel::massWeight(int ires, Energy moff,
					      Energy lower,Energy upper) {
  InvEnergy2 wgt = ZERO;
  if(lower>upper) {
    throw DecayPhaseSpaceError() << "DecayPhaseSpaceChannel::massWeight not allowed " 
				 << ires << "   " << _intpart[ires]->id() << "   " 
				 << moff/GeV << " " << lower/GeV << " " << upper/GeV << Exception::eventerror;
  } 
  // use a Breit-Wigner 
  if ( _jactype[ires] == 0 ) {
    double rhomin  = atanhelper_(ires,lower);
    double rhomax  = atanhelper_(ires,upper) - rhomin;
    if ( rhomax != 0.0 ) {
      Energy2 moff2=moff*moff-_intmass2[ires];
      wgt = _intmwidth[ires]/rhomax/(moff2*moff2+_intmwidth[ires]*_intmwidth[ires]);
    }
    else {
      wgt = 1./((sqr(upper)-sqr(lower))*sqr(sqr(moff)-_intmass2[ires])/
		(sqr(lower)-_intmass2[ires])/(sqr(upper)-_intmass2[ires]));
    }
  } 
  // power law
  else if(_jactype[ires]==1) {
    double rhomin = pow(sqr(lower/MeV),_intpower[ires]+1.);
    double rhomax = pow(sqr(upper/MeV),_intpower[ires]+1.)-rhomin;
    wgt = (_intpower[ires]+1.)/rhomax*pow(sqr(moff/MeV),_intpower[ires])
      /MeV/MeV;
  }
  else if(_jactype[ires]==2) {
    wgt = 1./Constants::pi/_intmwidth[ires];
  } 
  else {
    throw DecayPhaseSpaceError() << "Unknown type of Jacobian in " 
				 << "DecayPhaseSpaceChannel::massWeight"
				 << Exception::eventerror;
  } 
  return wgt;
}
Energy DecayPhaseSpaceChannel::generateMass(int ires,Energy lower,Energy upper) {
  static const Energy eps=1e-9*MeV;
  if(lower<eps) lower=eps;
  Energy mass=ZERO;
  if(lower>upper) throw DecayPhaseSpaceError() << "DecayPhaseSpaceChannel::generateMass"
					       << " not allowed" 
					       << Exception::eventerror;
  if(abs(lower-upper)/(lower+upper)>2e-10) {
    lower +=1e-10*(lower+upper);
    upper -=1e-10*(lower+upper);
  }
  else 
    return 0.5*(lower+upper);
  // use a Breit-Wigner
  if(_jactype[ires]==0) {
    if(_intmwidth[ires]!=ZERO) {
      Energy2 lower2 = sqr(lower);
      Energy2 upper2 = sqr(upper);

      double rhomin = atan2((lower2 - _intmass2[ires]),_intmwidth[ires]);
      double rhomax = atan2((upper2 - _intmass2[ires]),_intmwidth[ires])-rhomin;
      double rho = rhomin+rhomax*UseRandom::rnd();
      Energy2 mass2 = max(lower2,min(upper2,_intmass2[ires]+_intmwidth[ires]*tan(rho)));
      if(mass2<ZERO) mass2 = ZERO;
      mass = sqrt(mass2);
    }
    else {
      mass = sqrt(_intmass2[ires]+
		  (sqr(lower)-_intmass2[ires])*(sqr(upper)-_intmass2[ires])/
		  (sqr(lower)-_intmass2[ires]-UseRandom::rnd()*(sqr(lower)-sqr(upper))));
    }
  }
  // use a power-law
  else if(_jactype[ires]==1) {
    double rhomin = pow(sqr(lower/MeV),_intpower[ires]+1.);
    double rhomax = pow(sqr(upper/MeV),_intpower[ires]+1.)-rhomin;
    double rho = rhomin+rhomax*UseRandom::rnd();
    mass = pow(rho,0.5/(_intpower[ires]+1.))*MeV;
  }
  else if(_jactype[ires]==2) {
    mass = _intmass[ires];
  } 
  else {
    throw DecayPhaseSpaceError() << "Unknown type of Jacobian in " 
				 << "DecayPhaseSpaceChannel::generateMass" 
				 << Exception::eventerror;
  }
  if(mass<lower+1e-10*(lower+upper))      mass=lower+1e-10*(lower+upper);
  else if(mass>upper-1e-10*(lower+upper)) mass=upper-1e-10*(lower+upper);
  return mass;
}

void DecayPhaseSpaceChannel::twoBodyDecay(const Lorentz5Momentum & p,
					  const Energy m1, const Energy m2,
					  Lorentz5Momentum & p1,
					  Lorentz5Momentum & p2 ) {
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
    throw DecayPhaseSpaceError() << "Two body decay cannot proceed "
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
