// -*- C++ -*-
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
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Utilities/Kinematics.h"
#include <ThePEG/Helicity/SpinInfo.h>

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::SpinInfo;

DecayPhaseSpaceChannel::DecayPhaseSpaceChannel(DecayPhaseSpaceModePtr in){_mode=in;}

DecayPhaseSpaceChannel::~DecayPhaseSpaceChannel() {}
  
void DecayPhaseSpaceChannel::persistentOutput(PersistentOStream & os) const 
{
  os << _intpart << _jactype << _intmass << _intwidth 
     << _intmass2 << _intmwidth << _intpower 
     << _intdau1 << _intdau2 << _intext << _mode;
}
  
void DecayPhaseSpaceChannel::persistentInput(PersistentIStream & is, int) 
{
  is >> _intpart >> _jactype >> _intmass >> _intwidth 
     >> _intmass2 >> _intmwidth >> _intpower
     >> _intdau1 >> _intdau2 >> _intext >> _mode;
}
  
  
ClassDescription<DecayPhaseSpaceChannel> DecayPhaseSpaceChannel::initDecayPhaseSpaceChannel;
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
    ("The \\classname{DecayPhaseSpaceChannel} class defines a channel"
     " for the multichannel integration of the phase space for a decay.");
  
}

// generate the momenta of the external particles
vector<Lorentz5Momentum> 
DecayPhaseSpaceChannel::generateMomenta(const Lorentz5Momentum & pin,
					vector<Energy> massext)
{
  // integers for loops
  unsigned int ix,iy,idau[2],iz;
  double ctheta,phi;
  // storage of the momenta of the external particles
  vector<Lorentz5Momentum> pexternal;
  // and the internal particles
  vector<Lorentz5Momentum> pinter; 
  // copy the momentum of the incoming particle
  pexternal.push_back(pin); pinter.push_back(pin);
  pexternal.resize(_mode->numberofParticles());
  pinter.resize(_intpart.size());
  // masses of the intermediate particles
  vector<Energy> massint; massint.resize(_intpart.size());massint[0]=pin[5];
  // generate all the decays in the chain
  Energy lower,upper,lowerb[2];
  for(ix=0;ix<_intpart.size();++ix)
    {
      idau[0] = abs(_intdau1[ix]);
      idau[1] = abs(_intdau2[ix]);
      // if both decay products off-shell
      if(_intdau1[ix]<0&&_intdau2[ix]<0)
	{
	  // lower limits on the masses of the two resonances
	  for(iy=0;iy<2;++iy)
	    {
	      lowerb[iy]=0.;
	      for(iz=0;iz<_intext[idau[iy]].size();++iz)
		{lowerb[iy]+=massext[_intext[idau[iy]][iz]];}
	    }
	  // randomize the order
	  if(CurrentGenerator::current().rnd()<0.5)
	    {
	      // mass of the first resonance
	      upper = massint[ix]-lowerb[1];
	      lower = lowerb[0];
	      massint[idau[0]]=generateMass(idau[0],lower,upper);
	      // mass of the second resonance
	      upper = massint[ix]-massint[idau[0]];
	      lower = lowerb[1];
	      massint[idau[1]]=generateMass(idau[1],lower,upper);
	    }
	  else
	    {
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
	  Kinematics::generateAngles(ctheta,phi);
	  Kinematics::twoBodyDecay(pinter[ix],massint[idau[0]],massint[idau[1]], 
				   ctheta,phi,pinter[idau[0]],pinter[idau[1]]); 
	}
      // only first off-shell
      else if(_intdau1[ix]<0)
	{
	  // compute the limits of integration
	  upper = massint[ix]-massext[idau[1]];
	  lower = 0.;
	  for(iy=0;iy<_intext[idau[0]].size();++iy)
	    {lower+=massext[_intext[idau[0]][iy]];}
	  massint[idau[0]]=generateMass(idau[0],lower,upper);
	  // generate the momenta of the decay products
	  Kinematics::generateAngles(ctheta,phi);
	  Kinematics::twoBodyDecay(pinter[ix],massint[idau[0]],massext[idau[1]], 
				   ctheta,phi,pinter[idau[0]],pexternal[idau[1]]);
	}
      // only second off-shell
      else if(_intdau2[ix]<0)
	{
	  // compute the limits of integration
	  upper = massint[ix]-massext[idau[0]];
	  lower = 0;
	  for(iy=0;iy<_intext[idau[1]].size();++iy)
	    {lower+=massext[_intext[idau[1]][iy]];}
	  massint[idau[1]]=generateMass(idau[1],lower,upper);
	  // generate the momenta of the decay products
	  Kinematics::generateAngles(ctheta,phi);
	  Kinematics::twoBodyDecay(pinter[ix],massext[idau[0]],massint[idau[1]], 
				   ctheta,phi,pexternal[idau[0]],pinter[idau[1]]); 
	}
      // both on-shell
      else
	{
	  // generate the momenta of the decay products
	  Kinematics::generateAngles(ctheta,phi);
	  Kinematics::twoBodyDecay(pinter[ix],massext[idau[0]],massext[idau[1]], 
				   ctheta,phi,pexternal[idau[0]],pexternal[idau[1]]);
	}
    }
  // return the external momenta
  return pexternal;
}

// generate the weight for this channel given a phase space configuration
double DecayPhaseSpaceChannel::generateWeight(const vector<Lorentz5Momentum> & output)
{
  // integers for loops
  unsigned int ix,iy,idau[2],iz;
  // include the prefactor due to the weight of the channel
  double wgt=1.;
  // work out the masses of the intermediate particles
  vector<Energy2> intmass2; vector<Energy> intmass;
  Lorentz5Momentum pinter;
  for(ix=0;ix<_intpart.size();++ix)
    {
      pinter=output[_intext[ix][0]];
      for(iz=1;iz<_intext[ix].size();++iz){pinter+=output[_intext[ix][iz]];}
      pinter.rescaleMass();
      intmass.push_back(pinter.mass());
      intmass2.push_back(intmass[ix]*intmass[ix]);
    }
  Energy2 scale(intmass2[0]);
  // calculate the terms for each of the decays
  Energy lower,upper,lowerb[2];
  for(ix=0;ix<_intpart.size();++ix)
    {
      idau[0] = abs(_intdau1[ix]);
      idau[1] = abs(_intdau2[ix]);
      // if both decay products off-shell
      Energy pcm;
      if(_intdau1[ix]<0&&_intdau2[ix]<0)
	{
	  // lower limits on the masses of the two resonances
	  for(iy=0;iy<2;++iy)
	    {
	      lowerb[iy]=0.;
	      for(iz=0;iz<_intext[idau[iy]].size();++iz)
		{lowerb[iy]+=output[_intext[idau[iy]][iz]][5];}
	    }
	  // randomize the order
	  if(CurrentGenerator::current().rnd()<0.5)
	    {
	      // contribution of first resonance
	      upper = intmass[ix]-lowerb[1];
	      lower = lowerb[0];
	      wgt *=scale*massWeight(idau[0],intmass[idau[0]],lower,upper);
	      // contribution of second resonance
	      upper = intmass[ix]-intmass[idau[0]];
	      lower = lowerb[1];
	      wgt *=scale*massWeight(idau[1],intmass[idau[1]],lower,upper);
	    }
	  else
	    {
	      upper = intmass[ix]-lowerb[0];
	      lower = lowerb[1];
	      wgt *=scale*massWeight(idau[1],intmass[idau[1]],lower,upper);
	      upper = intmass[ix]-intmass[idau[1]];
	      lower = lowerb[0];
	      wgt *=scale*massWeight(idau[0],intmass[idau[0]],lower,upper);
	    }
	  // factor for the kinematics
	  pcm = Kinematics::CMMomentum(intmass[ix],intmass[idau[0]],
				       intmass[idau[1]]);
	  wgt *= intmass[ix]*8.*pi*pi/pcm;
	}
      // only first off-shell
      else if(_intdau1[ix]<0)
	{
	  // compute the limits of integration
	  upper = intmass[ix]-output[idau[1]][5];
	  lower = 0.;
	  for(iy=0;iy<_intext[idau[0]].size();++iy)
	    {lower+=output[_intext[idau[0]][iy]][5];}
	  wgt *=scale*massWeight(idau[0],intmass[idau[0]],lower,upper);
	  pcm = Kinematics::CMMomentum(intmass[ix],intmass[idau[0]],
				       output[idau[1]][5]);
	  wgt *= intmass[ix]*8.*pi*pi/pcm;
	}
      // only second off-shell
      else if(_intdau2[ix]<0)
	{
	    // compute the limits of integration
	  upper = intmass[ix]-output[idau[0]][5]; 
	  lower = 0;
	  for(iy=0;iy<_intext[idau[1]].size();++iy)
	    {lower+=output[_intext[idau[1]][iy]][5];}
	  wgt *=scale*massWeight(idau[1],intmass[idau[1]],lower,upper);
	  pcm = Kinematics::CMMomentum(intmass[ix],intmass[idau[1]],
				       output[idau[0]][5]);
	  wgt *=intmass[ix]*8.*pi*pi/pcm;
	}
      // both on-shell
      else
	{
	  pcm = Kinematics::CMMomentum(intmass[ix],output[idau[1]][5],
				       output[idau[0]][5]);
	  wgt *=intmass[ix]*8.*pi*pi/pcm;
	}
    }
  // finally the overall factor
  wgt /= pi;
  // return the answer
  return wgt;
}

// output the information to a stream
ostream & operator<<(ostream & os, const DecayPhaseSpaceChannel & channel)
{
  // output of the external particles
  os << "Channel for the decay of " << channel._mode->externalParticles(0)->PDGName() 
     << " -> ";
  for(unsigned int ix=1;ix<channel._mode->numberofParticles();++ix)
    {os << channel._mode->externalParticles(ix)->PDGName() << " ";}
  os << endl;
  os << "Decay proceeds in following steps ";
  for(unsigned int ix=0;ix<channel._intpart.size();++ix)
    {
      os << channel._intpart[ix]->PDGName() << " -> ";
      if(channel._intdau1[ix]>0)
	{os << channel._mode->externalParticles(channel._intdau1[ix])->PDGName()  
	    << "(" << channel._intdau1[ix]<< ") ";}
      else
	{os << channel._intpart[-channel._intdau1[ix]]->PDGName() 
	    << "(" << channel._intdau1[ix]<< ") ";}
      if(channel._intdau2[ix]>0)
	{os << channel._mode->externalParticles(channel._intdau2[ix])->PDGName()  
	    << "(" <<channel._intdau2[ix] << ") ";}
      else
	{os << channel._intpart[-channel._intdau2[ix]]->PDGName() 
	    << "(" <<channel._intdau2[ix] << ") ";}
      os << endl;
    }
  return os;
}

// doinit method  
void DecayPhaseSpaceChannel::doinit() throw(InitException) {
  Interfaced::doinit();
  // check if the mode pointer exists
  if(!_mode){throw InitException() << "DecayPhaseSpaceChannel::doinit() the " 
				   << "channel must have a pointer to a decay mode "
				   << Exception::abortnow;}
  // masses and widths of the intermediate particles
  for(unsigned int ix=0;ix<_intpart.size();++ix)
    {
      _intmass.push_back(_intpart[ix]->mass());
      _intwidth.push_back(_intpart[ix]->width());
      _intmass2.push_back(_intpart[ix]->mass()*_intpart[ix]->mass());
      _intmwidth.push_back(_intpart[ix]->mass()*_intpart[ix]->width());
    }
  // external particles for each intermediate particle
  vector<int> temp;
  _intext.resize(_intpart.size());
  // loop over the intermediate particles
  int ix,iy;
  for(ix=_intpart.size()-1;ix>=0;--ix)
    {
      temp.resize(0);
      // add the first daughter
      if(_intdau1[ix]>=0)
	{temp.push_back(_intdau1[ix]);}
      else
	{
	  iy = -_intdau1[ix];
	  vector<int>::iterator istart=_intext[iy].begin();
	  vector<int>::iterator iend=_intext[iy].end();
	  for(;istart!=iend;++istart)
	    {temp.push_back(*istart);}
	}
      // add the second daughter
      if(_intdau2[ix]>=0)
	{temp.push_back(_intdau2[ix]);}
      else
	{
	  iy = -_intdau2[ix];
	  vector<int>::iterator istart=_intext[iy].begin();
	  vector<int>::iterator iend=_intext[iy].end();
	  for(;istart!=iend;++istart)
	    {temp.push_back(*istart);}
	}
      _intext[ix]=temp;
    }
}

// generate the final-state particles including the intermediate resonances
void DecayPhaseSpaceChannel::generateIntermediates(bool cc, const Particle & in,
						   ParticleVector & out)
{
  // integers for the loops
  unsigned int ix,iz;
  // create the particles
  // incoming particle
  ParticleVector external;
  external.push_back(const_ptr_cast<tPPtr>(&in));
  // outgoing
  for(ix=0;ix<out.size();++ix){external.push_back(out[ix]);}
  out.resize(0);
  // now create the intermediates
  ParticleVector resonance; resonance.push_back(external[0]);
  PPtr respart;
  tcPDPtr parttemp;
  Lorentz5Momentum pinter;
  for(ix=1;ix<_intpart.size();++ix)
    {
      pinter=external[_intext[ix][0]]->momentum();
      for(iz=1;iz<_intext[ix].size();++iz)
	{pinter+=external[_intext[ix][iz]]->momentum();}
      pinter.rescaleMass();
      if(cc&&_intpart[ix]->CC()){respart=new_ptr(Particle(_intpart[ix]->CC()));}
      else{respart=new_ptr(Particle(_intpart[ix]));}
      respart->set5Momentum(pinter);
      resonance.push_back(respart);
    }
  // set up the mother daughter relations
  for(ix=1;ix<_intpart.size();++ix)
    {
      if(_intdau1[ix]<0){resonance[ix]->addChild(resonance[-_intdau1[ix]]);}
      else{resonance[ix]->addChild(external[_intdau1[ix]]);}
      if(_intdau2[ix]<0){resonance[ix]->addChild(resonance[-_intdau2[ix]]);}
      else{resonance[ix]->addChild(external[_intdau2[ix]]);}
    }
  // construct the output with the particles in the first step
  if(_intdau1[0]>0){out.push_back(external[_intdau1[0]]);}
  else{out.push_back(resonance[-_intdau1[0]]);}
  if(_intdau2[0]>0){out.push_back(external[_intdau2[0]]);}
  else{out.push_back(resonance[-_intdau2[0]]);}
}

}

