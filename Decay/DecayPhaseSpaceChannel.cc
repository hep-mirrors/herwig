// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayPhaseSpaceChannel class.
//
// Author: Peter Richardson
// 

#include "DecayPhaseSpaceChannel.h"
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

DecayPhaseSpaceChannel::~DecayPhaseSpaceChannel() {}
  
void DecayPhaseSpaceChannel::persistentOutput(PersistentOStream & os) const 
{
  os << _intpart << _jactype << _intmass << _intwidth << _intpower 
     << _intdau1 << _intdau2 << _extpart << _extmass << _intext;
}
  
void DecayPhaseSpaceChannel::persistentInput(PersistentIStream & is, int) 
{
  is >> _intpart >> _jactype >> _intmass >> _intwidth >> _intpower
     >> _intdau1 >> _intdau2 >> _extpart >> _extmass >> _intext;
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
  
  static RefVector<DecayPhaseSpaceChannel,ParticleData> interfaceExternalParticles
    ("ExternalParticles",
     "The external particles in the decay chain.",
       &DecayPhaseSpaceChannel::_extpart, 0, false, false, true, false);
  
  static ClassDocumentation<DecayPhaseSpaceChannel> documentation
    ("The \\classname{DecayPhaseSpaceChannel} class defines a channel"
     " for the multichannel integration of the phase space for a decay.");
  
}

// generate the momenta of the external particles
vector<Lorentz5Momentum> 
DecayPhaseSpaceChannel::generateMomenta(const Lorentz5Momentum & pin)
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
  pexternal.resize(_extpart.size());
  pinter.resize(_intpart.size());
  // masses of the external particles
  vector<Energy> massext; massext.push_back(pin[5]);
  // masses of the intermediate particles
  vector<Energy> massint; massint.resize(_intpart.size());massint[0]=pin[5];
  // first generate the masses for the external particles
  for(ix=1;ix<_extpart.size();++ix)
    {massext.push_back(_extpart[ix]->generateMass());}
  // generate all the decays in the chain
  Energy lower,upper,lowerb[2];
  for(ix=0;ix<_intpart.size();++ix)
    {
      idau[0] = abs(int(_intdau1[ix]));
      idau[1] = abs(int(_intdau2[ix]));
      // if both decay products off-shell
      if(_intdau1[ix]<0&&_intdau2[ix]<0)
	{
	  // lower limits on the masses of the two resonances
	  for(iy=0;iy<2;++iy)
	    {
	      lowerb[iy]=0.;
	      for(iz=0;iz<_intext[idau[iy]].size();++iz)
		{lowerb[0]+=massext[_intext[idau[iy]][iz]];}
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
	      upper - massint[ix]-massint[idau[1]];
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
  // return the external momentm
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
  Energy2 mass2; Energy temp;
  for(ix=0;ix<_intpart.size();++ix)
    {
      mass2=0.;
      for(int iy=0;iy<4;++iy)
	{
	  temp=0.;
	  for(iz=0;iz<_intext[ix].size();++iz){temp+=output[_intext[ix][iz]][iy];}
	  if(iy<3){mass2-=temp*temp;}
	  else{mass2+=temp*temp;}
	}
      intmass2.push_back(mass2);intmass.push_back(sqrt(mass2));
    }
  // calculate the terms for each of the decays
  Energy lower,upper,lowerb[2];
  for(ix=0;ix<_intpart.size();++ix)
    {
      idau[0] = abs(int(_intdau1[ix]));
      idau[1] = abs(int(_intdau2[ix]));
      // if both decay products off-shell
      if(_intdau1[ix]<0&&_intdau2[ix]<0)
	{
	  // lower limits on the masses of the two resonances
	  for(iy=0;iy<2;++iy)
	    {
	      lowerb[iy]=0.;
	      for(iz=0;iz<_intext[idau[iy]].size();++iz)
		{lowerb[0]+=output[_intext[idau[iy]][iz]][5];}
	      // randomize the order
	      if(CurrentGenerator::current().rnd()<0.5)
		{
		  // contribution of first resonance
		  upper = intmass[ix]-lowerb[1];
		  lower = lowerb[0];
		  wgt *=massWeight(idau[0],intmass[idau[0]],lower,upper);
		  // contribution of second resonance
		  upper = intmass[ix]-intmass[idau[0]];
		  lower = lowerb[1];
		  wgt *=massWeight(idau[1],intmass[idau[1]],lower,upper);
		}
	      else
		{
		  upper = intmass[ix]-lowerb[0];
		  lower = lowerb[1];
		  wgt *=massWeight(idau[1],intmass[idau[1]],lower,upper);
		  upper - intmass[ix]-intmass[idau[1]];
		  lower = lowerb[0];
		  wgt *=massWeight(idau[0],intmass[idau[0]],lower,upper);
		}
	      // factor for the kinematics
	      temp = Kinematics::CMMomentum(intmass[ix],intmass[idau[0]],
					    intmass[idau[1]]);
	      wgt *= intmass[ix]*8.*pi*pi/temp;
	    }
	}
      // only first off-shell
      else if(_intdau1[ix]<0)
	{
	  // compute the limits of integration
	  upper = intmass[ix]-output[idau[1]][5];
	  lower = 0.;
	  for(iy=0;iy<_intext[idau[0]].size();++iy)
	    {lower+=output[_intext[idau[0]][iy]][5];}
	  wgt *=massWeight(idau[0],intmass[idau[0]],lower,upper);
	  temp = Kinematics::CMMomentum(intmass[ix],intmass[idau[0]],
					output[idau[1]][5]);
	  wgt *= intmass[ix]*8.*pi*pi/temp;
	}
      // only second off-shell
      else if(_intdau2[ix]<0)
	{
	    // compute the limits of integration
	  upper = intmass[ix]-output[idau[0]][5]; 
	  lower = 0;
	  for(iy=0;iy<_intext[idau[1]].size();++iy)
	    {lower+=output[_intext[idau[1]][iy]][5];}
	  wgt *=massWeight(idau[1],intmass[idau[1]],lower,upper);
	  temp = Kinematics::CMMomentum(intmass[ix],intmass[idau[1]],
					output[idau[0]][5]);
	  wgt *=intmass[ix]*8.*pi*pi/temp;
	}
      // both on-shell
      else
	{
	  temp = Kinematics::CMMomentum(intmass[ix],output[idau[1]][5],
					output[idau[0]][5]);
	  wgt *=intmass[ix]*8.*pi*pi/temp;
	}
    }
  // finally the overall factor
  wgt *= intmass[0]/pi;
  // return the answer
  return wgt;
}

// output the information to a stream
ostream & Herwig::operator<<(ostream & os, const DecayPhaseSpaceChannel & channel)
{
  // output of the external particles
  os << "Channel for the decay of " << channel._extpart[0]->PDGName() 
     << " -> ";
  for(unsigned int ix=1;ix<channel._extpart.size();++ix)
    {os << channel._extpart[ix]->PDGName() << " ";}
  os << endl;
  os << "Decay proceeds in following steps ";
  for(unsigned int ix=0;ix<channel._intpart.size();++ix)
    {
      os << channel._intpart[ix]->PDGName() << " -> ";
      if(channel._intdau1[ix]>0)
	{os << channel._extpart[int(channel._intdau1[ix])]->PDGName()  
	    << "(" << channel._intdau1[ix]<< ") ";}
      else
	{os << channel._intpart[-int(channel._intdau1[ix])]->PDGName() 
	    << "(" << channel._intdau1[ix]<< ") ";}
      if(channel._intdau2[ix]>0)
	{os << channel._extpart[int(channel._intdau2[ix])]->PDGName()  
	    << "(" <<channel._intdau2[ix] << ") ";}
      else
	{os << channel._intpart[-int(channel._intdau2[ix])]->PDGName() 
	    << "(" <<channel._intdau2[ix] << ") ";}
      os << endl;
    }
  return os;
}

// doinit method  
void DecayPhaseSpaceChannel::doinit() throw(InitException) {
  Interfaced::doinit();
  // masses and widths of the intermediate particles
  for(unsigned int ix=0;ix<_intpart.size();++ix)
    {
      _intmass.push_back(_intpart[ix]->mass());
      _intwidth.push_back(_intpart[ix]->width());
    }
  // masses of the external particles
  for(unsigned int ix=0;ix<_extpart.size();++ix)
    {
      _extmass.push_back(_extpart[ix]->mass());
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
	{temp.push_back(int(_intdau1[ix]));}
      else
	{
	  iy = -int(_intdau1[ix]);
	  vector<int>::iterator istart=_intext[iy].begin();
	  vector<int>::iterator iend=_intext[iy].end();
	  for(;istart!=iend;++istart)
	    {temp.push_back(*istart);}
	}
      // add the second daughter
      if(_intdau2[ix]>=0)
	{temp.push_back(int(_intdau2[ix]));}
      else
	{
	  iy = -int(_intdau2[ix]);
	  vector<int>::iterator istart=_intext[iy].begin();
	  vector<int>::iterator iend=_intext[iy].end();
	  for(;istart!=iend;++istart)
	    {temp.push_back(*istart);}
	}
      _intext[ix]=temp;
    }
}

// generate the final-state particles including the intermediate resonances
void DecayPhaseSpaceChannel::generateIntermediates(const Particle & in,
							  ParticleVector & out)
{
  // integers for loops
  unsigned int ix,iy,iz;
  // momenta fo the external particles
  vector<Lorentz5Momentum> output(out.size()+1);output[0]=in.momentum();
  vector<tcSpinPtr> outspin(out.size()+1);outspin[0]=in.spinInfo();
  // map the external particles onto the internal ones
  bool found; int id;
  vector<bool> done(_extpart.size(),false);done[0]=true;
  if(in.id()==_extpart[0]->id())
    {
      for(iy=0;iy<out.size();++iy)
	{
	  found = false;
	  id=out[iy]->id();
	  ix=0;
	  do
	    {
	      ++ix;
	      if(!done[ix]&&_extpart[ix]->id()==id)
		{
		  done[ix]=true;
		  output[ix]=out[iy]->momentum();
		  outspin[ix]=out[iy]->spinInfo();
		  found=true;
		}
	    }
	  while(!found&&ix<_extpart.size()-1);
	  if(!found)
	    {
	      cerr << "error in DecayPhaseSpaceChannel::generateIntermediates "
		   << "can't find particle" 
		   << id << endl;
	      for(unsigned int ix=0;ix<_extpart.size();++ix)
		{cerr << "particle in decayer " << ix 
		      << _extpart[ix]->id()<<endl;}
	      for(unsigned int ix=0;ix<out.size();++ix)
		{cerr << "external particle " << ix 
		      << out[ix]->id()<<endl;}
	    }
	}
    }
  else
    {
      tcPDPtr anti;
      for(iy=0;iy<out.size();++iy)
	{
	  found = false;
	  anti=(out[iy]->dataPtr())->CC();
	  if(anti){id=anti->id();}
	  else{id=out[iy]->id();}
	  ix=0;
	  do 
	    {
	      ++ix;
	      if(!done[ix]&&_extpart[ix]->id()==id)
		{
		  done[ix]=true;
		  output[ix]=out[iy]->momentum();
		  outspin[ix]=out[iy]->spinInfo();
		  found=true;
		}
	    }
	  while(!found&&ix<_extpart.size()-1);
	  if(!found)
	    {
	      cerr << "error in DecayPhaseSpaceChannel::generateIntermediates " 
		   << "can't find particle" 
		   << id << endl;
	      for(unsigned int ix=0;ix<_extpart.size();++ix)
		{cerr << "particle in decayer " << ix 
		      << _extpart[ix]->id()<<endl;}
	      for(unsigned int ix=0;ix<out.size();++ix)
		{cerr << "external particle " << ix 
		      << out[ix]->id()<<endl;}
	    }
	}
    }
  // create the particles
  bool anti=(in.id()!=_extpart[0]->id());
  PPtr respart;
  tcPDPtr parttemp;
  // final-state
  ParticleVector external;
  external.push_back(const_ptr_cast<tPPtr>(&in));
  for( ix=1;ix<_extpart.size();++ix)
    {
      if(anti&&_extpart[ix]->CC()){parttemp=_extpart[ix]->CC();}
      else{parttemp=_extpart[ix];}
      respart=new_ptr(Particle(parttemp));
      respart->set5Momentum(output[ix]);
      respart->spinInfo(const_ptr_cast<tSpinPtr>(outspin[ix]));
      external.push_back(respart);
    }
  out.resize(0);
  // intermediates
  ParticleVector resonance; resonance.push_back(const_ptr_cast<tPPtr>(&in));
  for(ix=1;ix<_intpart.size();++ix)
    {
      Lorentz5Momentum temp(0.,0.,0.,0.);
      for(iz=0;iz<_intext[ix].size();++iz){temp+=output[_intext[ix][iz]];}
      temp.rescaleMass();
      if(anti&&_intpart[ix]->CC()){parttemp=_intpart[ix]->CC();}
      else{parttemp=_intpart[ix];}
      respart=new_ptr(Particle(parttemp));
      respart->set5Momentum(temp);
      resonance.push_back(respart);
    }
  // set up the mother daughter relations
  for(ix=1;ix<_intpart.size();++ix)
    {
      if(_intdau1[ix]<0)
	{resonance[ix]->addChild(resonance[-int(_intdau1[ix])]);}
      else
	{resonance[ix]->addChild(external[ int(_intdau1[ix])]);}
      if(_intdau2[ix]<0)
	{resonance[ix]->addChild(resonance[-int(_intdau2[ix])]);}
      else
	{resonance[ix]->addChild(external[ int(_intdau2[ix])]);}
    }
  // construct the output with the particles in the first step
  if(_intdau1[0]>0){out.push_back(external[int(_intdau1[0])]);}
  else{out.push_back(resonance[-int(_intdau1[0])]);}
  if(_intdau2[0]>0){out.push_back(external[int(_intdau2[0])]);}
  else{out.push_back(resonance[-int(_intdau2[0])]);}
}

}

