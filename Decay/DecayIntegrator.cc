// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayIntegrator class.
//
// Author: Peter Richardson
// 

#include "DecayIntegrator.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Utilities/Timer.h"

namespace Herwig {
using namespace ThePEG;

DecayIntegrator::~DecayIntegrator() {}

bool DecayIntegrator::accept(const DecayMode & dm) const {return false;}

ParticleVector DecayIntegrator::decay(const DecayMode & dm,
				      const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}
  
void DecayIntegrator::persistentOutput(PersistentOStream & os) const {
  os << _channels <<_Initialize << _MaxWeight << _niter << _npoint << _ntry 
     << _channelon << _channelwgts;
  }
  
void DecayIntegrator::persistentInput(PersistentIStream & is, int) {
  is >> _channels >>_Initialize >> _MaxWeight >> _niter >> _npoint >> _ntry
     >> _channelon >> _channelwgts;
}
  
AbstractClassDescription<DecayIntegrator> DecayIntegrator::initDecayIntegrator;
// Definition of the static class description member.
  
void DecayIntegrator::Init() {
    
  static RefVector<DecayIntegrator,DecayPhaseSpaceChannel> interfaceChannels
    ("Channels",
     "The phase space integration channels.",
     &DecayIntegrator::_channels, 0, false, false, true, true);
  
  static ClassDocumentation<DecayIntegrator> documentation
    ("The \\classname{DecayIntegrator} class is a base decayer class "
     "including a multi-channel integrator.");
  
  static Switch<DecayIntegrator,bool> interfaceInitialize
    ("Initialize",
     "Initialization of the phase space calculation",
     &DecayIntegrator::_Initialize, false, false, false);

  static SwitchOption interfaceInitializeon
    (interfaceInitialize,
     "on",
     "At initialisation find max weight and optimise the integration",
     true);

  static SwitchOption interfaceInitializeoff
    (interfaceInitialize,
     "off",
     "Use the maximum weight and channel weights supplied for the integration",
     false);
  
  static Parameter<DecayIntegrator,int> interfaceIteration
    ("Iteration",
     "Number of iterations for the initialization of the phase space",
     &DecayIntegrator::_niter, 10, 0, 100,
     false, false, true);  
  
  static Parameter<DecayIntegrator,int> interfacePoints
    ("Points",
     "number of phase space points to generate in the initialisation.",
     &DecayIntegrator::_npoint, 10000, 100, 100000000,
     false, false, true);
  
  static Parameter<DecayIntegrator,int> interfaceNtry
    ("Ntry",
     "Number of attempts to generate the decay",
     &DecayIntegrator::_ntry, 500, 0, 100000,
     false, false, true);
  
}

// return the matrix element squared
double DecayIntegrator::me2(bool vertex,
			    const int imode, const int ichan, const Particle &,
			    const ParticleVector & decay) const {return 0.;}

// flat phase space generation and weight
double DecayIntegrator::flatPhaseSpace(const Particle & inpart,
				       ParticleVector & outpart) const
{
  double wgt=0.;
  // masses of the particles
  Energy inmass=inpart.mass();
  vector<Energy> mass;
  // momenta of the particles
  vector<Lorentz5Momentum> part;
  part.resize(outpart.size());
  // generate the mass of the outgoing particles
  for(unsigned int ix=0;ix<outpart.size();++ix)
    {mass.push_back((outpart[ix]->dataPtr())->generateMass());}
  // two body decay
  if(outpart.size()==2)
    {
      double ctheta,phi;
      Kinematics::generateAngles(ctheta,phi);
      Kinematics::twoBodyDecay(inpart.momentum(), mass[0], mass[1],
			       ctheta, phi,part[0],part[1]);
      wgt = Kinematics::CMMomentum(inmass,mass[0],mass[1])/8./pi/inmass/inmass;
      outpart[0]->setMomentum(part[0]);
      outpart[1]->setMomentum(part[1]);
    }
  else
    {cout << "only the two body decay is currently implemented" << endl;}
  return wgt;
}

// initialise the phase space
void DecayIntegrator::initializePhaseSpace(const unsigned int imode,
					   const PDVector & decay)
{
  if(!_Initialize){return;}
  // create a particle vector from the particle data one
  ThePEG::PPtr inpart=decay[0]->produceParticle();
  ParticleVector particles;
  for(unsigned int ix=1;ix<decay.size();++ix)
    {particles.push_back(decay[ix]->produceParticle());}
  // now if using flat phase space
  _MaxWeight[imode]=0.;
  double wsum=0.,wsqsum=0.;
  double totsum(0.),totsq(0.);
  if(_channels.size()==0)
    {
      double wgt;
      Energy m0;
      for(int ix=0;ix<_npoint;++ix)
	{
	  // set the mass of the decaying particle
	  m0 = (inpart->dataPtr())->generateMass();
	  inpart->set5Momentum(Lorentz5Momentum(0.0,0.0,0.0,m0,m0));
	  // generate the weight for this point
	  int ichan;
	  wgt = weight(imode,ichan,*inpart,particles);
	  if(wgt>_MaxWeight[imode]){_MaxWeight[imode]=wgt;}
	  wsum=wsum+wgt;
	  wsqsum=wsqsum+wgt*wgt;
	}
      wsum=wsum/_npoint;
      wsqsum=wsqsum/_npoint-wsum*wsum;
      if(wsqsum<0.){wsqsum=0.;}
      wsqsum=sqrt(wsqsum/_npoint);
      // ouptut the information on the initialisation
      cout << "Initialized the phase space for the decay " 
	   << decay[0]->PDGName() << " -> ";
      for(unsigned int ix=1,N=decay.size();ix<N;++ix)
	{cout << decay[ix]->PDGName() << " ";}
      cout << endl;
      cout << "The partial width is " << wsum << " +/- " << wsqsum << " MeV" << endl;
      cout << "The maximum weight is " << _MaxWeight[imode] << endl;
    }
  else
    {
      // ensure that the starting weights add up to one
      double temp=0.;
      for(unsigned int ix=0;ix<_channels.size();++ix)
	{if(_channelon[imode][ix]){temp+=_channelwgts[imode][ix];}}
      for(unsigned int ix=0;ix<_channels.size();++ix)
	{if(_channelon[imode][ix]){_channelwgts[imode][ix]/=temp;}}
      for(int iy=0;iy<_niter;++iy)
	{
	  // zero the maximum weight
	  _MaxWeight[imode]=0.;
	  vector<double> wsum(_channels.size(),0.),wsqsum(_channels.size(),0.);
	  vector<int> nchan(_channels.size(),0);
	  totsum = 0.; totsq = 0.;
	  for(int ix=0;ix<_npoint;++ix)
	    {
	      double wgt;
	      Energy m0;
	      m0 = (inpart->dataPtr())->generateMass();
	      inpart->set5Momentum(Lorentz5Momentum(0.0,0.0,0.0,m0,m0));
	      // generate the weight for this point
	      int ichan;
	      wgt = weight(imode,ichan,*inpart,particles);
	      if(wgt>_MaxWeight[imode]){_MaxWeight[imode]=wgt;}
	      wsum[ichan]=wsum[ichan]+wgt;
	      totsum+=wgt;
	      wsqsum[ichan]=wsqsum[ichan]+wgt*wgt;
	      totsq+=wgt*wgt;
	      ++nchan[ichan];
	    }
	  totsum=totsum/_npoint;
	  totsq=totsq/_npoint-totsum*totsum;
	  if(totsq<0.){totsq=0.;}
	  totsq=sqrt(totsq/_npoint);
	  cout << "The partial width is " << iy << " " 
	       << totsum << " +/- " << totsq << " MeV" << endl;
	  // compute the individual terms
	  double total(0.);
	  for(unsigned int ix=0;ix<_channels.size();++ix)
	    {
	      if(nchan[ix]!=0)
		{
		  wsum[ix]=wsum[ix]/nchan[ix];
		  wsqsum[ix]=wsqsum[ix]/nchan[ix]-wsum[ix]*wsum[ix];
		  if(wsqsum[ix]<0.){wsqsum[ix]=0.;}
		  wsqsum[ix]=sqrt(wsqsum[ix]/nchan[ix]);
		}
	      else
		{
		  wsum[ix]=0;
		  wsqsum[ix]=0;
		}
	      total+=sqrt(wsqsum[ix])*_channelwgts[imode][ix];
	    }
	  double temp;
	  for(unsigned int ix=0;ix<_channels.size();++ix)
	    {
	      temp=sqrt(wsqsum[ix])*_channelwgts[imode][ix]/total;
	      _channelwgts[imode][ix]=temp;
	    }
	}
      // ouptut the information on the initialisation
      cout << "Initialized the phase space for the decay " 
	   << decay[0]->PDGName() << " -> ";
      for(unsigned int ix=1,N=decay.size();ix<N;++ix)
	{cout << decay[ix]->PDGName() << " ";}
      cout << endl;
      cout << "The partial width is " << totsum << " +/- " << totsq << " MeV" << endl;
      cout << "The maximum weight is " << _MaxWeight[imode] << endl;
      cout << "The weights for the different phase space channels are " << endl;
      for(unsigned int ix=0,N=_channels.size();ix<N;++ix)
	{
	  cout << "Channel " << ix << " was switched ";
	  if(_channelon[imode][ix]){cout << "on";}
	  else{cout << "off";}
	  cout << " and had weight " << _channelwgts[imode][ix] << endl;
	}
    }
}
  
// generate a phase-space point using multichannel phase space
double DecayIntegrator::channelPhaseSpace(unsigned int imode,
					  int & ichan, const Particle & inpart,
					  ParticleVector & outpart) const
{
  // select the channel
  vector<Lorentz5Momentum> momenta;
  double wgt=rnd();
  // select a channel
  ichan=-1;
  do
    {
      ++ichan;
      if(_channelon[imode][ichan]){wgt-=_channelwgts[imode][ichan];}
    }
  while(ichan<int(_channels.size())&&wgt>0.);
  // generate the momenta
  momenta=_channels[ichan]->generateMomenta(inpart.momentum());
  // compute the denominator of the weight
  wgt=0.;
  for(unsigned int ix=0;ix<_channels.size();++ix)
    {
      if(_channelon[imode][ix])
	{
	  wgt+=_channelwgts[imode][ix]*_channels[ix]->generateWeight(momenta);
	}
    }
  // now we need to set the momenta of the particles
  PDVector extpart=_channels[ichan]->externalParticles();
  vector<bool> done(extpart.size(),false);done[0]=true;
  bool found; int id;
  if(inpart.id()==extpart[0]->id())
    {
      for(unsigned int iy=0;iy<outpart.size();++iy)
	{
	  found = false;
	  id=outpart[iy]->id();
	  unsigned int ix=0;
	  do
	    {
	      ++ix;
	      if(!done[ix]&&extpart[ix]->id()==id)
		{
		  done[ix]=true;
		  outpart[iy]->setMomentum(momenta[ix]);
		  found=true;
		}
	    }
	  while(!found&&ix<extpart.size()-1);
	  if(!found)
	    {
	      cerr << "error in DecayIntegrator::channelPhaseSpace can't find particle"
		   << id << " for channel" << ichan << endl;
	      cout << "testing the id" << inpart.id() << "   " << extpart[0]->id()<<endl;
	    }
	}
    }
  else
    {
      tcPDPtr anti;
      for(unsigned int iy=0;iy<outpart.size();++iy)
	{
	  found = false;
	  anti=(outpart[iy]->dataPtr())->CC();
	  if(anti){id=anti->id();}
	  else{id=outpart[iy]->id();}
	  unsigned int ix=0;
	  do 
	    {
	      ++ix;
	      if(!done[ix]&&extpart[ix]->id()==id)
		{
		  done[ix]=true;
		  outpart[iy]->setMomentum(momenta[ix]);
		  found=true;
		}
	    }
	  while(!found&&ix<extpart.size()-1);
	  if(!found)
	    {
	      cerr << "error in DecayIntegrator::channelPhaseSpace can't find particle"
		   << id << " for channel" << ichan << endl;
	      cout << "testing the idB " << inpart.id() << "   " << extpart[0]->id()<<endl;
	    }
	}
    }
  // return the weight
  return 1./wgt;
}

// generate the decay
void DecayIntegrator::generate(bool intermediates,
			       unsigned int imode, const Particle & inpart,
			       ParticleVector & particles) const
{
  // construct a new particle which is at rest
  tcSpinInfoPtr tempspin =dynamic_ptr_cast<tcSpinInfoPtr>(inpart.spinInfo());
  Particle inrest=inpart;
  Hep3Vector boostv = -inpart.momentum().boostVector();
  LorentzMomentum test = inpart.momentum();
  inrest.setMomentum(test.boost(boostv));
  int ncount=0; double wgt;
  do
    {
      int ichan;
      wgt=weight(imode,ichan,inrest,particles);
      ++ncount;
      if(wgt>_MaxWeight[imode]){_MaxWeight[imode]=wgt;}
    }
  while(_MaxWeight[imode]*rnd()>wgt&&ncount<_ntry);
  if(ncount>=_ntry){particles.resize(0); return;}
  // set up the vertex for spin correlations
  const_ptr_cast<tPPtr>(&inpart)->spinInfo(inrest.spinInfo());
  constructVertex(inpart,particles);
  // return if intermediate particles not required
  if(_channelwgts[imode].size()==0||!intermediates)
    {
      Hep3Vector boostv = inpart.momentum().boostVector();
      LorentzMomentum test;
      for(unsigned int ix=0;ix<particles.size();++ix)
	{
	  test=particles[ix]->momentum();
	  particles[ix]->setMomentum(test.boost(boostv));
	}
    }
  // find the intermediate particles
  else
    {
      // select the channel
      int ichan=selectChannel(imode,inpart,particles);
      Hep3Vector boostv = -inpart.momentum().boostVector();
      LorentzMomentum test;
      for(unsigned int ix=0;ix<particles.size();++ix)
	{
	  test=particles[ix]->momentum();
	  particles[ix]->setMomentum(test.boost(boostv));
	}
      // generate the particle vector
      _channels[ichan]->generateIntermediates(inpart,particles);
    }
  return;
}

// output info on the integrator
ostream & operator<<(ostream & os, const DecayIntegrator & decay)
{
  os << "The integrator has " << decay._channels.size() << " channels"  << endl;
 for(unsigned int ix=0;ix<decay._channels.size();++ix)
   {
     os << "Information on channel " << ix << endl;
     os << *(decay._channels[ix]);
   }
 return os;
}

}
