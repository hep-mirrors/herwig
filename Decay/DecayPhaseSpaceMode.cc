// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayPhaseSpaceMode class.
//

#include "DecayPhaseSpaceMode.h"
#include "Herwig++/PDT/GenericWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/RefVector.h"
#include "Herwig++/PDT/GenericMassGenerator.h"
#include "DecayIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Helicity/HelicityVertex.h"
#include "Herwig++/Helicity/Correlations/DecayVertex.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DecayPhaseSpaceMode.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig{
using namespace ThePEG;
using ThePEG::Helicity::VertexPtr;
using Herwig::Helicity::DecayVertex;
using Herwig::Helicity::DVertexPtr;
using ThePEG::Helicity::tcSpinfoPtr;
using ThePEG::Helicity::SpinfoPtr;

// default constructor
DecayPhaseSpaceMode::DecayPhaseSpaceMode() 
{
  _niter=10;
  _npoint=10000;
  _ntry=500;
  _partial=-1;
}

// constructor with decayer and particles
DecayPhaseSpaceMode::DecayPhaseSpaceMode(PDVector in,DecayIntegratorPtr intin) 
{
  _niter=10;
  _npoint=10000;
  _ntry=500;
  _extpart=in;
  _integrator=intin;
  _partial=-1;
}

// copy constructor
DecayPhaseSpaceMode::DecayPhaseSpaceMode(const DecayPhaseSpaceMode & x)
  : Interfaced(x),_integrator(x._integrator), _channels(x._channels),
    _channelwgts(x._channelwgts), _MaxWeight(x._MaxWeight),_niter(x._niter),
    _npoint(x._npoint), _ntry(x._ntry),_extpart(x._extpart), _partial(x._partial),
    _widthgen(x._widthgen), _massgen(x._massgen) {}

// destructor
DecayPhaseSpaceMode::~DecayPhaseSpaceMode() {}

void DecayPhaseSpaceMode::persistentOutput(PersistentOStream & os) const {
  os << _integrator << _channels << _channelwgts << _MaxWeight << _niter 
     << _npoint << _ntry << _extpart << _partial << _widthgen;
}

void DecayPhaseSpaceMode::persistentInput(PersistentIStream & is, int) {
  is >> _integrator >> _channels >> _channelwgts >> _MaxWeight >> _niter 
     >> _npoint >> _ntry >> _extpart >> _partial >> _widthgen;
}

ClassDescription<DecayPhaseSpaceMode> DecayPhaseSpaceMode::initDecayPhaseSpaceMode;
// Definition of the static class description member.

void DecayPhaseSpaceMode::Init() {

  static ClassDocumentation<DecayPhaseSpaceMode> documentation
    ("The \\classname{DecayPhaseSpaceMode} class contains a number of phase space"
     " channels for the integration of a particular decay mode");

  static RefVector<DecayPhaseSpaceMode,DecayPhaseSpaceChannel> interfaceChannels
    ("Channels",
     "The phase space integration channels.",
     &DecayPhaseSpaceMode::_channels, 0, false, false, true, true);
 
}

// flat phase space generation and weight
double DecayPhaseSpaceMode::flatPhaseSpace(bool cc, const Particle & inpart,
					   ParticleVector & outpart) const
{
  double wgt(1.);
  if(outpart.empty())
    {
      for(unsigned int ix=1;ix<_extpart.size();++ix)
	{
	  if(cc&&_extpart[ix]->CC())
	    {outpart.push_back((_extpart[ix]->CC())->produceParticle());}
	  else{outpart.push_back(_extpart[ix]->produceParticle());}
	}
    }
  // masses of the particles
  Energy inmass(inpart.mass());
  vector<Energy> mass = externalMasses(inmass,wgt);
  // momenta of the particles
  vector<Lorentz5Momentum> part(outpart.size());
  // two body decay
  if(outpart.size()==2)
    {
      double ctheta,phi;
      Kinematics::generateAngles(ctheta,phi);
      Kinematics::twoBodyDecay(inpart.momentum(), mass[1], mass[2],
			       ctheta, phi,part[0],part[1]);
      wgt *= Kinematics::CMMomentum(inmass,mass[1],mass[2])/8./pi/inmass;
      outpart[0]->set5Momentum(part[0]);
      outpart[1]->set5Momentum(part[1]);
    }
  else
    {throw DecayIntegratorError() << "DecayPhaseSpaceMode::flatPhaseSpace "
				   << "only two body modes currently implemented" 
				   << Exception::abortnow;}
  return wgt*inmass;
}

// initialise the phase space
void DecayPhaseSpaceMode::initializePhaseSpace(bool init)
{
  // ensure that the weights add up to one
  if(!_channels.empty())
    {
      double temp=0.;
      for(unsigned int ix=0;ix<_channels.size();++ix){temp+=_channelwgts[ix];}
      for(unsigned int ix=0;ix<_channels.size();++ix){_channelwgts[ix]/=temp;}
    }
  if(!init){return;}
  // create a particle vector from the particle data one
  ThePEG::PPtr inpart=_extpart[0]->produceParticle();
  ParticleVector particles;
  // now if using flat phase space
  _MaxWeight=0.;
  double wsum=0.,wsqsum=0.;
  double totsum(0.),totsq(0.);
  InvEnergy pre=1.;Energy prewid;
  if(_channels.empty())
    {
      double wgt;
      Energy m0;
      for(int ix=0;ix<_npoint;++ix)
	{
	  // set the mass of the decaying particle
	  m0 = (inpart->dataPtr())->generateMass();
	  inpart->set5Momentum(Lorentz5Momentum(0.0,0.0,0.0,m0,m0));
	  // compute the prefactor
	  if(_widthgen&&_partial>=0)
	    {prewid=_widthgen->partialWidth(_partial,inpart->mass());}
	  else
	    {prewid=(inpart->dataPtr()->width());}
	  if(prewid>0.){pre=1./prewid;}
	  else{pre=1./MeV;}
	  // generate the weight for this point
	  int ichan;
	  try{wgt = pre*weight(false,false,ichan,*inpart,particles);}
	  catch (Veto){wgt=0.;}
	  if(wgt>_MaxWeight){_MaxWeight=wgt;}
	  wsum=wsum+wgt;
	  wsqsum=wsqsum+wgt*wgt;
	}
      wsum=wsum/_npoint;
      wsqsum=wsqsum/_npoint-wsum*wsum;
      if(wsqsum<0.){wsqsum=0.;}
      wsqsum=sqrt(wsqsum/_npoint);
      Energy fact;
      if(_widthgen&&_partial>=0)
	{
	  fact=_widthgen->partialWidth(_partial,inpart->nominalMass());
	  CurrentGenerator::current().log() << "testing the prefactor A" << fact << endl;
	}
      else
	{
	  fact=(inpart->dataPtr()->width());
	  if(fact==0.){fact=MeV;}
	  CurrentGenerator::current().log()<< "testing the prefactor B" << fact << endl;
	}
      if(fact==0.){fact=1.;}
      // ouptut the information on the initialisation
      CurrentGenerator::current().log() << "Initialized the phase space for the decay " 
					<< _extpart[0]->PDGName() << " -> ";
      for(unsigned int ix=1,N=_extpart.size();ix<N;++ix)
	{CurrentGenerator::current().log() << _extpart[ix]->PDGName() << " ";}
      CurrentGenerator::current().log() << "  " << _partial << endl;
      CurrentGenerator::current().log() << "The partial width is " << wsum 
					<< " +/- " << wsqsum << " MeV" << endl;
      CurrentGenerator::current().log() << "The partial width is " 
					<< wsum*fact/6.58212E-22 
	   << " +/- " << wsqsum*fact/6.58212E-22<< " s-1" << endl;
      CurrentGenerator::current().log() << "The maximum weight is " 
					<< _MaxWeight << endl;
    }
  else
    {
      for(int iy=0;iy<_niter;++iy)
	{
	  // zero the maximum weight
	  _MaxWeight=0.;
	  vector<double> wsum(_channels.size(),0.),wsqsum(_channels.size(),0.);
	  vector<int> nchan(_channels.size(),0);
	  totsum = 0.; totsq = 0.;
	  double wgt; Energy m0;
	  for(int ix=0;ix<_npoint;++ix)
	    {
	      m0 = (inpart->dataPtr())->generateMass();
	      inpart->set5Momentum(Lorentz5Momentum(0.0,0.0,0.0,m0,m0));
	      // compute the prefactor
	      if(_widthgen&&_partial>=0)
		{prewid=_widthgen->partialWidth(_partial,inpart->mass());}
	      else
		{prewid=(inpart->dataPtr()->width());}
	      if(prewid>0){pre=1./prewid;}
	      else{pre=1./MeV;}
	      // generate the weight for this point
	      int ichan;
	      try {wgt = pre*weight(false,false,ichan,*inpart,particles);}
	      catch (Veto){wgt=0.;}
	      if(wgt>_MaxWeight){_MaxWeight=wgt;}
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
	  CurrentGenerator::current().log() << "The partial width is " << iy << " " 
					    << totsum << " +/- " << totsq 
					    << " MeV " << _MaxWeight << endl;
	  // compute the individual terms
	  double total(0.);
	  for(unsigned int ix=0;ix<_channels.size();++ix)
	    {
	      if(nchan[ix]!=0)
		{
		  wsum[ix]=wsum[ix]/nchan[ix];
		  wsqsum[ix]=wsqsum[ix]/nchan[ix];
		  //wsum[ix]*wsum[ix];
		  if(wsqsum[ix]<0.){wsqsum[ix]=0.;}
		  wsqsum[ix]=sqrt(wsqsum[ix]/nchan[ix]);
		}
	      else
		{
		  wsum[ix]=0;
		  wsqsum[ix]=0;
		}
	      total+=sqrt(wsqsum[ix])*_channelwgts[ix];
	    }
	  double temp;
	  for(unsigned int ix=0;ix<_channels.size();++ix)
	    {
	      temp=sqrt(wsqsum[ix])*_channelwgts[ix]/total;
	      _channelwgts[ix]=temp;
	    }
	}
      // ouptut the information on the initialisation
      Energy fact;
      if(_widthgen&&_partial>=0)
	{
	  fact=_widthgen->partialWidth(_partial,inpart->nominalMass());
	  CurrentGenerator::current().log()<< "testing the prefactor A" << fact << endl;
	}
      else
	{
	  fact=(inpart->dataPtr()->width());
	  if(fact==0.){fact=MeV;}
	  CurrentGenerator::current().log()<< "testing the prefactor B" << fact << endl;
	}
      if(fact==0.){fact=1.;}
      CurrentGenerator::current().log() << "Initialized the phase space for the decay " 
					<< _extpart[0]->PDGName() << " -> ";
      for(unsigned int ix=1,N=_extpart.size();ix<N;++ix)
	{CurrentGenerator::current().log() << _extpart[ix]->PDGName() << " ";}
      CurrentGenerator::current().log() << endl;
      CurrentGenerator::current().log() << "The partial width is " << totsum 
					<< " +/- " << totsq << " MeV" << endl;
      CurrentGenerator::current().log() << "The partial width is " 
					<< totsum*fact/6.58212E-22 
					<< " +/- " << totsq*fact/6.58212E-22 
					<< " s-1" << endl;
      CurrentGenerator::current().log() << "The maximum weight is " 
					<< _MaxWeight << endl;
      CurrentGenerator::current().log() << "The weights for the different phase" 
					<< " space channels are " << endl;
      for(unsigned int ix=0,N=_channels.size();ix<N;++ix)
	{CurrentGenerator::current().log() << "Channel " << ix 
					   << " had weight " << _channelwgts[ix] 
					   << endl;}
    }
}
// generate a phase-space point using multichannel phase space
Energy DecayPhaseSpaceMode::channelPhaseSpace(bool cc,
					      int & ichan, const Particle & inpart,
					      ParticleVector & outpart) const
{
  // select the channel
  vector<Lorentz5Momentum> momenta;
  double wgt(CurrentGenerator::current().rnd());
  // select a channel
  ichan=-1;
  do{++ichan;wgt-=_channelwgts[ichan];}
  while(ichan<int(_channels.size())&&wgt>0.);
  // generate the momenta
  if(ichan==int(_channels.size()))
    {throw DecayIntegratorError() << "DecayPhaseSpaceMode::channelPhaseSpace()"
				  << " failed to select a channel" 
				  << Exception::abortnow;}
  // generate the masses of the external particles
  double masswgt(1.);
  vector<Energy> mass(externalMasses(inpart.mass(),masswgt));
  momenta=_channels[ichan]->generateMomenta(inpart.momentum(),mass);
  // compute the denominator of the weight
  wgt=0.;
  unsigned int ix;
  for(ix=0;ix<_channels.size();++ix)
    {wgt+=_channelwgts[ix]*_channels[ix]->generateWeight(momenta);}
  // now we need to set the momenta of the particles
  // create the particles if they don't exist
  if(outpart.empty())
    {
      for(ix=1;ix<_extpart.size();++ix)
	{
	  if(cc&&_extpart[ix]->CC())
	    {outpart.push_back((_extpart[ix]->CC())->produceParticle());}
	  else{outpart.push_back(_extpart[ix]->produceParticle());}
	}
    }
  // set up the momenta
  for(ix=0;ix<outpart.size();++ix){outpart[ix]->set5Momentum(momenta[ix+1]);}
  // return the weight
  return inpart.mass()*masswgt/wgt;
}

// generate the decay
ParticleVector DecayPhaseSpaceMode::generate(bool intermediates,bool cc,
					     const Particle & inpart) const
{
  // compute the prefactor
  InvEnergy pre(1.);Energy prewid;
  if(_widthgen&&_partial>=0){prewid=_widthgen->partialWidth(_partial,inpart.mass());}
  else                      {prewid=(inpart.dataPtr()->width());}
  if(prewid>0.){pre=1./prewid;}
  else{pre=1./MeV;}
  // Particle vector for the output
  ParticleVector particles;
  // construct a new particle which is at rest
  Particle inrest(inpart);
  inrest.boost(-inpart.momentum().boostVector());
  int ncount(0),ichan; double wgt;
  unsigned int ix,iy;
  try {
    do
      {
	for(iy=0;iy<particles.size();++iy){particles[iy]->spinInfo(SpinfoPtr());}
	wgt=pre*weight(true,cc,ichan,inrest,particles);
	++ncount;
	if(wgt>_MaxWeight)
	  {
	    CurrentGenerator::current().log() << "Resetting max weight for decay " 
					      << inrest.PDGName() << " -> ";
	    for(ix=0;ix<particles.size();++ix)
	      {CurrentGenerator::current().log()  << "  " << particles[ix]->PDGName();}
	    CurrentGenerator::current().log() << "  " << _MaxWeight << "   " << wgt 
					      << "   " << inrest.mass() << endl;
	    _MaxWeight=wgt;
	  }
      }
    while(_MaxWeight*CurrentGenerator::current().rnd()>wgt&&ncount<_ntry);
    if(ncount>=_ntry)
      {
	CurrentGenerator::current().log() << "The decay " << inrest.PDGName() << " -> ";
	for(ix=0;ix<particles.size();++ix)
	  {CurrentGenerator::current().log()  << "  " << particles[ix]->PDGName();}
	CurrentGenerator::current().log() << "  " << _MaxWeight << " " << _ntry 
					  << " is too inefficient vetoing event " << endl;
	CurrentGenerator::current().log() << inrest.mass() 
					  << "  " << pre 
					  << "  " << wgt << endl;
	if(_widthgen&&_partial>=0)
	  {CurrentGenerator::current().log() << "Used running width " 
					     << _widthgen << "  " << _partial << endl;}
	else
	  {CurrentGenerator::current().log() << _widthgen << "  " << _partial << endl;}
	particles.resize(0);
	throw Veto();
	return particles;
      }
  }
  catch (Veto)
    {
      // restore the incoming particle to its original state
      Hep3Vector boostv(inpart.momentum().boostVector());
      inrest.boost(boostv);
      throw Veto();
    }
  // set up the vertex for spin correlations
  const_ptr_cast<tPPtr>(&inpart)->spinInfo(inrest.spinInfo());
  constructVertex(inpart,particles);
  // return if intermediate particles not required
  Hep3Vector boostv(inpart.momentum().boostVector());
  if(_channelwgts.empty()||!intermediates)
    {for(ix=0;ix<particles.size();++ix){particles[ix]->boost(boostv);}}
  // find the intermediate particles
  else
    {
      // select the channel
      int ichan(selectChannel(inpart,particles));
      for(ix=0;ix<particles.size();++ix)
	{particles[ix]->boost(boostv);}
      // generate the particle vector
      _channels[ichan]->generateIntermediates(cc,inpart,particles);
    }
  return particles;
}

// construct the vertex for spin corrections
void DecayPhaseSpaceMode::constructVertex(const Particle & inpart, 
					  const ParticleVector & decay) const
 {
   // construct the decay vertex
   VertexPtr vertex(new_ptr(DecayVertex()));
   DVertexPtr Dvertex(dynamic_ptr_cast<DVertexPtr>(vertex));
   // set the incoming particle for the decay vertex
   dynamic_ptr_cast<tcSpinfoPtr>(inpart.spinInfo())->setDecayVertex(vertex);
   for(unsigned int ix=0;ix<decay.size();++ix)
     {dynamic_ptr_cast<tcSpinfoPtr>
	 (decay[ix]->spinInfo())->setProductionVertex(vertex);}
   // set the matrix element
   Dvertex->ME().reset(_integrator->ME());
}

// output info on the mode
ostream & operator<<(ostream & os, const DecayPhaseSpaceMode & decay)
{
  os << "The mode has " << decay._channels.size() << " channels"  << endl;
  os << "This is a mode for the decay of " << decay._extpart[0]->PDGName() 
     << " to ";
  for(unsigned int iz=1,N=decay._extpart.size();iz<N;++iz)
    {os << decay._extpart[iz]->PDGName() << " ";}
  os << endl;
  for(unsigned int ix=0;ix<decay._channels.size();++ix)
    {
      os << "Information on channel " << ix << endl;
      os << *(decay._channels[ix]);
    }
  return os;
}

void DecayPhaseSpaceMode::setPartialWidth(int in){_partial=in;}
// return the phase space weight for a given point

Energy DecayPhaseSpaceMode::weight(bool vertex,bool cc,int & ichan,
				   const Particle & inpart,
				   ParticleVector & particles) const
{
  double mewgt(0.);Energy phwgt(0.);
  ichan=0;
  // generate the phase space point and get the weight
  if(_channels.size()==0){phwgt = flatPhaseSpace(cc,inpart,particles);} 
  else{phwgt = channelPhaseSpace(cc,ichan,inpart,particles);}
  // generate the matrix element
  mewgt = me2(vertex,-1,inpart,particles);
  //cout << "testing the partial width  " 
  //     << mewgt*phwgt << "  "
  //     << mewgt       << "  " 
  //     << phwgt       << endl;
  return mewgt*phwgt;
}

void DecayPhaseSpaceMode::doinitrun() {
  Interfaced::doinitrun();
  _massgen.resize(_extpart.size());
  if(_extpart[0]->widthGenerator())
    {
      _widthgen=
	dynamic_ptr_cast<cGenericWidthGeneratorPtr>(_extpart[0]->widthGenerator());
      const_ptr_cast<GenericWidthGeneratorPtr>(_widthgen)->initrun();
    }
  tcGenericWidthGeneratorPtr wtemp;
  for(unsigned int ix=0;ix<_extpart.size();++ix)
    {
      _massgen[ix]=
	dynamic_ptr_cast<cGenericMassGeneratorPtr>(_extpart[ix]->massGenerator());
      if(ix>0)
	{
	  wtemp=
	    dynamic_ptr_cast<tcGenericWidthGeneratorPtr>(_extpart[ix]->widthGenerator());
	  if(wtemp){const_ptr_cast<tGenericWidthGeneratorPtr>(wtemp)->initrun();}
	}
    }
}

// generate the masses of the external particles
vector<Energy> DecayPhaseSpaceMode::externalMasses(Energy inmass,double & wgt) const
{
  // CurrentGenerator::current().log() 
  //   << "testing generating the off-shell masses" << endl;
  vector<Energy> mass(1,inmass);
  vector<int> notdone;
  Energy mlow(0.);
  // set masses of stable particles and limits 
  for(unsigned int ix=1;ix<_extpart.size();++ix)
    {
      // get the mass of the particle if can't use weight
      if(!_massgen[ix])
	{mass.push_back(_extpart[ix]->generateMass());mlow+=mass[ix];}
      else
	{
	  mass.push_back(0.);
	  notdone.push_back(ix);
	  mlow+=_extpart[ix]->mass()-_extpart[ix]->widthLoCut();
	}
    }
  if(mlow>inmass){throw Veto();}
  // now we need to generate the masses for the particles we haven't
  unsigned int iloc;
  double wgttemp;
  Energy low=0.;
  for( ;!notdone.empty();)
    {
      iloc=long(CurrentGenerator::current().rnd()*(notdone.size()-1));
      low=_extpart[notdone[iloc]]->mass()-_extpart[notdone[iloc]]->widthLoCut();
      mlow-=low;
      mass[notdone[iloc]]=
	_massgen[notdone[iloc]]->mass(*_extpart[notdone[iloc]],wgttemp,low,inmass-mlow);
      wgt*=wgttemp;
      mlow+=mass[notdone[iloc]];
      notdone.erase(notdone.begin()+iloc);
    }
  return mass;
}

}
