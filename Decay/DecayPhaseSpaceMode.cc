// -*- C++ -*-
//
// DecayPhaseSpaceMode.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayPhaseSpaceMode class.
//

#include "DecayPhaseSpaceMode.h"
#include "Herwig++/PDT/GenericWidthGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/RefVector.h"
#include "DecayIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/Helicity/HelicityVertex.h"
#include "Herwig++/Decay/DecayVertex.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/SpinInfo.h"
#include "ThePEG/EventRecord/Event.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void DecayPhaseSpaceMode::persistentOutput(PersistentOStream & os) const {
  os << _integrator << _channels << _channelwgts << _maxweight << _niter 
     << _npoint << _ntry << _extpart << _partial << _widthgen << _massgen
     << _testOnShell;
}

void DecayPhaseSpaceMode::persistentInput(PersistentIStream & is, int) {
  is >> _integrator >> _channels >> _channelwgts >> _maxweight >> _niter 
     >> _npoint >> _ntry >> _extpart >> _partial >> _widthgen >> _massgen
     >> _testOnShell;
}

ClassDescription<DecayPhaseSpaceMode> DecayPhaseSpaceMode::initDecayPhaseSpaceMode;
// Definition of the static class description member.

void DecayPhaseSpaceMode::Init() {

  static ClassDocumentation<DecayPhaseSpaceMode> documentation
    ("The DecayPhaseSpaceMode class contains a number of phase space"
     " channels for the integration of a particular decay mode");

  static RefVector<DecayPhaseSpaceMode,DecayPhaseSpaceChannel> interfaceChannels
    ("Channels",
     "The phase space integration channels.",
     &DecayPhaseSpaceMode::_channels, 0, false, false, true, true);
 
}

// flat phase space generation and weight
Energy DecayPhaseSpaceMode::flatPhaseSpace(bool cc, const Particle & inpart,
					   ParticleVector & outpart) const {
  double wgt(1.);
  if(outpart.empty()) {
    outpart.reserve(_extpart.size()-1);
    for(unsigned int ix=1;ix<_extpart.size();++ix) {
      if(cc&&_extpart[ix]->CC()) {
	outpart.push_back((_extpart[ix]->CC())->produceParticle());
      }
      else {
	outpart.push_back(_extpart[ix]->produceParticle());
      }
    }
  }
  // masses of the particles
  Energy inmass(inpart.mass());
  vector<Energy> mass = externalMasses(inmass,wgt);
  // momenta of the particles
  vector<Lorentz5Momentum> part(outpart.size());
  // two body decay
  assert(outpart.size()==2);
  double ctheta,phi;
  Kinematics::generateAngles(ctheta,phi);
  Kinematics::twoBodyDecay(inpart.momentum(), mass[1], mass[2],
			   ctheta, phi,part[0],part[1]);
  wgt *= Kinematics::pstarTwoBodyDecay(inmass,mass[1],mass[2])/8./Constants::pi/inmass;
  outpart[0]->set5Momentum(part[0]);
  outpart[1]->set5Momentum(part[1]);
  return wgt*inmass;
}

// initialise the phase space
Energy DecayPhaseSpaceMode::initializePhaseSpace(bool init) {
  Energy output(ZERO);
  // ensure that the weights add up to one
  if(!_channels.empty()) {
    double temp=0.;
    for(unsigned int ix=0;ix<_channels.size();++ix) temp+=_channelwgts[ix];
    for(unsigned int ix=0;ix<_channels.size();++ix) _channelwgts[ix]/=temp;
  }
  if(!init) return ZERO;
  // create a particle vector from the particle data one
  ThePEG::PPtr inpart=_extpart[0]->produceParticle();
  ParticleVector particles;
  // now if using flat phase space
  _maxweight=0.;
  double totsum(0.),totsq(0.);
  InvEnergy pre=1./MeV;
  Energy prewid;
  if(_channels.empty()) {
    double wsum=0.,wsqsum=0.;
    Energy m0,mmin(ZERO);
    for(unsigned int ix=1;ix<_extpart.size();++ix) {
      mmin+=_extpart[ix]->massMin();
    }
    for(int ix=0;ix<_npoint;++ix) {
      int ichan;
      // set the mass of the decaying particle
      m0 = (inpart->dataPtr())->generateMass();
      double wgt=0.;
      if(m0>mmin) {
	inpart->set5Momentum(Lorentz5Momentum(m0));
	// compute the prefactor
	prewid = (_widthgen&&_partial>=0) ?
	  _widthgen->partialWidth(_partial,inpart->mass()) :
	  inpart->dataPtr()->width();
	pre = prewid>ZERO ? 1./prewid : 1./MeV;
	// generate the weight for this point
	try {
	  wgt = pre*weight(false,ichan,*inpart,particles,true);
	}
	catch (Veto) {
	  wgt=0.;
	}
      }
      if(wgt>_maxweight) _maxweight=wgt;
      wsum=wsum+wgt;
      wsqsum=wsqsum+wgt*wgt;
    }
    wsum=wsum/_npoint;
    wsqsum=wsqsum/_npoint-sqr(wsum);
    if(wsqsum<0.) wsqsum=0.;
    wsqsum=sqrt(wsqsum/_npoint);
    Energy fact = (_widthgen&&_partial>=0) ? 
      _widthgen->partialWidth(_partial,inpart->nominalMass()) :
      inpart->dataPtr()->width();
    if(fact==ZERO) fact=MeV;
    // factor for the weight with spin correlations
    _maxweight *= inpart->dataPtr()->iSpin()==1 ? 1.1 : 1.6;
    // ouptut the information on the initialisation
    CurrentGenerator::log() << "Initialized the phase space for the decay " 
			    << _extpart[0]->PDGName() << " -> ";
    for(unsigned int ix=1,N=_extpart.size();ix<N;++ix)
      CurrentGenerator::log() << _extpart[ix]->PDGName() << " ";
    CurrentGenerator::log() << "\n";
    if(fact!=MeV) CurrentGenerator::log() << "The branching ratio is " << wsum 
					  << " +/- " << wsqsum << "\n";
    CurrentGenerator::log() << "The partial width is " << wsum*fact/MeV 
			    << " +/- " << wsqsum*fact/MeV << " MeV\n";
    CurrentGenerator::log() << "The partial width is " 
			    << wsum*fact/6.58212E-22/MeV 
			    << " +/- " << wsqsum*fact/6.58212E-22/MeV<< " s-1\n";
    CurrentGenerator::log() << "The maximum weight is " 
			    << _maxweight << endl;
    output=wsum*fact;
  }
  else {
    for(int iy=0;iy<_niter;++iy) {
      // zero the maximum weight
      _maxweight=0.;
      vector<double> wsum(_channels.size(),0.),wsqsum(_channels.size(),0.);
      vector<int> nchan(_channels.size(),0);
      totsum = 0.; totsq = 0.;
      Energy m0,mmin(ZERO);
      for(unsigned int ix=1;ix<_extpart.size();++ix) {
	mmin+=_extpart[ix]->massMin();
      }
      for(int ix=0;ix<_npoint;++ix) {
	m0 = (inpart->dataPtr())->generateMass();
	double wgt=0.; 
	int ichan;
	if(m0>mmin) {
	  inpart->set5Momentum(Lorentz5Momentum(m0));
	  // compute the prefactor
	  prewid= (_widthgen&&_partial>=0) ? 
	    _widthgen->partialWidth(_partial,inpart->mass()) :
	    inpart->dataPtr()->width();
	  pre = prewid>ZERO ? 1./prewid : 1./MeV;
	  // generate the weight for this point
	  try {
	    wgt = pre*weight(false,ichan,*inpart,particles,true);
	  }
	  catch (Veto) {
	    wgt=0.;
	  }
	}
	if(wgt>_maxweight) _maxweight=wgt;
	wsum[ichan]=wsum[ichan]+wgt;
	totsum+=wgt;
	wsqsum[ichan]=wsqsum[ichan]+sqr(wgt);
	totsq+=wgt*wgt;
	++nchan[ichan];
      }
      totsum=totsum/_npoint;
      totsq=totsq/_npoint-sqr(totsum);
      if(totsq<0.) totsq=0.;
      totsq=sqrt(totsq/_npoint);
      CurrentGenerator::log() << "The branching ratio is " << iy << " " 
			      << totsum << " +/- " << totsq 
			      << _maxweight << "\n";
      // compute the individual terms
      double total(0.);
      for(unsigned int ix=0;ix<_channels.size();++ix) {
	if(nchan[ix]!=0) {
	  wsum[ix]=wsum[ix]/nchan[ix];
	  wsqsum[ix]=wsqsum[ix]/nchan[ix];
	  if(wsqsum[ix]<0.) wsqsum[ix]=0.;
	  wsqsum[ix]=sqrt(wsqsum[ix]/nchan[ix]);
	}
	else {
	  wsum[ix]=0;
	  wsqsum[ix]=0;
	}
	total+=sqrt(wsqsum[ix])*_channelwgts[ix];
      }
      double temp;
      for(unsigned int ix=0;ix<_channels.size();++ix) {
	temp=sqrt(wsqsum[ix])*_channelwgts[ix]/total;
	_channelwgts[ix]=temp;
      }
    }
    // factor for the weight with spin correlations
    _maxweight*= inpart->dataPtr()->iSpin()==1 ? 1.1 : 1.6;
    // ouptut the information on the initialisation
    Energy fact = (_widthgen&&_partial>=0) ? 
      _widthgen->partialWidth(_partial,inpart->nominalMass()) :
      inpart->dataPtr()->width();
    if(fact==ZERO) fact=MeV;
    CurrentGenerator::log() << "Initialized the phase space for the decay " 
			    << _extpart[0]->PDGName() << " -> ";
    for(unsigned int ix=1,N=_extpart.size();ix<N;++ix)
      CurrentGenerator::log() << _extpart[ix]->PDGName() << " ";
    CurrentGenerator::log() << "\n";
    if(fact!=MeV) CurrentGenerator::log() << "The branching ratio is " << totsum 
					  << " +/- " << totsq << "\n";
    CurrentGenerator::log() << "The partial width is " << totsum*fact/MeV 
			    << " +/- " << totsq*fact/MeV << " MeV\n";
    CurrentGenerator::log() << "The partial width is " 
			    << totsum*fact/6.58212E-22/MeV 
			    << " +/- " << totsq*fact/6.58212E-22/MeV 
			    << " s-1\n";
    CurrentGenerator::log() << "The maximum weight is " 
			    << _maxweight << "\n";
    CurrentGenerator::log() << "The weights for the different phase" 
			    << " space channels are \n";
    for(unsigned int ix=0,N=_channels.size();ix<N;++ix) {
      CurrentGenerator::log() << "Channel " << ix 
			      << " had weight " << _channelwgts[ix] 
			      << "\n";
    }
    CurrentGenerator::log() << flush;
    output=totsum*fact;
  }
  return output;
}
 
// generate a phase-space point using multichannel phase space
Energy DecayPhaseSpaceMode::channelPhaseSpace(bool cc,
					      int & ichan, const Particle & inpart,
					      ParticleVector & outpart) const {
  // select the channel
  vector<Lorentz5Momentum> momenta;
  double wgt(UseRandom::rnd());
  // select a channel
  ichan=-1;
  do{++ichan;wgt-=_channelwgts[ichan];}
  while(ichan<int(_channels.size())&&wgt>0.);
  // generate the momenta
  if(ichan==int(_channels.size())) {
    throw DecayIntegratorError() << "DecayPhaseSpaceMode::channelPhaseSpace()"
				 << " failed to select a channel" 
				 << Exception::abortnow;
  }
  // generate the masses of the external particles
  double masswgt(1.);
  vector<Energy> mass(externalMasses(inpart.mass(),masswgt));
  momenta=_channels[ichan]->generateMomenta(inpart.momentum(),mass);
  // compute the denominator of the weight
  wgt=0.;
  unsigned int ix;
  for(ix=0;ix<_channels.size();++ix) {
    wgt+=_channelwgts[ix]*_channels[ix]->generateWeight(momenta);
  }
  // now we need to set the momenta of the particles
  // create the particles if they don't exist
  if(outpart.empty()) {
    for(ix=1;ix<_extpart.size();++ix) {
      if(cc&&_extpart[ix]->CC()) {
	outpart.push_back((_extpart[ix]->CC())->produceParticle());
      }
      else {
	outpart.push_back(_extpart[ix]->produceParticle());
      }
    }
  }
  // set up the momenta
  for(ix=0;ix<outpart.size();++ix) outpart[ix]->set5Momentum(momenta[ix+1]);
  // return the weight
  return inpart.mass()*masswgt/wgt;
}

// generate the decay
ParticleVector DecayPhaseSpaceMode::generate(bool intermediates,bool cc,
					     const Particle & inpart) const {
  // compute the prefactor
  InvEnergy pre(1./MeV); 
  Energy prewid;
  if(_widthgen&&_partial>=0) prewid=_widthgen->partialWidth(_partial,inpart.mass());
  else                       prewid=(inpart.dataPtr()->width());
  pre = prewid>ZERO ? 1./prewid : 1./MeV;
  // Particle vector for the output
  ParticleVector particles;
  // construct a new particle which is at rest
  Particle inrest(inpart);
  inrest.boost(-inpart.momentum().boostVector());
  int ncount(0),ichan; double wgt(0.);
  unsigned int ix;
  try {
    do {
      wgt=pre*weight(cc,ichan,inrest,particles,ncount==0);
      ++ncount;
      if(wgt>_maxweight) {
	CurrentGenerator::log() << "Resetting max weight for decay " 
				<< inrest.PDGName() << " -> ";
	for(ix=0;ix<particles.size();++ix)
	  CurrentGenerator::log() << "  " << particles[ix]->PDGName();
	CurrentGenerator::log() << "  " << _maxweight << "  " << wgt 
				<< "  " << inrest.mass()/MeV << "\n";
	_maxweight=wgt;
      }
    }
    while(_maxweight*UseRandom::rnd()>wgt&&ncount<_ntry);
    if(ncount>=_ntry) {
      CurrentGenerator::log() << "The decay " << inrest.PDGName() << " -> ";
      for(ix=0;ix<particles.size();++ix) 
	CurrentGenerator::log()  << "  " << particles[ix]->PDGName();
      CurrentGenerator::log() << "  " << _maxweight << " " << _ntry 
			      << " is too inefficient for the particle "
			      << inpart << "vetoing the decay \n";
      particles.clear();
      throw Veto();
    }
  }
  catch (Veto) {
    // restore the incoming particle to its original state
    Boost boostv(inpart.momentum().boostVector());
    inrest.boost(boostv);
    throw Veto();
  }
  // set up the vertex for spin correlations
  me2(-1,inrest,particles,DecayIntegrator::Terminate);
  const_ptr_cast<tPPtr>(&inpart)->spinInfo(inrest.spinInfo());
  constructVertex(inpart,particles);
  // return if intermediate particles not required
  Boost boostv(inpart.momentum().boostVector());
  if(_channelwgts.empty()||!intermediates) {
    for(ix=0;ix<particles.size();++ix) particles[ix]->boost(boostv);
  }
  // find the intermediate particles
  else {
    // select the channel
    _ichannel = selectChannel(inpart,particles);
    for(ix=0;ix<particles.size();++ix) particles[ix]->boost(boostv);
    // generate the particle vector
    _channels[_ichannel]->generateIntermediates(cc,inpart,particles);
  }
  return particles;
}

// construct the vertex for spin corrections
void DecayPhaseSpaceMode::constructVertex(const Particle & inpart, 
					  const ParticleVector & decay) const {
  // construct the decay vertex
  VertexPtr vertex(new_ptr(DecayVertex()));
  DVertexPtr Dvertex(dynamic_ptr_cast<DVertexPtr>(vertex));
  // set the incoming particle for the decay vertex
  dynamic_ptr_cast<tcSpinfoPtr>(inpart.spinInfo())->setDecayVertex(vertex);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    dynamic_ptr_cast<tcSpinfoPtr>(decay[ix]->spinInfo())->setProductionVertex(vertex);
  }
  // set the matrix element
  Dvertex->ME().reset(_integrator->ME());
}

// output info on the mode
ostream & Herwig::operator<<(ostream & os, const DecayPhaseSpaceMode & decay) {
  os << "The mode has " << decay._channels.size() << " channels\n";
  os << "This is a mode for the decay of " << decay._extpart[0]->PDGName() 
     << " to ";
  for(unsigned int iz=1,N=decay._extpart.size();iz<N;++iz) {
    os << decay._extpart[iz]->PDGName() << " ";
  }
  os << "\n";
  for(unsigned int ix=0;ix<decay._channels.size();++ix) {
    os << "Information on channel " << ix << "\n";
    os << *(decay._channels[ix]);
  }
  return os;
}

void DecayPhaseSpaceMode::doinit() {
  // check the size of the weight vector is the same as the number of channels
  if(_channelwgts.size()!=numberChannels()) {
    throw Exception() << "Inconsistent number of channel weights and channels"
		      << " in DecayPhaseSpaceMode " << Exception::abortnow;
  }
  Interfaced::doinit();
  _massgen.resize(_extpart.size());
  if(_extpart[0]->widthGenerator()) {
    _widthgen=
      dynamic_ptr_cast<cGenericWidthGeneratorPtr>(_extpart[0]->widthGenerator());
    //const_ptr_cast<GenericWidthGeneratorPtr>(_widthgen)->init();
  }
  else {
    _widthgen=cGenericWidthGeneratorPtr();
  }
  tcGenericWidthGeneratorPtr wtemp;
  for(unsigned int ix=0;ix<_extpart.size();++ix) {
    assert(_extpart[ix]);
    _massgen[ix]= dynamic_ptr_cast<cGenericMassGeneratorPtr>(_extpart[ix]->massGenerator());
    if(ix>0) {
      wtemp=
	dynamic_ptr_cast<tcGenericWidthGeneratorPtr>(_extpart[ix]->widthGenerator());
      //if(wtemp) const_ptr_cast<tGenericWidthGeneratorPtr>(wtemp)->init();
    }
  }
  for(unsigned int ix=0;ix<_channels.size();++ix) {
    _channels[ix]->init();
  }
}
  
void DecayPhaseSpaceMode::doinitrun() {
  // check the size of the weight vector is the same as the number of channels
  if(_channelwgts.size()!=numberChannels()) {
    throw Exception() << "Inconsistent number of channel weights and channels"
		      << " in DecayPhaseSpaceMode " << Exception::abortnow;
  }
  for(unsigned int ix=0;ix<_channels.size();++ix) {
    _channels[ix]->initrun();
  }
  if(_widthgen) const_ptr_cast<GenericWidthGeneratorPtr>(_widthgen)->initrun();
  tcGenericWidthGeneratorPtr wtemp;
  for(unsigned int ix=1;ix<_extpart.size();++ix) {
    wtemp=
      dynamic_ptr_cast<tcGenericWidthGeneratorPtr>(_extpart[ix]->widthGenerator());
    if(wtemp) const_ptr_cast<tGenericWidthGeneratorPtr>(wtemp)->initrun();
  }
  Interfaced::doinitrun();
}

// generate the masses of the external particles
vector<Energy> DecayPhaseSpaceMode::externalMasses(Energy inmass,double & wgt) const {
  vector<Energy> mass(1,inmass);
  mass.reserve(_extpart.size());
  vector<int> notdone;
  Energy mlow(ZERO);
  // set masses of stable particles and limits 
  for(unsigned int ix=1;ix<_extpart.size();++ix) {
    // get the mass of the particle if can't use weight
    if(!_massgen[ix] || _extpart[ix]->stable()) {
      mass.push_back(_extpart[ix]->generateMass());
      mlow+=mass[ix];
    }
    else {
      mass.push_back(ZERO);
      notdone.push_back(ix);
      mlow+=_extpart[ix]->mass()-_extpart[ix]->widthLoCut();
    }
  }
  if(mlow>inmass) {
    CurrentGenerator::log() << "Decay mode " << _extpart[0]->PDGName() << " -> ";
    for(unsigned int ix=1;ix<_extpart.size();++ix) {
      CurrentGenerator::log() << _extpart[ix]->PDGName() << " ";
    }
    CurrentGenerator::log() << "is below threshold in DecayPhaseMode::externalMasses()"
			    << "the threshold is " << mlow/GeV 
			    << "GeV and the parent mass is " << inmass/GeV 
			    << " GeV\n";
    throw Veto();
  }
  // now we need to generate the masses for the particles we haven't
  unsigned int iloc;
  double wgttemp;
  Energy low=ZERO;
  for( ;!notdone.empty();) {
    iloc=long(UseRandom::rnd()*(notdone.size()-1));
    low=_extpart[notdone[iloc]]->mass()-_extpart[notdone[iloc]]->widthLoCut();
    mlow-=low;
    mass[notdone[iloc]]=
      _massgen[notdone[iloc]]->mass(wgttemp,*_extpart[notdone[iloc]],low,inmass-mlow);
    wgt*=wgttemp;
    mlow+=mass[notdone[iloc]];
    notdone.erase(notdone.begin()+iloc);
  }
  return mass;
}
