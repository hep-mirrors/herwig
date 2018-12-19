// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhaseSpaceMode class.
//

#include "PhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PhaseSpaceMode::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << channels_ << maxWeight_ << outgoing_ << outgoingCC_
     << partial_ << widthGen_ << massGen_ << testOnShell_
     << ounit(eMax_,GeV) << nRand_;
}

void PhaseSpaceMode::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> channels_ >> maxWeight_ >> outgoing_ >> outgoingCC_
     >> partial_ >> widthGen_ >> massGen_ >> testOnShell_
     >> iunit(eMax_,GeV) >> nRand_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PhaseSpaceMode,Base>
describeHerwigPhaseSpaceMode("Herwig::PhaseSpaceMode", "libHerwig.so");

void PhaseSpaceMode::Init() {

  static ClassDocumentation<PhaseSpaceMode> documentation
    ("The PhaseSpaceMode classs contains a number of phase space"
     " channels for the integration of a particular decay mode");

}

ParticleVector
PhaseSpaceMode::generateDecay(const Particle & inpart,
			      tcDecayIntegratorPtr decayer,
			      bool intermediates,bool cc) {
  decayer->ME(DecayMEPtr());
  eps_=decayer->eps_;
  // compute the prefactor
  Energy prewid = (widthGen_&&partial_>=0) ?
    widthGen_->partialWidth(partial_,inpart.mass()) :
    incoming_.first->width();
  InvEnergy pre = prewid>ZERO ? 1./prewid : 1./MeV;
  // boosts to/from rest
  Boost bv =-inpart.momentum().boostVector();
  double gammarest = inpart.momentum().e()/inpart.momentum().mass();
  LorentzRotation boostToRest( bv,gammarest);
  LorentzRotation boostFromRest(-bv,gammarest);
  // construct a new particle which is at rest
  Particle inrest(inpart);
  inrest.transform(boostToRest);
  double wgt(0.);
  vector<Lorentz5Momentum> momenta(outgoing_.size());
  int ichan;
  unsigned int ncount(0);
  try {
    do {
      // phase-space pieces of the weight
      fillStack();
      wgt = pre*weight(ichan,inrest.momentum(),momenta);
      // matrix element piece
      wgt *= decayer->me2(-1,inrest,!cc ? outgoing_ : outgoingCC_,
			  momenta,
			  ncount ==0 ? DecayIntegrator::Initialize : DecayIntegrator::Calculate);
      ++ncount;
      if(wgt>maxWeight_) {
	CurrentGenerator::log() << "Resetting max weight for decay " 
				<< inrest.PDGName() << " -> ";
	for(tcPDPtr part : outgoing_)
	  CurrentGenerator::log() << "  " << part->PDGName();
	CurrentGenerator::log() << "  " << maxWeight_ << "  " << wgt 
				<< "  " << inrest.mass()/MeV << "\n";
	maxWeight_=wgt;
      }
    }
    while(maxWeight_*UseRandom::rnd()>wgt&&ncount<decayer->nTry_);
  }
  catch (Veto) {
    // restore the incoming particle to its original state
    inrest.transform(boostFromRest);
    while(!rStack_.empty()) rStack_.pop();
    throw Veto();
  }
  if(ncount>=decayer->nTry_) {
    CurrentGenerator::log() << "The decay " << inrest.PDGName() << " -> ";
    for(tcPDPtr part : outgoing_)
      CurrentGenerator::log() << "  " << part->PDGName();
    CurrentGenerator::log() << "  " << maxWeight_ << " " << decayer->nTry_
			    << " is too inefficient for the particle "
			    << inpart << "vetoing the decay \n";
    momenta.clear();
    throw Veto();
  }
  // construct the particles
  ParticleVector output(outgoing_.size());
  for(unsigned int ix=0;ix<outgoing_.size();++ix)
    output[ix] = (!cc ? outgoing_[ix] : outgoingCC_[ix] )->produceParticle(momenta[ix]);
  // set up the spin correlations
  decayer->constructSpinInfo(inrest,output);
  const_ptr_cast<tPPtr>(&inpart)->spinInfo(inrest.spinInfo());
  constructVertex(inpart,output,decayer);
  // return if intermediate particles not required
  if(channels_.empty()||!intermediates) {
    for(tPPtr part : output) part->transform(boostFromRest);
  }
  // find the intermediate particles
  else {
    // select the channel
    vector<double> mewgts(channels_.size(),0.0);
    double total=0.;
    for(unsigned int ix=0,N=channels_.size();ix<N;++ix) {
      mewgts[ix]=decayer->me2(ix,inrest,!cc ? outgoing_ : outgoingCC_,
			      momenta,DecayIntegrator::Calculate);
      total+=mewgts[ix];
    }
    // randomly pick a channel
    total *= UseRandom::rnd();
    int iChannel = -1;
    do {
      ++iChannel;
      total-=mewgts[iChannel];
    }
    while(iChannel<int(channels_.size()) && total>0.);
    iChannel_ = iChannel;
    // apply boost
    for(tPPtr part : output) part->transform(boostFromRest);
    channels_[iChannel].generateIntermediates(cc,inpart,output);
  }
  decayer->ME(DecayMEPtr());
  // return particles;
  return output;
}

// flat phase space generation and weight
Energy PhaseSpaceMode::flatPhaseSpace(const Lorentz5Momentum & in,
				      vector<Lorentz5Momentum> & momenta,
				      bool onShell) const {
  double wgt(1.);
  // masses of the particles
  vector<Energy> mass = externalMasses(in.mass(),wgt,onShell);
  // two body decay
  double ctheta  = 2.*rStack_.top() - 1.;
  rStack_.pop();
  double phi  = Constants::twopi*rStack_.top();
  rStack_.pop();
  if(! Kinematics::twoBodyDecay(in, mass[0], mass[1],
				ctheta, phi, momenta[0],momenta[1])) 
    throw Exception() << "Incoming mass - Outgoing mass negative in "
 		      << "PhaseSpaceMode::flatPhaseSpace()"
		      << Exception::eventerror;
  wgt *= Kinematics::pstarTwoBodyDecay(in.mass(),mass[0],mass[1])
    /8./Constants::pi/in.mass();
  return wgt*in.mass();
}

// generate the masses of the external particles
vector<Energy> PhaseSpaceMode::externalMasses(Energy inmass,double & wgt,
					      bool onShell) const {
  vector<Energy> mass(outgoing_.size());
  vector<int> notdone;
  Energy mlow(ZERO);
  // set masses of stable particles and limits 
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    // get the mass of the particle if can't use weight
    if(onShell)
      mass[ix] = outgoing_[ix]->mass();
    else if(!massGen_[ix] || outgoing_[ix]->stable()) {
      mass[ix] = outgoing_[ix]->generateMass();
      mlow += mass[ix];
    }
    else {
      mass[ix] = ZERO;
      notdone.push_back(ix);
      mlow+=max(outgoing_[ix]->mass()-outgoing_[ix]->widthLoCut(),ZERO);
    }
  }
  if(mlow>inmass) throw Veto();
  // now we need to generate the masses for the particles we haven't
  for( ;!notdone.empty();) {
    unsigned int iloc=long(UseRandom::rnd()*(notdone.size()-1)); 
    Energy low = max(outgoing_[notdone[iloc]]->mass()-outgoing_[notdone[iloc]]->widthLoCut(),ZERO);
    mlow-=low;
    double wgttemp;
    mass[notdone[iloc]]= massGen_[notdone[iloc]]->mass(wgttemp,*outgoing_[notdone[iloc]],low,inmass-mlow,rStack_.top());
    rStack_.pop();
    assert(mass[notdone[iloc]]>=low&&mass[notdone[iloc]]<=inmass-mlow);
    wgt   *= wgttemp;
    mlow  += mass[notdone[iloc]];
    notdone.erase(notdone.begin()+iloc);
  }
  return mass;
}

// construct the vertex for spin corrections
void PhaseSpaceMode::constructVertex(const Particle & inpart, 
				     const ParticleVector & decay,
				     tcDecayIntegratorPtr decayer) const {
  // construct the decay vertex
  VertexPtr vertex(new_ptr(DecayVertex()));
  DVertexPtr Dvertex(dynamic_ptr_cast<DVertexPtr>(vertex));
  // set the incoming particle for the decay vertex
  (inpart.spinInfo())->decayVertex(vertex);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    decay[ix]->spinInfo()->productionVertex(vertex);
  }
  // set the matrix element
  Dvertex->ME(decayer->ME());
  decayer->ME(DecayMEPtr());
}

// initialise the phase space
Energy PhaseSpaceMode::initializePhaseSpace(bool init, tcDecayIntegratorPtr decayer,
					    bool onShell) {
  decayer->ME(DecayMEPtr());
  eps_=decayer->eps_;
  Energy output(ZERO);
  // ensure that the weights add up to one
  if(!channels_.empty()) {
    double temp=0.;
    for(const PhaseSpaceChannel & channel : channels_) temp+= channel.weight();
    for(PhaseSpaceChannel & channel : channels_) channel.weight(channel.weight()/temp);
  }
  if(!init) return ZERO;
  ThePEG::PPtr inpart=incoming_.first->produceParticle();
  // now if using flat phase space
  maxWeight_=0.;
  if(channels_.empty()) {
    vector<Lorentz5Momentum> momenta(outgoing_.size());
    double wsum=0.,wsqsum=0.;
    Energy m0,mmin(ZERO);
    for(tcPDPtr part : outgoing_) mmin += part->massMin();
    for(unsigned int ix=0;ix<decayer->nPoint_;++ix) {
      // set the mass of the decaying particle
      m0 = !onShell ? inpart->dataPtr()->generateMass() : inpart->dataPtr()->mass();
      double wgt=0.;
      if(m0<=mmin) continue;
      inpart->set5Momentum(Lorentz5Momentum(m0));
      // compute the prefactor
      Energy prewid = (widthGen_&&partial_>=0) ?
	widthGen_->partialWidth(partial_,inpart->mass()) :
	incoming_.first->width();
      InvEnergy pre = prewid>ZERO ? 1./prewid : 1./MeV;
      // generate the weight for this point
      try {
	int dummy;
	// phase-space piece
	fillStack();
	wgt = pre*weight(dummy,inpart->momentum(),momenta,onShell);
	// matrix element piece
	wgt *= decayer->me2(-1,*inpart,outgoing_,momenta,DecayIntegrator::Initialize);
      }
      catch (Veto) {
	while(!rStack_.empty()) rStack_.pop();
	wgt=0.;
      }
      if(wgt>maxWeight_) maxWeight_ = wgt;
      wsum   += wgt;
      wsqsum += sqr(wgt);
    }
    wsum /= decayer->nPoint_;
    wsqsum=wsqsum/decayer->nPoint_-sqr(wsum);
    if(wsqsum<0.) wsqsum=0.;
    wsqsum=sqrt(wsqsum/decayer->nPoint_);
    Energy fact = (widthGen_&&partial_>=0) ? 
      widthGen_->partialWidth(partial_,inpart->nominalMass()) :
      inpart->dataPtr()->width();
    if(fact==ZERO) fact=MeV;
     // factor for the weight with spin correlations
    maxWeight_ *= inpart->dataPtr()->iSpin()==1 ? 1.1 : 1.6;
    if ( Debug::level > 1 ) {
      // ouptut the information on the initialisation
      CurrentGenerator::log() << "Initialized the phase space for the decay " 
  			      << incoming_.first->PDGName() << " -> ";
      for(tPDPtr part : outgoing_)
	CurrentGenerator::log() << part->PDGName() << " ";
      CurrentGenerator::log() << "\n";
      if(fact!=MeV) CurrentGenerator::log() << "The branching ratio is " << wsum 
 					    << " +/- " << wsqsum << "\n";
      CurrentGenerator::log() << "The partial width is " << wsum*fact/MeV 
 			      << " +/- " << wsqsum*fact/MeV << " MeV\n";
      CurrentGenerator::log() << "The partial width is " 
 			      << wsum*fact/6.58212E-22/MeV 
 			      << " +/- " << wsqsum*fact/6.58212E-22/MeV<< " s-1\n";
      CurrentGenerator::log() << "The maximum weight is " 
 			      << maxWeight_ << endl;
    }
    output=wsum*fact;
  }
  else {
    vector<Lorentz5Momentum> momenta(outgoing_.size());
    double totsum(0.),totsq(0.);
    for(unsigned int iy=0;iy<decayer->nIter_;++iy) {
      // zero the maximum weight
      maxWeight_=0.;
      vector<double> wsum(channels_.size(),0.),wsqsum(channels_.size(),0.);
      vector<int> nchan(channels_.size(),0);
      totsum = 0.;
      totsq  = 0.;
      Energy mmin(ZERO);
      for(tcPDPtr part : outgoing_) mmin += part->massMin();
      for(unsigned int ix=0;ix<decayer->nPoint_;++ix) {
     	Energy m0 = !onShell ? incoming_.first->generateMass() : incoming_.first->mass();
     	double wgt=0.; 
     	int ichan(-1);
     	if(m0>mmin) {
	  inpart->set5Momentum(Lorentz5Momentum(m0));
	  // compute the prefactor
	  Energy prewid= (widthGen_&&partial_>=0) ? 
	    widthGen_->partialWidth(partial_,inpart->mass()) :
	    inpart->dataPtr()->width();
	  InvEnergy pre = prewid>ZERO ? 1./prewid : 1./MeV;
	  // generate the weight for this point
	  try {
	    fillStack();
	    wgt = pre*weight(ichan,inpart->momentum(),momenta,onShell);
	    // matrix element piece
	    wgt *= decayer->me2(-1,*inpart,outgoing_,momenta,DecayIntegrator::Initialize);
	  }
	  catch (Veto) {
	    wgt=0.;
	  }
     	}
     	if(wgt>maxWeight_) maxWeight_=wgt;
    	if(ichan>=0) {
    	  wsum[ichan]   += wgt;
    	  wsqsum[ichan] += sqr(wgt);
    	  ++nchan[ichan];
    	}
    	totsum+=wgt;
    	totsq+=wgt*wgt;
      }
      totsum=totsum/decayer->nPoint_;
      totsq=totsq/decayer->nPoint_-sqr(totsum);
      if(totsq<0.) totsq=0.;
      totsq=sqrt(totsq/decayer->nPoint_);
      if ( Debug::level > 1 )
    	CurrentGenerator::log() << "The branching ratio is " << iy << " " 
    				<< totsum << " +/- " << totsq 
    				<< maxWeight_ << "\n";
      // compute the individual terms
      double total(0.);
      for(unsigned int ix=0;ix<channels_.size();++ix) {
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
     	total+=sqrt(wsqsum[ix])*channels_[ix].weight();
      }
      if(total>0.) {
    	for(unsigned int ix=0;ix<channels_.size();++ix) {
	  channels_[ix].weight(sqrt(wsqsum[ix])*channels_[ix].weight()/total);
     	}
      }
    }
    // factor for the weight with spin correlations
    maxWeight_*= inpart->dataPtr()->iSpin()==1 ? 1.1 : 1.6;
    // output the information on the initialisation
    Energy fact = (widthGen_&&partial_>=0) ? 
      widthGen_->partialWidth(partial_,inpart->nominalMass()) :
      inpart->dataPtr()->width();
    output=totsum*fact;
    if(fact==ZERO) fact=MeV;
    if ( Debug::level > 1 ) {
      CurrentGenerator::log() << "Initialized the phase space for the decay " 
     			      << incoming_.first->PDGName() << " -> ";
      for(tcPDPtr part : outgoing_)
	CurrentGenerator::log() << part->PDGName() << " ";
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
     			      << maxWeight_ << "\n";
      CurrentGenerator::log() << "The weights for the different phase" 
			      << " space channels are \n";
      for(unsigned int ix=0,N=channels_.size();ix<N;++ix) {
    	CurrentGenerator::log() << "Channel " << ix 
    				<< " had weight " << channels_[ix].weight()
    				<< "\n";
      }
    }
    CurrentGenerator::log() << flush;
  }
  decayer->ME(DecayMEPtr());
  return output;
}
 
// generate a phase-space point using multichannel phase space
Energy PhaseSpaceMode::channelPhaseSpace(int & ichan, const Lorentz5Momentum & in, 
					 vector<Lorentz5Momentum> & momenta,
					 bool onShell) const {
  double wgt(rStack_.top());
  rStack_.pop();
  // select a channel
  ichan=-1;
  do {
    ++ichan;
    wgt-=channels_[ichan].weight();
  }
  while(ichan<int(channels_.size())&&wgt>0.);
  // generate the momenta
  if(ichan==int(channels_.size())) {
    throw Exception() << "PhaseSpaceMode::channelPhaseSpace()"
		      << " failed to select a channel" 
		      << Exception::abortnow;
  }
  // generate the masses of the external particles
  double masswgt(1.);
  vector<Energy> mass(externalMasses(in.mass(),masswgt,onShell));
  momenta=channels_[ichan].generateMomenta(in,mass);
  // compute the denominator of the weight
  wgt=0.;
  for(const PhaseSpaceChannel & channel : channels_) 
    wgt += channel.generateWeight(momenta);
  // return the weight
  return wgt!=0. ? in.mass()*masswgt/wgt : ZERO;
}

void PhaseSpaceMode::init() {
  // get mass and width generators
  outgoingCC_.clear();
  for(tPDPtr part : outgoing_) {
    outgoingCC_.push_back(part->CC() ? part->CC() : part);
  }
  massGen_.resize(outgoing_.size());
  widthGen_ = dynamic_ptr_cast<cGenericWidthGeneratorPtr>(incoming_.first->widthGenerator());
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    assert(outgoing_[ix]);
    massGen_[ix]= dynamic_ptr_cast<cGenericMassGeneratorPtr>(outgoing_[ix]->massGenerator());
  }
  // get max energy for decays
  if(!incoming_.second)
    eMax_ = testOnShell_ ? incoming_.first->mass() : incoming_.first->massMax();
  // work out how many random numbers we need
  nRand_ = 3*outgoing_.size()-4;
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    if(massGen_[ix]) ++nRand_;
  }
  if(channels_.empty()) return;
  ++nRand_;
  // ensure weights sum to one
  double sum(0.);
  for(const PhaseSpaceChannel & channel : channels_)
    sum+=channel.weight();
  for(PhaseSpaceChannel & channel : channels_)
    channel.weight(channel.weight()/sum);
}

void PhaseSpaceMode::initrun() {
  // update the mass and width generators
  if(incoming_.first->widthGenerator()!=widthGen_)
    widthGen_ = dynamic_ptr_cast<cGenericWidthGeneratorPtr>(incoming_.first->widthGenerator());
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    if(massGen_[ix]!=outgoing_[ix]->massGenerator())
      massGen_[ix] = dynamic_ptr_cast<cGenericMassGeneratorPtr>(outgoing_[ix]->massGenerator());
  }
  for(PhaseSpaceChannel & channel : channels_) channel.initrun(this);
  if(widthGen_) const_ptr_cast<GenericWidthGeneratorPtr>(widthGen_)->initrun();
  tcGenericWidthGeneratorPtr wtemp;
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    wtemp=
      dynamic_ptr_cast<tcGenericWidthGeneratorPtr>(outgoing_[ix]->widthGenerator());
    if(wtemp) const_ptr_cast<tGenericWidthGeneratorPtr>(wtemp)->initrun();
  }
  // ensure weights sum to one
  double sum(0.);
  for(const PhaseSpaceChannel & channel : channels_)
    sum+=channel.weight();
  for(PhaseSpaceChannel & channel : channels_)
    channel.weight(channel.weight()/sum);
  nRand_ = 3*outgoing_.size()-4;
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    if(massGen_[ix]) ++nRand_;
  }
  if(channels_.empty()) return;
  ++nRand_;
}
