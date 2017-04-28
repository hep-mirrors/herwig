// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EvtGenDecayer class.
//

#include "EvtGenDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

bool EvtGenDecayer::accept(const DecayMode &) const {
  return true;
}

ParticleVector EvtGenDecayer::decay(const DecayMode & dm,
				    const Particle & parent) const {
  unsigned int ntry=0;
  ParticleVector output;
  do {
    if(evtOpt_==0)
      output=evtgen_->decay(parent,false,dm);
    else if(evtOpt_==1)
      output=evtgen_->decay(parent, true,dm);
    else
      throw Exception() << "Unknown option in EvtGenDecayer::decay() " 
			<< Exception::runerror;
    ++ntry;
  }
  while(output.empty()&&ntry<10);
  if(output.empty()) {
      throw Exception() << "EvtGenDecayer::decay() failed to decay"
			<< parent
			<< Exception::eventerror;
  }
  if(check_) {
    Lorentz5Momentum ptotal=parent.momentum();
    int charge=parent.dataPtr()->iCharge();
    for(unsigned int ix=0;ix<output.size();++ix) {
      ptotal-=output[ix]->momentum();
      charge-=output[ix]->dataPtr()->iCharge();
      checkDecay(output[ix]);
    }
    if(abs(ptotal.x())>0.001*MeV||abs(ptotal.y())>0.001*MeV||
       abs(ptotal.z())>0.001*MeV||abs(ptotal.e())>0.001*MeV) {
      if(check_==1 || !rescale(parent,output)) {
	generator()->log() << "Decay of " << parent.PDGName()<<" -> ";
	for (auto const & dec : output)
	  generator()->log()<< dec->dataPtr()->PDGName()<<" ";
	generator()->log() << " violates momentum conservation in"
			   << " EvtGenDecayer::decay in event "
			   << generator()->currentEventNumber() << "\n";
      }
    }
    if(charge!=0) {
      generator()->log() << "Decay of " << parent.PDGName() 
			 << " violates charge conservation in"
			 << " EvtGenDecayer::decay in event "
			 << generator()->currentEventNumber() << "\n";
    }
  }
  return output;
}


IBPtr EvtGenDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr EvtGenDecayer::fullclone() const {
  return new_ptr(*this);
}

void EvtGenDecayer::persistentOutput(PersistentOStream & os) const {
  os << evtgen_ << check_ << evtOpt_;
}

void EvtGenDecayer::persistentInput(PersistentIStream & is, int) {
  is >> evtgen_ >> check_ >> evtOpt_;  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EvtGenDecayer,Decayer>
describeHerwigEvtGenDecayer("Herwig::EvtGenDecayer",
			    "HwEvtGenInterface.so");

void EvtGenDecayer::Init() {

  static ClassDocumentation<EvtGenDecayer> documentation
    ("The EvtGenDecayer class allows the EvtGen decay package to be used as"
     " a decayer inside Herwig");

  static Reference<EvtGenDecayer,EvtGenInterface> interfaceEvtGen
    ("EvtGen",
     "Pointer to the EvtGenInterface object which encapsulates the EvtGen decay package.",
     &EvtGenDecayer::evtgen_, false, false, true, false, false);

  static Switch<EvtGenDecayer,unsigned int> interfaceCheck
    ("Check",
     "Perform some basic checks of the decay",
     &EvtGenDecayer::check_, false, false, false);
  static SwitchOption interfaceCheckCheck
    (interfaceCheck,
     "Yes",
     "Perform the checks",
     1);
  static SwitchOption interfaceCheckCheckRescale
    (interfaceCheck,
     "Rescale",
     "Perform the checks and rescale if momentum violation",
     2);
  static SwitchOption interfaceCheckNoCheck
    (interfaceCheck,
     "No",
     "Don't perform the checks",
     0);

  static Switch<EvtGenDecayer,unsigned int> interfaceOption
    ("Option",
     "The way in which EvtGen is used.",
     &EvtGenDecayer::evtOpt_, 0, false, false);
  static SwitchOption interfaceOptionParent
    (interfaceOption,
     "Parent",
     "EvtGen decays the particle and returns the decay products to be decayed by"
     " Herwig.",
     0);
  static SwitchOption interfaceOptionAll
    (interfaceOption,
     "All",
     "EvtGen decays the particle and all the unstable particles produced in the decay.",
     1);
}

void EvtGenDecayer::checkDecay(PPtr in) const {
  Lorentz5Momentum ptotal=in->momentum();
  int charge = in->dataPtr()->iCharge();
  if(in->children().empty()) return;
  for(unsigned int ix=0;ix<in->children().size();++ix) {
    checkDecay(in->children()[ix]);
    ptotal-=in->children()[ix]->momentum();
    charge-=in->children()[ix]->dataPtr()->iCharge();
  }
  if(abs(ptotal.x())>0.001*MeV||abs(ptotal.y())>0.001*MeV||
     abs(ptotal.z())>0.001*MeV||abs(ptotal.e())>0.001*MeV) {
    if(check_==1 || !rescale(*in,in->children())) {
      generator()->log() 
	<< "SubDecay of " << in->PDGName() <<" -> ";
      for (auto const & dec : in->children())
	generator()->log()<< dec->dataPtr()->PDGName()<<" ";
      generator()->log() <<" violates momentum conservation"
			 << " in EvtGenDecayer::checkDecay in event "
			 << generator()->currentEventNumber() << "\n";
    }
  }
  if(charge!=0) generator()->log() << "Decay of " << in->PDGName() 
				   << " violates charge conservation in "
				   << "EvtGenDecayer::checkDecay in event "
				   << generator()->currentEventNumber() << "\n";
}
namespace {

// momentum test function and derivative for N-R method
pair<Energy,Energy> momTest(double lambda,
			    const vector<Energy2> & p2,
			    const vector<Energy2> & m2) {
  pair<Energy,Energy> output;
  for(unsigned int ix=0;ix<p2.size();++ix) {
    Energy en = sqrt(lambda*p2[ix]+m2[ix]);
    output.first  += en;
    output.second += 0.5*p2[ix]/en;
  }
  return output;
}

// rescaling boost from qtilde shower
LorentzRotation solveBoost(const Lorentz5Momentum & q, 
			   const Lorentz5Momentum & p ) {
  Energy modp = p.vect().mag();
  Energy modq = q.vect().mag();
  double betam = (p.e()*modp-q.e()*modq)/(sqr(modq)+sqr(modp)+p.mass2());
  assert( abs(betam)-1. < 0. );
  Boost beta = -betam*q.vect().unit();
  ThreeVector<Energy2> ax = p.vect().cross( q.vect() ); 
  double delta = p.vect().angle( q.vect() );
  LorentzRotation R;
  using Constants::pi;
  assert( beta.mag2() - 1. < 0. );
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta );
  } 
  else {
    R.boost( beta );
  } 
  return R;
}
  
}

bool EvtGenDecayer::rescale(const Particle & in,
			    const ParticleVector & children) const {
  LorentzRotation rot(-in.momentum().boostVector());
  unsigned int npi(0),ngamma(0);
  Lorentz5Momentum psum;
  for(unsigned int ix=0;ix<children.size();++ix) {
    long id = children[ix]->id();
    if(id==ParticleID::gamma)
      ++ngamma;
    else if(abs(id)==ParticleID::piplus || id==ParticleID::pi0)
      ++npi;
    psum+=children[ix]->momentum();
  }
  // two dodgy cases from EvtGen
  int flavour = (abs(in.id())%1000)/10;
  bool isBd = (abs(in.id())==ParticleID::Bplus || abs(in.id())==ParticleID::B0);
  // B -> 3pi
  bool bDecay = isBd && npi==3 && children.size()==3;
  bool photon = (flavour==44 || flavour==55 || isBd || abs(in.id())==ParticleID::B_s0) && ngamma ==1;
  // not a dodgy case return
  if(!bDecay && !photon) return false;
  // ensure in rest frame of system
  psum *=rot;
  LorentzRotation rotInv = rot.inverse();
  rot.boost(-psum.boostVector());
  // setup Newton-Raphson for rescaling
  Lorentz5Momentum psum2;
  vector<Lorentz5Momentum> pold;
  vector<Energy2> p2,m2;
  for(unsigned int ix=0;ix<children.size();++ix) {
    pold.push_back(rot*children[ix]->momentum());
    psum2+=pold.back();
    p2.push_back(pold.back().vect().mag2());
    m2.push_back(sqr(children[ix]->mass()));
  }
  // calculate the rescaling using N-R
  double lambda = 1.;
  pair<Energy,Energy> test = momTest(lambda,p2,m2);
  test.first -=in.mass();
  do {
    lambda -= test.first/test.second;
    test = momTest(lambda,p2,m2);
    test.first -=in.mass();
  }
  while(abs(test.first/MeV)>1e-8);
  lambda = sqrt(lambda);
  // apply the rescaleing boost
  Lorentz5Momentum ptest;
  for(unsigned int ix=0;ix<children.size();++ix) {
    Lorentz5Momentum pnew(lambda*pold[ix].x(),lambda*pold[ix].y(),
			  lambda*pold[ix].z(),sqrt(sqr(lambda)*p2[ix]+m2[ix]),
			  children[ix]->mass());
    children[ix]->deepTransform(rotInv*solveBoost(pnew,pold[ix])*rot);
  }
  return true;
}
