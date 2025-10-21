// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Gto3S1OctetSplitFn class.
//

#include "Gto3S1OctetSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Shower/QTilde/Kinematics/FS_QTildeShowerKinematics1to1.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

IBPtr Gto3S1OctetSplitFn::clone() const {
  return new_ptr(*this);
}

IBPtr Gto3S1OctetSplitFn::fullclone() const {
  return new_ptr(*this);
}

void Gto3S1OctetSplitFn::doinit() {
  SudakovFormFactor::doinit();
  int iq=4+state_;
  long pid = iq*110+3 + (n_-1)*100000;
  O8_ = params_->octetMEProduction<0>(state_,1,1,pid);
  m_ = getParticleData(4+state_)->mass();
}

void Gto3S1OctetSplitFn::persistentOutput(PersistentOStream & os) const {
  os << params_ << ounit(O8_,GeV*GeV2) << ounit(m_,GeV) << oenum(state_) << n_
     << enhancementFactor_;
}

void Gto3S1OctetSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> params_ >> iunit(O8_,GeV*GeV2) >> iunit(m_,GeV) >> ienum(state_) >> n_
      >> enhancementFactor_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Gto3S1OctetSplitFn,SudakovFormFactor>
  describeHerwigGto3S1OctetSplitFn("Herwig::Gto3S1OctetSplitFn",
				    "HwOniumParameters.so HwOniumShower.so");

void Gto3S1OctetSplitFn::Init() {

  static ClassDocumentation<Gto3S1OctetSplitFn> documentation
    ("The Gto3S1OctetSplitFn class does the splitting g-> 3S1 in an octet state");

  static Reference<Gto3S1OctetSplitFn,OniumParameters> interfaceParameters
    ("Parameters",
     "Quarkonium parameters",
     &Gto3S1OctetSplitFn::params_, false, false, true, false, false);
  
  static Switch<Gto3S1OctetSplitFn,OniumState> interfaceState
    ("State",
     "The type of onium state",
     &Gto3S1OctetSplitFn::state_, ccbar, false, false);
  static SwitchOption interfaceStateccbar
    (interfaceState,
     "ccbar",
     "Charmonium state",
     ccbar);
  static SwitchOption interfaceStatebbbar
    (interfaceState,
     "bbbar",
     "Bottomonium state",
     bbbar);
  
  static Parameter<Gto3S1OctetSplitFn,unsigned int> interfacePrincipalQuantumNumber
    ("PrincipalQuantumNumber",
     "The principle quantum number of the states",
     &Gto3S1OctetSplitFn::n_, 1, 1, 10,
     false, false, Interface::limited);

  static Parameter<Gto3S1OctetSplitFn,double> interfaceEnhancementFactor
    ("EnhancementFactor",
     "Factor by which to enhance the splitting, compenstated by weight events.",
     &Gto3S1OctetSplitFn::enhancementFactor_, 1., 0.0, 1e6,
     false, false, Interface::limited);
}

ShoKinPtr Gto3S1OctetSplitFn::generateNextTimeBranching(const Energy ,
							const IdList &ids,
							const RhoDMatrix &,
							double , double ) {
  int iq = (ids[1]->id()%1000)/100;
  Energy mq = getParticleData(iq)->mass();
  double prob = Constants::pi*alpha()->value(4.*sqr(mq))/24.*O8_/pow<3,1>(mq);
  if(prob*enhancementFactor_>1.) {
    throw Exception() << "g-> 3S1 branching probability too large in Gto3S1OctetSplitFn::generateNextTimeBranching"
		      << " for " << fullName() << Exception::runerror;
  }
  // access the shower handler and step enhancement
  double enhance = 1.;
  tShowerHandlerPtr ch;
  if(ShowerHandler::currentHandlerIsSet()) {
    ch = ShowerHandler::currentHandler();
    enhance =enhancementFactor_;
  }
  if(UseRandom::rnd()<=enhance*prob) {
    if(enhance!=1.) ch->reweight(ch->reweight()/enhance);
    return new_ptr(FS_QTildeShowerKinematics1to1(2.*ids[1]->mass(),this)); 
  }
  else {
    if(enhance!=1.) ch->reweight(ch->reweight()*(1.-prob)/(1.-enhance*prob));
    return ShoKinPtr();
  }
}
