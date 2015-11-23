// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnomalousWWHVertex class.
//

#include "AnomalousWWHVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

AnomalousWWHVertex::AnomalousWWHVertex() 
  : interactionType_(0), Lambda_(1000.*GeV),
    couplast_(0.), q2last_(ZERO), mw_(ZERO), zfact_(0.) {
  // particles
  addToList(24,-24,25);
  addToList(23,23,25);
  // calculate the kinematic invariants needed
  kinematics(true);
}

IBPtr AnomalousWWHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr AnomalousWWHVertex::fullclone() const {
  return new_ptr(*this);
}

void AnomalousWWHVertex::doinit() {
  // parameters
  mw_ = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  zfact_ = 1./(1.-generator()->standardModel()->sin2ThetaW());
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  // base class
  GeneralVVSVertex::doinit();
}

void AnomalousWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << interactionType_ << ounit(mw_,GeV) << zfact_ << ounit(Lambda_,GeV);
}

void AnomalousWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> interactionType_ >> iunit(mw_,GeV) >> zfact_ >> iunit(Lambda_,GeV);
}

ClassDescription<AnomalousWWHVertex> AnomalousWWHVertex::initAnomalousWWHVertex;
// Definition of the static class description member.

void AnomalousWWHVertex::Init() {

  static ClassDocumentation<AnomalousWWHVertex> documentation
    ("The AnomalousWWHVertex class implemenets the Higgs coupling to two electroweak"
     " vector bosons including the option of anomalous couplings.");

  static Switch<AnomalousWWHVertex,unsigned int> interfaceInteractionType
    ("InteractionType",
     "The type of interaction",
     &AnomalousWWHVertex::interactionType_, 0, false, false);
  static SwitchOption interfaceInteractionTypeSM
    (interfaceInteractionType,
     "SM",
     "Use the standard model form",
     0);
  static SwitchOption interfaceInteractionTypeCPOdd
    (interfaceInteractionType,
     "CPOdd",
     "Use a CP odd form",
     1);
  static SwitchOption interfaceInteractionTypeCPEven
    (interfaceInteractionType,
     "CPEven",
     "Use a CP even form",
     2);

  static Parameter<AnomalousWWHVertex,Energy> interfaceLambda
    ("Lambda",
     "The scale of new physics",
     &AnomalousWWHVertex::Lambda_, GeV, 1000.0*GeV, ZERO, 1000000.0*GeV,
     false, false, Interface::limited);

}

void AnomalousWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr) {
  int ibos=abs(a->id());
  a00(0.);
  a11(0.);
  a12(0.);
  a21(0.);
  a22(0.);
  aEp(0.);
  switch(interactionType_) {
  case 0:
    // first the overall normalisation
    if(q2!=q2last_) {
      couplast_ = weakCoupling(q2) * UnitRemoval::InvE * mw_;
      q2last_   = q2;
    }
    if(ibos==24)      norm(couplast_          );
    else if(ibos==23) norm(couplast_ * zfact_ );
    else
      throw HelicityConsistencyError() << "AnomalousWWHVertex::setCoupling "
				       << "Invalid particles in WWH Vertex" 
				       << Exception::runerror;
    a00(double(UnitRemoval::E2/invariant(1,2)));
    break;
  case 1:
    norm(double(UnitRemoval::E/Lambda_));
    aEp(1.);
    break;
  case 2:
    norm(double(UnitRemoval::E/Lambda_));
    a00( 1.);
    a21(-1.);
    break;
  default:
    throw HelicityConsistencyError() << "AnomalousWWHVertex::setCoupling "
				     << "Unknown type of interaction" 
				     << Exception::runerror;
  }
}
