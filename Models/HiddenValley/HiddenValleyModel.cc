// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiddenValleyModel class.
//

#include "HiddenValleyModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HiddenValleyModel::HiddenValleyModel() : groupType_(SU), Nc_(3), Nf_(1), 
					 qL_(-0.2), uR_(-0.2), dR_(0.6),
					 lL_(0.6), lR_(-0.2), gPrime_(1./7.),
					 qCharge_(1,-0.2)
{}

IBPtr HiddenValleyModel::clone() const {
  return new_ptr(*this);
}

IBPtr HiddenValleyModel::fullclone() const {
  return new_ptr(*this);
}

void HiddenValleyModel::persistentOutput(PersistentOStream & os) const {
  os << oenum(groupType_) << Nc_ << Nf_ << FFZPVertex_
     << qL_ << uR_ << dR_ << lL_ << lR_ << qCharge_ << gPrime_;
}

void HiddenValleyModel::persistentInput(PersistentIStream & is, int) {
  is >> ienum(groupType_) >> Nc_ >> Nf_ >> FFZPVertex_
     >> qL_ >> uR_ >> dR_ >> lL_ >> lR_ >> qCharge_ >> gPrime_;
}

ClassDescription<HiddenValleyModel> HiddenValleyModel::initHiddenValleyModel;
// Definition of the static class description member.

void HiddenValleyModel::Init() {

  static ClassDocumentation<HiddenValleyModel> documentation
    ("The HiddenValleyModel class is the main class for the Hidden Valley Model.");

  static Switch<HiddenValleyModel,GroupType> interfaceGroupType
    ("GroupType",
     "Type of the new unbroken group",
     &HiddenValleyModel::groupType_, SU, false, false);
  static SwitchOption interfaceGroupTypeSU
    (interfaceGroupType,
     "SU",
     "Use an SU(N) group",
     SU);

  static Parameter<HiddenValleyModel,unsigned int> interfaceGroupOrder
    ("GroupOrder",
     "The value of the number of 'colours' for the new symmetry group",
     &HiddenValleyModel::Nc_, 3, 1, 1000,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,unsigned int> interfaceNumberOfFermions
    ("NumberOfFermions",
     "The number of fermions charged under the new group",
     &HiddenValleyModel::Nf_, 1, 1, 100,
     false, false, Interface::limited);

  static Reference<HiddenValleyModel,AbstractFFVVertex> interfaceVertexFFZPrime
    ("Vertex/FFZPrime",
     "The vertex coupling the Zprime to the fermions",
     &HiddenValleyModel::FFZPVertex_, false, false, true, false, false);

  static ParVector<HiddenValleyModel,double> interfaceQuirkCharges
    ("QuirkCharges",
     "The charges under the new U(1) of the new quirks",
     &HiddenValleyModel::qCharge_, 1, -0.2, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,double> interfaceZPrimeQL
    ("ZPrime/QL",
     "The charge of the left-handed quarks under the new U(1)",
     &HiddenValleyModel::qL_, -0.2, 10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,double> interfaceZPrimeUR
    ("ZPrime/UR",
     "The charge of the right-handed up quarks under the new U(1)",
     &HiddenValleyModel::uR_, -0.2, 10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,double> interfaceZPrimeDR
    ("ZPrime/DR",
     "The charge of the right-handed down quarks under the new U(1)",
     &HiddenValleyModel::uR_, 0.6, 10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,double> interfaceZPrimeLL
    ("ZPrime/LL",
     "The charge of the left-handed leptons under the new U(1)",
     &HiddenValleyModel::lL_, 0.6, 10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,double> interfaceZPrimeLR
    ("ZPrime/LR",
     "The charge of the right-handed charged leptons under the new U(1)",
     &HiddenValleyModel::lR_, -0.2, 10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyModel,double> interfaceGPrime
    ("GPrime",
     "The new gauge coupling",
     &HiddenValleyModel::gPrime_, 1./7., 0.0, 10.0,
     false, false, Interface::limited);

}

void HiddenValleyModel::doinit() {
  addVertex(FFZPVertex_);
  StandardModel::doinit();
  FFZPVertex_->init();
  if(qCharge_.size()!=Nf_) 
    throw InitException() << "Number of new fermions and  number of new charges"
			  << "different in HiddenValleyModel::doinit()" 
			  << Exception::runerror;
  for(unsigned int ix=0;ix<qCharge_.size();++ix)
    cerr << qCharge_[ix] << "\n";
}
