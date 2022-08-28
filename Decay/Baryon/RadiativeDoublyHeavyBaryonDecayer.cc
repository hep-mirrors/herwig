// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeDoublyHeavyBaryonDecayer class.
//

#include "RadiativeDoublyHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeDoublyHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxWeight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) maxWeight_.push_back(mode(ix)->maxWeight());
      else         maxWeight_.push_back(1.);
    }
  }
}

void RadiativeDoublyHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the parameters are consistent
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size() ||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeDoublyHeavyBaryonDecayer::doinit()" 
			  << Exception::abortnow;
  // the decay modes
  tPDPtr photon = getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix]),photon};
    PhaseSpaceModePtr mode;
    if(in&&out[0]) {
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else {
      mode = PhaseSpaceModePtr();
    }
    addMode(mode);
  }
}

int RadiativeDoublyHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
						  const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id()),ibaryon;
  if(id1==ParticleID::gamma){ibaryon=id2;}
  else if(id2==ParticleID::gamma){ibaryon=id1;}
  else {return imode;}
  unsigned int ix(0);
  do {
    if(     id0== incoming_[ix]&&ibaryon== outgoing_[ix]) {
      imode=ix;
      cc=false;
    }
    else if(id0==-incoming_[ix]&&ibaryon==-outgoing_[ix]) {
      imode=ix;
      cc=true;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void RadiativeDoublyHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1./GeV2) << incoming_ << outgoing_ << maxWeight_;
}

void RadiativeDoublyHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1./GeV2) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeDoublyHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeDoublyHeavyBaryonDecayer("Herwig::RadiativeDoublyHeavyBaryonDecayer", "HwBaryonDecay.so");

void RadiativeDoublyHeavyBaryonDecayer::Init() {

  static ClassDocumentation<RadiativeDoublyHeavyBaryonDecayer> documentation
    ("The RadiativeDoublyHeavyBaryonDecayer class is designed for the radiative decays of"
     " doubly heavy baryons.");

  static Command<RadiativeDoublyHeavyBaryonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryons, M1(1/GeV2) coupling and max weight for a decay",
     &RadiativeDoublyHeavyBaryonDecayer::setUpDecayMode, false);

}

void RadiativeDoublyHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix] << " " << " " 
	   << coupling_[ix]*GeV2 << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void RadiativeDoublyHeavyBaryonDecayer::halfHalfVectorCoupling(int imode,Energy m0,Energy m1,
							       Energy,Complex& A1,
							       Complex& A2,Complex& B1,
							       Complex& B2) const {
  useMe();
  A1 = -0.5*sqr(m0+m1)*coupling_[imode];
  A2 =      sqr(m0+m1)*coupling_[imode];
  B1 = 0.;
  B2 = 0.;
}

void RadiativeDoublyHeavyBaryonDecayer::threeHalfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A1,Complex& A2,
							      Complex& A3,Complex& B1,
							      Complex& B2,
							      Complex& B3) const {
  useMe();
  A1 = 0.;
  A2 = 0.;
  A3 = 0.;
  B1 = 0.5*sqr(m0+m1)*coupling_[imode];
  B2 =- m0*   (m0+m1)*coupling_[imode];
  B3 =-    sqr(m0+m1)*coupling_[imode];
}

string RadiativeDoublyHeavyBaryonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2 or 3/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "Outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Outgoing particle with id " + std::to_string(out) + "does not have spin 1/2";
  // M1 coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy2 g = stof(stype)/GeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
