// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeHyperonDecayer class.
//

#include "RadiativeHyperonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeHyperonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

void RadiativeHyperonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(incomingB_.size());
  if(isize!=outgoingB_.size()||isize!=A_.size()||
     isize!=B_.size()        ||isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeHyperonDecayer::doinit()" 
			  << Exception::runerror;
  // set up the decay modes
  tPDPtr photon = getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    tPDPtr in = getParticleData(incomingB_[ix]);
    tPDVector out={getParticleData(outgoingB_[ix]),photon};
    double wgtmax = maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

int RadiativeHyperonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing pa4rticles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0=parent->id();
  int id1(children[0]->id());
  int id2(children[1]->id());
  if(id1==ParticleID::gamma) swap(id1,id2);
  if(id2!=ParticleID::gamma) return imode;
  unsigned int ix(0);
  do {
    if(id0==incomingB_[ix] && id1==outgoingB_[ix]) imode=ix;
    else if(id0==-incomingB_[ix] && id1==-outgoingB_[ix]) imode=ix;
    ++ix;
  }
  while(ix<incomingB_.size()&&imode<0);
  // charge conjugate
  cc=id0<0;
  // return the answer
  return imode;
}

IBPtr RadiativeHyperonDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr RadiativeHyperonDecayer::fullclone() const {
  return new_ptr(*this);
}

void RadiativeHyperonDecayer::persistentOutput(PersistentOStream & os) const {
  os << incomingB_ << outgoingB_ << ounit(A_,1./GeV) << ounit(B_,1./GeV) 
     << maxweight_;
}

void RadiativeHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incomingB_ >> outgoingB_ >> iunit(A_,1./GeV) >> iunit(B_,1./GeV) 
     >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeHyperonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeHyperonDecayer("Herwig::RadiativeHyperonDecayer", "HwBaryonDecay.so");

void RadiativeHyperonDecayer::Init() {

  static ClassDocumentation<RadiativeHyperonDecayer> documentation
    ("The RadiativeHyperonDecayer class performs the radiative decays of hyperons.",
     "The radiative hyperons decays were simulated using the RadiativeHyperonDecayer"
     " class which implements the results of \\cite{Borasoy:1999nt}.",
     "\\bibitem{Borasoy:1999nt}\n"
     "B.~Borasoy and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 59} (1999) 054019 [arXiv:hep-ph/9902431].\n"
     "%%CITATION = PHRVA,D59,054019;%%\n");

  static Command<RadiativeHyperonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryon , A and B couplings (1/GeV) and max weight for a decay",
     &RadiativeHyperonDecayer::setUpDecayMode, false);
  
  static Deleted<RadiativeHyperonDecayer> interfaceMaxWeight
    ("MaxWeight",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceIncomingBaryon
    ("IncomingBaryon",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceOutgoingBaryon
    ("OutgoingBaryon",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceCouplingA
    ("CouplingA",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceCouplingB
    ("CouplingB",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

}

void RadiativeHyperonDecayer::halfHalfVectorCoupling(int imode, Energy m0, Energy m1,
						     Energy,Complex& A1,Complex& A2,
						     Complex& B1,Complex& B2) const {
  useMe();
  A1 = A_[imode]*(m0+m1);
  B1 = B_[imode]*(m0-m1);
  A2 = 2.*A_[imode]*(m0+m1);
  B2 = 2.*B_[imode]*(m0+m1);
}

void RadiativeHyperonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incomingB_[ix] << " " << outgoingB_[ix] << " "
	   << A_[ix]*GeV << " " << B_[ix]*GeV << " "
	   << maxweight_[ix] << "\n";
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}


string RadiativeHyperonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out;
  out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 1/2";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  A_.push_back(g/GeV);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  g = stof(stype);
  B_.push_back(g/GeV);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incomingB_.push_back(in);
  outgoingB_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
