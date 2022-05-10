// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeHeavyBaryonDecayer class.
//

#include "RadiativeHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxWeight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) maxWeight_.push_back(mode(ix)->maxWeight());
      else         maxWeight_.push_back(1.);
    }
  }
}

void RadiativeHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the parameters are consistent
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size() ||isize!=maxWeight_.size()||isize!=M1Coupling_.size()||
     isize!=E1Coupling_.size()||isize!=modeType_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeHeavyBaryonDecayer::doinit()" 
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

int RadiativeHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					    const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
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

void RadiativeHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(M1Coupling_,1./GeV) << ounit(E1Coupling_,1./GeV2) 
     << incoming_ << outgoing_ << modeType_ << maxWeight_;
}

void RadiativeHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(M1Coupling_,1./GeV) >> iunit(E1Coupling_,1./GeV2) 
     >> incoming_ >> outgoing_ >> modeType_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeHeavyBaryonDecayer("Herwig::RadiativeHeavyBaryonDecayer", "HwBaryonDecay.so");

void RadiativeHeavyBaryonDecayer::Init() {

  static ClassDocumentation<RadiativeHeavyBaryonDecayer> documentation
    ("The RadiativeHeavyBaryonDecayer class is designed for the radiative decays of"
     " heavy baryons.",
     "The radiative decays of the heavy baryons were simulated using the results of"
     "\\cite{Ivanov:1999bk,Ivanov:1998wj}.",
     "\\bibitem{Ivanov:1999bk}\n"
     "M.~A.~Ivanov, J.~G.~Korner, V.~E.~Lyubovitskij and A.~G.~Rusetsky,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 094002\n"
     "[arXiv:hep-ph/9904421].\n"
     "%%CITATION = PHRVA,D60,094002;%%\n"
     "\\bibitem{Ivanov:1998wj}\n"
     "M.~A.~Ivanov, J.~G.~Korner and V.~E.~Lyubovitskij,\n"
     "Phys.\\ Lett.\\  B {\\bf 448} (1999) 143 [arXiv:hep-ph/9811370].\n"
     "%%CITATION = PHLTA,B448,143;%%\n");

  static Command<RadiativeHeavyBaryonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryons, type of decay, M1(1/MeV) and E1(1/MeV2) couplings and max weight for a decay",
     &RadiativeHeavyBaryonDecayer::setUpDecayMode, false);

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceOutgoingB
    ("OutgoingB","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceModeType
    ("ModeType","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<RadiativeHeavyBaryonDecayer> interfaceM1Coupling
    ("M1Coupling","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceE1Coupling
    ("E1Coupling","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");
}

void RadiativeHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix] << " " << " " << modeType_[ix]
	   << M1Coupling_[ix]*MeV << " " << E1Coupling_[ix]*MeV2 << " "
	   << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void RadiativeHeavyBaryonDecayer::halfHalfVectorCoupling(int imode,Energy m0,Energy m1,
							 Energy,Complex& A1,
							 Complex& A2,Complex& B1,
							 Complex& B2) const {
  useMe();
  if(modeType_[imode]==0) {
    B1=-0.5*(sqr(m0) - sqr(m1))*E1Coupling_[imode];
    B2=-    (m0+m1)*(m0+m1)*E1Coupling_[imode];
    A1=0.;
    A2=0.;
  }
  else if(modeType_[imode]==1) {
    A1=-2.*(m0+m1)*M1Coupling_[imode];
    A2= 4.*(m0+m1)*M1Coupling_[imode];
    B1=0.;
    B2=0.;
  }
  else {
    throw Exception() << "Unknown type of mode " << modeType_[imode] 
				 << " in RadiativeHeavyBaryonDecayer::"
				 << "halfHalfVectorCoupling()" << Exception::abortnow;
  }
}

void RadiativeHeavyBaryonDecayer::threeHalfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A1,Complex& A2,
							      Complex& A3,Complex& B1,
							      Complex& B2,
							      Complex& B3) const {
  useMe();
  if(modeType_[imode]==0) {
    A1=-0.5*(m0*m0 - m1*m1)*E1Coupling_[imode];
    A3=-    (m0+m1)*(m0+m1)*E1Coupling_[imode];    
    A2=0.;
    B1=0.;
    B2=0.;
    B3=0.;
  }
  else if(modeType_[imode]==1) {
    B1=-(m0+m1)*M1Coupling_[imode];
    B2= (m0+m1)*M1Coupling_[imode];
    B3=0.;
    A1=0.;
    A2=0.;
    A3=0.;
  }
  else {
    throw Exception() << "Unknown type of mode " << modeType_[imode] 
				 << " in RadiativeHeavyBaryonDecayer::"
				 << "threeHalfHalfVectorCoupling()" 
				 << Exception::abortnow;
  }
}

string RadiativeHeavyBaryonDecayer::setUpDecayMode(string arg) {
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
  // type of mode
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int itype = stoi(stype);
  // M1 coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy m1 = stof(stype)/MeV;
  // E1 coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy2 e1 = stof(stype)/MeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  modeType_.push_back(itype);
  M1Coupling_.push_back(m1);
  E1Coupling_.push_back(e1);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
