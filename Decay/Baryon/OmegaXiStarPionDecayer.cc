// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaXiStarPionDecayer class.
//

#include "OmegaXiStarPionDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

OmegaXiStarPionDecayer::OmegaXiStarPionDecayer()  {
  // the ids of the particles
  idin_  = 3334;
  idout_ = 3324;
  // the couplings from the paper
  Acomm_ =  20.91e-8;
  AP_    =  -9.20e-8;
  AS_    =  -6.32e-8;
  BP_    =  230.1e-8;
  BS_    = -100.8e-8;
  // maximum weight for the decay
  wgtmax_=0.0032;
  // intermediates
  generateIntermediates(false);
}

void OmegaXiStarPionDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the phase space
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  vector<double> wgt(0);
  extpart[0]=getParticleData(idin_);
  extpart[1]=getParticleData(idout_);
  extpart[2]=getParticleData(-211);
  mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
  addMode(mode,wgtmax_,wgt);
}

void OmegaXiStarPionDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) wgtmax_=mode(0)->maxWeight();
}

int OmegaXiStarPionDecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  if(id0==idin_) {
    if((id1==idout_&&id2==-211)||
       (id2==idout_&&id1==-211)) imode=0;
  }
  else if(id0==-idin_) {
    if((id1==-idout_&&id2==211)||
       (id2==-idout_&&id1==211)) imode=0;
  }
  // charge conjugation
  cc=id0<0;
  // return the answer
  return imode;
}

void OmegaXiStarPionDecayer::persistentOutput(PersistentOStream & os) const {
  os << Acomm_ << AP_ << AS_ << BP_ << BS_ << idin_ << idout_ <<  wgtmax_;
}

void OmegaXiStarPionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> Acomm_ >> AP_ >> AS_ >> BP_ >> BS_ >> idin_ >> idout_ >>  wgtmax_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OmegaXiStarPionDecayer,Baryon1MesonDecayerBase>
describeHerwigOmegaXiStarPionDecayer("Herwig::OmegaXiStarPionDecayer", "HwBaryonDecay.so");

void OmegaXiStarPionDecayer::Init() {
  
  static ClassDocumentation<OmegaXiStarPionDecayer> documentation
    ("The OmegaXiStarPionDecayer class performs the weak decay"
     " of the Omega to Xi*0 and pi-",
     "The decay of the $\\Omega^-$ to $\\Xi^{*0}\\pi^-$ was simulated"
     " using the model of \\cite{Duplancic:2004dy}.",
     "\\bibitem{Duplancic:2004dy}\n"
     "G.~Duplancic, H.~Pasagic and J.~Trampetic,\n"
     "Phys.\\ Rev.\\  D {\\bf 70} (2004) 077506 [arXiv:hep-ph/0405162].\n"
     "%%CITATION = PHRVA,D70,077506;%%\n");

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAcomm
    ("Acomm",
     "The Acomm coupling for the decay",
     &OmegaXiStarPionDecayer::Acomm_, 20.91e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAP
    ("AP",
     "The A_P coupling for the decay",
     &OmegaXiStarPionDecayer::AP_, -9.20e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAS
    ("AS",
     "The A_S coupling for the decay",
     &OmegaXiStarPionDecayer::AS_, -6.32e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceBP
    ("BP",
     "The B_P coupling for the decay",
     &OmegaXiStarPionDecayer::BP_, 230.1e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceBS
    ("BS",
     "The B_S coupling for the decay",
     &OmegaXiStarPionDecayer::BS_, -100.8e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the decay",
     &OmegaXiStarPionDecayer::wgtmax_, 0.0032, 0., 100.,
     false, false, false);

  static Parameter<OmegaXiStarPionDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDF code for the incoming baryon",
     &OmegaXiStarPionDecayer::idin_, 3334, 0, 1000000,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,int> interfaceOutgoing
    ("Outgoing",
     "The PDF code for the outgoing baryon",
     &OmegaXiStarPionDecayer::idout_, 3324, 0, 1000000,
     false, false, true);
}

// couplings for spin-3/2 to spin-3/2 spin-0
void OmegaXiStarPionDecayer::
threeHalfThreeHalfScalarCoupling(int,Energy,Energy,Energy,
				 Complex&A1,Complex&A2,Complex&B1,Complex&B2) const {
  useMe();
  A2=0.;
  B2=0.;
  A1=Acomm_+AP_+AS_;
  B1=BP_+BS_;
}

void OmegaXiStarPionDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Acomm " << Acomm_ << "\n";
  output << "newdef " << name() << ":AP " << AP_ << "\n";
  output << "newdef " << name() << ":AS " << AS_ << "\n";
  output << "newdef " << name() << ":BP " << BP_ << "\n";
  output << "newdef " << name() << ":BS " << BS_ << "\n";
  output << "newdef " << name() << ":MaximumWeight " << wgtmax_ << "\n";
  output << "newdef " << name() << ":Incoming " << idin_ << "\n";
  output << "newdef " << name() << ":Outgoing " << idout_ << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
