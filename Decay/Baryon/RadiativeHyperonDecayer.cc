// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeHyperonDecayer class.
//

#include "RadiativeHyperonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
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
  tPDVector extpart(3);
  extpart[2]=getParticleData(ParticleID::gamma);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    extpart[0]=getParticleData(incomingB_[ix]);
    extpart[1]=getParticleData(outgoingB_[ix]);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax = maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
  }
}

RadiativeHyperonDecayer::RadiativeHyperonDecayer() {
  // Sigma+ -> p+
  incomingB_.push_back(3222); outgoingB_.push_back(2212);
  A_.push_back(-1.81e-7/GeV); B_.push_back(0.47e-7/GeV);
  maxweight_.push_back(1.);
  // Xi- -> Sigma- 
  incomingB_.push_back(3312); outgoingB_.push_back(3112);
  A_.push_back(0.08e-7/GeV); B_.push_back(0.15e-7/GeV);
  maxweight_.push_back(1.);
  // Sigma0 -> n
  incomingB_.push_back(3212); outgoingB_.push_back(2112);
  A_.push_back(-0.02e-7/GeV); B_.push_back(-0.45e-7/GeV);
  maxweight_.push_back(1.);
  // Lambda0 -> n
  incomingB_.push_back(3122); outgoingB_.push_back(2112);
  A_.push_back(-0.52e-7/GeV); B_.push_back(-0.05e-7/GeV);
  maxweight_.push_back(1.);
  // Xi0  -> Sigma0
  incomingB_.push_back(3322); outgoingB_.push_back(3212);
  A_.push_back( 0.05e-7/GeV); B_.push_back( 0.70e-7/GeV);
  maxweight_.push_back(1.);
  // Xi0  -> Lambda0 
  incomingB_.push_back(3322); outgoingB_.push_back(3122);
  A_.push_back(-0.34e-7/GeV); B_.push_back(-0.08e-7/GeV);
  maxweight_.push_back(1.);
  // initial size of vectors
  initsize_ = A_.size();
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
     << maxweight_ << initsize_;
}

void RadiativeHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incomingB_ >> outgoingB_ >> iunit(A_,1./GeV) >> iunit(B_,1./GeV) 
     >> maxweight_ >> initsize_;
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

  static ParVector<RadiativeHyperonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &RadiativeHyperonDecayer::maxweight_,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<RadiativeHyperonDecayer,long> interfaceIncomingBaryon
    ("IncomingBaryon",
     "The PDG code for the incoming baryon.",
     &RadiativeHyperonDecayer::incomingB_,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<RadiativeHyperonDecayer,long> interfaceOutgoingBaryon
    ("OutgoingBaryon",
     "The PDG code for the outgoing baryon.",
     &RadiativeHyperonDecayer::outgoingB_,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<RadiativeHyperonDecayer,InvEnergy> interfaceCouplingA
    ("CouplingA",
     "The A coupling",
     &RadiativeHyperonDecayer::A_, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, Interface::limited);

  static ParVector<RadiativeHyperonDecayer,InvEnergy> interfaceCouplingB
    ("CouplingB",
     "The B coupling",
     &RadiativeHyperonDecayer::B_, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, Interface::limited);

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
    if(ix<initsize_) {
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << maxweight_[ix] << "\n";
      output << "newdef " << name() << ":IncomingBaryon " << ix << " " 
	     << incomingB_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingBaryon " << ix << " " 
	     << outgoingB_[ix] << "\n";
      output << "newdef " << name() << ":CouplingA " << ix << " " 
	     << A_[ix]*GeV << "\n";
      output << "newdef " << name() << ":CouplingB " << ix << " " 
	     << B_[ix]*GeV << "\n";
    }
    else {
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << maxweight_[ix] << "\n";
      output << "insert " << name() << ":IncomingBaryon " << ix << " " 
	     << incomingB_[ix] << "\n";
      output << "insert " << name() << ":OutgoingBaryon " << ix << " " 
	     << outgoingB_[ix] << "\n";
      output << "insert " << name() << ":CouplingA " << ix << " " 
	     << A_[ix]*GeV << "\n";
      output << "insert " << name() << ":CouplingB " << ix << " " 
	     << B_[ix]*GeV << "\n";
    }
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

