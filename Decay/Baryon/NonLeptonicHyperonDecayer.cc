// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonLeptonicHyperonDecayer class.
//

#include "NonLeptonicHyperonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void NonLeptonicHyperonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

void NonLeptonicHyperonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size() || isize!=a_.size()||
     isize!=b_.size()        || isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "NonLeptonicHyperonDecayer::doinit()" 
			  << Exception::runerror;
  // set up the decay modes
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first ),
    		     getParticleData(outgoing_[ix].second)};
    double wgtmax = maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    PhaseSpaceModePtr mode(new_ptr(PhaseSpaceMode(in,out,wgtmax)));
    addMode(mode);
    // test of the asummetries
    // Energy e1 = (sqr(in->mass())+sqr(out[0]->mass())-
    // 		 sqr(out[1]->mass()))/2./in->mass();
    // double btemp = b_[ix]*sqrt((e1-out[0]->mass())/(e1+out[0]->mass()));
    // double alpha = -2.*(a_[ix]*btemp)/(sqr(a_[ix])+sqr(btemp));
    // generator()->log() << "Asymmetry parameter for " << in->PDGName() << "->"
    // 		       << out[0]->PDGName() << "," << out[1]->PDGName()
    // 		       << " = " << alpha << "\n";
  }
}

NonLeptonicHyperonDecayer::NonLeptonicHyperonDecayer() {
  // intermediates
  generateIntermediates(false);
}

int NonLeptonicHyperonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing pa4rticles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0=parent->id();
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do {
    if(id0==incoming_[ix]) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) imode=ix;
    }
    else if(id0==-incoming_[ix]) {
      if((id1==-outgoing_[ix].first&&id2==-outgoing_[ix].second)||
	 (id2==-outgoing_[ix].first&&id1==-outgoing_[ix].second)) imode=ix;
      if(((id1==-outgoing_[ix].first&&id2==outgoing_[ix].second)||
	  (id2==-outgoing_[ix].first&&id1==outgoing_[ix].second))&&
	 (outgoing_[ix].second==111||outgoing_[ix].second==221||outgoing_[ix].second==331||
	  outgoing_[ix].second==223||outgoing_[ix].second==333)) imode=ix;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  // charge conjugate
  cc=id0<0;
  // return the answer
  return imode;
}

void NonLeptonicHyperonDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << a_ << b_ << maxweight_;
}

void NonLeptonicHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> a_ >> b_ >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NonLeptonicHyperonDecayer,Baryon1MesonDecayerBase>
describeHerwigNonLeptonicHyperonDecayer("Herwig::NonLeptonicHyperonDecayer", "HwBaryonDecay.so");

void NonLeptonicHyperonDecayer::Init() {

  static ClassDocumentation<NonLeptonicHyperonDecayer> documentation
    ("The NonLeptonicHyperonDecayer class performs the non-leptonic"
     " weak decay of the hyperons.",
     "The non-leptonic hyperon decays were simulated using the "
     "NonLeptonicHyperonDecayer class which implements the model of"
     "\\cite{Borasoy:1999md}",
     "\\bibitem{Borasoy:1999md}\n"
     "B.~Borasoy and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 59} (1999) 094025 [arXiv:hep-ph/9902351].\n"
     "%%CITATION = PHRVA,D59,094025;%%\n");

  static Command<NonLeptonicHyperonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryon & meson, A and B couplings and max weight for a decay",
     &NonLeptonicHyperonDecayer::setUpDecayMode, false);

  static Deleted<NonLeptonicHyperonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceIncomingBaryon
    ("IncomingBaryon","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceOutgoingBaryon
    ("OutgoingBaryon","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceOutgoingMeson
    ("OutgoingMeson","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceA
    ("CouplingA","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceB
    ("CouplingB","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");
}

// couplings for spin-1/2 to spin-1/2 spin-0
void NonLeptonicHyperonDecayer::halfHalfScalarCoupling(int imode,
						       Energy,Energy,Energy,
						       Complex& A,Complex& B) const {
  useMe();
  A=a_[imode];
  B=b_[imode];
}

void NonLeptonicHyperonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix].first << " " << outgoing_[ix].second
	   << " " << a_[ix] << " " << b_[ix] << " " << maxweight_[ix] << "\n";
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string NonLeptonicHyperonDecayer::setUpDecayMode(string arg) {
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
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  a_.push_back(g);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  g = stof(stype);
  b_.push_back(g);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
