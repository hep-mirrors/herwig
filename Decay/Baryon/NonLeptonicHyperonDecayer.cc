// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonLeptonicHyperonDecayer class.
//

#include "NonLeptonicHyperonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void NonLeptonicHyperonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

void NonLeptonicHyperonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(_incomingB.size());
  if(isize!=_outgoingB.size()||isize!=_outgoingM.size()||isize!=_a.size()||
     isize!=_b.size()        ||isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in "
			  << "NonLeptonicHyperonDecayer::doinit()" 
			  << Exception::runerror;
  // set up the decay modes
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double wgtmax;
  vector<double> wgt(0);
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    extpart[0]=getParticleData(_incomingB[ix]);
    extpart[1]=getParticleData(_outgoingB[ix]);
    extpart[2]=getParticleData(_outgoingM[ix]);
    mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    wgtmax = _maxweight.size()>numberModes() ? 
      _maxweight[numberModes()] : 1.;
    addMode(mode,wgtmax,wgt);
    // test of the asummetries
//     Energy e1 = (sqr(extpart[0]->mass())+sqr(extpart[1]->mass())-
// 		 sqr(extpart[2]->mass()))/2./extpart[0]->mass();
//     double btemp = _b[ix]*sqrt((e1-extpart[1]->mass())/(e1+extpart[1]->mass()));
//     double alpha = -2.*(_a[ix]*btemp)/(sqr(_a[ix])+sqr(btemp));
//     generator()->log() << "Asymmetry parameter for " << extpart[0]->PDGName() << "->"
// 		       << extpart[1]->PDGName() << "," << extpart[2]->PDGName()
// 		       << " = " << alpha << "\n";
  }
}

NonLeptonicHyperonDecayer::NonLeptonicHyperonDecayer() {
  // lambda -> p pi-
  _incomingB.push_back(3122);_outgoingB.push_back(2212);_outgoingM.push_back(-211);
  _a.push_back(3.25e-7);_b.push_back(-23.4e-7);
  _maxweight.push_back(1.4);
  // lambda -> n p0
  _incomingB.push_back(3122);_outgoingB.push_back(2112);_outgoingM.push_back(111);
  _a.push_back(-2.30e-7);_b.push_back(16.5e-7);
  _maxweight.push_back(0.7);
  // xi-    -> lambda pi-
  _incomingB.push_back(3312);_outgoingB.push_back(3122);_outgoingM.push_back(-211);
  _a.push_back(-4.51e-7);_b.push_back(-14.8e-7);
  _maxweight.push_back(2.0);
  // xi0    -> lambda pi0
  _incomingB.push_back(3322);_outgoingB.push_back(3122);_outgoingM.push_back(111);
  _a.push_back(3.19e-7);_b.push_back(10.5e-7);
  _maxweight.push_back(2.0);
  // sigma+ -> p pi0
  _incomingB.push_back(3222);_outgoingB.push_back(2212);_outgoingM.push_back(111);
  _a.push_back(-2.93e-7);_b.push_back(-32.5e-7);
  _maxweight.push_back(1.2);
  // sigma- -> n pi-
  _incomingB.push_back(3112);_outgoingB.push_back(2112);_outgoingM.push_back(-211);
  _a.push_back(4.27e-7);_b.push_back(1.52e-7);
  _maxweight.push_back(2.);
  // sigma+ -> n pi+
  _incomingB.push_back(3222);_outgoingB.push_back(2112);_outgoingM.push_back(211);
  _a.push_back(0.13e-7);_b.push_back(-44.4e-7);
  _maxweight.push_back(1.1);
  // initial size of the vectors
  _initsize=_incomingB.size();
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
    if(id0==_incomingB[ix]) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) imode=ix;
    }
    else if(id0==-_incomingB[ix]) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) imode=ix;
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) imode=ix;
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  // charge conjugate
  cc=id0<0;
  // return the answer
  return imode;
}

void NonLeptonicHyperonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incomingB << _outgoingB << _outgoingM << _a << _b << _maxweight << _initsize;
}

void NonLeptonicHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incomingB >> _outgoingB >> _outgoingM >> _a >> _b >> _maxweight >> _initsize;
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

  static ParVector<NonLeptonicHyperonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &NonLeptonicHyperonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,long> interfaceIncomingBaryon
    ("IncomingBaryon",
     "The PDG code for the incoming baryon.",
     &NonLeptonicHyperonDecayer::_incomingB,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,long> interfaceOutgoingBaryon
    ("OutgoingBaryon",
     "The PDG code for the outgoing baryon.",
     &NonLeptonicHyperonDecayer::_outgoingB,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,long> interfaceOutgoingMeson
    ("OutgoingMeson",
     "The PDG code for the outgoing meson.",
     &NonLeptonicHyperonDecayer::_outgoingM,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,double> interfaceA
    ("CouplingA",
     "The A coupling for the decay",
     &NonLeptonicHyperonDecayer::_a,
     0, 0, 0, -1e-5, 1e-5, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,double> interfaceB
    ("CouplingB",
     "The B coupling for the decay",
     &NonLeptonicHyperonDecayer::_b,
     0, 0, 0, -1e-5, 1e-5, false, false, true);
}

// couplings for spin-1/2 to spin-1/2 spin-0
void NonLeptonicHyperonDecayer::halfHalfScalarCoupling(int imode,
						       Energy,Energy,Energy,
						       Complex& A,Complex& B) const {
  useMe();
  A=_a[imode];
  B=_b[imode];
}

void NonLeptonicHyperonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
      output << "newdef " << name() << ":IncomingBaryon " << ix << " " 
	     << _incomingB[ix] << "\n";
      output << "newdef " << name() << ":OutgoingBaryon " << ix << " " 
	     << _outgoingB[ix] << "\n";
      output << "newdef " << name() << ":OutgoingMeson " << ix << " " 
	     << _outgoingM[ix] << "\n";
      output << "newdef " << name() << ":CouplingA " << ix << " " 
	     << _a[ix] << "\n";
      output << "newdef " << name() << ":CouplingB " << ix << " " 
	     << _b[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
      output << "insert " << name() << ":IncomingBaryon " << ix << " " 
	     << _incomingB[ix] << "\n";
      output << "insert " << name() << ":OutgoingBaryon " << ix << " " 
	     << _outgoingB[ix] << "\n";
      output << "insert " << name() << ":OutgoingMeson " << ix << " " 
	     << _outgoingM[ix] << "\n";
      output << "insert " << name() << ":CouplingA " << ix << " " 
	     << _a[ix] << "\n";
      output << "insert " << name() << ":CouplingB " << ix << " " 
	     << _b[ix] << "\n";
    }
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
