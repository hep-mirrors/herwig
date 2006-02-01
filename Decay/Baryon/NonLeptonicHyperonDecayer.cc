// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonLeptonicHyperonDecayer class.
//

#include "NonLeptonicHyperonDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "NonLeptonicHyperonDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

NonLeptonicHyperonDecayer::NonLeptonicHyperonDecayer() 
{
  // lambda -> p pi-
  _incomingB.push_back(3122);_outgoingB.push_back(2212);_outgoingM.push_back(-211);
  _A.push_back(3.25e-7);_B.push_back(23.4e-7);
  _maxweight.push_back(1.4);
  // lambda -> n p0
  _incomingB.push_back(3122);_outgoingB.push_back(2112);_outgoingM.push_back(111);
  _A.push_back(-2.30e-7);_B.push_back(-16.5e-7);
  _maxweight.push_back(0.7);
  // xi-    -> lambda pi-
  _incomingB.push_back(3312);_outgoingB.push_back(3122);_outgoingM.push_back(-211);
  _A.push_back(-4.51e-7);_B.push_back(14.8e-7);
  _maxweight.push_back(2.0);
  // xi0    -> lambda pi0
  _incomingB.push_back(3322);_outgoingB.push_back(3122);_outgoingM.push_back(111);
  _A.push_back(3.18e-7);_B.push_back(-10.5e-7);
  _maxweight.push_back(2.0);
  // sigma+ -> p pi0
  _incomingB.push_back(3222);_outgoingB.push_back(2212);_outgoingM.push_back(111);
  _A.push_back(-2.93e-7);_B.push_back(32.5e-7);
  _maxweight.push_back(1.2);
  // sigma- -> n pi-
  _incomingB.push_back(3112);_outgoingB.push_back(2112);_outgoingM.push_back(-211);
  _A.push_back(4.27e-7);_B.push_back(-1.52e-7);
  _maxweight.push_back(2.);
  // sigma+ -> n pi+
  _incomingB.push_back(3222);_outgoingB.push_back(2112);_outgoingM.push_back(211);
  _A.push_back(0.13e-7);_B.push_back(44.4e-7);
  _maxweight.push_back(1.1);
  // initial size of the vectors
  _initsize=_incomingB.size();
  // intermediates
  generateIntermediates(false);
}

NonLeptonicHyperonDecayer::~NonLeptonicHyperonDecayer() {}

int NonLeptonicHyperonDecayer::modeNumber(bool & cc,const DecayMode & dm) const
{
  int imode(-1);
  // must be two outgoing particles
  if(dm.products().size()!=2){return imode;}
  // ids of the particles
  int id0=dm.parent()->id();
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1=(**pit).id();++pit;
  int id2=(**pit).id();
  unsigned int ix(0);
  do
    {
      if(id0==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;}
	}
      else if(id0==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;}
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
  os << _incomingB << _outgoingB << _outgoingM << _A << _B << _maxweight << _initsize;
}

void NonLeptonicHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incomingB >> _outgoingB >> _outgoingM >> _A >> _B >> _maxweight >> _initsize;
}

ClassDescription<NonLeptonicHyperonDecayer> NonLeptonicHyperonDecayer::initNonLeptonicHyperonDecayer;
// Definition of the static class description member.

void NonLeptonicHyperonDecayer::Init() {

  static ClassDocumentation<NonLeptonicHyperonDecayer> documentation
    ("The NonLeptonicHyperonDecayer class performs the non-leptonic"
     " weak decay of the hyperons.");

  static ParVector<NonLeptonicHyperonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &NonLeptonicHyperonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,int> interfaceIncomingBaryon
    ("IncomingBaryon",
     "The PDG code for the incoming baryon.",
     &NonLeptonicHyperonDecayer::_incomingB,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,int> interfaceOutgoingBaryon
    ("OutgoingBaryon",
     "The PDG code for the outgoing baryon.",
     &NonLeptonicHyperonDecayer::_outgoingB,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,int> interfaceOutgoingMeson
    ("OutgoingMeson",
     "The PDG code for the outgoing meson.",
     &NonLeptonicHyperonDecayer::_outgoingM,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,double> interfaceA
    ("CouplingA",
     "The A coupling for the decay",
     &NonLeptonicHyperonDecayer::_A,
     0, 0, 0, -1e-5, 1e-5, false, false, true);

  static ParVector<NonLeptonicHyperonDecayer,double> interfaceB
    ("CouplingB",
     "The B coupling for the decay",
     &NonLeptonicHyperonDecayer::_B,
     0, 0, 0, -1e-5, 1e-5, false, false, true);
}

// couplings for spin-1/2 to spin-1/2 spin-0
void NonLeptonicHyperonDecayer::halfHalfScalarCoupling(int imode,
						       Energy m0,Energy m1,Energy m2,
						       Complex& A,Complex& B) const
{A=_A[imode];B=_B[imode];}

  void NonLeptonicHyperonDecayer::dataBaseOutput(ofstream & output,bool header) const
{
  if(header){output << "update decayers set parameters=\"";}
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incomingB.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	  output << "set " << fullName() << ":IncomingBaryon " << ix << " " 
		 << _incomingB[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingBaryon " << ix << " " 
		 << _outgoingB[ix] << "\n";
	  output << "set " << fullName() << ":OutgoingMeson " << ix << " " 
		 << _outgoingM[ix] << "\n";
	  output << "set " << fullName() << ":CouplingA " << ix << " " 
		 << _A[ix] << "\n";
	  output << "set " << fullName() << ":CouplingB " << ix << " " 
		 << _B[ix] << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":MaxWeight " << ix << " " 
		 << _maxweight[ix] << "\n";
	  output << "insert " << fullName() << ":IncomingBaryon " << ix << " " 
		 << _incomingB[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingBaryon " << ix << " " 
		 << _outgoingB[ix] << "\n";
	  output << "insert " << fullName() << ":OutgoingMeson " << ix << " " 
		 << _outgoingM[ix] << "\n";
	  output << "insert " << fullName() << ":CouplingA " << ix << " " 
		 << _A[ix] << "\n";
	  output << "insert " << fullName() << ":CouplingB " << ix << " " 
		 << _B[ix] << "\n";
	}
    }
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
}
