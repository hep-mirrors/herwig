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

NonLeptonicHyperonDecayer::~NonLeptonicHyperonDecayer() {}

bool NonLeptonicHyperonDecayer::accept(const DecayMode & dm) const {
  // is this mode allowed
  bool allowed(false);
  // must be two outgoing particles
  if(dm.products().size()!=2){return allowed;}
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
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){allowed=true;}
	}
      else if(id0==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){allowed=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){allowed=true;}
	}
      ++ix;
    }
  while(ix<_incomingB.size()&&!allowed);
  return allowed;
}

ParticleVector NonLeptonicHyperonDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(-1),id(parent.id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);bool cc(false);
  do 
    {
      if(id==_incomingB[ix])
	{
	  if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	     (id2==_outgoingB[ix]&&id1==_outgoingM[ix])){imode=ix;cc=false;}
	}
      else if(id==-_incomingB[ix])
	{
	  if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	     (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])){imode=ix;cc=true;}
	  if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	      (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	     (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	      _outgoingM[ix]==223||_outgoingM[ix]==333)){imode=ix;cc=true;}
	}
      ++ix;
    }
  while(ix<_incomingB.size()&&imode<0);
  // generate the decay
  return generate(false,cc,imode,parent);
}


void NonLeptonicHyperonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incomingB << _outgoingB << _outgoingM << _A << _B << _maxweight << _initsize;
}

void NonLeptonicHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incomingB >> _outgoingB >> _outgoingM >> _A >> _B >> _maxweight << _initsize;
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
