// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarFormFactor class.
//

#include "ScalarFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ScalarFormFactor::~ScalarFormFactor() {}

void ScalarFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _incomingid << _outgoingid << _outgoingJ << _spectator << _inquark << _outquark
     << _numbermodes;
}

void ScalarFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _incomingid >> _outgoingid >> _outgoingJ >> _spectator >> _inquark >> _outquark
     >> _numbermodes;
}

AbstractClassDescription<ScalarFormFactor> ScalarFormFactor::initScalarFormFactor;
// Definition of the static class description member.

void ScalarFormFactor::Init() {

  static ClassDocumentation<ScalarFormFactor> documentation
    ("The \\classname{ScalarFormFactor} class is the base class for"
     " the implementation of the form factors for weak decays of scalar mesons.");

  static ParVector<ScalarFormFactor,int> interfaceIncoming
    ("Incoming",
     "The id of the incoming mesons.",
     &ScalarFormFactor::_incomingid,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<ScalarFormFactor,int> interfaceOutgoing
    ("Outgoing",
     "The id of the outgoing mesons.",
     &ScalarFormFactor::_outgoingid,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<ScalarFormFactor,int> interfaceSpin
    ("Spin",
     "The spin of the outgoing mesons to decide which form factors to use.",
     &ScalarFormFactor::_outgoingJ,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<ScalarFormFactor,int> interfaceSpectator
    ("Spectator",
     "The PDG code for the spectator quark.",
     &ScalarFormFactor::_spectator,
     0, 0, 0, -6, 6, false, false, true);

  static ParVector<ScalarFormFactor,int> interfaceInQuark
    ("InQuark",
     "The PDG code for the decaying quark.",
     &ScalarFormFactor::_inquark,
     0, 0, 0, -6, 6, false, false, true);

  static ParVector<ScalarFormFactor,int> interfaceOutQuark
    ("OutQuark",
     "The PDG code for the quark produced in the decay.",
     &ScalarFormFactor::_outquark,
     0, 0, 0, -6, 6, false, false, true);
}

// form-factor for scalar to scalar
void ScalarFormFactor::ScalarScalarFormFactor(Energy2,unsigned int,int,int,Energy,Energy,
					      Complex &,Complex &) const
{
  throw Exception() << "Error in ScalarFormFactor::ScalarScalarFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to vector
void ScalarFormFactor::ScalarVectorFormFactor(Energy2,unsigned int,int,int,Energy,Energy,
					      Complex &,Complex &,
					      Complex &,Complex &) const
{
  throw Exception() << "Error in ScalarFormFactor::ScalarVectorFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to tensor
void ScalarFormFactor::ScalarTensorFormFactor(Energy2,unsigned int,int,int,Energy,Energy,
					      Complex &,Complex &,
					      Complex &,Complex &) const
{
  throw Exception() << "Error in ScalarFormFactor::ScalarTensorFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to scalar (sigma)
void ScalarFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,
						   int id1,Energy m0, Energy m1,
						   Complex & fT) const
{
  throw Exception() << "Error in ScalarFormFactor::ScalarScalarSigmaFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}
// form-factor for scalar to vector (sigma)
void ScalarFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
						   Energy m0, Energy m1, Complex & T1,
						   Complex & T2, Complex & T3) const
{
  throw Exception() << "Error in ScalarFormFactor::ScalarVectorSigmaFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}
void  ScalarFormFactor::dataBaseOutput(ofstream & output)
{
  for(unsigned int ix=0;ix<_incomingid.size();++ix)
    {
      if(ix<_numbermodes)
	{
	  output << "set " << fullName() << ":Incoming "  << ix << " " 
		 << _incomingid[ix] << endl;
	  output << "set " << fullName() << ":Outgoing "  << ix << " " 
		 << _outgoingid[ix] << endl;
	  output << "set " << fullName() << ":Spin "      << ix << " " 
		 << _outgoingJ[ix] << endl;
	  output << "set " << fullName() << ":Spectator " << ix << " " 
		 << _spectator[ix] << endl;
	  output << "set " << fullName() << ":InQuark "   << ix << " " 
		 << _inquark[ix] << endl;
	  output << "set " << fullName() << ":OutQuark "  << ix << " " 
		 << _outquark[ix]<< endl;
	}
      else
	{
	  output << "insert " << fullName() << ":Incoming "  << ix << " " 
		 << _incomingid[ix] << endl;
	  output << "insert " << fullName() << ":Outgoing "  << ix << " " 
		 << _outgoingid[ix] << endl;
	  output << "insert " << fullName() << ":Spin "      << ix << " " 
		 << _outgoingJ[ix] << endl;
	  output << "insert " << fullName() << ":Spectator " << ix << " "
		 << _spectator[ix] << endl;
	  output << "insert " << fullName() << ":InQuark "   << ix << " " 
		 << _inquark[ix] << endl;
	  output << "insert " << fullName() << ":OutQuark "  << ix << " " 
		 << _outquark[ix]<< endl;
	}
    }
}
}
