// -*- C++ -*-
//
// ScalarFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarFormFactor class.
//

#include "ScalarFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void ScalarFormFactor::doinit() {
  Interfaced::doinit();
  // check the consistency of the parameters
  unsigned int isize=_incomingid.size();
  if(isize!=_outgoingid.size() || isize!=_outgoingJ.size()||
     isize!=_spectator.size()  || isize!=_inquark.size()||
     isize!=_outquark.size())
    throw InitException() << "Inconsistent parameters in ScalarFormFactor::doinit() " 
			  << Exception::abortnow;
}

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
    ("The ScalarFormFactor class is the base class for"
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
					      Complex &,Complex &) const {
  throw Exception() << "Error in ScalarFormFactor::ScalarScalarFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to vector
void ScalarFormFactor::ScalarVectorFormFactor(Energy2,unsigned int,int,int,Energy,Energy,
					      Complex &,Complex &,
					      Complex &,Complex &) const {
  throw Exception() << "Error in ScalarFormFactor::ScalarVectorFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to tensor
void ScalarFormFactor::ScalarTensorFormFactor(Energy2,unsigned int,int,int,Energy,Energy,
					      complex<InvEnergy2> & ,
					      Complex &, complex<InvEnergy2> &,
					      complex<InvEnergy2> &) const {
  throw Exception() << "Error in ScalarFormFactor::ScalarTensorFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to scalar (sigma)
void ScalarFormFactor::ScalarScalarSigmaFormFactor(Energy2,unsigned int,int,
						   int,Energy, Energy,
						   Complex &) const {
  throw Exception() << "Error in ScalarFormFactor::ScalarScalarSigmaFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// form-factor for scalar to vector (sigma)
void ScalarFormFactor::ScalarVectorSigmaFormFactor(Energy2,unsigned int,int,int,
						   Energy, Energy, Complex &,
						   Complex &, Complex &) const {
  throw Exception() << "Error in ScalarFormFactor::ScalarVectorSigmaFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

void  ScalarFormFactor::dataBaseOutput(ofstream & output,bool header,
				       bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ScalarFormFactor " << name() << " \n";
  for(unsigned int ix=0;ix<_incomingid.size();++ix) {
    if(ix<_numbermodes) {
      output << "newdef " << name() << ":Incoming "  << ix << " " 
	     << _incomingid[ix] << "\n";
      output << "newdef " << name() << ":Outgoing "  << ix << " " 
	     << _outgoingid[ix] << "\n";
      output << "newdef " << name() << ":Spin "      << ix << " " 
	     << _outgoingJ[ix] << "\n";
      output << "newdef " << name() << ":Spectator " << ix << " " 
	     << _spectator[ix] << "\n";
      output << "newdef " << name() << ":InQuark "   << ix << " " 
	     << _inquark[ix] << "\n";
      output << "newdef " << name() << ":OutQuark "  << ix << " " 
	     << _outquark[ix]<< "\n";
    }
    else {
      output << "insert " << name() << ":Incoming "  << ix << " " 
	     << _incomingid[ix] << "\n";
      output << "insert " << name() << ":Outgoing "  << ix << " " 
	     << _outgoingid[ix] << "\n";
      output << "insert " << name() << ":Spin "      << ix << " " 
	     << _outgoingJ[ix] << "\n";
      output << "insert " << name() << ":Spectator " << ix << " "
	     << _spectator[ix] << "\n";
      output << "insert " << name() << ":InQuark "   << ix << " " 
	     << _inquark[ix] << "\n";
      output << "insert " << name() << ":OutQuark "  << ix << " " 
	     << _outquark[ix]<< "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
