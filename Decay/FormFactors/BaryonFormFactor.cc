// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonFormFactor class.
//

#include "BaryonFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

void BaryonFormFactor::doinit() {
  Interfaced::doinit();
  // check the consistency of the parameters
  unsigned int isize(_incomingid.size());
  if(isize!=_outgoingid.size() || isize!=_incomingJ.size()||isize!=_outgoingJ.size()||
     isize!=_spectator1.size()|| isize!=_spectator2.size()|| isize!=_inquark.size()||
     isize!=_outquark.size())
    {throw InitException() << "Inconsistent parameters in BaryonFormFactor::doinit() " 
			   << Exception::abortnow;}
}

void BaryonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _incomingid << _outgoingid << _incomingJ << _outgoingJ << _spectator1 
     << _spectator2 << _inquark << _outquark << _numbermodes;
}

void BaryonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _incomingid >> _outgoingid >> _incomingJ >> _outgoingJ >> _spectator1 
     >> _spectator2 >> _inquark >> _outquark >> _numbermodes;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<BaryonFormFactor,Interfaced>
describeHerwigBaryonFormFactor("Herwig::BaryonFormFactor", "Herwig.so");

void BaryonFormFactor::Init() {

  static ClassDocumentation<BaryonFormFactor> documentation
    ("The BaryonFormFactor class is the base class for"
     " the implementation of the form factors for weak decays of baryon.");

  static ParVector<BaryonFormFactor,int> interfaceIncoming
    ("Incoming",
     "The id of the incoming baryons.",
     &BaryonFormFactor::_incomingid,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceOutgoing
    ("Outgoing",
     "The id of the outgoing baryons.",
     &BaryonFormFactor::_outgoingid,
     0, 0, 0, 0, 1000000, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceInSpin
    ("InSpin",
     "The spin of the incoming baryons to decide which form factors to use.",
     &BaryonFormFactor::_incomingJ,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceOutSpin
    ("OutSpin",
     "The spin of the outgoing baryons to decide which form factors to use.",
     &BaryonFormFactor::_outgoingJ,
     0, 0, 0, 0, 4, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceSpectator1
    ("Spectator1",
     "The first specator quark.",
     &BaryonFormFactor::_spectator1,
     0, 0, 0, 0, 5, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceSpectator2
    ("Spectator2",
     "The second specator quark.",
     &BaryonFormFactor::_spectator2,
     0, 0, 0, 0, 5, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceInQuark
    ("InQuark",
     "The PDG code for the decaying quark.",
     &BaryonFormFactor::_inquark,
     0, 0, 0, 0, 5, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceOutQuark
    ("OutQuark",
     "The PDG code for the quark produced in the decay.",
     &BaryonFormFactor::_outquark,
     0, 0, 0, 0, 5, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void BaryonFormFactor::SpinHalfSpinHalfFormFactor(Energy2,int,int,int,Energy,Energy,
						  Complex &,Complex &,Complex &,
						  Complex &,Complex &,Complex &) {
  throw Exception() << "Error in BaryonFormFactor::SpinHalfSpinHalfFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}
void BaryonFormFactor::SpinHalfSpinThreeHalfFormFactor(Energy2,int,int,int,Energy,
						       Energy,Complex &,Complex &,
						       Complex &,Complex &,Complex &,
						       Complex &,Complex &,Complex &) {
  throw Exception() << "Error in BaryonFormFactor::SpinHalfSpinThreeHalfFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// output the information for the database
void BaryonFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BaryonFormFactor " << name() << " \n";
  for(unsigned int ix=0;ix<_incomingid.size();++ix) {
    if(ix<_numbermodes) {
      output << "newdef " << name() << ":Incoming "  << ix << " " 
	     << _incomingid[ix] << endl;
      output << "newdef " << name() << ":Outgoing "  << ix << " " 
	     << _outgoingid[ix] << endl;
      output << "newdef " << name() << ":InSpin "      << ix << " " 
	     << _incomingJ[ix] << endl;
      output << "newdef " << name() << ":OutSpin "      << ix << " " 
	     << _outgoingJ[ix] << endl;
      output << "newdef " << name() << ":Spectator1 " << ix << " " 
	     << _spectator1[ix] << endl;
      output << "newdef " << name() << ":Spectator2 " << ix << " " 
	     << _spectator2[ix] << endl;
      output << "newdef " << name() << ":InQuark "   << ix << " " 
	     << _inquark[ix] << endl;
      output << "newdef " << name() << ":OutQuark "  << ix << " " 
	     << _outquark[ix]<< endl;
    }
    else {
      output << "insert " << name() << ":Incoming "  << ix << " " 
	     << _incomingid[ix] << endl;
      output << "insert " << name() << ":Outgoing "  << ix << " " 
	     << _outgoingid[ix] << endl;
      output << "insert " << name() << ":InSpin "      << ix << " " 
	     << _incomingJ[ix] << endl;
      output << "insert " << name() << ":OutSpin "      << ix << " " 
	     << _outgoingJ[ix] << endl;
      output << "insert " << name() << ":Spectator1 " << ix << " " 
	     << _spectator1[ix] << endl;
      output << "insert " << name() << ":Spectator2 " << ix << " " 
	     << _spectator2[ix] << endl;
      output << "insert " << name() << ":InQuark "   << ix << " " 
	     << _inquark[ix] << endl;
      output << "insert " << name() << ":OutQuark "  << ix << " " 
	     << _outquark[ix]<< endl;
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

}
