// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonFormFactor class.
//

#include "BaryonFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

BaryonFormFactor::~BaryonFormFactor() {}

void BaryonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _incomingid << _outgoingid << _incomingJ << _outgoingJ << _spectator1 
     << _spectator2 << _inquark << _outquark;
}

void BaryonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _incomingid >> _outgoingid >> _incomingJ >> _outgoingJ >> _spectator1 
     >> _spectator2 >> _inquark >> _outquark;
}

AbstractClassDescription<BaryonFormFactor> BaryonFormFactor::initBaryonFormFactor;
// Definition of the static class description member.

void BaryonFormFactor::Init() {

  static ClassDocumentation<BaryonFormFactor> documentation
    ("The \\classname{BaryonFormFactor} class is the base class for"
     " the implementation of the form factors for weak decays of baryon.");

  static ParVector<BaryonFormFactor,int> interfaceIncoming
    ("Incoming",
     "The id of the incoming baryons.",
     &BaryonFormFactor::_incomingid,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceOutgoing
    ("Outgoing",
     "The id of the outgoing baryons.",
     &BaryonFormFactor::_outgoingid,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceInSpin
    ("InSpin",
     "The spin of the incoming baryons to decide which form factors to use.",
     &BaryonFormFactor::_incomingJ,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceOutSpin
    ("OutSpin",
     "The spin of the outgoing baryons to decide which form factors to use.",
     &BaryonFormFactor::_outgoingJ,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceSpectator1
    ("Spectator1",
     "The first specator quark.",
     &BaryonFormFactor::_spectator1,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceSpectator2
    ("Spectator2",
     "The second specator quark.",
     &BaryonFormFactor::_spectator2,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceInQuark
    ("InQuark",
     "The PDG code for the decaying quark.",
     &BaryonFormFactor::_inquark,
     0, 0, 0, 0, 2, false, false, true);

  static ParVector<BaryonFormFactor,int> interfaceOutQuark
    ("OutQuark",
     "The PDG code for the quark produced in the decay.",
     &BaryonFormFactor::_outquark,
     0, 0, 0, 0, 2, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void BaryonFormFactor::SpinHalfSpinHalfFormFactor(Energy2,int,int,int,Energy,Energy,
						  Complex &,Complex &,Complex &,
						  Complex &,Complex &,Complex &)
{
  throw Exception() << "Error in BaryonFormFactor::SpinHalfSpinHalfFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}
void BaryonFormFactor::SpinHalfSpinThreeHalfFormFactor(Energy2,int,int,int,Energy,
						       Energy,Complex &,Complex &,
						       Complex &,Complex &,Complex &,
						       Complex &,Complex &,Complex &)
{
  throw Exception() << "Error in BaryonFormFactor::SpinHalfSpinThreeHalfFormFactor"
		    << " not implemented"
		    << Exception::abortnow;
}

// output the information for the database
void BaryonFormFactor::dataBaseOutput(ofstream&)
{}
}
