// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MECorrection class.
//

#include "MECorrection.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h" 
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;


MECorrection::~MECorrection() {}


void MECorrection::persistentOutput(PersistentOStream & os) const {
  os << _correctionMode
     << _hardProcessME
     << _hardProcessPlusJetME
     << _decayProcessME
     << _decayProcessPlusJetME;
}

void MECorrection::persistentInput(PersistentIStream & is, int) {
  is >> _correctionMode
     >> _hardProcessME
     >> _hardProcessPlusJetME
     >> _decayProcessME
     >> _decayProcessPlusJetME;
}


AbstractClassDescription<MECorrection> MECorrection::initMECorrection;
// Definition of the static class description member.


void MECorrection::Init() {

  static ClassDocumentation<MECorrection> documentation
    ("This is the abstract class from which any matrix element correction object ",
     "should inherit from.");

  static Switch<MECorrection, int> interfaceMECorrectionMode
    ("OnOffMECorrectionMode",
     "Choice of the on-off matrix element ccrrections switch mode",
     &MECorrection::_correctionMode, 1, false, false);
  static SwitchOption interfaceMECorrectionMode0
    (interfaceMECorrectionMode,"MECorrection-OFF",
     "ME correction is OFF", 0);
  static SwitchOption interfaceMECorrectionMode1
    (interfaceMECorrectionMode,"MECorrection-ON",
     "ME correction is ON", 1);

  static Reference<MECorrection,MEBase> 
    interfaceHardProcessME("HardProcessME", 
                           "A reference to the Hard Process Matrix Element object", 
                           &Herwig::MECorrection::_hardProcessME,
                           false, false, true, true);
  static Reference<MECorrection,MEBase> 
    interfaceHardProcessPlusJetME("HardProcessPlusJetME", 
				  "A reference to the Hard Process + Jet Matrix Element object", 
				  &Herwig::MECorrection::_hardProcessPlusJetME,
				  false, false, true, true);
  static Reference<MECorrection,Decayer> 
    interfaceDecayProcessME("DecayProcessME", 
			    "A reference to the Decay Process Matrix Element object", 
			    &Herwig::MECorrection::_decayProcessME,
			    false, false, true, true);
  static Reference<MECorrection,Decayer> 
    interfaceDecayProcessPlusJetME("DecayProcessPlusJetME", 
				   "A reference to the Decay Process + Jet Matrix Element object", 
				   &Herwig::MECorrection::_decayProcessPlusJetME,
				   false, false, true, true);

}


void MECorrection::hardMECorrection() throw(Veto, Stop, Exception) {}


