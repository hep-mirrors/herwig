// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MECorrections class.
//

#include "MECorrections.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;


MECorrections::~MECorrections() {}


void MECorrections::persistentOutput(PersistentOStream & os) const {
  os << _vecMECorrection
     << _MECorrectionsMode
     << _composeMECorrectionsMode;
}


void MECorrections::persistentInput(PersistentIStream & is, int) {
  is >> _vecMECorrection
     >>	_MECorrectionsMode
     >>	_composeMECorrectionsMode;
}


ClassDescription<MECorrections> MECorrections::initMECorrections;
// Definition of the static class description member.


void MECorrections::Init() {

  static ClassDocumentation<MECorrections> documentation
    ("This class is responsible for managing the matrix element corrections.");

  static RefVector<MECorrections,MECorrection> interfaceVecMECorrection
    ("VecMECorrection",
     "The collection of ME corrections. ",
     &MECorrections::_vecMECorrection, 0, false, false, true, false);

  static Switch<MECorrections, int> interfaceMECorrectionsMode
    ("OnOffMECorrectionMode",
     "Choice of the on-off matrix element ccrrections switch mode",
     &MECorrections::_MECorrectionsMode, 0, false, false);
  static SwitchOption interfaceMECorrectionsMode0
    (interfaceMECorrectionsMode,"MECorrections-OFF","ME Corrections are OFF", 0);
  static SwitchOption interfaceMECorrectionsMode1
    (interfaceMECorrectionsMode,"MECorrections-ON","ME Corrections are ON", 1);

  static Switch<MECorrections, int> interfaceComposeMECorrectionsMode
    ("OnOffComposeMECorrectionMode",
     "Choice of the on-off composition of matrix element ccrrections switch mode",
     &MECorrections::_composeMECorrectionsMode, 0, false, false);
  static SwitchOption interfaceComposeMECorrectionsMode0
    (interfaceComposeMECorrectionsMode,
     "composeMECorrections-OFF","Composition of ME Corrections is OFF", 0);
  static SwitchOption interfaceComposeMECorrectionsMode1
    (interfaceComposeMECorrectionsMode,
     "composeMECorrections-ON","Composition of ME Corrections is ON", 1);

}

//----------------------------------------------------------------------------

void MECorrections::initializeRun() {

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower ) {
    for ( vector<MECorrectionPtr>::const_iterator cit = _vecMECorrection.begin();
	  cit != _vecMECorrection.end(); ++cit ) {
      if ( ! ( (   (*cit)->hardProcessME()   &&    (*cit)->hardProcessPlusJetME()  && 
		 ! (*cit)->decayProcessME()  &&  ! (*cit)->decayProcessPlusJetME() ) 
	       || 
	       ( ! (*cit)->hardProcessME()   &&  ! (*cit)->hardProcessPlusJetME()  && 
		   (*cit)->decayProcessME()  &&    (*cit)->decayProcessPlusJetME() ) ) ) { 
	generator()->logWarning( Exception("MECorrections::initializeRun "
					   "***Some inconsistency in a MECorrection*** ", 
					   Exception::warning) );
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "         ===>" << endl 
	                     << "\t hardProcessME()         : " 
			     << ( (*cit)->hardProcessME() ? "1" : "0" ) << endl
	                     << "\t hardProcessPlusJetME()  : " 
			     << ( (*cit)->hardProcessPlusJetME() ? "1" : "0" ) << endl
	                     << "\t decayProcessME()        : " 
			     << ( (*cit)->decayProcessME() ? "1" : "0" ) << endl
	                     << "\t decayProcessPlusJetME() : " 
			     << ( (*cit)->decayProcessPlusJetME() ? "1" : "0" ) << endl
			     << endl;
	}
      }
    }
  }

  int index = 0;  
  for ( vector<MECorrectionPtr>::const_iterator cit = _vecMECorrection.begin();
	cit != _vecMECorrection.end(); ++cit, ++index ) {
    if ( (*cit)->hardProcessME() ) {
      _mapHardProcesses.insert( pair<tMEPtr,int>( (*cit)->hardProcessME(), index ) );
      _mapHardPlusJetProcesses.insert( pair<tMEPtr,int>( (*cit)->hardProcessPlusJetME(), 
							 index ) );
    } else if ( (*cit)->decayProcessME() ) {
      _mapDecayProcesses.insert( pair<tDecayerPtr,int>( (*cit)->decayProcessME(), index ) );
      _mapDecayPlusJetProcesses.insert( pair<tDecayerPtr,int>( (*cit)->decayProcessPlusJetME(),
							       index ) );
    }
  }
}


tMECorrectionPtr MECorrections::getMECorrection( const tMEPtr hardProcess ) const {
  tMECorrectionPtr theMECorrection = tMECorrectionPtr();
  int index = -1;
  if ( _mapHardProcesses.find( hardProcess ) != _mapHardProcesses.end() ) {
    index = _mapHardProcesses.find( hardProcess )->second;
  } else if ( _mapHardPlusJetProcesses.find( hardProcess ) != 
	      _mapHardPlusJetProcesses.end() ) {
    index = _mapHardPlusJetProcesses.find( hardProcess )->second;    
  }
  if ( index >= 0  &&  (_vecMECorrection[index])->isMECorrectionON() ) {
    theMECorrection = _vecMECorrection[index];
  }
  return theMECorrection;
}


tMECorrectionPtr MECorrections::getMECorrection( const tDecayerPtr decayProcess ) const {
  tMECorrectionPtr theMECorrection = tMECorrectionPtr();
  int index = -1;
  if ( _mapDecayProcesses.find( decayProcess ) != _mapDecayProcesses.end() ) {
    index = _mapDecayProcesses.find( decayProcess )->second;
  } else if ( _mapDecayPlusJetProcesses.find( decayProcess ) != 
	      _mapDecayPlusJetProcesses.end() ) {
    index = _mapDecayPlusJetProcesses.find( decayProcess )->second;    
  }
  if ( index >= 0  &&  (_vecMECorrection[index])->isMECorrectionON() ) {
    theMECorrection = _vecMECorrection[index];
  }
  return theMECorrection;
}


