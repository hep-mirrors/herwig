// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingGenerator class.
//

#include "SplittingGenerator.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
#include "Pythia7/Interface/Switch.h"
#include "Pythia7/Interface/Reference.h"
#include "Pythia7/PDT/EnumParticles.h"
#include "ShowerParticle.h"
#include "Pythia7/Handlers/PartialCollisionHandler.h"
#include "Pythia7/Repository/FullEventGenerator.h"
#include "ShowerAlpha.h"
#include "SplitFun.h"
#include "SplitFun1to2.h"
#include "IS_QtoQGSplitFun.h"
#include "FS_QtoQGSplitFun.h"
#include "IS_QtoQGammaSplitFun.h"
#include "FS_QtoQGammaSplitFun.h"
#include "QtoQGSudakovFormFactor.h"
#include "QtoQGammaSudakovFormFactor.h"
#include "IS_GtoGGSplitFun.h"
#include "FS_GtoGGSplitFun.h"
#include "GtoGGSudakovFormFactor.h"
#include "IS_GtoQQbarSplitFun.h"
#include "FS_GtoQQbarSplitFun.h"
#include "GtoQQbarSudakovFormFactor.h"
#include "IS_QtildaShowerKinematics1to2.h"
#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ShowerConfig.h"

using namespace Herwig;


SplittingGenerator::~SplittingGenerator() {}


void SplittingGenerator::persistentOutput(PersistentOStream & os) const {
  os << _QCDinteractionMode
     << _QEDinteractionMode
     << _EWKinteractionMode
     << _ISR_Mode
     << _ISR_QCDMode
     << _ISR_QEDMode
     << _ISR_EWKMode
     << _FSR_Mode
     << _FSR_QCDMode
     << _FSR_QEDMode
     << _FSR_EWKMode
     << _UtoUGsplittingMode
     << _DtoDGsplittingMode
     << _StoSGsplittingMode
     << _CtoCGsplittingMode
     << _BtoBGsplittingMode
     << _TtoTGsplittingMode
     << _GtoGGsplittingMode
     << _GtoUUbarsplittingMode
     << _GtoDDbarsplittingMode
     << _GtoSSbarsplittingMode
     << _GtoCCbarsplittingMode
     << _GtoBBbarsplittingMode
     << _GtoTTbarsplittingMode
     << _UtoUGammasplittingMode
     << _DtoDGammasplittingMode
     << _StoSGammasplittingMode
     << _CtoCGammasplittingMode
     << _BtoBGammasplittingMode
     << _TtoTGammasplittingMode
     << _pointerIS_ShowerAlphaQCD
     << _pointerFS_ShowerAlphaQCD  
     << _pointerIS_ShowerAlphaQED
     << _pointerFS_ShowerAlphaQED  
     << _pointerShowerConstrainer;  
}


void SplittingGenerator::persistentInput(PersistentIStream & is, int) {
  is >>	_QCDinteractionMode
     >>	_QEDinteractionMode
     >>	_EWKinteractionMode
     >>	_ISR_Mode
     >>	_ISR_QCDMode
     >>	_ISR_QEDMode
     >>	_ISR_EWKMode
     >>	_FSR_Mode
     >>	_FSR_QCDMode
     >>	_FSR_QEDMode
     >>	_FSR_EWKMode
     >>	_UtoUGsplittingMode
     >>	_DtoDGsplittingMode
     >>	_StoSGsplittingMode
     >>	_CtoCGsplittingMode
     >>	_BtoBGsplittingMode
     >>	_TtoTGsplittingMode
     >>	_GtoGGsplittingMode
     >>	_GtoUUbarsplittingMode
     >>	_GtoDDbarsplittingMode
     >>	_GtoSSbarsplittingMode
     >>	_GtoCCbarsplittingMode
     >>	_GtoBBbarsplittingMode
     >>	_GtoTTbarsplittingMode
     >> _UtoUGammasplittingMode
     >> _DtoDGammasplittingMode
     >> _StoSGammasplittingMode
     >> _CtoCGammasplittingMode
     >> _BtoBGammasplittingMode
     >> _TtoTGammasplittingMode
     >> _pointerIS_ShowerAlphaQCD
     >> _pointerFS_ShowerAlphaQCD  
     >> _pointerIS_ShowerAlphaQED
     >> _pointerFS_ShowerAlphaQED  
     >> _pointerShowerConstrainer;  
}


ClassDescription<SplittingGenerator> SplittingGenerator::initSplittingGenerator;
// Definition of the static class description member.

void SplittingGenerator::Init() {

  static ClassDocumentation<SplittingGenerator> documentation
    ("There class is responsible for initializing the Sudakov form factors ",
     "and generating splittings.");

  static Switch<SplittingGenerator, int> interfaceQCDinteractionMode
    ("OnOffQCDinteractionMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_QCDinteractionMode, 1, false, false);
  static SwitchOption interfaceQCDinteractionMode0
    (interfaceQCDinteractionMode,"QCDinteraction-OFF","QCD interaction is OFF", 0);
  static SwitchOption interfaceQCDinteractionMode1
    (interfaceQCDinteractionMode,"QCDinteraction-ON","QCD interaction is ON", 1);

  static Switch<SplittingGenerator, int> interfaceQEDinteractionMode
    ("OnOffQEDinteractionMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_QEDinteractionMode, 0, false, false);
  static SwitchOption interfaceQEDinteractionMode0
    (interfaceQEDinteractionMode,"QEDinteraction-OFF","QED interaction is OFF", 0);
  static SwitchOption interfaceQEDinteractionMode1
    (interfaceQEDinteractionMode,"QEDinteraction-ON","QED interaction is ON", 1);

  static Switch<SplittingGenerator, int> interfaceEWKinteractionMode
    ("OnOffEWKinteractionMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_EWKinteractionMode, 0, false, false);
  static SwitchOption interfaceEWKinteractionMode0
    (interfaceEWKinteractionMode,"EWKinteraction-OFF","EWK interaction is OFF", 0);
  static SwitchOption interfaceEWKinteractionMode1
    (interfaceEWKinteractionMode,"EWKinteraction-ON","EWK interaction is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISRMode
    ("OnOffISRMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_ISR_Mode, 0, false, false);
  static SwitchOption interfaceISRMode0
    (interfaceISRMode,"ISR-OFF","ISR (Initial State Radiation) is OFF", 0);
  static SwitchOption interfaceISRMode1
    (interfaceISRMode,"ISR-ON","ISR (Initial State Radiation) is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISR_QCDMode
    ("OnOffISR_QCDMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_ISR_QCDMode, 1, false, false);
  static SwitchOption interfaceISR_QCDMode0
    (interfaceISR_QCDMode,"ISR_QCD-OFF","QCD ISR is OFF", 0);
  static SwitchOption interfaceISR_QCDMode1
    (interfaceISR_QCDMode,"ISR_QCD-ON","QCD ISR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISR_QEDMode
    ("OnOffISR_QEDMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_ISR_QEDMode, 1, false, false);
  static SwitchOption interfaceISR_QEDMode0
    (interfaceISR_QEDMode,"ISR_QED-OFF","QED ISR is OFF", 0);
  static SwitchOption interfaceISR_QEDMode1
    (interfaceISR_QEDMode,"ISR_QED-ON","QED ISR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceISR_EWKMode
    ("OnOffISR_EWKMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_ISR_EWKMode, 1, false, false);
  static SwitchOption interfaceISR_EWKMode0
    (interfaceISR_EWKMode,"ISR_EWK-OFF","EWK ISR is OFF", 0);
  static SwitchOption interfaceISR_EWKMode1
    (interfaceISR_EWKMode,"ISR_EWK-ON","EWK ISR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSRMode
    ("OnOffFSRMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_FSR_Mode, 1, false, false);
  static SwitchOption interfaceFSRMode0
    (interfaceFSRMode,"FSR-OFF","FSR (Final State Radiation) is OFF", 0);
  static SwitchOption interfaceFSRMode1
    (interfaceFSRMode,"FSR-ON","FSR (Final State Radiation) is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSR_QCDMode
    ("OnOffFSR_QCDMode",
     "Choice of the on-off QCD interaction switch mode",
     &SplittingGenerator::_FSR_QCDMode, 1, false, false);
  static SwitchOption interfaceFSR_QCDMode0
    (interfaceFSR_QCDMode,"FSR_QCD-OFF","QCD FSR is OFF", 0);
  static SwitchOption interfaceFSR_QCDMode1
    (interfaceFSR_QCDMode,"FSR_QCD-ON","QCD FSR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSR_QEDMode
    ("OnOffFSR_QEDMode",
     "Choice of the on-off QED interaction switch mode",
     &SplittingGenerator::_FSR_QEDMode, 1, false, false);
  static SwitchOption interfaceFSR_QEDMode0
    (interfaceFSR_QEDMode,"FSR_QED-OFF","QED FSR is OFF", 0);
  static SwitchOption interfaceFSR_QEDMode1
    (interfaceFSR_QEDMode,"FSR_QED-ON","QED FSR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceFSR_EWKMode
    ("OnOffFSR_EWKMode",
     "Choice of the on-off EWK interaction switch mode",
     &SplittingGenerator::_FSR_EWKMode, 1, false, false);
  static SwitchOption interfaceFSR_EWKMode0
    (interfaceFSR_EWKMode,"FSR_EWK-OFF","EWK FSR is OFF", 0);
  static SwitchOption interfaceFSR_EWKMode1
    (interfaceFSR_EWKMode,"FSR_EWK-ON","EWK FSR is ON", 1);

  static Switch<SplittingGenerator, int> interfaceUtoUGsplittingMode
    ("OnOffUtoUGsplittingMode",
     "Choice of the on-off  U -> U G  splitting switch mode",
     &SplittingGenerator::_UtoUGsplittingMode, 1, false, false);
  static SwitchOption interfaceUtoUGsplittingMode0
    (interfaceUtoUGsplittingMode,"UtoUG-OFF","U -> U G splitting is OFF", 0);
  static SwitchOption interfaceUtoUGsplittingMode1
    (interfaceUtoUGsplittingMode,"UtoUG-ON","U -> U G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceDtoDGsplittingMode
    ("OnOffDtoDGsplittingMode",
     "Choice of the on-off  D -> D G  splitting switch mode",
     &SplittingGenerator::_DtoDGsplittingMode, 1, false, false);
  static SwitchOption interfaceDtoDGsplittingMode0
    (interfaceDtoDGsplittingMode,"DtoDG-OFF","D -> D G splitting is OFF", 0);
  static SwitchOption interfaceDtoDGsplittingMode1
    (interfaceDtoDGsplittingMode,"DtoDG-ON","D -> D G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceStoSGsplittingMode
    ("OnOffStoSGsplittingMode",
     "Choice of the on-off  S -> S G  splitting switch mode",
     &SplittingGenerator::_StoSGsplittingMode, 1, false, false);
  static SwitchOption interfaceStoSGsplittingMode0
    (interfaceStoSGsplittingMode,"StoSG-OFF","S -> S G splitting is OFF", 0);
  static SwitchOption interfaceStoSGsplittingMode1
    (interfaceStoSGsplittingMode,"StoSG-ON","S -> S G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceCtoCGsplittingMode
    ("OnOffCtoCGsplittingMode",
     "Choice of the on-off  C -> C G  splitting switch mode",
     &SplittingGenerator::_CtoCGsplittingMode, 1, false, false);
  static SwitchOption interfaceCtoCGsplittingMode0
    (interfaceCtoCGsplittingMode,"CtoCG-OFF","C -> C G splitting is OFF", 0);
  static SwitchOption interfaceCtoCGsplittingMode1
    (interfaceCtoCGsplittingMode,"CtoCG-ON","C -> C G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceBtoBGsplittingMode
    ("OnOffBtoBGsplittingMode",
     "Choice of the on-off  B -> B G  splitting switch mode",
     &SplittingGenerator::_BtoBGsplittingMode, 1, false, false);
  static SwitchOption interfaceBtoBGsplittingMode0
    (interfaceBtoBGsplittingMode,"BtoBG-OFF","B -> B G splitting is OFF", 0);
  static SwitchOption interfaceBtoBGsplittingMode1
    (interfaceBtoBGsplittingMode,"BtoBG-ON","B -> B G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceTtoTGsplittingMode
    ("OnOffTtoTGsplittingMode",
     "Choice of the on-off  T -> T G  splitting switch mode",
     &SplittingGenerator::_TtoTGsplittingMode, 1, false, false);
  static SwitchOption interfaceTtoTGsplittingMode0
    (interfaceTtoTGsplittingMode,"TtoTG-OFF","T -> T G splitting is OFF", 0);
  static SwitchOption interfaceTtoTGsplittingMode1
    (interfaceTtoTGsplittingMode,"TtoTG-ON","T -> T G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoGGsplittingMode
    ("OnOffGtoGGsplittingMode",
     "Choice of the on-off  G -> G G  splitting switch mode",
     &SplittingGenerator::_GtoGGsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoGGsplittingMode0
    (interfaceGtoGGsplittingMode,"GtoGG-OFF","G -> G G splitting is OFF", 0);
  static SwitchOption interfaceGtoGGsplittingMode1
    (interfaceGtoGGsplittingMode,"GtoGG-ON","G -> G G splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoUUbarsplittingMode
    ("OnOffGtoUUbarsplittingMode",
     "Choice of the on-off  G -> U Ubar  splitting switch mode",
     &SplittingGenerator::_GtoUUbarsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoUUbarsplittingMode0
    (interfaceGtoUUbarsplittingMode,"GtoUUbar-OFF","G -> U Ubar splitting is OFF", 0);
  static SwitchOption interfaceGtoUUbarsplittingMode1
    (interfaceGtoUUbarsplittingMode,"GtoUUbar-ON","G -> U Ubar splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoDDbarsplittingMode
    ("OnOffGtoDDbarsplittingMode",
     "Choice of the on-off  G -> D Dbar  splitting switch mode",
     &SplittingGenerator::_GtoDDbarsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoDDbarsplittingMode0
    (interfaceGtoDDbarsplittingMode,"GtoDDbar-OFF","G -> D Dbar splitting is OFF", 0);
  static SwitchOption interfaceGtoDDbarsplittingMode1
    (interfaceGtoDDbarsplittingMode,"GtoDDbar-ON","G -> D Dbar splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoSSbarsplittingMode
    ("OnOffGtoSSbarsplittingMode",
     "Choice of the on-off  G -> S Sbar  splitting switch mode",
     &SplittingGenerator::_GtoSSbarsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoSSbarsplittingMode0
    (interfaceGtoSSbarsplittingMode,"GtoSSbar-OFF","G -> S Sbar splitting is OFF", 0);
  static SwitchOption interfaceGtoSSbarsplittingMode1
    (interfaceGtoSSbarsplittingMode,"GtoSSbar-ON","G -> S Sbar splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoCCbarsplittingMode
    ("OnOffGtoCCbarsplittingMode",
     "Choice of the on-off  G -> C Cbar  splitting switch mode",
     &SplittingGenerator::_GtoCCbarsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoCCbarsplittingMode0
    (interfaceGtoCCbarsplittingMode,"GtoCCbar-OFF","G -> C Cbar splitting is OFF", 0);
  static SwitchOption interfaceGtoCCbarsplittingMode1
    (interfaceGtoCCbarsplittingMode,"GtoCCbar-ON","G -> C Cbar splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoBBbarsplittingMode
    ("OnOffGtoBBbarsplittingMode",
     "Choice of the on-off  G -> B Bbar  splitting switch mode",
     &SplittingGenerator::_GtoBBbarsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoBBbarsplittingMode0
    (interfaceGtoBBbarsplittingMode,"GtoBBbar-OFF","G -> B Bbar splitting is OFF", 0);
  static SwitchOption interfaceGtoBBbarsplittingMode1
    (interfaceGtoBBbarsplittingMode,"GtoBBbar-ON","G -> B Bbar splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceGtoTTbarsplittingMode
    ("OnOffGtoTTbarsplittingMode",
     "Choice of the on-off  G -> T Tbar  splitting switch mode",
     &SplittingGenerator::_GtoTTbarsplittingMode, 1, false, false);
  static SwitchOption interfaceGtoTTbarsplittingMode0
    (interfaceGtoTTbarsplittingMode,"GtoTTbar-OFF","G -> T Tbar splitting is OFF", 0);
  static SwitchOption interfaceGtoTTbarsplittingMode1
    (interfaceGtoTTbarsplittingMode,"GtoTTbar-ON","G -> T Tbar splitting is ON", 1);

  // this is new:
  static Switch<SplittingGenerator, int> interfaceUtoUGammasplittingMode
    ("OnOffUtoUGammasplittingMode",
     "Choice of the on-off  U -> U Gamma  splitting switch mode",
     &SplittingGenerator::_UtoUGammasplittingMode, 1, false, false);
  static SwitchOption interfaceUtoUGammasplittingMode0
    (interfaceUtoUGammasplittingMode,"UtoUGamma-OFF","U -> U Gamma splitting is OFF", 0);
  static SwitchOption interfaceUtoUGammasplittingMode1
    (interfaceUtoUGammasplittingMode,"UtoUGamma-ON","U -> U Gamma splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceDtoDGammasplittingMode
    ("OnOffDtoDGammasplittingMode",
     "Choice of the on-off  D -> D Gamma  splitting switch mode",
     &SplittingGenerator::_DtoDGammasplittingMode, 1, false, false);
  static SwitchOption interfaceDtoDGammasplittingMode0
    (interfaceDtoDGammasplittingMode,"DtoDGamma-OFF","D -> D Gamma splitting is OFF", 0);
  static SwitchOption interfaceDtoDGammasplittingMode1
    (interfaceDtoDGammasplittingMode,"DtoDGamma-ON","D -> D Gamma splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceStoSGammasplittingMode
    ("OnOffStoSGammasplittingMode",
     "Choice of the on-off  S -> S Gamma  splitting switch mode",
     &SplittingGenerator::_StoSGammasplittingMode, 1, false, false);
  static SwitchOption interfaceStoSGammasplittingMode0
    (interfaceStoSGammasplittingMode,"StoSGamma-OFF","S -> S Gamma splitting is OFF", 0);
  static SwitchOption interfaceStoSGammasplittingMode1
    (interfaceStoSGammasplittingMode,"StoSGamma-ON","S -> S Gamma splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceCtoCGammasplittingMode
    ("OnOffCtoCGammasplittingMode",
     "Choice of the on-off  C -> C Gamma  splitting switch mode",
     &SplittingGenerator::_CtoCGammasplittingMode, 1, false, false);
  static SwitchOption interfaceCtoCGammasplittingMode0
    (interfaceCtoCGammasplittingMode,"CtoCGamma-OFF","C -> C Gamma splitting is OFF", 0);
  static SwitchOption interfaceCtoCGammasplittingMode1
    (interfaceCtoCGammasplittingMode,"CtoCGamma-ON","C -> C Gamma splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceBtoBGammasplittingMode
    ("OnOffBtoBGammasplittingMode",
     "Choice of the on-off  B -> B Gamma  splitting switch mode",
     &SplittingGenerator::_BtoBGammasplittingMode, 1, false, false);
  static SwitchOption interfaceBtoBGammasplittingMode0
    (interfaceBtoBGammasplittingMode,"BtoBGamma-OFF","B -> B Gamma splitting is OFF", 0);
  static SwitchOption interfaceBtoBGammasplittingMode1
    (interfaceBtoBGammasplittingMode,"BtoBGamma-ON","B -> B Gamma splitting is ON", 1);

  static Switch<SplittingGenerator, int> interfaceTtoTGammasplittingMode
    ("OnOffTtoTGammasplittingMode",
     "Choice of the on-off  T -> T Gamma  splitting switch mode",
     &SplittingGenerator::_TtoTGammasplittingMode, 1, false, false);
  static SwitchOption interfaceTtoTGammasplittingMode0
    (interfaceTtoTGammasplittingMode,"TtoTGamma-OFF","T -> T Gamma splitting is OFF", 0);
  static SwitchOption interfaceTtoTGammasplittingMode1
    (interfaceTtoTGammasplittingMode,"TtoTGamma-ON","T -> T Gamma splitting is ON", 1);

  static Reference<SplittingGenerator,ShowerAlpha> interfaceIS_ShowerAlphaQED
    ("IS_ShowerAlphaQED", "A reference to the IS_ShowerAlphaQED object", 
     &Herwig::SplittingGenerator::_pointerIS_ShowerAlphaQED, false, false, true, false);
  
  static Reference<SplittingGenerator,ShowerAlpha> interfaceFS_ShowerAlphaQED
    ("FS_ShowerAlphaQED", "A reference to the FS_ShowerAlphaQED object", 
     &Herwig::SplittingGenerator::_pointerFS_ShowerAlphaQED, false, false, true, false);
  // up to here.

  static Reference<SplittingGenerator,ShowerAlpha> interfaceIS_ShowerAlphaQCD
    ("IS_ShowerAlphaQCD", "A reference to the IS_ShowerAlphaQCD object", 
     &Herwig::SplittingGenerator::_pointerIS_ShowerAlphaQCD, false, false, true, false);
  
  static Reference<SplittingGenerator,ShowerAlpha> interfaceFS_ShowerAlphaQCD
    ("FS_ShowerAlphaQCD", "A reference to the FS_ShowerAlphaQCD object", 
     &Herwig::SplittingGenerator::_pointerFS_ShowerAlphaQCD, false, false, true, false);
  
  static Reference<SplittingGenerator,ShowerConstrainer> interfaceShowerConstrainer
    ("ShowerConstrainer", "A reference to the ShowerConstrainer object", 
     &Herwig::SplittingGenerator::_pointerShowerConstrainer, false, false, true, false);
  
}


//--------------------------------------------------------------------------------

bool SplittingGenerator::isInteractionON(const ShowerIndex::InteractionType interaction) const {
  int mode = 0;
  switch ( interaction ) {
  case ShowerIndex::QCD : mode = _QCDinteractionMode; break; 
  case ShowerIndex::QED : mode = _QEDinteractionMode; break; 
  case ShowerIndex::EWK : mode = _EWKinteractionMode; break; 
  }
  return mode;
}

bool SplittingGenerator::isISRadiationON(const ShowerIndex::InteractionType interaction) const {
  int mode = 0;
  if ( isInteractionON(interaction) && isISRadiationON() ) { 
    switch ( interaction ) {
    case ShowerIndex::QCD : mode = _ISR_QCDMode; break; 
    case ShowerIndex::QED : mode = _ISR_QEDMode; break; 
    case ShowerIndex::EWK : mode = _ISR_EWKMode; break; 
    }
  }
  return mode;
}  

bool SplittingGenerator::isFSRadiationON(const ShowerIndex::InteractionType interaction) const {
  int mode = 0;
  if ( isInteractionON(interaction) && isFSRadiationON() ) { 
    switch ( interaction ) {
    case ShowerIndex::QCD : mode = _FSR_QCDMode; break; 
    case ShowerIndex::QED : mode = _FSR_QEDMode; break; 
    case ShowerIndex::EWK : mode = _FSR_EWKMode; break; 
    }
  }
  return mode;
}


pair<ShoKinPtr, tSudakovFormFactorPtr> SplittingGenerator::chooseForwardBranching
( tPartCollHdlPtr ch, ShowerParticle & particle, const bool reverseAngularOrder ) const {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "SplittingGenerator::chooseForwardBranching(): full _____________________________" << endl; 
  }
  
  Energy newQ = Energy();
  tSudakovFormFactorPtr sudakov = tSudakovFormFactorPtr();
  
  // First, find the eventual branching, corresponding to the highest scale.
  for (int i = 0; i < ShowerIndex::NumInteractionTypes; ++i ) {

    ShowerIndex index;
    index.id = abs( particle.data().id() ); 
    index.interaction = ShowerIndex::int2Interaction( i );
    index.timeFlag = ShowerIndex::FS; 
    if ( _multimapSudakov.find( index ) != _multimapSudakov.end() ) {
      for ( CollecIndexSudakov::const_iterator cit =  
	      _multimapSudakov.lower_bound( index ); 
	    cit != _multimapSudakov.upper_bound( index ); ++cit ) {
	tSudakovFormFactorPtr candidateSudakov = cit->second;
	Energy candidateNewQ = 0*MeV;
	if ( candidateSudakov ) {
	  // *** ACHTUNG *** some preliminary angular ordering...
	  // the z-factor should be put in here by hand.
	  candidateNewQ = candidateSudakov->
	    generateNextBranching( ch, particle.evolutionScales()[i], reverseAngularOrder );

	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	    generator()->log() << "  trying "
			       << candidateSudakov->splitFun()->idEmitter() << " -> ("
			       << (dynamic_ptr_cast<tSplitFun1to2Ptr>(candidateSudakov->splitFun()))->
	      idFirstProduct() << ", "
			       << (dynamic_ptr_cast<tSplitFun1to2Ptr>(candidateSudakov->splitFun()))->
	      idSecondProduct() << ") in "
			       << candidateSudakov->splitFun()->interactionType() << "; "
			       << "(Q0 -> Q) = (" 
			       << particle.evolutionScales()[index.interaction]/MeV
			       << " -> " 
			       << (candidateNewQ/MeV > 0 ? candidateNewQ/MeV : 0)
			       << ") MeV"  << endl;
	  }
	  
	  if ( ( candidateNewQ > newQ  &&  ! reverseAngularOrder ) ||
	       ( candidateNewQ < newQ  &&  reverseAngularOrder ) ) {
	    newQ = candidateNewQ;
	    sudakov = candidateSudakov;
	  } 
	}
      } 
    }
  }
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    generator()->log() << "  " << particle.data().PDGName();
      //		       << " [" << particle.number() << "]: ";
    if( sudakov ) 
      generator()->log() << sudakov->splitFun()->idEmitter() << " -> ("
			 << (dynamic_ptr_cast<tSplitFun1to2Ptr>(sudakov->splitFun()))->
	idFirstProduct() << ", "
			 << (dynamic_ptr_cast<tSplitFun1to2Ptr>(sudakov->splitFun()))->
	idSecondProduct() << ") in "
			 << sudakov->splitFun()->interactionType() << "; "
			 << "(Q0 -> Q) = (" 
			 << particle.evolutionScales()[sudakov->splitFun()->interactionType()]/MeV 
			 << " -> " 
			 << (newQ/MeV > 0 ? newQ/MeV : 0)
			 << ") MeV"  << endl;
    else 
      generator()->log() << " won't branch." << endl; 
  }

  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  if ( newQ && sudakov ) {
    if ( sudakov->splitFun() ) {

      // For the time being we are considering only 1->2 branching
      tSplitFun1to2Ptr splitFun = 
	dynamic_ptr_cast< tSplitFun1to2Ptr >( sudakov->splitFun() );
      if ( splitFun ) {	  
	Lorentz5Momentum p, n, ppartner, pcm;
	double th;
	Energy nnorm;
	if ( particle.isFromHardSubprocess() ) {
	  p = particle.momentum();
	  //***LOOKHERE*** We choose  n  naively for the time being.  
	  // Lorentz5Momentum n = Lorentz5Momentum( 0.0, - p.vect() ); 
	  // for test purposes: n points into the negative z-direction:
	  ppartner = particle.partners()[sudakov->splitFun()->interactionType()]
	    ->momentum();
	  pcm = p; 
	  pcm.boost( (p + ppartner).findBoostToCM() );	  
	  th = (pi);
	  nnorm = pcm.vect().mag(); 
	  nnorm*=1.;
	  Vector3 nv = cos(th)*pcm.vect().unit() 
	    + sin(th)*pcm.vect().orthogonal().unit(); 
	  n = Lorentz5Momentum( 0.0, nnorm*nv ); 
	  n.boost( -(p + ppartner).findBoostToCM() );
	  // n = Lorentz5Momentum( 2000.0, 3000.0, -1000.0, sqrt(14.)*1000.0); 
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	    generator()->log() << "  chosen normalization = " << nnorm 
			       << ", theta = " << th << endl 
			       << "  with p_partner = " << ppartner << endl
			       << "  and    p in cm = " << pcm << endl;
	  }
	} else {
	  p = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	                     ->showerKinematics()->getBasis()[0];
	  n = dynamic_ptr_cast<ShowerParticlePtr>(particle.parents()[0])
	                     ->showerKinematics()->getBasis()[1];
	} 

	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	  generator()->log() << "  create ShowerKinematics with " 
			     << endl 
			     << "  p = " << p << " and n = " << n << endl;
	}
	
	Ptr< FS_QtildaShowerKinematics1to2 >::pointer showerKin = 
	  new_ptr (FS_QtildaShowerKinematics1to2( p, n ) );

        showerKin->qtilde( newQ );
	showerKin->z( sudakov->z() );
	showerKin->phi( sudakov->phi() );

// 	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
// 	  generator()->log() << "SplittingGenerator::chooseForwardBranching(): end full _________________________" << endl; 
// 	}

	return pair<ShoKinPtr,tSudakovFormFactorPtr>( showerKin, sudakov );
      }
    }
  }
  
//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     generator()->log() << "SplittingGenerator::chooseForwardBranching(): end full _________________________" << endl; 
//   }
  
  return pair<ShoKinPtr,tSudakovFormFactorPtr>( ShoKinPtr(), tSudakovFormFactorPtr() );

}
 

pair<ShoKinPtr, tSudakovFormFactorPtr> SplittingGenerator::chooseBackwardBranching
( tPartCollHdlPtr ch, ShowerParticle & particle ) const {

  Energy newQ = Energy();
  tSudakovFormFactorPtr sudakov = tSudakovFormFactorPtr();
  
  // First, find the eventual branching, corresponding to the highest scale.
  for (int i = 0; i < ShowerIndex::NumInteractionTypes; ++i ) {
    ShowerIndex index;
    index.id = abs( particle.data().id() ); 
    index.interaction = ShowerIndex::int2Interaction( i );
    index.timeFlag = ShowerIndex::IS; 
    if ( _multimapSudakov.find( index ) != _multimapSudakov.end() ) {
      for ( CollecIndexSudakov::const_iterator cit =  
	      _multimapSudakov.lower_bound( index ); 
	    cit != _multimapSudakov.upper_bound( index ); ++cit ) {
	tSudakovFormFactorPtr candidateSudakov = cit->second;
	Energy candidateNewQ = Energy();
	if ( candidateSudakov ) {
	  candidateNewQ = candidateSudakov->generateNextBranching( ch, particle.evolutionScales()[i] );
	  if ( candidateNewQ > newQ ) {
	    newQ = candidateNewQ;
	    sudakov = candidateSudakov;
	  } 
	}
      } 
    }
  }

  // Then, if a branching has been selected, create the proper 
  // ShowerKinematics object which contains the kinematics information
  // about such branching. Notice that the cases 1->2 and 1->3
  // branching should be treated separately.
  if ( newQ && sudakov ) {
    if ( sudakov->splitFun() ) {

      // For the time being we are considering only 1->2 branching
      tSplitFun1to2Ptr splitFun = 
	dynamic_ptr_cast< tSplitFun1to2Ptr >( sudakov->splitFun() );
      if ( splitFun ) {	  

        //***LOOKHERE*** Do something similar as in chooseForwardBranching
        //               but use IS_QtildaShowerKinematics1to2 instead.

      }
    }
  }

  return pair<ShoKinPtr,tSudakovFormFactorPtr>( ShoKinPtr(), tSudakovFormFactorPtr() );

}


void SplittingGenerator::generateBranchingKinematics 
( tPartCollHdlPtr ch, ShowerParticle & particle,
  tShoKinPtr showerKin, const tSudakovFormFactorPtr sudakov ) const {

  //***LOOKHERE*** Complete the kinematics of the branching by filling
  //               the eventual missing bits of the ShowerKinematics
  //               object created, and already at least partially filled,
  //               by  chooseForwardBranching  or  chooseBackwardBranching.
  //               Notice that this part could remain empty if such 
  //               ShowerKinematics object is already completely filled.

}
 

void SplittingGenerator::initializeRun() {

  tFEGPtr fulleg = dynamic_ptr_cast<tFEGPtr>( generator() );
  Energy maxCMEnergy = Energy();
  if (fulleg) maxCMEnergy = fulleg->maximumCMEnergy(); 

  ShowerIndex index;
  SplitFunPtr splitFun;
  SudakovFormFactorPtr sudakov;
  Energy minQValue = Energy();
  Energy maxQValue = _pointerShowerConstrainer->convertMassScaleToQScale( maxCMEnergy );  

  //===========
  //=== QCD ===
  //===========
  index.interaction = ShowerIndex::QCD;
  minQValue = _pointerShowerConstrainer->cutoffQScale( index.interaction );
  if ( isInteractionON( index.interaction ) ) {

    index.id = ParticleID::d; //--- D ---
    //---  D -> D + G  ---
    if ( isDtoDGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) ); 
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }

    index.id = ParticleID::u; //--- U ---
    //---  U -> U + G  ---
    if ( isUtoUGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );   
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) ); 
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }

    index.id = ParticleID::s; //--- S ---
    //---  S -> S + G  ---
    if ( isStoSGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) ); 
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }

    index.id = ParticleID::c; //--- C ---
    //---  C -> C + G  ---
    if ( isCtoCGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }

    index.id = ParticleID::b; //--- B ---
    //---  B -> B + G  ---
    if ( isBtoBGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }

    index.id = ParticleID::t; //--- T ---
    //---  T -> T + G  ---
    if ( isTtoTGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }

    index.id = ParticleID::g; //--- GLUON ---
    //---  G -> G + G  ---
    if ( isGtoGGsplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoGGSplitFun() );
	sudakov  = new_ptr( GtoGGSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoGGSplitFun() );     
	sudakov  = new_ptr( GtoGGSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						    minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
    //---  G -> U + Ubar  ---
    if ( isGtoUUbarsplittingON() ) { 
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoQQbarSplitFun( ParticleID::u,
						 getParticleData( ParticleID::u )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoQQbarSplitFun( ParticleID::u,
						 getParticleData( ParticleID::u )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }    
    //---  G -> D + Dbar  ---
    if ( isGtoDDbarsplittingON() ) { 
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoQQbarSplitFun( ParticleID::d,
						 getParticleData( ParticleID::d )->mass() ) );     
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoQQbarSplitFun( ParticleID::d,
						 getParticleData( ParticleID::d )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
    //---  G -> S + Sbar  ---
    if ( isGtoSSbarsplittingON() ) { 
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoQQbarSplitFun( ParticleID::s,
						 getParticleData( ParticleID::s )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoQQbarSplitFun( ParticleID::s,
						 getParticleData( ParticleID::s )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
    //---  G -> C + Cbar  ---
    if ( isGtoCCbarsplittingON() ) { 
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoQQbarSplitFun( ParticleID::c,
						 getParticleData( ParticleID::c )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoQQbarSplitFun( ParticleID::c,
						 getParticleData( ParticleID::c )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
    //---  G -> B + Bbar  ---
    if ( isGtoBBbarsplittingON() ) { 
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoQQbarSplitFun( ParticleID::b,
						 getParticleData( ParticleID::b )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoQQbarSplitFun( ParticleID::b,
						 getParticleData( ParticleID::b )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }    
    //---  G -> T + Tbar  ---
    if ( isGtoTTbarsplittingON() ) { 
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_GtoQQbarSplitFun( ParticleID::t,
						 getParticleData( ParticleID::t )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_GtoQQbarSplitFun( ParticleID::t,
						 getParticleData( ParticleID::t )->mass() ) );
	sudakov  = new_ptr( GtoQQbarSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQCD,
						       minQValue, maxQValue ) );
        sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
    
  } // === end QCD ===

  //===========
  //=== QED ===
  //===========
  index.interaction = ShowerIndex::QED;
  minQValue = _pointerShowerConstrainer->cutoffQScale( index.interaction );
  if ( isInteractionON( index.interaction ) ) {

    index.id = ParticleID::u; //--- U ---
    //---  U -> U + Gamma  ---
    if ( isUtoUGammasplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
 
    index.id = ParticleID::d; //--- D ---
    //---  D -> D + Gamma  ---
    if ( isDtoDGammasplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
 
    index.id = ParticleID::s; //--- S ---
    //---  S -> S + Gamma  ---
    if ( isStoSGammasplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
 
    index.id = ParticleID::c; //--- C ---
    //---  C -> C + Gamma  ---
    if ( isCtoCGammasplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
 
    index.id = ParticleID::b; //--- B ---
    //---  B -> B + Gamma  ---
    if ( isBtoBGammasplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
 
    index.id = ParticleID::t; //--- T ---
    //---  T -> T + Gamma  ---
    if ( isTtoTGammasplittingON() ) {
      if ( isISRadiationON( index.interaction ) ) {  // Initial State Radiation
	index.timeFlag = ShowerIndex::IS;
	splitFun = new_ptr( IS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerIS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
      if ( isFSRadiationON( index.interaction ) ) {  // Final State Radiation
	index.timeFlag = ShowerIndex::FS;
	splitFun = new_ptr( FS_QtoQGammaSplitFun( index.id, getParticleData( index.id )->mass() ) );
	sudakov  = new_ptr( QtoQGammaSudakovFormFactor( splitFun, _pointerFS_ShowerAlphaQED, minQValue, maxQValue ) );
	sudakov->setupLookupTables();
	_multimapSudakov.insert( pair<ShowerIndex,SudakovFormFactorPtr>( index, sudakov ) );
      }
    }
 
  } // === end QED ===
  
  //===========
  //=== EWK ===
  //===========
  index.interaction = ShowerIndex::EWK;
  minQValue = _pointerShowerConstrainer->cutoffQScale( index.interaction );
  if ( isInteractionON( index.interaction ) ) {

    //...

  } // === end EWK ===

  // Debugging
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    debuggingInfo();
  }

}


void SplittingGenerator::debuggingInfo() {

  generator()->log() << "SplittingGenerator::debuggingInfo() begin ______________________________________" << endl; 
  generator()->log() << "  no of initialized Sudakov FF's = " 
		     << _multimapSudakov.size() << endl
		     << "  id\t" << "int\t" << "timeFlag" << endl;
  for ( CollecIndexSudakov::const_iterator cit = _multimapSudakov.begin();
  	cit != _multimapSudakov.end(); ++cit ) {
    generator()->log() << "  " << cit->first.id 
		       << "\t" << cit->first.interaction 
		       << "\t" << cit->first.timeFlag 
		       << endl;
  }

  //  generator()->log() << "SplittingGenerator::debuggingInfo() end ________________________________________" << endl; 
}
