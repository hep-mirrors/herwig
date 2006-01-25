// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerVariables class.
//

#include "ShowerVariables.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerVariables.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ShowerVariables::~ShowerVariables() {}

void ShowerVariables::persistentOutput(PersistentOStream & os) const {
  os << _multiScaleShowerMode 
     << _decayBeforeShowerMode
     << _cutoffQCDMassScale
     << _cutoffQEDMassScale
     << _cutoffEWKMassScale
     << _kinCutoffScale
     << _meCorrMode
     << _qqgPSMode
     << _particlesDecayBeforeShower;
}

void ShowerVariables::persistentInput(PersistentIStream & is, int) {
  is >> _multiScaleShowerMode 
     >> _decayBeforeShowerMode 
     >> _cutoffQCDMassScale
     >> _cutoffQEDMassScale
     >> _cutoffEWKMassScale
     >> _kinCutoffScale
     >> _meCorrMode
     >> _qqgPSMode
     >> _particlesDecayBeforeShower;
}

Energy ShowerVariables::HUGEMASS = 1.0e+20 * GeV;  // more then the Plank scale!

ClassDescription<ShowerVariables> ShowerVariables::initShowerVariables;
// Definition of the static class description member.

void ShowerVariables::Init() {

  static ClassDocumentation<ShowerVariables> documentation
    ("This class is responsible for keeping all of the constraints on the showering.");

  static Switch<ShowerVariables, int> interfaceMultiScaleShowerMode
    ("OnOffMultiScaleShowerMode",
     "Choice of the multi-scale shower mode",
     &ShowerVariables::_multiScaleShowerMode, 1, false, false);
  static SwitchOption interfaceMultiScaleShowerMode0
    (interfaceMultiScaleShowerMode,"MultiScaleShower-OFF","multi-scale shower is OFF", 0);
  static SwitchOption interfaceMultiScaleShowerMode1
    (interfaceMultiScaleShowerMode,"MultiScaleShower-ON","multi-scale shower is ON", 1);
  
  static Switch<ShowerVariables, int> ifaceMECorrMode
    ("MECorrMode",
     "Choice of the ME Correction Mode",
     &ShowerVariables::_meCorrMode, 1, false, false);
  static SwitchOption off
    (ifaceMECorrMode,"MEC-off","MECorrections off", 0);
  static SwitchOption on
    (ifaceMECorrMode,"MEC-on","hard+soft on", 1);
  static SwitchOption hard
    (ifaceMECorrMode,"MEC-hard","only hard on", 2);
  static SwitchOption soft
    (ifaceMECorrMode,"MEC-soft","only soft on", 3);

  static Switch<ShowerVariables, int> ifaceqqgPSMode
    ("qqgPSMode",
     "Choice of initial conditions, tested for qqg only",
     &ShowerVariables::_qqgPSMode, 0, false, false);
  static SwitchOption symm
    (ifaceqqgPSMode,"PS-symm",
     "most symmetric choice of initial conditions (default)", 0);
  static SwitchOption asy
    (ifaceqqgPSMode,"PS-asy",
     "most asymmetric choice of initial conditions, quark larger Q0", 1);
  static SwitchOption rnd
    (ifaceqqgPSMode,"PS-rnd",
     "asymmetric, large Q0 assigned randomly", 2);

  static Switch<ShowerVariables, int> interfaceDecayBeforeShowerMode
    ("OnOffDecayBeforeShowerMode",
     "Choice of the decay before shower mode",
     &ShowerVariables::_decayBeforeShowerMode, 0, false, false);
  static SwitchOption interfaceDecayBeforeShowerMode0
    (interfaceDecayBeforeShowerMode,"DecayBeforeShower-OFF","decay before shower is OFF", 0);
  static SwitchOption interfaceDecayBeforeShowerMode1
    (interfaceDecayBeforeShowerMode,"DecayBeforeShower-ON","decay before shower is ON", 1);

  static Parameter<ShowerVariables,Energy>
    interfaceCutoffQCD ("CutoffQCDMassScale",
			"low energy cutoff mass scale for QCD radiation  (unit [GeV])",
			&ShowerVariables::_cutoffQCDMassScale, GeV, 
			0.0*GeV, 0.0*GeV, 10.0*GeV,false,false,false);
  static Parameter<ShowerVariables,Energy>
    interfaceCutoffQED ("CutoffQEDMassScale",
			"low energy cutoff mass scale for QED radiation  (unit [GeV])",
			&ShowerVariables::_cutoffQEDMassScale, GeV, 
			0.0005*GeV, 0.0*GeV, 10.0*GeV,false,false,false);
  static Parameter<ShowerVariables,Energy>
    interfaceCutoffEWK ("CutoffEWKMassScale",
			"low energy cutoff mass scale for EWK radiation  (unit [GeV])",
			&ShowerVariables::_cutoffEWKMassScale, GeV, 
			91.0*GeV, 0.0*GeV, 1000.0*GeV,false,false,false);

  static Parameter<ShowerVariables,Energy>
    interfaceKinScale ("cutoffKinScale",
		       "kinematic cutoff scale for the parton shower phase space (unit [GeV])",
		       &ShowerVariables::_kinCutoffScale, GeV, 
		       0.75*GeV, 0.001*GeV, 10.0*GeV,false,false,false);

}

//-------------------------------------------------------------------------------

void ShowerVariables::reset() {
  _stopShowerAtMassScale = _vetoBelowPtScale = Energy(); 
  _vetoAbovePtScale = HUGEMASS;
}

void ShowerVariables::initialize() {

  //***LOOKHERE*** Here is the list of particles that should decay before
  //               showering, in the case the switch for decay before shower
  //               is 1 (ON).   Insert only positive id.
  //
  _particlesDecayBeforeShower.insert( 1000001 ); //  SUSY_d_L 
  _particlesDecayBeforeShower.insert( 1000002 ); //  SUSY_u_L 
  _particlesDecayBeforeShower.insert( 1000003 ); //  SUSY_s_L 
  _particlesDecayBeforeShower.insert( 1000004 ); //  SUSY_c_L 
  _particlesDecayBeforeShower.insert( 1000005 ); //  SUSY_b_1 
  _particlesDecayBeforeShower.insert( 1000006 ); //  SUSY_t_1 
  _particlesDecayBeforeShower.insert( 1000011 ); //  SUSY_e_Lminus 
  _particlesDecayBeforeShower.insert( 1000012 ); //  SUSY_nu_eL 
  _particlesDecayBeforeShower.insert( 1000013 ); //  SUSY_mu_Lminus 
  _particlesDecayBeforeShower.insert( 1000014 ); //  SUSY_nu_muL 
  _particlesDecayBeforeShower.insert( 1000015 ); //  SUSY_tau_1minus 
  _particlesDecayBeforeShower.insert( 1000016 ); //  SUSY_nu_tauL 
  _particlesDecayBeforeShower.insert( 1000021 ); //  SUSY_g 
  _particlesDecayBeforeShower.insert( 1000022 ); //  SUSY_chi_10 
  _particlesDecayBeforeShower.insert( 1000023 ); //  SUSY_chi_20 
  _particlesDecayBeforeShower.insert( 1000024 ); //  SUSY_chi_1plus 
  _particlesDecayBeforeShower.insert( 1000025 ); //  SUSY_chi_30 
  _particlesDecayBeforeShower.insert( 1000035 ); //  SUSY_chi_40 
  _particlesDecayBeforeShower.insert( 1000037 ); //  SUSY_chi_2plus 
  _particlesDecayBeforeShower.insert( 1000039 ); //  SUSY_gravitino 
  _particlesDecayBeforeShower.insert( 2000001 ); //  SUSY_d_R 
  _particlesDecayBeforeShower.insert( 2000002 ); //  SUSY_u_R 
  _particlesDecayBeforeShower.insert( 2000003 ); //  SUSY_s_R 
  _particlesDecayBeforeShower.insert( 2000004 ); //  SUSY_c_R 
  _particlesDecayBeforeShower.insert( 2000005 ); //  SUSY_b_2 
  _particlesDecayBeforeShower.insert( 2000006 ); //  SUSY_t_2 
  _particlesDecayBeforeShower.insert( 2000011 ); //  SUSY_e_Rminus 
  _particlesDecayBeforeShower.insert( 2000012 ); //  SUSY_nu_eR 
  _particlesDecayBeforeShower.insert( 2000013 ); //  SUSY_mu_Rminus 
  _particlesDecayBeforeShower.insert( 2000014 ); //  SUSY_nu_muR 
  _particlesDecayBeforeShower.insert( 2000015 ); //  SUSY_tau_2minus 
  _particlesDecayBeforeShower.insert( 2000016 ); //  SUSY_nu_tauR 
  //***endLOOKHERE*** 

}


Energy ShowerVariables::cutoffMassScale(const ShowerIndex::InteractionType interaction) const {
  Energy cutoff = Energy();
  switch ( interaction ) {
  case ShowerIndex::QCD : cutoff = _cutoffQCDMassScale; break; 
  case ShowerIndex::QED : cutoff = _cutoffQEDMassScale; break; 
  case ShowerIndex::EWK : cutoff = _cutoffEWKMassScale; break; 
  }
  return cutoff;
}


Energy ShowerVariables::cutoffQScale(const ShowerIndex::InteractionType interaction) const {
  return convertMassScaleToQScale( cutoffMassScale( interaction ) );
}





