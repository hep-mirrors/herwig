// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerVariables class.
//

#include "ShowerVariables.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerVariables.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

const Energy ShowerVariables::HUGEMASS = 1.0e+20 * GeV;  // more then the Plank scale!

ShowerVariables::ShowerVariables() :
  _cutoffQCDMassScale( 1.0*GeV ),
  _cutoffQEDMassScale( 0.51*MeV ),
  _cutoffEWKMassScale( 91.0*GeV ),
  _kinCutoffScale( .75*GeV ),
  _meCorrMode(1),
  _stopShowerAtMassScale( Energy() ),
  _vetoAbovePtScale( HUGEMASS ), 
  _vetoBelowPtScale( Energy() ),
  _a(0.3), _b(2.3), _c(0.3*GeV),
  _initialenhance(1.),_finalenhance(1.),
  _finalFinalConditions(0),
  _initialFinalDecayConditions(0),
  _useMEForT2(false)
{
  _inputparticlesDecayInShower.push_back( 6 ); //  top
  _inputparticlesDecayInShower.push_back( 1000001 ); //  SUSY_d_L 
  _inputparticlesDecayInShower.push_back( 1000002 ); //  SUSY_u_L 
  _inputparticlesDecayInShower.push_back( 1000003 ); //  SUSY_s_L 
  _inputparticlesDecayInShower.push_back( 1000004 ); //  SUSY_c_L 
  _inputparticlesDecayInShower.push_back( 1000005 ); //  SUSY_b_1 
  _inputparticlesDecayInShower.push_back( 1000006 ); //  SUSY_t_1 
  _inputparticlesDecayInShower.push_back( 1000011 ); //  SUSY_e_Lminus 
  _inputparticlesDecayInShower.push_back( 1000012 ); //  SUSY_nu_eL 
  _inputparticlesDecayInShower.push_back( 1000013 ); //  SUSY_mu_Lminus 
  _inputparticlesDecayInShower.push_back( 1000014 ); //  SUSY_nu_muL 
  _inputparticlesDecayInShower.push_back( 1000015 ); //  SUSY_tau_1minus 
  _inputparticlesDecayInShower.push_back( 1000016 ); //  SUSY_nu_tauL 
  _inputparticlesDecayInShower.push_back( 1000021 ); //  SUSY_g 
  _inputparticlesDecayInShower.push_back( 1000022 ); //  SUSY_chi_10 
  _inputparticlesDecayInShower.push_back( 1000023 ); //  SUSY_chi_20 
  _inputparticlesDecayInShower.push_back( 1000024 ); //  SUSY_chi_1plus 
  _inputparticlesDecayInShower.push_back( 1000025 ); //  SUSY_chi_30 
  _inputparticlesDecayInShower.push_back( 1000035 ); //  SUSY_chi_40 
  _inputparticlesDecayInShower.push_back( 1000037 ); //  SUSY_chi_2plus 
  _inputparticlesDecayInShower.push_back( 1000039 ); //  SUSY_gravitino 
  _inputparticlesDecayInShower.push_back( 2000001 ); //  SUSY_d_R 
  _inputparticlesDecayInShower.push_back( 2000002 ); //  SUSY_u_R 
  _inputparticlesDecayInShower.push_back( 2000003 ); //  SUSY_s_R 
  _inputparticlesDecayInShower.push_back( 2000004 ); //  SUSY_c_R 
  _inputparticlesDecayInShower.push_back( 2000005 ); //  SUSY_b_2 
  _inputparticlesDecayInShower.push_back( 2000006 ); //  SUSY_t_2 
  _inputparticlesDecayInShower.push_back( 2000011 ); //  SUSY_e_Rminus 
  _inputparticlesDecayInShower.push_back( 2000012 ); //  SUSY_nu_eR 
  _inputparticlesDecayInShower.push_back( 2000013 ); //  SUSY_mu_Rminus 
  _inputparticlesDecayInShower.push_back( 2000014 ); //  SUSY_nu_muR 
  _inputparticlesDecayInShower.push_back( 2000015 ); //  SUSY_tau_2minus 
  _inputparticlesDecayInShower.push_back( 2000016 ); //  SUSY_nu_tauR 
  _inputparticlesDecayInShower.push_back( 25      ); //  h0
  _inputparticlesDecayInShower.push_back( 35      ); //  H0
  _inputparticlesDecayInShower.push_back( 36      ); //  A0
  _inputparticlesDecayInShower.push_back( 37      ); //  H+
}

void ShowerVariables::persistentOutput(PersistentOStream & os) const {
  os << _cutoffQCDMassScale
     << _cutoffQEDMassScale
     << _cutoffEWKMassScale
     << _kinCutoffScale
     << _meCorrMode
     << _inputparticlesDecayInShower
     << _particlesDecayInShower << _a << _b << _c
     << _globalParameters
     << _finalFinalConditions
     << _initialFinalDecayConditions << _useMEForT2;
}

void ShowerVariables::persistentInput(PersistentIStream & is, int) {
  is >> _cutoffQCDMassScale
     >> _cutoffQEDMassScale
     >> _cutoffEWKMassScale
     >> _kinCutoffScale
     >> _meCorrMode
     >> _inputparticlesDecayInShower
     >> _particlesDecayInShower >> _a >> _b >> _c
     >> _globalParameters
     >> _finalFinalConditions
     >> _initialFinalDecayConditions >> _useMEForT2;
}

ClassDescription<ShowerVariables> ShowerVariables::initShowerVariables;
// Definition of the static class description member.

void ShowerVariables::Init() {

  static ClassDocumentation<ShowerVariables> documentation
    ("This class is responsible for keeping all of the constraints on the showering.");

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
		       "kinematic cutoff scale for the parton shower phase"
		       " space (unit [GeV])",
		       &ShowerVariables::_kinCutoffScale, GeV, 
		       0.75*GeV, 0.001*GeV, 10.0*GeV,false,false,false);

  static Parameter<ShowerVariables,double> interfaceaParameter
    ("aParameter",
     "The a parameter for the kinematic cut-off",
     &ShowerVariables::_a, 0.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ShowerVariables,double> interfacebParameter
    ("bParameter",
     "The b parameter for the kinematic cut-off",
     &ShowerVariables::_b, 2.3, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<ShowerVariables,Energy> interfacecParameter
    ("cParameter",
     "The c parameter for the kinematic cut-off",
     &ShowerVariables::_c, GeV, 0.3*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);


  static Reference<ShowerVariables,GlobalParameters> interfaceGlobalParameters
    ("GlobalParameters",
     "Pointer to the GlobalParameters object",
     &ShowerVariables::_globalParameters, false, false, true, false, false);

  static ParVector<ShowerVariables,long> interfaceDecayInShower
    ("DecayInShower",
     "PDG codes of the particles to be decayed in the shower",
     &ShowerVariables::_inputparticlesDecayInShower, -1, 0l, -10000000l, 10000000l,
     false, false, Interface::limited);

  static Switch<ShowerVariables,unsigned int> interfaceFinalFinalConditions
    ("FinalFinalConditions",
     "The initial conditions for the shower of a final-final colour connection",
     &ShowerVariables::_finalFinalConditions, 0, false, false);
  static SwitchOption interfaceFinalFinalConditionsSymmetric
    (interfaceFinalFinalConditions,
     "Symmetric",
     "The symmetric choice",
     0);
  static SwitchOption interfaceFinalFinalConditionsColoured
    (interfaceFinalFinalConditions,
     "Coloured",
     "Maximal radiation from the coloured particle",
     1);
  static SwitchOption interfaceFinalFinalConditionsAntiColoured
    (interfaceFinalFinalConditions,
     "AntiColoured",
     "Maximal emission from the anticoloured particle",
     2);
  static SwitchOption interfaceFinalFinalConditionsRandom
    (interfaceFinalFinalConditions,
     "Random",
     "Randomly selected maximal emission from one of the particles",
     3);

  static Switch<ShowerVariables,unsigned int> interfaceInitialFinalDecayConditions
    ("InitialFinalDecayConditions",
     "The initial conditions for the shower of an initial-final"
     " decay colour connection.",
     &ShowerVariables::_initialFinalDecayConditions, 0, false, false);
  static SwitchOption interfaceInitialFinalDecayConditionsSymmetric
    (interfaceInitialFinalDecayConditions,
     "Symmetric",
     "The symmetric choice",
     0);
  static SwitchOption interfaceInitialFinalDecayConditionsMaximal
    (interfaceInitialFinalDecayConditions,
     "Maximal",
     "Maximal radiation from the decay product",
     1);
  static SwitchOption interfaceInitialFinalDecayConditionsSmooth
    (interfaceInitialFinalDecayConditions,
     "Smooth",
     "Smooth matching in the soft limit",
     2);

  static Switch<ShowerVariables,bool> interfaceUseMEForT2
    ("UseMEForT2",
     "Use the matrix element correction, if available to fill the T2"
     " region for the decay shower and don't fill using the shower",
     &ShowerVariables::_useMEForT2, false, false, false);
  static SwitchOption interfaceUseMEForT2Shower
    (interfaceUseMEForT2,
     "Shower",
     "Use the shower to fill the T2 region",
     false);
  static SwitchOption interfaceUseMEForT2ME
    (interfaceUseMEForT2,
     "ME",
     "Use the Matrix element to fill the T2 region",
     true);

}

