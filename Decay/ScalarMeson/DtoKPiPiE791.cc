// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiE791 class.
//

#include "DtoKPiPiE791.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

DtoKPiPiE791::DtoKPiPiE791() {
  // which model to use
  _imodel=1;
  _aANR     = 1.00     ; _phiANR     =    0.;
  _aAK892   = 0.39     ; _phiAK892   =   54.;
  _aAK14300 = 0.58*GeV2; _phiAK14300 =   54.;
  _aAK14302 = 0.07/GeV2; _phiAK14302 =   33.;
  _aAK1680  = 0.19     ; _phiAK1680  =   66.;
  _aBNR     = 2.72     ; _phiBNR     = - 49.;
  _aBK892   = 1.00     ; _phiBK892   =    0.;
  _aBK14300 = 1.54*GeV2; _phiBK14300 =    6.;
  _aBK14302 = 0.21/GeV2; _phiBK14302 = -  3.;
  _aBK1680  = 0.56     ; _phiBK1680  =   36.;
  _aCNR     = 1.03     ; _phiCNR     = - 11.;
  _aCkappa  = 1.97*GeV2; _phiCkappa  =  187.;
  _aCK892   = 1.00     ; _phiCK892   =    0.;
  _aCK14300 = 1.01*GeV2; _phiCK14300 =   48.;
  _aCK14302 = 0.20/GeV2; _phiCK14302 = - 54.;
  _aCK1680  = 0.45     ; _phiCK1680  =   28.;
  // radii for the different models
 _rD0A = 5.0/GeV; _rresA = 1.5/GeV;
 _rD0B = 0.8/GeV; _rresB = 1.8/GeV;
 _rD0C = 5.0/GeV; _rresC = 1.6/GeV;
 // use local values for masses and widths
 _localparameters=true;
 // masses and widths
 _mkappa  =  797. *MeV; _wkappa  = 410. *MeV;
 _mK892   =  896.1*MeV; _wK892   =  50.7*MeV;
 _mK1430A = 1412. *MeV; _wK1430A = 294. *MeV;
 _mK1430B = 1416. *MeV; _wK1430B = 250. *MeV;
 _mK1430C = 1459. *MeV; _wK1430C = 175. *MeV;
 _mK14302 = 1432.4*MeV; _wK14302 = 109. *MeV;
 _mK1680  = 1717. *MeV; _wK1680  = 322. *MeV;
}

void DtoKPiPiE791::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // amplitudes for model A
  double fact = pi/180.;
  if(_imodel==1) {
    _cNR     = _aANR    *Complex(cos(_phiANR    *fact),sin(_phiANR    *fact));
    _ckappa  = 0.;
    _cK892   = _aAK892  *Complex(cos(_phiAK892  *fact),sin(_phiAK892  *fact));
    _cK14300 = _aAK14300*Complex(cos(_phiAK14300*fact),sin(_phiAK14300*fact));
    _cK14302 = _aAK14302*Complex(cos(_phiAK14302*fact),sin(_phiAK14302*fact));
    _cK1680  = _aAK1680 *Complex(cos(_phiAK1680 *fact),sin(_phiAK1680 *fact));
    _mK1430 = _mK1430A;
    _wK1430 = _wK1430A;
  }
  // amplitudes for model B
  else if(_imodel==2) {
    _cNR     = _aBNR    *Complex(cos(_phiBNR   *fact),sin(_phiBNR   *fact));
    _ckappa  = 0.;
    _cK892   = _aBK892  *Complex(cos(_phiBK892  *fact),sin(_phiBK892  *fact));
    _cK14300 = _aBK14300*Complex(cos(_phiBK14300*fact),sin(_phiBK14300*fact));
    _cK14302 = _aBK14302*Complex(cos(_phiBK14302*fact),sin(_phiBK14302*fact));
    _cK1680  = _aBK1680 *Complex(cos(_phiBK1680 *fact),sin(_phiBK1680 *fact));
    _mK1430 = _mK1430B;
    _wK1430 = _wK1430B;
  }
  // amplitudes for model C
  else if(_imodel==3) {
    _cNR     = _aCNR    *Complex(cos(_phiCNR    *fact),sin(_phiCNR    *fact));
    _ckappa  = _aCkappa *Complex(cos(_phiCkappa *fact),sin(_phiCkappa *fact));
    _cK892   = _aCK892  *Complex(cos(_phiCK892  *fact),sin(_phiCK892  *fact));
    _cK14300 = _aCK14300*Complex(cos(_phiCK14300*fact),sin(_phiCK14300*fact));
    _cK14302 = _aCK14302*Complex(cos(_phiCK14302*fact),sin(_phiCK14302*fact));
    _cK1680  = _aCK1680 *Complex(cos(_phiCK1680 *fact),sin(_phiCK1680 *fact));
    _mK1430 = _mK1430C;
    _wK1430 = _wK1430C;
  }
  // intermediate resonances
  tPDPtr kappa  = getParticleData(-9000311);
  tPDPtr k14300 = getParticleData(ParticleID::Kstar_0bar0);
  tPDPtr k14302 = getParticleData(ParticleID::Kstar_2bar0);
  tPDPtr k1680  = getParticleData(-30313);
  tPDPtr k892   = getParticleData(ParticleID::Kstarbar0);
  // external particles
  PDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::Dplus);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piplus);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  unsigned int ix=0;
  if(kappa) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(kappa,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(kappa,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k892) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k892,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k14300) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k14300,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k14300,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k14302) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k14302,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k14302,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k1680) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1680,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1680,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(ix!=_weights.size()) _weights=vector<double>(ix,1./float(ix));
  addMode(mode,_maxwgt,_weights);
  // reset the resonance parameters in the integration if needed
  if(_localparameters) {
    resetIntermediate(kappa  ,_mkappa ,_wkappa );
    resetIntermediate(k892   ,_mK892  ,_wK892  );
    resetIntermediate(k14300 ,_mK1430 ,_wK1430 );
    resetIntermediate(k14302 ,_mK14302,_wK14302);
    resetIntermediate(k1680  ,_mK1680 ,_wK1680 );
  }
  // get values from the ParticleData objects if needed
  else {
    _mkappa  = kappa  ->mass();
    _mK892   = k892   ->mass();
    _mK1430A = k14300 ->mass();
    _mK1430B = k14300 ->mass();
    _mK1430C = k14300 ->mass();
    _mK14302 = k14302 ->mass();
    _mK1680  = k1680  ->mass();
    _wkappa  = kappa  ->width();
    _wK892   = k892   ->width();
    _wK1430A = k14300 ->width();
    _wK1430B = k14300 ->width();
    _wK1430C = k14300 ->width();
    _wK14302 = k14302 ->width();
    _wK1680  = k1680  ->width();
  }
  cerr << "testing non-resonant " << _cNR << "\n";
  cerr << "testing kappa        " << _ckappa/GeV2<< "\n";
  cerr << "testing k*           " << _cK892<< "\n";
  cerr << "testing k0           " << _cK14300/GeV2<< "\n";
  cerr << "testing k2           " << _cK14302*GeV2<< "\n";
  cerr << "testing k1680        " << _cK1680<< "\n";
}

ClassDescription<DtoKPiPiE791> DtoKPiPiE791::initDtoKPiPiE791;
// Definition of the static class description member.

void DtoKPiPiE791::Init() {

  static ClassDocumentation<DtoKPiPiE791> documentation
    ("There is no documentation for the DtoKPiPiE791 class");


  static Switch<DtoKPiPiE791,unsigned int> interfaceModel
    ("Model",
     "Which of the different fits for E791 should be used",
     &DtoKPiPiE791::_imodel, 3, false, false);
  static SwitchOption interfaceModelModelA
    (interfaceModel,
     "ModelA",
     "Use Model A",
     1);
  static SwitchOption interfaceModelModelB
    (interfaceModel,
     "ModelB",
     "Use model B",
     2);
  static SwitchOption interfaceModelModelC
    (interfaceModel,
     "ModelC",
     "Use Model C",
     3);

  static Parameter<DtoKPiPiE791,double> interfaceModelAAmplitudeNR
    ("ModelAAmplitudeNR",
     "Amplitude for the non-resonant component in model A",
     &DtoKPiPiE791::_aANR, 1., 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAPhaseNR
    ("ModelAPhaseNR",
     "Phase for the non-resonant component in model A",
     &DtoKPiPiE791::_phiANR, 0., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAAmplitudeK892
    ("ModelAAmplitudeK892",
     "Amplitude for the K*(892) component in model A",
     &DtoKPiPiE791::_aAK892, 0.39, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAPhaseK892
    ("ModelAPhaseK892",
     "Phase for the K*(892) component in model A",
     &DtoKPiPiE791::_phiAK892, 54., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAAmplitudeK_01430
    ("ModelAAmplitudeK_01430",
     "Amplitude for the K_0*(1430) component in model A",
     &DtoKPiPiE791::_aAK14300, 0.58, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAPhaseK_01430
    ("ModelAPhaseK_01430",
     "Phase for the K_0*(1430) component in model A",
     &DtoKPiPiE791::_phiAK14300, 54., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAAmplitudeK_21430
    ("ModelAAmplitudeK_21430",
     "Amplitude for the K_2*(1430) component in model A",
     &DtoKPiPiE791::_aAK14302, 0.07, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAPhaseK_21430
    ("ModelAPhaseK_21430",
     "Phase for the K_2*(1430) component in model A",
     &DtoKPiPiE791::_phiAK14302, 33., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAAmplitudeK1680
    ("ModelAAmplitudeK1680",
     "Amplitude for the K*(1680) component in model A",
     &DtoKPiPiE791::_aAK1680, 0.19, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelAPhaseK1680
    ("ModelAPhaseK1680",
     "Phase for the K*(1680) component in model A",
     &DtoKPiPiE791::_phiAK1680, 66., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBAmplitudeNR
    ("ModelBAmplitudeNR",
     "Amplitude for the non-resonant component in model B",
     &DtoKPiPiE791::_aBNR, 2.72, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBPhaseNR
    ("ModelBPhaseNR",
     "Phase for the non-resonant component in model B",
     &DtoKPiPiE791::_phiBNR, -49., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBAmplitudeK892
    ("ModelBAmplitudeK892",
     "Amplitude for the K*(892) component in model B",
     &DtoKPiPiE791::_aBK892, 1., 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBPhaseK892
    ("ModelBPhaseK892",
     "Phase for the K*(892) component in model B",
     &DtoKPiPiE791::_phiBK892, 0., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBAmplitudeK_01430
    ("ModelBAmplitudeK_01430",
     "Amplitude for the K_0*(1430) component in model B",
     &DtoKPiPiE791::_aBK14300, 1.54, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBPhaseK_01430
    ("ModelBPhaseK_01430",
     "Phase for the K_0*(1430) component in model B",
     &DtoKPiPiE791::_phiBK14300, 6., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBAmplitudeK_21430
    ("ModelBAmplitudeK_21430",
     "Amplitude for the K_2*(1430) component in model B",
     &DtoKPiPiE791::_aBK14302, 0.21, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBPhaseK_21430
    ("ModelBPhaseK_21430",
     "Phase for the K_2*(1430) component in model B",
     &DtoKPiPiE791::_phiBK14302, -3., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBAmplitudeK1680
    ("ModelBAmplitudeK1680",
     "Amplitude for the K*(1680) component in model B",
     &DtoKPiPiE791::_aBK1680, 0.56, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelBPhaseK1680
    ("ModelBPhaseK1680",
     "Phase for the K*(1680) component in model B",
     &DtoKPiPiE791::_phiBK1680, 36., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCAmplitudeNR
    ("ModelCAmplitudeNR",
     "Amplitude for the non-resonant component in model C",
     &DtoKPiPiE791::_aCNR, 1.03, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCPhaseNR
    ("ModelCPhaseNR",
     "Phase for the non-resonant component in model C",
     &DtoKPiPiE791::_phiCNR, -11., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCAmplitudeKappa
    ("ModelCAmplitudeKappa",
     "Amplitude for the kappa component in model C",
     &DtoKPiPiE791::_aCkappa, 1.97, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCPhaseKappa
    ("ModelCPhaseKappa",
     "Phase for the kappa component in model C",
     &DtoKPiPiE791::_phiCkappa, 187, -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCAmplitudeK892
    ("ModelCAmplitudeK892",
     "Amplitude for the K*(892) component in model C",
     &DtoKPiPiE791::_aCK892, 1., 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCPhaseK892
    ("ModelCPhaseK892",
     "Phase for the K*(892) component in model C",
     &DtoKPiPiE791::_phiCK892, 0., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCAmplitudeK_01430
    ("ModelCAmplitudeK_01430",
     "Amplitude for the K_0*(1430) component in model C",
     &DtoKPiPiE791::_aCK14300, 1.01, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCPhaseK_01430
    ("ModelCPhaseK_01430",
     "Phase for the K_0*(1430) component in model C",
     &DtoKPiPiE791::_phiCK14300, 48., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCAmplitudeK_21430
    ("ModelCAmplitudeK_21430",
     "Amplitude for the K_2*(1430) component in model C",
     &DtoKPiPiE791::_aCK14302, 0.20, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCPhaseK_21430
    ("ModelCPhaseK_21430",
     "Phase for the K_2*(1430) component in model C",
     &DtoKPiPiE791::_phiCK14302, -54., -360.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCAmplitudeK1680
    ("ModelCAmplitudeK1680",
     "Amplitude for the K*(1680) component in model C",
     &DtoKPiPiE791::_aCK1680, 0.45, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiE791,double> interfaceModelCPhaseK1680
    ("ModelCPhaseK1680",
     "Phase for the K*(1680) component in model C",
     &DtoKPiPiE791::_phiCK1680, 28., -360.0, 360.0,
     false, false, Interface::limited);


//   /**
//    *  Parameters for the Blatt-Weisskopf form-factors
//    */
//   //@{
//   /**
//    *  Radial size for the \f$D^0\f$ for model A
//    */
//   InvEnergy _rD0A;

//   /**
//    *  Radial size for the light resonances for model A
//    */
//   InvEnergy _rresA;

//   /**
//    *  Radial size for the \f$D^0\f$ for model B
//    */
//   InvEnergy _rD0B;

//   /**
//    *  Radial size for the light resonances for model B
//    */
//   InvEnergy _rresB;

//   /**
//    *  Radial size for the \f$D^0\f$ for model C
//    */
//   InvEnergy _rD0C;

//   /**
//    *  Radial size for the light resonances for model C
//    */
//   InvEnergy _rresC;
//   //@}

//   /**
//    *  Masses and Widths of the resonances
//    */
//   //@{
//   /**
//    *  Use local values for the masses and widths
//    */
//   bool _localparameters;

//   /**
//    *  Mass of \f$\kappa\f$
//    */
//   Energy _mkappa;

//   /**
//    *  Width of \f$\kappa\f$
//    */
//   Energy _wkappa;

//   /**
//    *  Mass of \f$K^*(892)\f$
//    */
//   Energy _mK892;

//   /**
//    *  Width of \f$K^*(892)\f$
//    */
//   Energy _wK892;

//   /**
//    *  Mass of \f$K^*_0(1430)\f$ in model A
//    */
//   Energy _mK1430A;

//   /**
//    *  Width of \f$K^*_0(1430)\f$ in model A
//    */
//   Energy _wK1430A;

//   /**
//    *  Mass of \f$K^*_0(1430)\f$ in model B
//    */
//   Energy _mK1430B;

//   /**
//    *  Width of \f$K^*_0(1430)\f$ in model B
//    */
//   Energy _wK1430B;

//   /**
//    *  Mass of \f$K^*_0(1430)\f$ in model C
//    */
//   Energy _mK1430C;

//   /**
//    *  Width of \f$K^*_0(1430)\f$ in model C
//    */
//   Energy _wK1430C;

//   /**
//    *  Mass of \f$K^*_2(1430)\f$ 
//    */
//   Energy _mK14302;

//   /**
//    *  Width of \f$K^*_2(1430)\f$
//    */
//   Energy _wK14302;

//   /**
//    *  Mass of \f$K^*(1680)\f$
//    */
//   Energy _mK1680;

//   /**
//    *  Width of \f$K^*(1680)\f$
//    */
//   Energy _wK1680;
//   //@}
}


void DtoKPiPiE791::persistentOutput(PersistentOStream & os) const {
  os << _imodel << _aANR << _phiANR << _aAK892 << _phiAK892 << _aAK14300 
     << _phiAK14300 << _aAK14302 << _phiAK14302 << _aAK1680 << _phiAK1680 
     << _aBNR << _phiBNR << _aBK892 << _phiBK892 << _aBK14300 << _phiBK14300 
     << _aBK14302 << _phiBK14302 << _aBK1680 << _phiBK1680 << _aCNR << _phiCNR 
     << _aCkappa << _phiCkappa << _aCK892 << _phiCK892 << _aCK14300 << _phiCK14300 
     << _aCK14302 << _phiCK14302 << _aCK1680 << _phiCK1680 << _ckappa << _cNR << _cK892 
     << _cK14300 << _cK14302 << _cK1680 << _rD0A << _rresA << _rD0B << _rresB
     << _rD0C << _rresC << _localparameters << _mK1430 << _wK1430
     << _mkappa << _wkappa << _mK892 << _wK892 << _mK1430A << _wK1430A << _mK1430B
     << _wK1430B << _mK1430C << _wK1430C << _mK14302 << _wK14302 << _mK1680 << _wK1680
     << _maxwgt << _weights;
}

void DtoKPiPiE791::persistentInput(PersistentIStream & is, int) {
  is >> _imodel >> _aANR >> _phiANR >> _aAK892 >> _phiAK892 >> _aAK14300 
     >> _phiAK14300 >> _aAK14302 >> _phiAK14302 >> _aAK1680 >> _phiAK1680 
     >> _aBNR >> _phiBNR >> _aBK892 >> _phiBK892 >> _aBK14300 >> _phiBK14300 
     >> _aBK14302 >> _phiBK14302 >> _aBK1680 >> _phiBK1680 >> _aCNR >> _phiCNR 
     >> _aCkappa >> _phiCkappa >> _aCK892 >> _phiCK892 >> _aCK14300 >> _phiCK14300 
     >> _aCK14302 >> _phiCK14302 >> _aCK1680 >> _phiCK1680 >> _ckappa >> _cNR >> _cK892 
     >> _cK14300 >> _cK14302 >> _cK1680 >> _rD0A >> _rresA >> _rD0B >> _rresB
     >> _rD0C >> _rresC >> _localparameters >> _mK1430 >> _wK1430
     >> _mkappa >> _wkappa >> _mK892 >> _wK892 >> _mK1430A >> _wK1430A >> _mK1430B
     >> _wK1430B >> _mK1430C >> _wK1430C >> _mK14302 >> _wK14302 >> _mK1680 >> _wK1680
     >> _maxwgt >> _weights;
}

int DtoKPiPiE791::modeNumber(bool & cc,const DecayMode & dm) const {
  int id0(dm.parent()->id());
  // incoming particle must be D+/-
  if(abs(id0)!=ParticleID::Dplus) return -1;
  cc = id0==ParticleID::Dminus;
  // must be three decay products
  if(dm.products().size()!=3) return -1;
  ParticleMSet::const_iterator pit = dm.products().begin();
  int isign = cc ? -1 : 1;
  unsigned int npip(0),nkm(0);
  for( ;pit!=dm.products().end();++pit) {
    id0=isign*(**pit).id();
    if(     id0==ParticleID::piplus)  ++npip;
    else if(id0==ParticleID::Kminus)   ++nkm;
  }
  return (npip==2&&nkm==1) ? 0 : -1;
}

double DtoKPiPiE791::me2(bool vertex, const int ichan,
			    const Particle & inpart,
			    const ParticleVector & decay) const {
  useMe();
  // wavefunnction for the decaying particle
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // wavefunctions for the outgoing particles
  for(unsigned int ix=0;ix<3;++ix) {
    PPtr mytemp = decay[ix]; 
    ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }
  // compute the angles and momenta we need
  Energy p2[2],p3[2];
  double ctheta[2];
  for(unsigned int ix=0;ix<2;++ix) {
    Lorentz5Momentum p[3]={decay[0]->momentum(),decay[1]->momentum(),
			   decay[2]->momentum()};
    if(ix==1) swap(p[1],p[2]);
    Lorentz5Momentum pres=p[0]+p[1];
    Hep3Vector boost=-pres.boostVector();
    for(unsigned int iy=0;iy<3;++iy) p[iy].boost(boost);
    p2[ix]=p[1].vect().mag();
    p3[ix]=p[2].vect().mag();
    ctheta[ix]=-p[1].vect().cosTheta(p[2].vect());
    //cerr << "testing p " << p2[ix] << " " << p3[ix] << "\n";
  }
//   cerr << "testing angles A " << ctheta[0] << " " << ctheta[1] << "\n";
//   Lorentz5Momentum pres1(decay[0]->momentum()+decay[1]->momentum());
//   pres1.rescaleMass();
//   cerr << "testing angle B " << decayAngle(inpart.momentum(),pres1,decay[0]->momentum())
//        << "\n";
//   Lorentz5Momentum pres2(decay[0]->momentum()+decay[2]->momentum());
//   pres2.rescaleMass();
//   cerr << "testing angle C " << decayAngle(inpart.momentum(),pres2,decay[0]->momentum())
//        << "\n";



  // masses of the particles and intermediates
  Energy mA  = decay[0]->mass();
  Energy mB  = decay[1]->mass();
  Energy mC  = decay[2]->mass();
  Energy mD  = inpart.mass();
  Energy mAB = (decay[0]->momentum()+decay[1]->momentum()).m();
  Energy mAC = (decay[0]->momentum()+decay[2]->momentum()).m();
  Complex amp;
  if(ichan<0) {
    amp=
      +_cNR
      +_ckappa/GeV2 *amplitude(0,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],_mkappa ,_wkappa )
      +_ckappa/GeV2 *amplitude(0,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],_mkappa ,_wkappa )
      +_cK892       *amplitude(1,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],_mK892  ,_wK892  )
      -_cK892       *amplitude(1,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],_mK892  ,_wK892  )
      +_cK14300/GeV2/2.5*amplitude(0,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],_mK1430 ,_wK1430 )
      +_cK14300/GeV2/2.5*amplitude(0,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],_mK1430 ,_wK1430 )
      +_cK14302*GeV2*8.5*amplitude(2,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],_mK14302,_wK14302)
      +_cK14302*GeV2*8.5*amplitude(2,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],_mK14302,_wK14302)
      +_cK1680*5.9*amplitude(1,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],_mK1680 ,_wK1680 )
      -_cK1680*5.9*amplitude(1,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],_mK1680 ,_wK1680)
      ;
  }
  else if(ichan==0) {
    amp=_ckappa/GeV2 *amplitude(0,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],
				_mkappa ,_wkappa );
  }
  else if(ichan==1) {
    amp=_ckappa/GeV2 *amplitude(0,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],
				_mkappa ,_wkappa );
  }
  else if(ichan==2) {
    amp=_cK892       *amplitude(1,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],
				_mK892  ,_wK892  );
  }
  else if(ichan==3) {
    amp=_cK892       *amplitude(1,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],
				_mK892  ,_wK892  );
  }
  else if(ichan==4) {
    amp=_cK14300/GeV2*amplitude(0,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],
				_mK1430 ,_wK1430 );
  }
  else if(ichan==5) {
    amp=_cK14300/GeV2*amplitude(0,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],
				_mK1430 ,_wK1430 );
  }
  else if(ichan==6) {
    amp=_cK14302*GeV2*amplitude(2,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],
				_mK14302,_wK14302);
  }
  else if(ichan==7) {
    amp=_cK14302*GeV2*amplitude(2,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1],
				_mK14302,_wK14302);
  }
  else if(ichan==8) {
    amp=_cK1680      *amplitude(1,mD,mA,mB,mC,mAB,ctheta[0],p2[0],p3[0],
				_mK1680 ,_wK1680 );
  }
  else if(ichan==9) {
    amp=_cK1680      *amplitude(1,mD,mA,mC,mB,mAC,ctheta[1],p2[1],p3[1]
				,_mK1680 ,_wK1680);;
  }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newME(0,0,0,0)=amp;
  ME(newME);
  return real(amp*conj(amp));
}

void DtoKPiPiE791::dataBaseOutput(ofstream & output, bool header) const {
}
