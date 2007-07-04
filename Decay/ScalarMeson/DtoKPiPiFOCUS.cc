// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiFOCUS class.
//

#include "DtoKPiPiFOCUS.h"
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

DtoKPiPiFOCUS::DtoKPiPiFOCUS() {
  // the option
  _imodel=0;
  // amplitudes and phases for the isobar model
  _a1NR     = 1.47     ; _phi1NR     = 325.;
  _a1kappa  = 1.28*GeV2; _phi1kappa  = 199.;
  _a1K892   = 1.00     ; _phi1K892   =   0.;
  _a1K1410  = 0.12     ; _phi1K1410  = 350.;
  _a1K1680  = 0.36     ; _phi1K1680  =   3.;
  _a1K14300 = 1.13*GeV2; _phi1K14300 =  36.;
  _a1K14302 = 0.17/GeV2; _phi1K14302 = 319.;
  // amplitudes and phases for the K-matrix model
  _a2K892   = 1.00      ;_phi2K892   =   0.;
  _a2K1410  = 0.188     ;_phi2K1410  = 293.;
  _a2K1680  = 0.373     ;_phi2K1680  =   1.;
  _a2K14302 = 0.169/GeV2;_phi2K14302 = 296.;
  // Masses and Widths for the intermediate resonances
  _localparameters=true;
  _mkappa  =  856. *MeV; _wkappa  = 464. *MeV;
  _mK892   =  896.0*MeV; _wK892   =  50.3*MeV;
  _mK1410  = 1414. *MeV; _wK1410  = 232. *MeV;
  _mK1680  = 1717. *MeV; _wK1680  = 322. *MeV;
  _mK14300 = 1461. *MeV; _wK14300 = 177. *MeV;
  _mK14302 = 1432.4*MeV; _wK14302 = 109. *MeV;
  // radial sizes
  _rD0  = 5.0/GeV;
  _rres = 1.5/GeV;
  // parameters for the K-matrix
  _g1 = 0.31072*GeV; _g2 = -0.02323*GeV;
  _c110 = 0.79299; _c111 =-0.15099 ; _c112 = 0.00811   ;
  _c120 = 0.15040; _c121 =-0.038266; _c122 = 0.0022596 ;
  _c220 = 0.17054; _c221 =-0.0219  ; _c222 = 0.00085655;
  _d110 =-0.22147; _d111 = 0.026637; _d112 =-0.00092057;
  _s0half      = 0.23*GeV2;
  _s0threehalf = 0.27*GeV2;
  _s1 = 1.7919*GeV2;
  // parameters for the production vector
  _beta  = 3.389*GeV;
  _theta = 286.;
  _c10 =  1.655; _c20 = 17.182; _c30 =  0.734;
  _c11 =  0.780; _c21 =  0.   ; _c31 =  0.   ;
  _c12 = -0.954; _c22 =  0.   ; _c32 =  0.   ;
  _gamma1 = 304.;
  _gamma2 = 126.;
  _gamma3 = 211.;
  // masses of the stable particles (reset in doinit)
  _mK=0.*MeV;
  _mpi=0.*MeV;
  _metap=0.*MeV;
  // integration parameters
  _maxwgt = 1.;
}

void DtoKPiPiFOCUS::persistentOutput(PersistentOStream & os) const {
  os << _imodel << _a1NR << _phi1NR << ounit(_a1kappa,GeV2) << _phi1kappa 
     << _a1K892 << _phi1K892 << _a1K1410 << _phi1K1410 << _a1K1680 
     << _phi1K1680 << ounit(_a1K14300,GeV2) << _phi1K14300 
     << ounit(_a1K14302,1./GeV2) << _phi1K14302 << _a2K892 << _phi2K892 
     << _a2K1410 << _phi2K1410 << _a2K1680 << _phi2K1680 
     << ounit(_a2K14302,1./GeV2) << _phi2K14302 << ounit(_g1,GeV) << ounit(_g2,GeV) 
     << _c110 << _c111 << _c112 << _c120 << _c121 << _c122 << _c220 << _c221 
     << _c222 << _d110 << _d111 << _d112 << ounit(_s0half,GeV2) << ounit(_s0threehalf,GeV2)
     << ounit(_s1,GeV2) << ounit(_beta,GeV) << _theta << _c10 << _c20 << _c30 
     << _c11 << _c21 << _c31 << _c12 << _c22 << _c32 << _gamma1 << _gamma2 
     << _gamma3 << ounit(_mK,GeV) << ounit(_mpi,GeV) << ounit(_metap,GeV) 
     << _cNR << ounit(_ckappa,GeV2) << _cK892 << _cK1410 << _cK1680 << ounit(_cK14300,GeV2)
     << ounit(_cK14302,1./GeV2) << _localparameters << ounit(_mkappa,GeV) 
     << ounit(_wkappa,GeV) << ounit(_mK892,GeV) << ounit(_wK892,GeV) 
     << ounit(_mK1410,GeV) << ounit(_wK1410,GeV) << ounit(_mK1680,GeV) 
     << ounit(_wK1680,GeV) << ounit(_mK14300,GeV) << ounit(_wK14300,GeV) 
     << ounit(_mK14302,GeV) << ounit(_wK14302,GeV) << ounit(_rD0,1./GeV) 
     << ounit(_rres,1./GeV) << _maxwgt << _weights;
}

void DtoKPiPiFOCUS::persistentInput(PersistentIStream & is, int) {
  is >> _imodel >> _a1NR >> _phi1NR >> iunit(_a1kappa,GeV2) >> _phi1kappa 
     >> _a1K892 >> _phi1K892 >> _a1K1410 >> _phi1K1410 >> _a1K1680 
     >> _phi1K1680 >> iunit(_a1K14300,GeV2) >> _phi1K14300 
     >> iunit(_a1K14302,1./GeV2) >> _phi1K14302 >> _a2K892 >> _phi2K892 
     >> _a2K1410 >> _phi2K1410 >> _a2K1680 >> _phi2K1680 
     >> iunit(_a2K14302,1./GeV2) >> _phi2K14302 >> iunit(_g1,GeV) >> iunit(_g2,GeV) 
     >> _c110 >> _c111 >> _c112 >> _c120 >> _c121 >> _c122 >> _c220 >> _c221 
     >> _c222 >> _d110 >> _d111 >> _d112 >> iunit(_s0half,GeV2) >> iunit(_s0threehalf,GeV2)
     >> iunit(_s1,GeV2) >> iunit(_beta,GeV) >> _theta >> _c10 >> _c20 >> _c30 
     >> _c11 >> _c21 >> _c31 >> _c12 >> _c22 >> _c32 >> _gamma1 >> _gamma2 
     >> _gamma3 >> iunit(_mK,GeV) >> iunit(_mpi,GeV) >> iunit(_metap,GeV) 
     >> _cNR >> iunit(_ckappa,GeV2) >> _cK892 >> _cK1410 >> _cK1680 >> iunit(_cK14300,GeV2)
     >> iunit(_cK14302,1./GeV2) >> _localparameters >> iunit(_mkappa,GeV) 
     >> iunit(_wkappa,GeV) >> iunit(_mK892,GeV) >> iunit(_wK892,GeV) 
     >> iunit(_mK1410,GeV) >> iunit(_wK1410,GeV) >> iunit(_mK1680,GeV) 
     >> iunit(_wK1680,GeV) >> iunit(_mK14300,GeV) >> iunit(_wK14300,GeV) 
     >> iunit(_mK14302,GeV) >> iunit(_wK14302,GeV) >> iunit(_rD0,1./GeV) 
     >> iunit(_rres,1./GeV) >> _maxwgt >> _weights;
}

ClassDescription<DtoKPiPiFOCUS> DtoKPiPiFOCUS::initDtoKPiPiFOCUS;
// Definition of the static class description member.

void DtoKPiPiFOCUS::Init() {

  static ClassDocumentation<DtoKPiPiFOCUS> documentation
    ("The DtoKPiPiFOCUS class implements the model of FOCUS for D+ -> K-pi+pi+",
     "The FOCUS fit of \\cite{Collaboration:2007se} was used for the decay "
     "$D^+\\to K^-\\pi^+\\pi^+$",
     "\\bibitem{Collaboration:2007se} FOCUS~Collaboration and M.~R.~Pennington,"
     "arXiv:0705.2248 [hep-ex].");

  static Switch<DtoKPiPiFOCUS,unsigned int> interfaceModel
    ("Model",
     "Which model to use",
     &DtoKPiPiFOCUS::_imodel, 1, false, false);
  static SwitchOption interfaceModelIsobar
    (interfaceModel,
     "Isobar",
     "Use the isobar model",
     0);
  static SwitchOption interfaceModelKMatrix
    (interfaceModel,
     "KMatrix",
     "Use the K-matrix model",
     1);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarNonResonantAmplitude
    ("IsobarNonResonantAmplitude",
     "The amplitude of the non-resonant component for the isobar model",
     &DtoKPiPiFOCUS::_a1NR, 1.47, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarNonResonantPhase
    ("IsobarNonResonantPhase",
     "The phase of the non-resonant component for the isobar model",
     &DtoKPiPiFOCUS::_phi1NR, 325., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy2> interfaceIsobarKappaAmplitude
    ("IsobarKappaAmplitude",
     "The amplitude for the kappa component for the isobar model",
     &DtoKPiPiFOCUS::_a1kappa, GeV2, 1.28*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarKappaPhase
    ("IsobarKappaPhase",
     "The phase for the kappa component for the isobar model",
     &DtoKPiPiFOCUS::_phi1kappa, 199., 0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK892Amplitude
    ("IsobarK892Amplitude",
     "The amplitude for the K*(892) component for the isobar model",
     &DtoKPiPiFOCUS::_a1K892, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK892Phase
    ("IsobarK892Phase",
     "The phase for the K*(892) component for the isobar model",
     &DtoKPiPiFOCUS::_phi1K892, 0.0, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK1410Amplitude
    ("IsobarK1410Amplitude",
     "The amplitude for the K*(1410) component for the isobar model",
     &DtoKPiPiFOCUS::_a1K1410, 0.12, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK1410Phase
    ("IsobarK1410Phase",
     "The phase for the K*(1410) component for the isobar model",
     &DtoKPiPiFOCUS::_phi1K1410, 350., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK1680Amplitude
    ("IsobarK1680Amplitude",
     "The amplitude for the K*(1680) component for the isobar model",
     &DtoKPiPiFOCUS::_a1K1680, 0.36, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK1680Phase
    ("IsobarK1680Phase",
     "The phase for the K*(1680) component for the isobar model",
     &DtoKPiPiFOCUS::_phi1K1680, 3., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy2> interfaceIsobarK14300Amplitude
    ("IsobarK14300Amplitude",
     "The amplitude for the K(1430)0 component for the isobar model",
     &DtoKPiPiFOCUS::_a1kappa, GeV2, 1.13*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK14300Phase
    ("IsobarK14300Phase",
     "The phase for the K(1430)0 component for the isobar model",
     &DtoKPiPiFOCUS::_phi1kappa, 36., 0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,InvEnergy2> interfaceIsobarK14302Amplitude
    ("IsobarK14302Amplitude",
     "The amplitude ofr the K(1430)2 component for the isobar model",
     &DtoKPiPiFOCUS::_a1K14302, 1./GeV2, 0.17/GeV2, 0./GeV2, 10.0/GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceIsobarK14302Phase
    ("IsobarK14302Phase",
     "The phase ofr the K(1430)2 component for the isobar model",
     &DtoKPiPiFOCUS::_phi1K14302, 319., 0., 360.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK892Amplitude
    ("KMatrixK892Amplitude",
     "The amplitude for the K*(892) component for the isobar model",
     &DtoKPiPiFOCUS::_a2K892, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK892Phase
    ("KMatrixK892Phase",
     "The phase for the K*(892) component for the isobar model",
     &DtoKPiPiFOCUS::_phi2K892, 0.0, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK1410Amplitude
    ("KMatrixK1410Amplitude",
     "The amplitude for the K*(1410) component for the isobar model",
     &DtoKPiPiFOCUS::_a2K1410, 0.188, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK1410Phase
    ("KMatrixK1410Phase",
     "The phase for the K*(1410) component for the isobar model",
     &DtoKPiPiFOCUS::_phi2K1410, 293., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK1680Amplitude
    ("KMatrixK1680Amplitude",
     "The amplitude for the K*(1680) component for the isobar model",
     &DtoKPiPiFOCUS::_a2K1680, 0.373, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK1680Phase
    ("KMatrixK1680Phase",
     "The phase for the K*(1680) component for the isobar model",
     &DtoKPiPiFOCUS::_phi2K1680, 1., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,InvEnergy2> interfaceKMatrixK14302Amplitude
    ("KMatrixK14302Amplitude",
     "The amplitude ofr the K(1430)2 component for the isobar model",
     &DtoKPiPiFOCUS::_a2K14302, 1./GeV2, 0.169/GeV2, 0./GeV2, 10.0/GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceKMatrixK14302Phase
    ("KMatrixK14302Phase",
     "The phase ofr the K(1430)2 component for the isobar model",
     &DtoKPiPiFOCUS::_phi2K14302, 319., 0., 296.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceg1
    ("g1",
     "The g1 coupling for the K-matrix",
     &DtoKPiPiFOCUS::_g1, GeV, 0.31072*GeV, -10.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceg2
    ("g2",
     "The g2 coupling for the K-matrix",
     &DtoKPiPiFOCUS::_g2, GeV, -0.02323*GeV, -10.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec110
    ("c110",
     "The c110 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c110, 0.79299, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec111
    ("c111",
     "The c111 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c111, -0.15099, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec112
    ("c112",
     "The c112 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c112, 0.00811, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec120
    ("c120",
     "The c120 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c120, 0.15040, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec121
    ("c121",
     "The c121 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c121, -0.038266, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec122
    ("c122",
     "The c122 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c122, 0.0022596, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec220
    ("c220",
     "The c220 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c220, 0.17054, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec221
    ("c221",
     "The c221 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c221, -0.0219, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec222
    ("c222",
     "The c222 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_c222, 0.00085655, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaced110
    ("d110",
     "The d110 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_d110, -0.22147, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaced111
    ("d111",
     "The d111 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_d111, 0.026637, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaced112
    ("d112",
     "The d112 coefficent for the K-matrix",
     &DtoKPiPiFOCUS::_d112, -0.00092057, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy2> interfaces0half
    ("s0half",
     "The position of the Adler zero in the I=1/2 channel",
     &DtoKPiPiFOCUS::_s0half, GeV2, 0.23*GeV2, -10.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy2> interfaces0threehalf
    ("s0threehalf",
     "The position of the Adler zero in the I=3/2 channel",
     &DtoKPiPiFOCUS::_s0threehalf, GeV2, 0.27*GeV2, -10.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy2> interfaces1
    ("s1",
     "The position of the pole in the I=1/2 channel",
     &DtoKPiPiFOCUS::_s1, GeV2, 1.7919*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceBeta
    ("Beta",
     "The magnitude of the production coupling for the pole in the K-matrix",
     &DtoKPiPiFOCUS::_beta, GeV, 3.389*GeV, -10.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceTheta
    ("Theta",
     "The phase of the production coupling for the pole in the K-matrix",
     &DtoKPiPiFOCUS::_theta, 286., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec10
    ("c10",
     "The c10 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c10, 1.655, -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec11
    ("c11",
     "The c11 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c11, 0.780, -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec12
    ("c12",
     "The c12 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c12, -0.954, -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec20
    ("c20",
     "The c20 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c20, 17.182, -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec21
    ("c21",
     "The c21 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c21, 0., -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec22
    ("c22",
     "The c22 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c22, 0., -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec30
    ("c30",
     "The c30 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c30, 0.734, -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec31
    ("c31",
     "The c31 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c31, 0., -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfacec32
    ("c32",
     "The c32 coefficient for the production vector",
     &DtoKPiPiFOCUS::_c32, 0., -50., 50.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceGamma1
    ("Gamma1",
     "The phase gamma1 for the production vector",
     &DtoKPiPiFOCUS::_gamma1, 304., 0.0, 360.0,
     true, true, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceGamma2
    ("Gamma2",
     "The phase gamma2 for the production vector",
     &DtoKPiPiFOCUS::_gamma2, 126., 0.0, 360.0,
     true, true, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceGamma3
    ("Gamma3",
     "The phase gamma3 for the production vector",
     &DtoKPiPiFOCUS::_gamma3, 211., 0.0, 360.0,
     true, true, Interface::limited);

  static Switch<DtoKPiPiFOCUS,bool> interfaceLocalParameters
    ("LocalParameters",
     "Whether to use local values for the masses and widths or"
     " those from the ParticleData objects",
     &DtoKPiPiFOCUS::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the ParticleData objects",
     false);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKappaMass
    ("KappaMass",
     "The mass of the kappa meson",
     &DtoKPiPiFOCUS::_mkappa, MeV, 856.0  *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKstar892Mass
    ("Kstar892Mass",
     "The mass of the K*(892) meson",
     &DtoKPiPiFOCUS::_mK892, MeV, 896.0  *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKstar1410Mass
    ("Kstar1410Mass",
     "The mass of the K*(1410) meson",
     &DtoKPiPiFOCUS::_mK1410, MeV, 1414   *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKstar1680Mass
    ("Kstar1680Mass",
     "The mass of the K*(1680) meson",
     &DtoKPiPiFOCUS::_mK1680, MeV, 1717   *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceK_01430Mass
    ("K_01430Mass",
     "The mass of the K_0(1430) meson",
     &DtoKPiPiFOCUS::_mK14300, MeV, 1461   *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceK_21430Mass
    ("K_21430Mass",
     "The mass of the K_2(1430) meson",
     &DtoKPiPiFOCUS::_mK14302, MeV, 1432.4 *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKappaWidth
    ("KappaWidth",
     "The width of the kappa meson",
     &DtoKPiPiFOCUS::_wkappa, MeV, 464.0  *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKstar892Width
    ("Kstar892Width",
     "The width of the K*(892) meson",
     &DtoKPiPiFOCUS::_wK892, MeV, 50.3  *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKstar1410Width
    ("Kstar1410Width",
     "The width of the K*(1410) meson",
     &DtoKPiPiFOCUS::_wK1410, MeV, 232.   *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceKstar1680Width
    ("Kstar1680Width",
     "The width of the K*(1680) meson",
     &DtoKPiPiFOCUS::_wK1680, MeV, 322.   *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceK_01430Width
    ("K_01430Width",
     "The width of the K_0(1430) meson",
     &DtoKPiPiFOCUS::_wK14300, MeV, 177.   *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,Energy> interfaceK_21430Width
    ("K_21430Width",
     "The width of the K_2(1430) meson",
     &DtoKPiPiFOCUS::_wK14302, MeV, 109. *MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,InvEnergy> interfaceDRadius
    ("DRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &DtoKPiPiFOCUS::_rD0, 1./GeV, 5./GeV, 0./GeV, 10./GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,InvEnergy> interfaceResonanceRadius
    ("ResonanceRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the"
     "intermediate resonances",
     &DtoKPiPiFOCUS::_rres, 1./GeV, 1.5/GeV, 0./GeV, 10./GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiFOCUS,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the unweighting of the decay",
     &DtoKPiPiFOCUS::_maxwgt, 1.0, 0.0, 1000.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiFOCUS,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &DtoKPiPiFOCUS::_weights, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
}

void DtoKPiPiFOCUS::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // complex amplitudes
  double fact = Constants::pi/180.;
  if(_imodel==0) {
    _cNR     = _a1NR    *Complex(cos(fact*_phi1NR    ),sin(fact*_phi1NR    ));
    _ckappa  = _a1kappa *Complex(cos(fact*_phi1kappa ),sin(fact*_phi1kappa ));
    _cK892   = _a1K892  *Complex(cos(fact*_phi1K892  ),sin(fact*_phi1K892  ));
    _cK1410  = _a1K1410 *Complex(cos(fact*_phi1K1410 ),sin(fact*_phi1K1410 ));
    _cK1680  = _a1K1680 *Complex(cos(fact*_phi1K1680 ),sin(fact*_phi1K1680 ));
    _cK14300 = _a1K14300*Complex(cos(fact*_phi1K14300),sin(fact*_phi1K14300));
    _cK14302 = _a1K14302*Complex(cos(fact*_phi1K14302),sin(fact*_phi1K14302));
  }
  else {
    _cNR     = 0.;
    _ckappa  = 0.*GeV2;
    _cK892   = _a2K892  *Complex(cos(fact*_phi2K892  ),sin(fact*_phi2K892  ));
    _cK1410  = _a2K1410 *Complex(cos(fact*_phi2K1410 ),sin(fact*_phi2K1410 ));
    _cK1680  = _a2K1680 *Complex(cos(fact*_phi2K1680 ),sin(fact*_phi2K1680 ));
    _cK14300 = 0.*GeV2;
    _cK14302 = _a2K14302*Complex(cos(fact*_phi2K14302),sin(fact*_phi2K14302));
  }
  // intermediate resonances
  tPDPtr kappa  = getParticleData(-9000311);
  tPDPtr k14300 = getParticleData(ParticleID::Kstar_0bar0);
  tPDPtr k14302 = getParticleData(ParticleID::Kstar_2bar0);
  tPDPtr k1410  = getParticleData(-100313);
  tPDPtr k1680  = getParticleData(-30313);
  tPDPtr k892   = getParticleData(ParticleID::Kstarbar0);
  // external particles
  PDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::Dplus);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piplus);
  // the integration channels
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
  if(k1410) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1410,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1410,0,0., 1,3);
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
  // add the mode
  if(ix!=_weights.size()) _weights=vector<double>(ix,1./float(ix));
  addMode(mode,_maxwgt,_weights);
  // masses of the intermediate resonances
  if(_localparameters) {
    resetIntermediate(kappa  ,_mkappa ,_wkappa );
    resetIntermediate(k892   ,_mK892  ,_wK892  );
    resetIntermediate(k14300 ,_mK14300,_wK14300);
    resetIntermediate(k14302 ,_mK14302,_wK14302);
    resetIntermediate(k1410  ,_mK1410 ,_wK1410 );
    resetIntermediate(k1680  ,_mK1680 ,_wK1680 );
  }
  else {
    _mkappa  = kappa  ->mass();
    _mK892   = k892   ->mass();
    _mK14300 = k14300 ->mass();
    _mK14302 = k14302 ->mass();
    _mK1410  = k1410  ->mass();
    _mK1680  = k1680  ->mass();
    _wkappa  = kappa  ->width();
    _wK892   = k892   ->width();
    _wK14300 = k14300 ->width();
    _wK14302 = k14302 ->width();
    _wK1410  = k1410  ->width();
    _wK1680  = k1680  ->width();
  }
  // masses of the stable particles
  _mK    = getParticleData(ParticleID::Kminus  )->mass();
  _mpi   = getParticleData(ParticleID::piplus  )->mass();
  _metap = getParticleData(ParticleID::etaprime)->mass();
  cerr << "set limits y 0 2\n";
  for(Energy2 s = sqr(_mK+_mpi);s<sqr(1.75*GeV);s+=0.01*GeV2) {
    Complex test = F(s);
    cerr << sqrt(s/GeV2) << " " << sqrt(real(test*conj(test))) << "\n";
  }
  cerr << "join \n new frame \n";
  cerr << "set limits y -200 200\n";
  for(Energy2 s = sqr(_mK+_mpi);s<sqr(1.75*GeV);s+=0.01*GeV2) {
    Complex test = F(s);
    cerr << sqrt(s/GeV2) << " " << atan2(test.imag(),test.real())/Constants::pi*180 << "\n";
  }
  cerr << "join\n";
}

int DtoKPiPiFOCUS::modeNumber(bool & cc,const DecayMode & dm) const {
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

double DtoKPiPiFOCUS::me2(bool vertex, const int, const Particle & inpart,
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
  // masses of the particles
  Energy mD  = inpart.mass();
  Energy mA  = decay[0]->mass();
  Energy mB  = decay[1]->mass();
  Energy mC  = decay[2]->mass();
  // resonances and angles
  double ct1,ct2;
  Energy2 E1,E2;
  Lorentz5Momentum pres1=decay[0]->momentum()+decay[1]->momentum();
  Lorentz5Momentum pres2=decay[0]->momentum()+decay[2]->momentum();
  pres1.rescaleMass();
  pres2.rescaleMass();
  decayAngle(inpart.momentum(),pres1,decay[1]->momentum(),ct1,E1);
  decayAngle(inpart.momentum(),pres2,decay[2]->momentum(),ct2,E2);
  Complex amp(0.);
  // isobar model
  if(_imodel==0) {
    // normalised to individual components
//     amp = 
//       _cNR*1.033
//       +_ckappa*0.458/GeV2*amplitude(0,mD,mA,mB,mC,pres1.mass(),_mkappa ,_wkappa ,E1,ct1)
//       +_ckappa*0.458/GeV2*amplitude(0,mD,mA,mC,mB,pres2.mass(),_mkappa ,_wkappa ,E2,ct2)
//       +_cK892            *amplitude(1,mD,mA,mB,mC,pres1.mass(),_mK892  ,_wK892  ,E1,ct1)
//       +_cK892            *amplitude(1,mD,mA,mC,mB,pres2.mass(),_mK892  ,_wK892  ,E2,ct2)
//       +_cK1410*3.352     *amplitude(1,mD,mA,mB,mC,pres1.mass(),_mK1410 ,_wK1410 ,E1,ct1)
//       +_cK1410*3.352     *amplitude(1,mD,mA,mC,mB,pres2.mass(),_mK1410 ,_wK1410 ,E2,ct2)
//       +_cK1680*7.773      *amplitude(1,mD,mA,mB,mC,pres1.mass(),_mK1680 ,_wK1680 ,E1,ct1)
//       +_cK1680*7.773      *amplitude(1,mD,mA,mC,mB,pres2.mass(),_mK1680 ,_wK1680 ,E2,ct2)
//       +_cK14300*0.379/GeV2*amplitude(0,mD,mA,mB,mC,pres1.mass(),_mK14300,_wK14300,E1,ct1)
//       +_cK14300*0.379/GeV2*amplitude(0,mD,mA,mC,mB,pres2.mass(),_mK14300,_wK14300,E2,ct2)
//       +_cK14302*33.32*GeV2*amplitude(2,mD,mA,mB,mC,pres1.mass(),_mK14302,_wK14302,E1,ct1)
//       +_cK14302*33.32*GeV2*amplitude(2,mD,mA,mC,mB,pres2.mass(),_mK14302,_wK14302,E2,ct2)
//       ;
// normalised to sum
    amp = 
      _cNR*1.067
      +Complex(_ckappa *0.407/GeV2)*amplitude(0,mD,mA,mB,mC,pres1.mass(),_mkappa ,_wkappa ,E1,ct1)
      +Complex(_ckappa *0.407/GeV2)*amplitude(0,mD,mA,mC,mB,pres2.mass(),_mkappa ,_wkappa ,E2,ct2)
      +_cK892             *amplitude(1,mD,mA,mB,mC,pres1.mass(),_mK892  ,_wK892  ,E1,ct1)
      +_cK892             *amplitude(1,mD,mA,mC,mB,pres2.mass(),_mK892  ,_wK892  ,E2,ct2)
      +_cK1410 *3.036     *amplitude(1,mD,mA,mB,mC,pres1.mass(),_mK1410 ,_wK1410 ,E1,ct1)
      +_cK1410 *3.036     *amplitude(1,mD,mA,mC,mB,pres2.mass(),_mK1410 ,_wK1410 ,E2,ct2)
      +_cK1680 *6.413     *amplitude(1,mD,mA,mB,mC,pres1.mass(),_mK1680 ,_wK1680 ,E1,ct1)
      +_cK1680 *6.413     *amplitude(1,mD,mA,mC,mB,pres2.mass(),_mK1680 ,_wK1680 ,E2,ct2)
      +Complex(_cK14300*0.351/GeV2)*amplitude(0,mD,mA,mB,mC,pres1.mass(),_mK14300,_wK14300,E1,ct1)
      +Complex(_cK14300*0.351/GeV2)*amplitude(0,mD,mA,mC,mB,pres2.mass(),_mK14300,_wK14300,E2,ct2)
      +Complex(_cK14302*31.20*GeV2)*amplitude(2,mD,mA,mB,mC,pres1.mass(),_mK14302,_wK14302,E1,ct1)
      +Complex(_cK14302*31.20*GeV2)*amplitude(2,mD,mA,mC,mB,pres2.mass(),_mK14302,_wK14302,E2,ct2)
      ;
  }
  // K-matrix model
  else {
    cerr << "testing k-matrix not implmented in toKPiPiFOCUS::me2()\n";
    exit(1);
  }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newME(0,0,0,0)=amp;
  ME(newME);
  return real(amp*conj(amp));
}

void DtoKPiPiFOCUS::dataBaseOutput(ofstream &, bool) const {}

Complex DtoKPiPiFOCUS::F(Energy2 s) {
  double rhoa[2]={(1.-sqr(_mpi  +_mK)/s)*(1.-sqr(_mpi  -_mK)/s),
		  (1.-sqr(_metap+_mK)/s)*(1.-sqr(_metap-_mK)/s)};
  Complex rho[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(rhoa[ix]>=0.) rho[ix] =            sqrt( rhoa[ix]) ;
    else             rho[ix] = Complex(0.,sqrt(-rhoa[ix]));
  }
  Energy2 snorm=sqr(_mpi)+sqr(_mK);
  Complex Khalf[2][2];
  double stilde = s/snorm-1.;
  Khalf[0][0] = (s-_s0half)/snorm*(sqr(_g1)/(_s1-s)+_c110+_c111*stilde
				   +_c112*sqr(stilde));
  Khalf[1][1] = (s-_s0half)/snorm*(sqr(_g2)/(_s1-s)+_c220+_c221*stilde
				   +_c222*sqr(stilde));
  Khalf[0][1] = (s-_s0half)/snorm*(_g1*_g2 /(_s1-s)+_c120+_c121*stilde
				   +_c122*sqr(stilde));
  Khalf[1][0] = Khalf[0][1];
  Complex K3half = double((s-_s0threehalf)/snorm*(_d110+_d111*stilde+_d112*sqr(stilde)));
  // multiply by i rho and add unit matrix
  Complex ii(0.,1.);
  Khalf[0][0] = 1.-ii*Khalf[0][0]*rho[0];
  Khalf[0][1] =   -ii*Khalf[0][1]*rho[1];
  Khalf[1][0] =   -ii*Khalf[1][0]*rho[0];
  Khalf[1][1] = 1.-ii*Khalf[1][1]*rho[1];
  K3half = 1.-ii*K3half*rho[0];
  // calculate inverse
  Complex Khalfinv[2][2],det;
  det = Khalf[0][0]*Khalf[1][1]-Khalf[0][1]*Khalf[1][0];
  Khalfinv[0][0] =  Khalf[1][1]/det;
  Khalfinv[0][1] = -Khalf[0][1]/det;
  Khalfinv[1][0] = -Khalf[1][0]/det;
  Khalfinv[1][1] =  Khalf[0][0]/det;
  Complex K3halfinv=1./K3half;
  double shat = (s-2.*GeV2)/GeV2;
  double conv(Constants::pi/180.);
  // calculate the production piece
  Complex Phalf[2],P3half,etheta(cos(_theta*conv),sin(_theta*conv));
  Phalf[0] = _beta*_g1*etheta/(_s1-s)+(_c10+_c11*shat+_c12*sqr(shat))*
    Complex(cos(_gamma1*conv),sin(_gamma1*conv));
  Phalf[1] = _beta*_g2*etheta/(_s1-s)+(_c20+_c21*shat+_c22*sqr(shat))*
    Complex(cos(_gamma2*conv),sin(_gamma2*conv));
  P3half   =                          (_c30+_c31*shat+_c32*sqr(shat))*
    Complex(cos(_gamma3*conv),sin(_gamma3*conv));
  // calculate the answer
  Complex output(0.);
  for(unsigned int ix=0;ix<2;++ix) {
    output += Khalfinv[0][ix]*Phalf[ix];
  }
  output += K3halfinv*P3half;
  return output;
}
