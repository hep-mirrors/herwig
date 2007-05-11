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

using namespace Herwig;

DtoKPiPiE791::DtoKPiPiE791() {
  // which model to use
  _imodel=3;
  _aANR     = 1.00; _phiANR    =    0.;
  _aAK892   = 0.39; _phiAK892  =   54.;
  _aAK14300 = 0.58; _phiA14300 =   54.;
  _aAK14302 = 0.07; _phiA14302 =   33.;
  _aAK1680  = 0.19; _phiAK1680 =   66.;
  _aBNR     = 2.72; _phiBNR    = - 49.;
  _aBK892   = 1.00; _phiBK892  =    0.;
  _aBK14300 = 1.54; _phiB14300 =    6.;
  _aBK14302 = 0.21; _phiB14302 = -  3.;
  _aBK1680  = 0.56; _phiBK1680 =   36.;
  _aCNR     = 1.03; _phiCNR    = - 11.;
  _aCkappa  = 1.97; _phiCkappa =  187.;
  _aCK892   = 1.00; _phiCK892  =    0.;
  _aCK14300 = 1.01; _phiC14300 =   48.;
  _aCK14302 = 0.20; _phiC14302 = - 54.;
  _aCK1680  = 0.45; _phiCK1680 =   28.;
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

//   /**
//    *  Amplitudes and Phases for the different components
//    */
//   //@{
//   /**
//    *  Amplitude for the non-resonant component in model A
//    */
//   double _aANR;

//   /**
//    *  Phase for the non-resonant component in model A
//    */
//   double _phiANR;

//   /**
//    *  Amplitude for the \f$K^*(892)^-\f$ component in model A
//    */
//   double _aAK892;

//   /**
//    *  Phase for the \f$K^*(892)^-\f$ component in model A
//    */
//   double _phiAK892;

//   /**
//    *  Amplitude for the \f$K^*_0(1430)^-\f$ component in model A
//    */
//   Energy2 _aAK14300;

//   /**
//    *  Phase for the \f$K^*_0(1430)^-\f$ component in model A
//    */
//   double _phiA14300;

//   /**
//    *  Amplitude for the \f$K^*_2(1430)^-\f$ component in model A
//    */
//   InvEnergy2 _aAK14302;

//   /**
//    *  Phase for the \f$K^*_2(1430)^-\f$ component in model A
//    */
//   double _phiA14302;

//   /**
//    *  Amplitude for the \f$K^*(1680)^-\f$ component in model A
//    */
//   double _aAK1680;

//   /**
//    *  Phase for the \f$K^*(1680)^-\f$ component in model A
//    */
//   double _phiAK1680;

//   /**
//    *  Amplitude for the non-resonant component in model B
//    */
//   double _aBNR;

//   /**
//    *  Phase for the non-resonant component in model B
//    */
//   double _phiBNR;

//   /**
//    *  Amplitude for the \f$K^*(892)^-\f$ component in model B
//    */
//   double _aBK892;

//   /**
//    *  Phase for the \f$K^*(892)^-\f$ component in model B
//    */
//   double _phiBK892;

//   /**
//    *  Amplitude for the \f$K^*_0(1430)^-\f$ component in model B
//    */
//   Energy2 _aBK14300;

//   /**
//    *  Phase for the \f$K^*_0(1430)^-\f$ component in model B
//    */
//   double _phiB14300;

//   /**
//    *  Amplitude for the \f$K^*_2(1430)^-\f$ component in model B
//    */
//   InvEnergy2 _aBK14302;

//   /**
//    *  Phase for the \f$K^*_2(1430)^-\f$ component in model B
//    */
//   double _phiB14302;

//   /**
//    *  Amplitude for the \f$K^*(1680)^-\f$ component in model B
//    */
//   double _aBK1680;

//   /**
//    *  Phase for the \f$K^*(1680)^-\f$ component in model B
//    */
//   double _phiBK1680;

//   /**
//    *  Amplitude for the non-resonant component in model C
//    */
//   double _aCNR;

//   /**
//    *  Phase for the non-resonant component in model C
//    */
//   double _phiCNR;

//   /**
//    *  Amplitude for the \f$\kappa\f$ component in model C
//    */
//   Energy2 _aCkappa;

//   /**
//    *  Phase for the \f$\kappa\f$ component in model C
//    */
//   Energy2 _phiCkappa;

//   /**
//    *  Amplitude for the \f$K^*(892)^-\f$ component in model C
//    */
//   double _aCK892;

//   /**
//    *  Phase for the \f$K^*(892)^-\f$ component in model C
//    */
//   double _phiCK892;

//   /**
//    *  Amplitude for the \f$K^*_0(1430)^-\f$ component in model C
//    */
//   Energy2 _aCK14300;

//   /**
//    *  Phase for the \f$K^*_0(1430)^-\f$ component in model C
//    */
//   double _phiC14300;

//   /**
//    *  Amplitude for the \f$K^*_2(1430)^-\f$ component in model C
//    */
//   InvEnergy2 _aCK14302;

//   /**
//    *  Phase for the \f$K^*_2(1430)^-\f$ component in model C
//    */
//   double _phiC14302;

//   /**
//    *  Amplitude for the \f$K^*(1680)^-\f$ component in model C
//    */
//   double _aCK1680;

//   /**
//    *  Phase for the \f$K^*(1680)^-\f$ component in model C
//    */
//   double _phiCK1680;

//   /**
//    *  Complex amplitude for the \f$\kappa\f$ component in 
//    */
//   Energy2 _ckappa;

//   /**
//    *  Complex amplitude for the \f$K^*(892)^-\f$ component in 
//    */
//   double _cK892;

//   /**
//    *  Complex amplitude for the \f$K^*_0(1430)^-\f$ component in 
//    */
//   Energy2 _cK14300;

//   /**
//    *  Complex amplitude for the \f$K^*_2(1430)^-\f$ component in 
//    */
//   InvEnergy2 _cK14302;

//   /**
//    *  Complex amplitude for the \f$K^*(1680)^-\f$ component in 
//    */
//   double _cK1680;
//   //@}


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
     << _phiA14300 << _aAK14302 << _phiA14302 << _aAK1680 << _phiAK1680 
     << _aBNR << _phiBNR << _aBK892 << _phiBK892 << _aBK14300 << _phiB14300 
     << _aBK14302 << _phiB14302 << _aBK1680 << _phiBK1680 << _aCNR << _phiCNR 
     << _aCkappa << _phiCkappa << _aCK892 << _phiCK892 << _aCK14300 << _phiC14300 
     << _aCK14302 << _phiC14302 << _aCK1680 << _phiCK1680 << _ckappa << _cK892 
     << _cK14300 << _cK14302 << _cK1680 << _rD0A << _rresA << _rD0B << _rresB
     << _rD0C << _rresC << _localparameters 
     << _mkappa << _wkappa << _mK892 << _wK892 << _mK1430A << _wK1430A << _mK1430B
     << _wK1430B << _mK1430C << _wK1430C << _mK14302 << _wK14302 << _mK1680 << _wK1680;
}

void DtoKPiPiE791::persistentInput(PersistentIStream & is, int) {
  is >> _imodel >> _aANR >> _phiANR >> _aAK892 >> _phiAK892 >> _aAK14300 
     >> _phiA14300 >> _aAK14302 >> _phiA14302 >> _aAK1680 >> _phiAK1680 
     >> _aBNR >> _phiBNR >> _aBK892 >> _phiBK892 >> _aBK14300 >> _phiB14300 
     >> _aBK14302 >> _phiB14302 >> _aBK1680 >> _phiBK1680 >> _aCNR >> _phiCNR 
     >> _aCkappa >> _phiCkappa >> _aCK892 >> _phiCK892 >> _aCK14300 >> _phiC14300 
     >> _aCK14302 >> _phiC14302 >> _aCK1680 >> _phiCK1680 >> _ckappa >> _cK892 
     >> _cK14300 >> _cK14302 >> _cK1680 >> _rD0A >> _rresA >> _rD0B >> _rresB
     >> _rD0C >> _rresC >> _localparameters 
     >> _mkappa >> _wkappa >> _mK892 >> _wK892 >> _mK1430A >> _wK1430A >> _mK1430B
     >> _wK1430B >> _mK1430C >> _wK1430C >> _mK14302 >> _wK14302 >> _mK1680 >> _wK1680;
}

int DtoKPiPiE791::modeNumber(bool & cc,const DecayMode & dm) const {
  return -1;
}

double DtoKPiPiE791::me2(bool vertex, const int ichan,
			 const Particle & inpart,
			 const ParticleVector & decay) const {
  return 0.;
}

void DtoKPiPiE791::dataBaseOutput(ofstream & output,
					       bool header) const {
}

