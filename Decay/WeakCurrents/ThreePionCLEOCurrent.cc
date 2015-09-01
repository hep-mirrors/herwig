// -*- C++ -*-
//
// ThreePionCLEOCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionCLEOCurrent class.
//

#include "ThreePionCLEOCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"

namespace {
  inline Energy  timesGeV (double x) { return x * GeV; }
  inline Energy2 timesGeV2(double x) { return x * GeV2; }
}

using namespace Herwig;

DescribeClass<ThreePionCLEOCurrent,ThreeMesonCurrentBase>
describeHerwigThreePionCLEOCurrent("Herwig::ThreePionCLEOCurrent",
				   "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(ThreePionCLEOCurrent,Energy,Energy2)

ThreePionCLEOCurrent::ThreePionCLEOCurrent() {
  // local particle properties
  _localparameters=true;
  // rho masses and widths
  if(_rhomass.size()==0) {
    _rhomass.push_back(0.7743*GeV);_rhowidth.push_back(0.1491*GeV);
    _rhomass.push_back(1.370*GeV);_rhowidth.push_back(0.386*GeV);
  }
  // f_2 mass and width
  _f2mass=1.275*GeV;_f2width=0.185*GeV;
  // f_0(1370) mass and width
  _f0mass=1.186*GeV;_f0width=0.350*GeV;
  // sigma mass and width
  _sigmamass = 0.860*GeV;_sigmawidth =0.880*GeV;
  // a1 mass and width
  _a1mass = 1.331*GeV;_a1width=0.814*GeV;
  // parameters for the K K* contribution to the a_1 running width
  _mKstar = 894*MeV;
  _mK=496*MeV;
  _gammk=3.32;
  // pion decay constant
  _fpi = 130.7*MeV/sqrt(2.);
  // couplings and phases for the different channels
  // p-wave rho and rho prime
  using Constants::pi;
  if(_rhomagP.size()==0) {
    _rhomagP.push_back(1.)  ;_rhophaseP.push_back(0.);
    _rhomagP.push_back(0.12);_rhophaseP.push_back(0.99*pi);
  }
  // d-wave rho and rho prime
  if(_rhomagD.size()==0) {
    _rhomagD.push_back(0.37/GeV2);_rhophaseD.push_back(-0.15*pi);
    _rhomagD.push_back(0.87/GeV2);_rhophaseD.push_back( 0.53*pi);
  }
  // f_2
  _f2mag=0.71/GeV2;_f2phase=0.56*pi;_f2coup=ZERO;
  // sigma
  _sigmamag=2.10;_sigmaphase=0.23*pi;_sigmacoup=0.;
  // f_0
  _f0mag=0.77;_f0phase=-0.54*pi;_f0coup=0.;
  // initialize the a_1 width
  _initializea1=false;
  _a1opt=true;
  double  a1q2in[200]={0      ,15788.6,31577.3,47365.9,63154.6,78943.2,
		       94731.9,110521 ,126309 ,142098 ,157886 ,173675 ,
		       189464 ,205252 ,221041 ,236830 ,252618 ,268407 ,     
		       284196 ,299984 ,315773 ,331562 ,347350 ,363139 ,     
		       378927 ,394716 ,410505 ,426293 ,442082 ,457871 ,     
		       473659 ,489448 ,505237 ,521025 ,536814 ,552603 ,      
		       568391 ,584180 ,599969 ,615757 ,631546 ,647334 ,      
		       663123 ,678912 ,694700 ,710489 ,726278 ,742066 ,      
		       757855 ,773644 ,789432 ,805221 ,821010 ,836798 ,      
		       852587 ,868375 ,884164 ,899953 ,915741 ,931530 ,     
		       947319 ,963107 ,978896 ,994685 ,     
		       1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06, 
		       1.07363e+06,1.08942e+06,1.10521e+06,1.12099e+06, 
		       1.13678e+06,1.15257e+06,1.16836e+06,1.18415e+06, 
		       1.19994e+06,1.21573e+06,1.23151e+06,1.24730e+06, 
		       1.26309e+06,1.27888e+06,1.29467e+06,1.31046e+06, 
		       1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06, 
		       1.38940e+06,1.40519e+06,1.42098e+06,1.43677e+06, 
		       1.45256e+06,1.46834e+06,1.48413e+06,1.49992e+06, 
		       1.51571e+06,1.53150e+06,1.54729e+06,1.56308e+06, 
		       1.57886e+06,1.59465e+06,1.61044e+06,1.62623e+06, 
		       1.64202e+06,1.65781e+06,1.67360e+06,1.68939e+06, 
		       1.70517e+06,1.72096e+06,1.73675e+06,1.75254e+06, 
		       1.76833e+06,1.78412e+06,1.79991e+06,1.81569e+06, 
		       1.83148e+06,1.84727e+06,1.86306e+06,1.87885e+06, 
		       1.89464e+06,1.91043e+06,1.92621e+06,1.94200e+06, 
		       1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06, 
		       2.02095e+06,2.03674e+06,2.05252e+06,2.06831e+06, 
		       2.08410e+06,2.09989e+06,2.11568e+06,2.13147e+06, 
		       2.14726e+06,2.16304e+06,2.17883e+06,2.19462e+06, 
		       2.21041e+06,2.22620e+06,2.24199e+06,2.25778e+06, 
		       2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06, 
		       2.33672e+06,2.35251e+06,2.36830e+06,2.38409e+06, 
		       2.39987e+06,2.41566e+06,2.43145e+06,2.44724e+06, 
		       2.46303e+06,2.47882e+06,2.49461e+06,2.51039e+06, 
		       2.52618e+06,2.54197e+06,2.55776e+06,2.57355e+06, 
		       2.58934e+06,2.60513e+06,2.62092e+06,2.63670e+06, 
		       2.65249e+06,2.66828e+06,2.68407e+06,2.69986e+06, 
		       2.71565e+06,2.73144e+06,2.74722e+06,2.76301e+06, 
		       2.77880e+06,2.79459e+06,2.81038e+06,2.82617e+06, 
		       2.84196e+06,2.85774e+06,2.87353e+06,2.88932e+06, 
		       2.90511e+06,2.92090e+06,2.93669e+06,2.95248e+06, 
		       2.96827e+06,2.98405e+06,2.99984e+06,3.01563e+06, 
		       3.03142e+06,3.04721e+06,3.06300e+06,3.07879e+06, 
		       3.09457e+06,3.11036e+06,3.12615e+06,3.14194e+06};
  double a1widthin[200]={0,0,0,0,0,0,0,0,
			 0,0,0,0.00021256,0.0107225,0.0554708,0.150142,0.303848,
			 0.522655,0.81121,1.1736,1.61381,2.13606,2.74499,3.44583,4.24454,
			 5.14795,6.16391,7.3014,8.57079,9.98398,11.5547,13.2987,15.2344,
			 17.3827,19.7683,22.4195,25.3695,28.6568,32.3264,36.4311,41.0322,
			 46.201,52.0203,58.5847,66.0011,74.3871,83.8666,94.5615,106.578,
			 119.989,134.807,150.968,168.315,186.615,205.576,224.893,244.28,
			 263.499,282.364,300.748,318.569,335.781,352.367,368.327,383.677,
			 398.438,412.638,426.306,439.472,452.167,464.421,476.263,487.719,
			 498.815,509.576,520.024,530.179,540.063,549.693,559.621,568.26,
			 577.229,586.005,594.604,603.035,611.314,619.447,627.446,635.321,
			 643.082,650.736,658.288,665.75,673.127,680.427,687.659,694.82,
			 701.926,708.977,715.983,722.944,729.862,736.752,743.619,750.452,
			 757.271,764.076,770.874,777.658,784.444,791.233,798.027,804.838,
			 811.649,818.485,825.342,832.224,839.139,846.082,853.059,860.079,
			 867.143,874.248,881.409,919.527,945.28,965.514,983.228,999.471,
			 1014.69,1029.15,1043.05,1056.49,1069.57,1082.36,1094.88,1107.2,
			 1120.89,1131.4,1143.33,1155.15,1166.92,1178.61,1190.27,1201.92,
			 1213.55,1225.18,1236.81,1250.06,1260.16,1271.86,1283.64,1295.46,
			 1307.36,1319.3,1331.34,1343.45,1355.64,1367.93,1380.31,1392.77,
			 1405.35,1418.03,1430.83,1443.75,1457.17,1469.94,1483.22,1496.64,
			 1510.18,1523.86,1537.67,1551.64,1565.72,1579.99,1594.38,1608.92,
			 1623.63,1642.08,1653.51,1668.69,1684.03,1699.53,1715.21,1731.04,
			 1747.05,1763.23,1779.59,1796.12,1812.83,1829.72,1846.79,1864.04,
			 1881.49,1899.11,1916.93,1934.93,1953.13,1971.52,1990.12,2008.89};
  if(_a1runwidth.empty()) {
    vector<double> tmp1(a1widthin,a1widthin+200);
    std::transform(tmp1.begin(), tmp1.end(),
		   back_inserter(_a1runwidth),
		   timesGeV);
    
    vector<double> tmp2(a1q2in,a1q2in+200);
    _a1runq2.clear();
    std::transform(tmp2.begin(), tmp2.end(),
		   back_inserter(_a1runq2),
		   timesGeV2);
  }
  // zero parameters which will be calculated later to avoid problems
  _pf2cc=ZERO; 
  _pf200=ZERO;
  _pf0cc=ZERO;
  _pf000=ZERO;
  _psigmacc=ZERO;
  _psigma00=ZERO; 
  _mpi0=ZERO;
  _mpic=ZERO;
  _fact=ZERO;
  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void ThreePionCLEOCurrent::doinit() {
  ThreeMesonCurrentBase::doinit();
  // pointers to the particles we need
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  // the different rho resonances
  tPDPtr rhom[3];
  rhom[0] = getParticleData(-213);
  rhom[1] = getParticleData(-100213);
  rhom[2] = getParticleData(-30213);
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // the f_2
  tPDPtr f2=getParticleData(225);
  // the f_0
  tPDPtr f0=getParticleData(10221);
  if(_localparameters) {
    // make sure the rho array has enough masses
    if(_rhomass.size()<3) {
      for(unsigned int ix=_rhomass.size();ix<3;++ix) {
	_rhomass.push_back(rhom[ix]->mass());
	_rhowidth.push_back(rhom[ix]->width());
      }
    }
  }
  // set the local variables if needed
  else {
    // masses and widths for the particles
    _rhomass.resize(3);_rhowidth.resize(3);
    for(unsigned int ix=0;ix<3;++ix)
      {_rhomass[ix]=rhom[ix]->mass();_rhowidth[ix]=rhom[ix]->width();}
    _f2mass=f2->mass();_f2width=f2->width();
    _f0mass=f0->mass();_f0width=f0->width();
    _sigmamass=sigma->mass();_sigmawidth=sigma->width();
    _a1mass=a1m->mass();_a1width=a1m->width();
    _mKstar=getParticleData(ParticleID::Kstarplus)->mass();
    _mK =getParticleData(ParticleID::Kplus)->mass();
  }
  // parameters for the breit-wigners
  _mpic=getParticleData(ParticleID::piplus)->mass();
  _mpi0=getParticleData(ParticleID::pi0)->mass();
  // momenta of the decay products for on-shell particles
  _psigmacc=Kinematics::pstarTwoBodyDecay(_sigmamass,_mpic,_mpic);
  _psigma00=Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi0,_mpi0);
  _pf2cc=Kinematics::pstarTwoBodyDecay(_f2mass,_mpic,_mpic);
  _pf200=Kinematics::pstarTwoBodyDecay(_f2mass,_mpi0,_mpi0); 
  _pf0cc=Kinematics::pstarTwoBodyDecay(_f0mass,_mpic,_mpic);
  _pf000=Kinematics::pstarTwoBodyDecay(_f0mass,_mpi0,_mpi0); 
  _prhocc.resize(3);_prhoc0.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    _prhocc[ix]=Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpic,_mpic);
    _prhoc0[ix]=Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpic,_mpi0);
  }
  // couplings for the different modes
  Complex ii(0.,1.);
  _rhocoupP.resize(_rhomagP.size());
  for(unsigned int ix=0;ix<_rhomagP.size();++ix)
    {_rhocoupP[ix]=_rhomagP[ix]*(cos(_rhophaseP[ix])+ii*sin(_rhophaseP[ix]));}
  _rhocoupD.resize(_rhomagD.size());
  for(unsigned int ix=0;ix<_rhomagD.size();++ix)
    {_rhocoupD[ix]=_rhomagD[ix]*(cos(_rhophaseD[ix])+ii*sin(_rhophaseD[ix]));}
  _f0coup=_f0mag*(cos(_f0phase)+ii*sin(_f0phase));
  _f2coup=_f2mag*(cos(_f2phase)+ii*sin(_f2phase));
  _sigmacoup=_sigmamag*(cos(_sigmaphase)+ii*sin(_sigmaphase));
  // overall coupling
  _fact = 2.*sqrt(2.)/_fpi/3.;
  // initialise the a_1 running width calculation
  inita1Width(-1);
}

void ThreePionCLEOCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) << ounit(_prhocc,GeV) << ounit(_prhoc0,GeV) 
     << ounit(_f2mass,GeV) << ounit(_f2width,GeV) << ounit(_pf2cc,GeV) 
     << ounit(_pf200,GeV) << ounit(_f0mass,GeV) << ounit(_f0width,GeV) 
     << ounit(_pf0cc,GeV) << ounit(_pf000,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_sigmawidth,GeV)
     << ounit(_psigmacc,GeV) << ounit(_psigma00,GeV) << ounit(_mpi0,GeV) << ounit(_mpic,GeV) 
     << ounit(_fpi,GeV) << ounit(_fact,1/GeV) 
     << _rhomagP << _rhophaseP 
     << _rhocoupP << ounit(_rhomagD,1/GeV2) << _rhophaseD 
     << ounit(_rhocoupD,1/GeV2) <<ounit(_f2mag,1/GeV2) << _f2phase << ounit(_f2coup ,1/GeV2)
     << _f0mag << _f0phase << _f0coup << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters 
     << ounit(_a1mass,GeV) << ounit(_a1width,GeV) << ounit(_a1runwidth,GeV) 
     << ounit(_a1runq2,GeV2) <<  _initializea1
     << ounit(_mKstar,GeV) << ounit(_mK,GeV) << _gammk << _a1opt 
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV) << _a1runinter;
}

void ThreePionCLEOCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) >> iunit(_prhocc,GeV) >> iunit(_prhoc0,GeV) 
     >> iunit(_f2mass,GeV) >> iunit(_f2width,GeV) >> iunit(_pf2cc,GeV)
     >> iunit(_pf200,GeV) >> iunit(_f0mass,GeV) >> iunit(_f0width,GeV) 
     >> iunit(_pf0cc,GeV) >> iunit(_pf000,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_sigmawidth,GeV)
     >> iunit(_psigmacc,GeV) >> iunit(_psigma00,GeV) >> iunit(_mpi0,GeV) >> iunit(_mpic,GeV) 
     >> iunit(_fpi,GeV) >> iunit(_fact,1/GeV) 
     >> _rhomagP >> _rhophaseP 
     >> _rhocoupP >> iunit(_rhomagD,1/GeV2) >> _rhophaseD >> iunit(_rhocoupD,1/GeV2) 
     >> iunit(_f2mag,1/GeV2) >> _f2phase >> iunit(_f2coup,1/GeV2) 
     >> _f0mag >> _f0phase >> _f0coup >> _sigmamag >> _sigmaphase >> _sigmacoup 
     >> _localparameters 
     >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV) >>  iunit(_a1runwidth,GeV) 
     >> iunit(_a1runq2,GeV2) >>  _initializea1
     >> iunit(_mKstar,GeV) >> iunit(_mK,GeV) >> _gammk >> _a1opt 
     >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV) >> _a1runinter;
}

void ThreePionCLEOCurrent::Init() {

  static ClassDocumentation<ThreePionCLEOCurrent> documentation
    ("The ThreePionCLEOCurrent class performs the decay of the"
     " tau to three pions using the currents from CLEO",
     "The decay of tau to three pions is modelled using the currents from "
     "\\cite{Asner:1999kj}.",
     "  %\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

  static ParVector<ThreePionCLEOCurrent,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &ThreePionCLEOCurrent::_rhomass,
     MeV, 0, ZERO, -10000*MeV, 10000*MeV, false, false, true);

  static ParVector<ThreePionCLEOCurrent,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &ThreePionCLEOCurrent::_rhowidth,
     MeV, 0, ZERO, -10000*MeV, 10000*MeV, false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &ThreePionCLEOCurrent::_f2mass, GeV, 1.275*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &ThreePionCLEOCurrent::_f2width, GeV, 0.185*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &ThreePionCLEOCurrent::_f0mass, GeV, 1.186*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &ThreePionCLEOCurrent::_f0width, GeV, 0.350*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &ThreePionCLEOCurrent::_sigmamass, GeV, 0.860*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &ThreePionCLEOCurrent::_sigmawidth, GeV, 0.880*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacea1Mass
    ("a1Mass",
     "The mass of the a_1 meson",
     &ThreePionCLEOCurrent::_a1mass, GeV, 1.331*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfacea1Width
    ("a1Width",
     "The width of the a_1 meson",
     &ThreePionCLEOCurrent::_a1width, GeV, 0.814*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceKaonMass
    ("KaonMass",
     "The mass of the kaon",
     &ThreePionCLEOCurrent::_mK, GeV, 0.496*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceKStarMass
    ("KStarMass",
     "The mass of the k* meson",
     &ThreePionCLEOCurrent::_mKstar, GeV, 0.894*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfaceKaonCoupling
    ("KaonCoupling",
     "The relative coupling for the kaon in the a_1 running width",
     &ThreePionCLEOCurrent::_gammk, 3.32, 0.0, 10.0,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &ThreePionCLEOCurrent::_fpi, MeV, 130.7*MeV/sqrt(2.), ZERO, 500.0*MeV,
     false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &ThreePionCLEOCurrent::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &ThreePionCLEOCurrent::_rhophaseP,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);

  static ParVector<ThreePionCLEOCurrent,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &ThreePionCLEOCurrent::_rhomagD,
     1/MeV2, 0, ZERO, ZERO, 10000/MeV2, false, false, true);

  static ParVector<ThreePionCLEOCurrent,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &ThreePionCLEOCurrent::_rhophaseD,
     0, 0, 0, -Constants::twopi, Constants::twopi, false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &ThreePionCLEOCurrent::_f0phase, 0.54*Constants::pi, -Constants::twopi, Constants::twopi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &ThreePionCLEOCurrent::_f2phase, 0.56*Constants::pi,-Constants::twopi, Constants::twopi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &ThreePionCLEOCurrent::_sigmaphase, 0.23*Constants::pi, -Constants::twopi, Constants::twopi,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &ThreePionCLEOCurrent::_f0mag, 0.77, 0.0, 10,
     false, false, true);


  static Parameter<ThreePionCLEOCurrent,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &ThreePionCLEOCurrent::_f2mag, 1./GeV2, 0.71/GeV2, ZERO, 10./GeV2,
     false, false, true);

  static Parameter<ThreePionCLEOCurrent,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &ThreePionCLEOCurrent::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static Switch<ThreePionCLEOCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &ThreePionCLEOCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particleData objects",
     false);

  static ParVector<ThreePionCLEOCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreePionCLEOCurrent::_a1runwidth,
     MeV, 0, ZERO, ZERO, 10000000*MeV, false, false, true);
  
  static ParVector<ThreePionCLEOCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreePionCLEOCurrent::_a1runq2,
     MeV2, 0, ZERO, ZERO, 10000000*MeV2, false, false, true);

  static Switch<ThreePionCLEOCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreePionCLEOCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Yes",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "No",
     "Use the default values",
     false);
  
  static Switch<ThreePionCLEOCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &ThreePionCLEOCurrent::_a1opt, true, false, false);
  static SwitchOption interfacea1WidthOptionLocal
    (interfacea1WidthOption,
     "Local",
     "Use a calculation of the running width based on the parameters as"
     " interpolation table.",
     true);
  static SwitchOption interfacea1WidthOptionParam
    (interfacea1WidthOption,
     "Kuhn",
     "Use the parameterization of Kuhn and Santamaria for default parameters."
     " This should only be used for testing vs TAUOLA",
     false);
}

// initialisation of the a_1 width
// (iopt=-1 initialises, iopt=0 starts the interpolation)
void ThreePionCLEOCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return;
    // parameters for the table of values
    Energy2 step=sqr(_maxmass)/200.;
    // function to be integrated to give the matrix element
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho=getParticleData(ParticleID::rhoplus)->mass();
    Energy wrho=getParticleData(ParticleID::rhoplus)->width();
    vector<Energy> inmass;inmass.push_back(mrho);inmass.push_back(mrho);
    vector<Energy> inwidth;inwidth.push_back(wrho);inwidth.push_back(wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<ThreePionCLEOCurrent>
      widthgenN(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi0,_mpi0,_mpic);
    ThreeBodyAllOnCalculator<ThreePionCLEOCurrent>
      widthgenC(inweights,intype,inmass,inwidth,inpow,*this,1,_mpic,_mpic,_mpic);
    // normalisation constant to give physical width if on shell
    double a1const = _a1width/(widthgenN.partialWidth(sqr(_a1mass))+
			       widthgenC.partialWidth(sqr(_a1mass)));
    // loop to give the values
    _a1runq2.clear();_a1runwidth.clear();
    for(Energy2 moff2=ZERO; moff2<=sqr(_maxmass); moff2+=step) {
      Energy moff=sqrt(moff2);
      _a1runq2.push_back(moff2);
      Energy charged=a1const*widthgenC.partialWidth(moff2);
      Energy neutral=a1const*widthgenN.partialWidth(moff2);
      Energy kaon = moff<=_mK+_mKstar ? ZERO : 2.870*_gammk*_gammk/8./Constants::pi*
	Kinematics::pstarTwoBodyDecay(moff,_mK,_mKstar)/moff2*GeV2;
      Energy total = charged + neutral + kaon;
      _a1runwidth.push_back(total);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}

// modes handled by this class
bool ThreePionCLEOCurrent::acceptMode(int imode) const {
  return imode>=0&&imode<=1;
}

ThreePionCLEOCurrent::FormFactors ThreePionCLEOCurrent::
  calculateFormFactors(const int ichan, const int imode,
		      Energy2 q2, Energy2 s1, Energy2 s2, Energy2 s3) const {
  Complex F1=0., F2=0., F3=0.;//F4=0.;F5=0.;
  // calculate the form factors without the a_1 piece
  CLEOFormFactor(imode,ichan,q2,s1,s2,s3,F1,F2,F3);
  // change sign of the f_2 term
  F2=-F2;
  // multiply by the a_1 factor
  complex<InvEnergy> a1fact = a1BreitWigner(q2) * _fact;
  return FormFactors(F1*a1fact, F2*a1fact, F3*a1fact);
}

void ThreePionCLEOCurrent::CLEOFormFactor(int imode,int ichan,
					  Energy2 q2,Energy2 s1, Energy2 s2, Energy2 s3,
					  Complex & F1, Complex & F2, 
					  Complex & F3) const {
  useMe();
  if(imode==0) {
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bws1,sigbws1,f2bws1,f0bws2,sigbws2,f2bws2;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,0);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,0);
    }
    f0bws1  =f0BreitWigner(s1,0);
    sigbws1 =sigmaBreitWigner(s1,0);
    f2bws1  =f2BreitWigner(s1,0);
    f0bws2  =f0BreitWigner(s2,0);
    sigbws2 =sigmaBreitWigner(s2,0);
    f2bws2  =f2BreitWigner(s2,0);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1-=_rhocoupP[ix]*rhos1bw[ix];
	F2-=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=1./3.*(s1-s3);
      Energy2 Dfact2=1./3.*(s2-s3);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1-=Dfact1*_rhocoupD[ix]*rhos2bw[ix];
	F2-=Dfact2*_rhocoupD[ix]*rhos1bw[ix];
	F3-=_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]);
      }
      // the scalar terms
      F1-=2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2-=2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3+=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	  +2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      complex<Energy2> sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1+=_f2coup*(0.5*(s3-s2)*f2bws1-sfact2);
      F2+=_f2coup*(0.5*(s3-s1)*f2bws2-sfact1);
      F3+=_f2coup*(-sfact1+sfact2);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      Energy2 Dfact2=1./3.*(s2-s3);
      if(ires<_rhocoupP.size()) F1-=_rhocoupP[ires]*rhos1bw[ires];
      if(ires<_rhocoupD.size()) {
	F2-=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	F3-=_rhocoupD[ires]*Dfact2*rhos1bw[ires];
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      Energy2 Dfact1=1./3.*(s1-s3);
      if(ires<_rhocoupP.size()) {
	F2-=_rhocoupP[ires]*rhos2bw[ires];
      }
      if(ires<_rhocoupD.size()) {
	F1-=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	F3+=_rhocoupD[ires]*Dfact1*rhos2bw[ires];
      }
    }
    else if(ichan==6) {
      F2-=2./3.*_sigmacoup*sigbws1;
      F3-=2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==7) {
      F1-=2./3.*_sigmacoup*sigbws2;
      F3+=2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==8) {
      complex<Energy2> sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      F1+=_f2coup*0.5*(s3-s2)*f2bws1;
      F2-=_f2coup*sfact1;
      F3-=_f2coup*sfact1;
    }
    else if(ichan==9) {
      complex<Energy2> sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1-=_f2coup*sfact2;
      F2+=_f2coup*0.5*(s3-s1)*f2bws2;
      F3+=_f2coup*sfact2;
    }
    else if(ichan==10) {
      F2-=2./3.*_f0coup*f0bws1;
      F3-=2./3.*_f0coup*f0bws1;
    }
    else if(ichan==11) {
      F1-=2./3.*_f0coup*f0bws2;
      F3+=2./3.*_f0coup*f0bws2;
    }
  }
  // calculate the pi0 pi0 pi+ factor
  else if(imode==1) {
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
    }
    f0bw  =f0BreitWigner(s3,1);
    sigbw =sigmaBreitWigner(s3,1);
    f2bw  =f2BreitWigner(s3,1);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Dfact1*_rhocoupD[ix]*rhos2bw[ix];
	F2+=Dfact2*_rhocoupD[ix]*rhos1bw[ix];
	F3+=_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]);
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;F2+=Dfact3;
      F3-=0.5*_f2coup*(s1-s2)*f2bw;
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()){F1+=_rhocoupP[ires]*rhos1bw[ires];}
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F2+=Dfact2*_rhocoupD[ires]*rhos1bw[ires];
	F3+=_rhocoupD[ires]*Dfact2*rhos1bw[ires];
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()){F2+=_rhocoupP[ires]*rhos2bw[ires];}
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F1+=Dfact1*_rhocoupD[ires]*rhos2bw[ires];
	F3-=_rhocoupD[ires]*Dfact1*rhos2bw[ires];
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;F2+=Dfact3;
      F3-=0.5*_f2coup*(s1-s2)*f2bw;
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  else {
    throw DecayIntegratorError() << "ThreePionCLEOCurrent Unknown Decay" << imode
				 << Exception::abortnow;
  }
  // identical particle factors
  double fact = 1./sqrt(2.);
  F1*=fact;F2*=fact;F3*=fact;
}

// complete the construction of the decay mode for integration
bool ThreePionCLEOCurrent::createMode(int icharge, unsigned int imode,
				      DecayPhaseSpaceModePtr mode,
				      unsigned int iloc,unsigned int ires,
				      DecayPhaseSpaceChannelPtr phase,Energy upp) {
  if(!acceptMode(imode)){return false;}
  int iq(0),ia(0);
  tPDVector extpart=particles(1,imode,iq,ia);
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  _maxmass=max(_maxmass,upp);
  // pointers to the particles we need
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  // the different rho resonances
  tPDPtr rhom[3];
  if(icharge==-3) {
    rhom[0] = getParticleData(-213);
    rhom[1] = getParticleData(-100213);
    rhom[2] = getParticleData(-30213);
  }
  else if(icharge==3) {
    rhom[0] = getParticleData(213);
    rhom[1] = getParticleData(100213);
    rhom[2] = getParticleData(30213);
  }
  else {
    return false;
  }
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // the f_2
  tPDPtr f2=getParticleData(225);
  // the f_0
  tPDPtr f0=getParticleData(10221);
  // set up the integration channels
  DecayPhaseSpaceChannelPtr newchannel;
  if(imode==0) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(!rho0[ix]) continue;
      // the neutral rho channels
      // first channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(rho0[ix],0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(rho0[ix],0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);      
    }
    // the sigma channels
    if(sigma) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(sigma,0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(sigma,0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);
    }
    // the f_2 channels
    if(f2) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(f2,0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(f2,0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);
    }
    if(f0) {
      // the f_0 channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(f0,0,0.0,iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // interchanged channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(f0,0,0.0,iloc,iloc+2);
      mode->addChannel(newchannel);
    }
  }
  else {
    for(unsigned int ix=0;ix<3;++ix) {
      if(!rhom[ix]) continue;
      // first rho+ channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(rhom[ix],0,0.0,iloc+2,iloc+1);
      mode->addChannel(newchannel);
      // second rho+ channel
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(rhom[ix],0,0.0,iloc+2,iloc);
      mode->addChannel(newchannel);
    }
    // the sigma channel
    if(sigma) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(sigma,0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
    //  the f_2  channel
    if(f2) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(f2,0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
    // the f_0 channel
    if(f0) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(a1m,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(f0,0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
  }
  if(_localparameters) {
    for(unsigned int iy=0;iy<_rhomass.size();++iy) {
      if(rho0[iy]) mode->resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
      if(rhom[iy]) mode->resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
    }
    if(sigma) mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
    if(f2)    mode->resetIntermediate(f2,_f2mass,_f2width);
    if(f0)    mode->resetIntermediate(f0,_f0mass,_f0width);
    mode->resetIntermediate(a1m,_a1mass,_a1width);
  }
  return true;
}

void ThreePionCLEOCurrent::dataBaseOutput(ofstream & output,bool header,
					  bool create) const {
  if(header){output << "update decayers set parameters=\"";}
  if(create) {
    output << "create Herwig::ThreePionCLEOCurrent " << name() 
	   << " HwWeakCurrents.so\n";
  }
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoMasses " << ix 
	     << " " << _rhomass[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":RhoMasses " << ix 
	     << " " << _rhomass[ix]/MeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoWidths " << ix 
	     << " " << _rhowidth[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":RhoWidths " << ix 
	     << " " << _rhowidth[ix]/MeV << "\n";
    }
  }
  output << "newdef " << name() << ":f_2Mass " << _f2mass/GeV << "\n";
  output << "newdef " << name() << ":f_2Width " << _f2width/GeV << "\n";
  output << "newdef " << name() << ":f_0Mass " << _f0mass/GeV << "\n";
  output << "newdef " << name() << ":f_0Width " << _f0width/GeV << "\n";
  output << "newdef " << name() << ":sigmaMass " << _sigmamass/GeV << "\n";
  output << "newdef " << name() << ":sigmaWidth " << _sigmawidth/GeV << "\n";
  output << "newdef " << name() << ":a1Mass " << _a1mass/GeV << "\n";
  output << "newdef " << name() << ":a1Width " <<_a1width /GeV << "\n";
  output << "newdef " << name() << ":KaonMass " << _mK/GeV << "\n";
  output << "newdef " << name() << ":KStarMass " << _mKstar/GeV << "\n";
  output << "newdef " << name() << ":KaonCoupling " << _gammk << "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":a1WidthOption " << _a1opt << "\n";
  for(unsigned int ix=0;ix<_rhomagP.size();++ix) {
      if(ix<2) {
	output << "newdef    " << name() << ":RhoPWaveMagnitude " << ix 
	       << " " << _rhomagP[ix] << "\n";
      }
      else {
	output << "insert " << name() << ":RhoPWaveMagnitude " << ix 
	       << " " << _rhomagP[ix] << "\n";
      }
  }
  for(unsigned int ix=0;ix<_rhophaseP.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoPWavePhase " << ix 
	     << " " << _rhophaseP[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":RhoPWavePhase " << ix 
	     << " " << _rhophaseP[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rhomagD.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoDWaveMagnitude " << ix 
	     << " " << _rhomagD[ix]*MeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":RhoDWaveMagnitude " << ix 
	     << " " << _rhomagD[ix]*MeV2 << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rhophaseD.size();++ix) {
    if(ix<2) {
      output << "newdef    " << name() << ":RhoDWavePhase " << ix 
	     << " " << _rhophaseD[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":RhoDWavePhase " << ix 
	     << " " << _rhophaseD[ix] << "\n";
    }
  }
  output << "newdef " << name() << ":f0Phase " << _f0phase << "\n";
  output << "newdef " << name() << ":f2Phase " <<_f2phase  << "\n";
  output << "newdef " << name() << ":sigmaPhase " <<_sigmaphase  << "\n";
  output << "newdef " << name() << ":f0Magnitude " << _f0mag << "\n";
  output << "newdef " << name() << ":f2Magnitude " << _f2mag*GeV2 << "\n";
  output << "newdef " << name() << ":sigmaMagnitude " <<_sigmamag  << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Initializea1 " <<_initializea1  << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    if(ix<200) {
      output << "newdef    " << name() << ":a1RunningWidth " << ix 
	     << " " << _a1runwidth[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":a1RunningWidth " << ix 
	     << " " << _a1runwidth[ix]/MeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    if(ix<200) {
      output << "newdef    " << name() << ":a1RunningQ2 " << ix 
	     << " " << _a1runq2[ix]/MeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":a1RunningQ2 " << ix 
	     << " " << _a1runq2[ix]/MeV2 << "\n";
    }
  }
  ThreeMesonCurrentBase::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void ThreePionCLEOCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  ThreeMesonCurrentBase::doinitrun();
}

void ThreePionCLEOCurrent::doupdate() {
  ThreeMesonCurrentBase::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

Energy ThreePionCLEOCurrent::a1width(Energy2 q2) const {
  Energy output;
  if(_a1opt) output=(*_a1runinter)(q2);
  else {
    double gam(0.);
    if(q2<0.1753*GeV2) {
      gam =0.;
    }
    else if(q2<0.823*GeV2) {
      double p=q2/GeV2-0.1753;
      gam = 5.80900*p*sqr(p)*(1.-3.00980*p+4.57920*sqr(p));
    }
    else {
      double p=q2/GeV2;
      gam = -13.91400+27.67900*p-13.39300*sqr(p)
	+3.19240*sqr(p)*p-0.10487*sqr(sqr(p));
    }
    if(q2<0.1676*GeV2) {
      gam+=0.;
    }
    else if(q2<0.823*GeV2) {
      double p=q2/GeV2-0.1676;
      gam+= 6.28450*p*sqr(p)*(1.-2.95950*p+4.33550*sqr(p));
    }
    else {
      double p=q2/GeV2;
      gam+= -15.41100+32.08800*p-17.66600*sqr(p)
	+4.93550*sqr(p)*p-0.37498*sqr(sqr(p));
    }
    Energy mkst=0.894*GeV,mk=0.496*GeV;
    Energy2 mk1sq=sqr(mkst+mk), mk2sq=sqr(mkst-mk);
    double c3pi=sqr(0.2384),ckst=sqr(4.7621)*c3pi;
    gam*=c3pi;
    if(q2>mk1sq) gam+=0.5*ckst*sqrt((q2-mk1sq)*(q2-mk2sq))/q2;
    gam = gam*_a1width*_a1mass/GeV2/1.331/0.814/1.0252088;
    output = gam*GeV2/sqrt(q2);
  }
  return output;
}

double 
ThreePionCLEOCurrent::threeBodyMatrixElement(const int iopt, const Energy2 q2,
					     const Energy2 s3, const Energy2 s2, 
					     const Energy2 s1, const Energy,
					     const Energy, const Energy) const {
  Energy p1[5],p2[5],p3[5];
  Energy2 p1sq, p2sq, p3sq;
  Energy q=sqrt(q2);
  Energy2 mpi2c=_mpic*_mpic;
  Energy2 mpi20=_mpi0*_mpi0;
  // construct the momenta for the 2 neutral 1 charged mode
  Complex F1,F2,F3;
  if(iopt==0) {
    // construct the momenta of the decay products
    p1[0] = 0.5*(q2+mpi20-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-mpi20);
    p2[0] = 0.5*(q2+mpi20-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-mpi20);
    p3[0] = 0.5*(q2+mpi2c-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-mpi2c);
    // take momentum of 1 parallel to z axis
    p1[1]=ZERO;p1[2]=ZERO;p1[3]=p1[4];
    // construct 2 
    double cos2 = 0.5*(p1sq+p2sq-p3sq-2.*mpi20+mpi2c)/p1[4]/p2[4];
    p2[1] = p2[4]*sqrt(1.-cos2*cos2); p2[2]=ZERO; p2[3]=-p2[4]*cos2;
    // construct 3
    double cos3 = 0.5*(p1sq-p2sq+p3sq-mpi2c)/p1[4]/p3[4];
    p3[1] =-p3[4]*sqrt(1.-cos3*cos3); p3[2]=ZERO; p3[3]=-p3[4]*cos3; 
    // calculate the form factors
    CLEOFormFactor(1,-1,q2,s1,s2,s3,F1,F2,F3);
  }
  // construct the momenta for the 3 charged mode 
  else {
    // construct the momenta of the decay products
    p1[0] = 0.5*(q2+mpi2c-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-mpi2c);
    p2[0] = 0.5*(q2+mpi2c-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-mpi2c);
    p3[0] = 0.5*(q2+mpi2c-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-mpi2c);
    // take momentum of 1 parallel to z axis
    p1[1]=ZERO;p1[2]=ZERO;p1[3]=p1[4];
    // construct 2 
    double cos2 = 0.5*(p1sq+p2sq-p3sq-mpi2c)/p1[4]/p2[4];
    p2[1] = p2[4]*sqrt(1.-cos2*cos2); p2[2]=ZERO; p2[3]=-p2[4]*cos2;
    // construct 3
    double cos3 = 0.5*(p1sq-p2sq+p3sq-mpi2c)/p1[4]/p3[4];
    p3[1] =-p3[4]*sqrt(1.-cos3*cos3); p3[2]=ZERO; p3[3]=-p3[4]*cos3; 
    // calculate the form factors
    CLEOFormFactor(0,-1,q2,s1,s2,s3,F1,F2,F3);
  }
  // construct a vector with the current
  complex<Energy> current[4];
  for(unsigned int ix=0;ix<4;++ix)
    current[ix] = F1*(p2[ix]-p3[ix])-F2*(p3[ix]-p1[ix])+F3*(p1[ix]-p2[ix]);
  complex<Energy2> dot1=current[0]*conj(current[0]);
  for(unsigned int ix=1;ix<4;++ix) dot1-=current[ix]*conj(current[ix]);
  complex<Energy2> dot2=current[0]*q;
  return(-dot1+dot2*conj(dot2)/q2).real() / sqr(_rhomass[0]);
}
