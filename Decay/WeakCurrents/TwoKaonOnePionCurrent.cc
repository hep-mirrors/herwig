// -*- C++ -*-
//
// TwoKaonOnePionCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoKaonOnePionCurrent class.
//

#include "TwoKaonOnePionCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

DescribeClass<TwoKaonOnePionCurrent,WeakCurrent>
describeHerwigTwoKaonOnePionCurrent("Herwig::TwoKaonOnePionCurrent",
				    "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(TwoKaonOnePionCurrent,Energy,Energy2)

IBPtr TwoKaonOnePionCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr TwoKaonOnePionCurrent::fullclone() const {
  return new_ptr(*this);
}

TwoKaonOnePionCurrent::TwoKaonOnePionCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(12);
  // rho parameters
  // use local values
  _rhoparameters=true;
  // rho parameters for axial-vector pieces
  _rho1wgts.push_back( 1.0  );   _rho1wgts.push_back(-0.145);
  _rho1wgts.push_back(0.);
  _rho1mass.push_back(0.773*GeV);_rho1mass.push_back(1.370*GeV);
  _rho1mass.push_back(1.750*GeV);
  _rho1width.push_back(0.145*GeV);_rho1width.push_back(0.510*GeV);
  _rho1width.push_back(0.120*GeV);
  // rho parameters for vector pieces
  _rho2wgts.push_back( 1.0  );   _rho2wgts.push_back(-0.25 );
  _rho2wgts.push_back(-0.038);
  _rho2mass.push_back(0.773*GeV);_rho2mass.push_back(1.500*GeV);
  _rho2mass.push_back(1.750*GeV);
  _rho2width.push_back(0.145*GeV);_rho2width.push_back(0.220*GeV);
  _rho2width.push_back(0.120*GeV);
  // K* parameters
  _kstarparameters=true;
  // K* parameters for the axial-vector pieces
  _kstar1wgts.push_back( 1.0  );   _kstar1wgts.push_back(-0.135);
  _kstar1wgts.push_back(0.);
  _kstar1mass.push_back(0.892*GeV);_kstar1mass.push_back(1.412*GeV);
  _kstar1mass.push_back(1.714*GeV);
  _kstar1width.push_back(0.050*GeV);_kstar1width.push_back(0.227*GeV);
  _kstar1width.push_back(0.323*GeV);
  // K* parameters for vector pieces
  _kstar2wgts.push_back( 1.0  );   _kstar2wgts.push_back(-0.25 );
  _kstar2wgts.push_back(-0.038);
  _kstar2mass.push_back(0.892*GeV);_kstar2mass.push_back(1.412*GeV);
  _kstar2mass.push_back(1.714*GeV);
  _kstar2width.push_back(0.050*GeV);_kstar2width.push_back(0.227*GeV);
  _kstar2width.push_back(0.323*GeV);
  // a_1 parameters
  _a1parameters = true;
  _initializea1 = false;
  _a1opt        = true;
  _a1mass  = 1.251*GeV;
  _a1width = 0.475*GeV;
  
  _a1runwidth = {0*GeV, 0*GeV, 0*GeV, 0*GeV, 0*GeV,
		 0*GeV, 0*GeV, 0*GeV, 0*GeV, 0*GeV,
		 0*GeV, 0*GeV, 1.47729e-06*GeV, 1.19209e-05*GeV, 3.884e-05*GeV,
		 8.83255e-05*GeV, 0.00016561*GeV, 0.000275439*GeV, 0.000422332*GeV,
		 0.000610773*GeV, 0.000845357*GeV, 0.00113092*GeV, 0.00147264*GeV,
		 0.00187616*GeV, 0.0023477*GeV, 0.00289413*GeV, 0.00352315*GeV,
		 0.00424342*GeV, 0.0050647*GeV, 0.00599808*GeV, 0.00705616*GeV,
		 0.00825335*GeV, 0.0096062*GeV, 0.0111337*GeV, 0.0128579*GeV,
		 0.0148041*GeV, 0.017002*GeV, 0.0194858*GeV, 0.0222956*GeV,
		 0.0254781*GeV, 0.0290874*GeV, 0.0331862*GeV, 0.0378467*GeV,
		 0.0431501*GeV, 0.0491862*GeV, 0.0560496*GeV, 0.0638341*GeV,
		 0.0726215*GeV, 0.0824662*GeV, 0.0933765*GeV, 0.105297*GeV,
		 0.118103*GeV, 0.131602*GeV, 0.145564*GeV, 0.159749*GeV,
		 0.173938*GeV, 0.18795*GeV, 0.201649*GeV, 0.214943*GeV,
		 0.227773*GeV, 0.240109*GeV, 0.25194*GeV, 0.263268*GeV,
		 0.274104*GeV, 0.284466*GeV, 0.294372*GeV, 0.303845*GeV,
		 0.312905*GeV, 0.321576*GeV, 0.329878*GeV, 0.337832*GeV,
		 0.345456*GeV, 0.35277*GeV, 0.35979*GeV, 0.366532*GeV,
		 0.373012*GeV, 0.379243*GeV, 0.38524*GeV, 0.391014*GeV,
		 0.396577*GeV, 0.401939*GeV, 0.407111*GeV, 0.412102*GeV,
		 0.416923*GeV, 0.421577*GeV, 0.426078*GeV, 0.430427*GeV,
		 0.434636*GeV, 0.43871*GeV, 0.442654*GeV, 0.446475*GeV,
		 0.450177*GeV, 0.453765*GeV, 0.457245*GeV, 0.460621*GeV,
		 0.463899*GeV, 0.467077*GeV, 0.470164*GeV, 0.473162*GeV,
		 0.476076*GeV, 0.478909*GeV, 0.481658*GeV, 0.484333*GeV,
		 0.486934*GeV, 0.489465*GeV, 0.491926*GeV, 0.494321*GeV,
		 0.496651*GeV, 0.49892*GeV, 0.501128*GeV, 0.503277*GeV,
		 0.505371*GeV, 0.507409*GeV, 0.509395*GeV, 0.511328*GeV,
		 0.513212*GeV, 0.515047*GeV, 0.516846*GeV, 0.518624*GeV,
		 0.520285*GeV, 0.52194*GeV, 0.523553*GeV, 0.525124*GeV,
		 0.526646*GeV, 0.52814*GeV, 0.529638*GeV, 0.531016*GeV,
		 0.532401*GeV, 0.533751*GeV, 0.535069*GeV, 0.536354*GeV,
		 0.537608*GeV, 0.538831*GeV, 0.540039*GeV, 0.541194*GeV,
		 0.542327*GeV, 0.543438*GeV, 0.544522*GeV, 0.545582*GeV,
		 0.546616*GeV, 0.54764*GeV, 0.548615*GeV, 0.549581*GeV,
		 0.550525*GeV, 0.551449*GeV, 0.552351*GeV, 0.55324*GeV,
		 0.554101*GeV, 0.554944*GeV, 0.555772*GeV, 0.556583*GeV,
		 0.557373*GeV, 0.558155*GeV, 0.558917*GeV, 0.559664*GeV,
		 0.560396*GeV, 0.561114*GeV, 0.561849*GeV, 0.562508*GeV,
		 0.563186*GeV, 0.563851*GeV, 0.564503*GeV, 0.565145*GeV,
		 0.565774*GeV, 0.566394*GeV, 0.567001*GeV, 0.567595*GeV,
		 0.568182*GeV, 0.56876*GeV, 0.56933*GeV, 0.569886*GeV,
		 0.570433*GeV, 0.570976*GeV, 0.571504*GeV, 0.572027*GeV,
		 0.572542*GeV, 0.573114*GeV, 0.573548*GeV, 0.574108*GeV,
		 0.574524*GeV, 0.575002*GeV, 0.575473*GeV, 0.575937*GeV,
		 0.576394*GeV, 0.576845*GeV, 0.57729*GeV, 0.57773*GeV,
		 0.578173*GeV, 0.5786*GeV, 0.579013*GeV, 0.579431*GeV,
		 0.579834*GeV, 0.580246*GeV, 0.580649*GeV, 0.581045*GeV,
		 0.581437*GeV, 0.581827*GeV, 0.582208*GeV, 0.582586*GeV, 0.582959*GeV};
  _a1runq2 = {  0*GeV2       , 0.0158678*GeV2, 0.0317356*GeV2, 0.0476034*GeV2, 0.0634712*GeV2,
		0.079339*GeV2, 0.0952068*GeV2,  0.111075*GeV2, 0.126942*GeV2, 0.14281*GeV2,
	        0.158678*GeV2, 0.174546*GeV2, 0.190414*GeV2, 0.206281*GeV2, 0.222149*GeV2,
		0.238017*GeV2, 0.253885*GeV2, 0.269753*GeV2, 0.285621*GeV2, 0.301488*GeV2,
		0.317356*GeV2, 0.333224*GeV2, 0.349092*GeV2, 0.36496*GeV2, 0.380827*GeV2,
		0.396695*GeV2, 0.412563*GeV2, 0.428431*GeV2, 0.444299*GeV2, 0.460166*GeV2,
		0.476034*GeV2, 0.491902*GeV2, 0.50777*GeV2, 0.523638*GeV2, 0.539505*GeV2,
		0.555373*GeV2, 0.571241*GeV2, 0.587109*GeV2, 0.602977*GeV2, 0.618844*GeV2,
		0.634712*GeV2, 0.65058*GeV2, 0.666448*GeV2, 0.682316*GeV2, 0.698183*GeV2,
		0.714051*GeV2, 0.729919*GeV2, 0.745787*GeV2, 0.761655*GeV2, 0.777523*GeV2,
		0.79339*GeV2, 0.809258*GeV2, 0.825126*GeV2, 0.840994*GeV2, 0.856862*GeV2,
		0.872729*GeV2, 0.888597*GeV2, 0.904465*GeV2, 0.920333*GeV2, 0.936201*GeV2,
		0.952068*GeV2, 0.967936*GeV2, 0.983804*GeV2, 0.999672*GeV2, 1.01554*GeV2,
		1.03141*GeV2, 1.04728*GeV2, 1.06314*GeV2, 1.07901*GeV2, 1.09488*GeV2,
		1.11075*GeV2, 1.12661*GeV2, 1.14248*GeV2, 1.15835*GeV2, 1.17422*GeV2,
		1.19009*GeV2, 1.20595*GeV2, 1.22182*GeV2, 1.23769*GeV2, 1.25356*GeV2,
		1.26942*GeV2, 1.28529*GeV2, 1.30116*GeV2, 1.31703*GeV2, 1.3329*GeV2,
		1.34876*GeV2, 1.36463*GeV2, 1.3805*GeV2, 1.39637*GeV2, 1.41223*GeV2,
		1.4281*GeV2, 1.44397*GeV2, 1.45984*GeV2, 1.47571*GeV2, 1.49157*GeV2,
		1.50744*GeV2, 1.52331*GeV2, 1.53918*GeV2, 1.55505*GeV2, 1.57091*GeV2,
		1.58678*GeV2, 1.60265*GeV2, 1.61852*GeV2, 1.63438*GeV2, 1.65025*GeV2,
		1.66612*GeV2, 1.68199*GeV2, 1.69786*GeV2, 1.71372*GeV2, 1.72959*GeV2,
		1.74546*GeV2, 1.76133*GeV2, 1.77719*GeV2, 1.79306*GeV2, 1.80893*GeV2,
		1.8248*GeV2, 1.84067*GeV2, 1.85653*GeV2, 1.8724*GeV2, 1.88827*GeV2,
		1.90414*GeV2, 1.92*GeV2, 1.93587*GeV2, 1.95174*GeV2, 1.96761*GeV2,
		1.98348*GeV2, 1.99934*GeV2, 2.01521*GeV2, 2.03108*GeV2, 2.04695*GeV2,
		2.06281*GeV2, 2.07868*GeV2, 2.09455*GeV2, 2.11042*GeV2, 2.12629*GeV2,
		2.14215*GeV2, 2.15802*GeV2, 2.17389*GeV2, 2.18976*GeV2, 2.20563*GeV2,
		2.22149*GeV2, 2.23736*GeV2, 2.25323*GeV2, 2.2691*GeV2, 2.28496*GeV2,
		2.30083*GeV2, 2.3167*GeV2, 2.33257*GeV2, 2.34844*GeV2, 2.3643*GeV2,
		2.38017*GeV2, 2.39604*GeV2, 2.41191*GeV2, 2.42777*GeV2, 2.44364*GeV2,
		2.45951*GeV2, 2.47538*GeV2, 2.49125*GeV2, 2.50711*GeV2, 2.52298*GeV2,
		2.53885*GeV2, 2.55472*GeV2, 2.57058*GeV2, 2.58645*GeV2, 2.60232*GeV2,
		2.61819*GeV2, 2.63406*GeV2, 2.64992*GeV2, 2.66579*GeV2, 2.68166*GeV2,
		2.69753*GeV2, 2.71339*GeV2, 2.72926*GeV2, 2.74513*GeV2, 2.761*GeV2,
		2.77687*GeV2, 2.79273*GeV2, 2.8086*GeV2, 2.82447*GeV2, 2.84034*GeV2,
		2.85621*GeV2, 2.87207*GeV2, 2.88794*GeV2, 2.90381*GeV2, 2.91968*GeV2,
		2.93554*GeV2, 2.95141*GeV2, 2.96728*GeV2, 2.98315*GeV2, 2.99902*GeV2,
		3.01488*GeV2, 3.03075*GeV2, 3.04662*GeV2, 3.06249*GeV2, 3.07835*GeV2,
		3.09422*GeV2, 3.11009*GeV2, 3.12596*GeV2, 3.14183*GeV2, 3.15769*GeV2};
  // K_1 parameters
  _k1parameters=true;
  _k1mass.push_back(1.270*GeV);_k1width.push_back(0.090*GeV);
  _k1mass.push_back(1.402*GeV);_k1width.push_back(0.174*GeV);
  _k1wgta.push_back(0.33);_k1wgta.push_back(1.);
  _k1wgtb.push_back(1.00);_k1wgtb.push_back(0.);
  // parameters for the T_omega function
  _omegaopt=true;
  _epsomega=0.05;
  _omegamass  = 0.782*GeV;
  _omegawidth = 0.00843*GeV;
  _phimass    = 1.020*GeV;
  _phiwidth   = 0.00443*GeV;
  _omegaKstarwgt=1./sqrt(2.);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi=ZERO;_mK=ZERO;
  _maxmass=ZERO;
  _maxcalc=ZERO;
}


void TwoKaonOnePionCurrent::persistentOutput(PersistentOStream & os) const {
  os << _a1runinter
     << _rho1wgts << ounit(_rho1mass,GeV) << ounit(_rho1width,GeV) 
     << _rho2wgts << ounit(_rho2mass,GeV) << ounit(_rho2width,GeV)
     << _kstar1wgts << ounit(_kstar1mass,GeV) << ounit(_kstar1width,GeV) 
     << _kstar2wgts << ounit(_kstar2mass,GeV) << ounit(_kstar2width,GeV) 
     << ounit(_a1mass,GeV) << ounit(_a1width,GeV) << ounit(_k1mass,GeV) 
     << ounit(_k1width,GeV) << _k1wgta << _k1wgtb 
     << ounit(_a1runwidth,GeV) << ounit(_a1runq2,GeV2) << _epsomega 
     << ounit(_omegamass,GeV) << ounit(_omegawidth,GeV) 
     << ounit(_phimass,GeV) << ounit(_phiwidth,GeV) << _omegaKstarwgt 
     << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV) 
     << _initializea1 << _rhoparameters << _kstarparameters << _a1parameters 
     << _k1parameters << _a1opt << _omegaopt 
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV);
}

void TwoKaonOnePionCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _a1runinter
     >> _rho1wgts >> iunit(_rho1mass,GeV) >> iunit(_rho1width,GeV) 
     >> _rho2wgts >> iunit(_rho2mass,GeV) >> iunit(_rho2width,GeV) 
     >> _kstar1wgts >> iunit(_kstar1mass,GeV) >> iunit(_kstar1width,GeV) 
     >> _kstar2wgts >> iunit(_kstar2mass,GeV) >> iunit(_kstar2width,GeV) 
     >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV) >> iunit(_k1mass,GeV) 
     >> iunit(_k1width,GeV) >> _k1wgta >> _k1wgtb 
     >> iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) >> _epsomega 
     >> iunit(_omegamass,GeV) >> iunit(_omegawidth,GeV) 
     >> iunit(_phimass,GeV) >> iunit(_phiwidth,GeV) >> _omegaKstarwgt 
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV) 
     >> _initializea1 >> _rhoparameters >> _kstarparameters >> _a1parameters 
     >> _k1parameters >> _a1opt >> _omegaopt 
     >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV);
}


void TwoKaonOnePionCurrent::Init() {

  static ClassDocumentation<TwoKaonOnePionCurrent> documentation
    ("The TwoKaonOnePionCurrent class implements the model of "
     "Z. Phys.  C 69 (1996) 243 [arXiv:hep-ph/9503474]"
     " for the weak current with three "
     "mesons, at least one of which is a kaon",
     "The TwoKaonOnePionCurrent class implements the model of "
     "\\cite{Finkemeier:1995sr} for the weak current with three "
     "mesons, at least one of which is a kaon.",
     "\\bibitem{Finkemeier:1995sr}\n"
     "M.~Finkemeier and E.~Mirkes,\n"
     "Z.\\ Phys.\\  C {\\bf 69} (1996) 243 [arXiv:hep-ph/9503474].\n"
     " %%CITATION = ZEPYA,C69,243;%%\n"

);

  static Switch<TwoKaonOnePionCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &TwoKaonOnePionCurrent::_initializea1, false, false, false);
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

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &TwoKaonOnePionCurrent::_a1width, GeV, 0.599*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &TwoKaonOnePionCurrent::_a1mass, GeV, 1.251*GeV, ZERO, 10.0*GeV,
     false, false, false);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &TwoKaonOnePionCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoAxialMasses
    ("RhoAxialMasses",
     "The masses for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoAxialWidths
    ("RhoAxialWidths",
     "The widths for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoVectorMasses
    ("RhoVectorMasses",
     "The masses for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceRhoVectorWidths
    ("RhoVectorWidths",
     "The widths for the rho resonances if used local values",
     &TwoKaonOnePionCurrent::_rho2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceKstarAxialMasses
    ("KstarAxialMasses",
     "The masses for the Kstar resonances if used local values",
     &TwoKaonOnePionCurrent::_kstar1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceKstarAxialWidths
    ("KstarAxialWidths",
     "The widths for the Kstar resonances if used local values",
     &TwoKaonOnePionCurrent::_kstar1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceKstarVectorMasses
    ("KstarVectorMasses",
     "The masses for the Kstar resonances if used local values",
     &TwoKaonOnePionCurrent::_kstar2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceKstarVectorWidths
    ("KstarVectorWidths",
     "The widths for the Kstar resonances if used local values",
     &TwoKaonOnePionCurrent::_kstar2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,double> interfaceAxialRhoWeight
    ("AxialRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_rho1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,double> interfaceAxialKStarWeight
    ("AxialKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_kstar1wgts,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<TwoKaonOnePionCurrent,double> interfaceVectorRhoWeight
    ("VectorRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_rho2wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<TwoKaonOnePionCurrent,double> interfaceVectorKStarWeight
    ("VectorKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &TwoKaonOnePionCurrent::_kstar2wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values of the rho meson masses and widths",
     &TwoKaonOnePionCurrent::_rhoparameters, true, false, false);
  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the masses and widths from the particle data objects",
     false);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values of the rho meson masses and widths",
     &TwoKaonOnePionCurrent::_kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the masses and widths from the particle data objects",
     false);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfacea1Parameters
    ("a1Parameters",
     "Use local values of the rho meson masses and widths",
     &TwoKaonOnePionCurrent::_a1parameters, true, false, false);
  static SwitchOption interfacea1Parameterstrue
    (interfacea1Parameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfacea1ParametersParticleData
    (interfacea1Parameters,
     "ParticleData",
     "Use the masses and widths from the particle data objects",
     false);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfaceK1Parameters
    ("K1Parameters",
     "Use local values of the rho meson masses and widths",
     &TwoKaonOnePionCurrent::_k1parameters, true, false, false);
  static SwitchOption interfaceK1Parameterstrue
    (interfaceK1Parameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceK1ParametersParticleData
    (interfaceK1Parameters,
     "ParticleData",
     "Use the masses and widths from the particle data objects",
     false);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &TwoKaonOnePionCurrent::_a1opt, true, false, false);
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

  static ParVector<TwoKaonOnePionCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &TwoKaonOnePionCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &TwoKaonOnePionCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceK1Masses
    ("K1Masses",
     "Masses of the K_1 mesons",
     &TwoKaonOnePionCurrent::_k1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<TwoKaonOnePionCurrent,Energy> interfaceK1Widths
    ("K1Widths",
     "Widths of the K_1 mesons",
     &TwoKaonOnePionCurrent::_k1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<TwoKaonOnePionCurrent,double> interfaceK1WeightKStarPi
    ("K1WeightKStarPi",
     "The relative weights for the K_1 resonances in the K* pi final-state",
     &TwoKaonOnePionCurrent::_k1wgta, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);

  static ParVector<TwoKaonOnePionCurrent,double> interfaceK1WeightRhoK
    ("K1WeightRhoK",
     "The relative weights for the K_1 resonances in the rho K final-state",
     &TwoKaonOnePionCurrent::_k1wgtb, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,double> interfaceEpsOmega
    ("EpsOmega",
     "The omega-phi mixing ",
     &TwoKaonOnePionCurrent::_epsomega, 0.05, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &TwoKaonOnePionCurrent::_omegamass, GeV, 0.782*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &TwoKaonOnePionCurrent::_omegawidth, GeV, 0.00843*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfacePhiMass
    ("PhiMass",
     "The mass of the phi meson",
     &TwoKaonOnePionCurrent::_phimass, GeV, 1.020*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,Energy> interfacePhiWidth
    ("PhiWidth",
     "The width of the phi meson",
     &TwoKaonOnePionCurrent::_phiwidth, GeV, 0.00443*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<TwoKaonOnePionCurrent,double> interfaceOmegaKStarWeight
    ("OmegaKStarWeight",
     "The relative weight of the omega-phi and K* terms",
     &TwoKaonOnePionCurrent::_omegaKstarwgt, 1./sqrt(2.), 0.0, 100.0,
     false, false, Interface::limited);
  
  static Switch<TwoKaonOnePionCurrent,bool> interfaceOmegaParameters
    ("OmegaParameters",
     "Use local values of the omega/phi meson masses and widths",
     &TwoKaonOnePionCurrent::_omegaopt, true, false, false);
  static SwitchOption interfaceOmegaParameterstrue
    (interfaceOmegaParameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceOmegaParametersParticleData
    (interfaceOmegaParameters,
     "ParticleData",
     "Use the masses and widths from the particle data objects",
     false);

}

void TwoKaonOnePionCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return; 
    // parameters for the table of values
    Energy2 step(sqr(_maxmass)/199.);
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho(getParticleData(ParticleID::rhoplus)->mass()),
      wrho(getParticleData(ParticleID::rhoplus)->width());
    vector<Energy> inmass(2,mrho),inwidth(2,wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<TwoKaonOnePionCurrent> 
      widthgen(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi,_mpi,_mpi);
    // normalisation constant to give physical width if on shell
    double a1const(_a1width/(widthgen.partialWidth(sqr(_a1mass))));
    // loop to give the values
    _a1runq2.clear();_a1runwidth.clear();
    for(Energy2 moff2 = ZERO; moff2<=sqr(_maxmass); moff2+=step) {
      _a1runwidth.push_back(widthgen.partialWidth(moff2)*a1const);
      _a1runq2.push_back(moff2);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}
  
// modes handled by this class
bool TwoKaonOnePionCurrent::acceptMode(int imode) const { 
  return imode>=2&&imode<=11&&imode!=8;
}



// complete the construction of the decay mode for integration
bool TwoKaonOnePionCurrent::createMode(int icharge, tcPDPtr resonance,
				       IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  if(abs(icharge)!=3) return false;
  int iq(0),ia(0);
  if(!acceptMode(imode)) return false;
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1    = getParticleData(ParticleID::a_1minus);
  tPDPtr k1[2] = {getParticleData(ParticleID::K_1minus),
		  getParticleData(ParticleID::Kstar_1minus)};
  _maxmass=max(_maxmass,upp);
  // the rho0 resonances
  tPDPtr rho0[3]  ={getParticleData( 113),getParticleData( 100113),
		    getParticleData( 30113)};
  // the charged rho resonances
  tPDPtr rhoc[3]  ={getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the K*0 resonances
  tPDPtr Kstar0[3]={getParticleData( 313),getParticleData( 100313),
		    getParticleData( 30313)};
  // the charged K* resonances
  tPDPtr Kstarc[3]={getParticleData(-323),getParticleData(-100323),
		    getParticleData(-30323)};
  if(icharge==3) {
    a1    = a1->CC();
    k1[0] = k1[0]->CC();
    k1[1] = k1[1]->CC();
    for(unsigned int ix=0;ix<3;++ix) {
      if(rhoc[ix]) rhoc[ix]=rhoc[ix]->CC();
      if(Kstar0[ix]) Kstar0[ix]=Kstar0[ix]->CC();
      if(Kstarc[ix]) Kstarc[ix]=Kstarc[ix]->CC();
    }
  }
  if(imode==2) {
    // channels for K- pi- K+
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstar0[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
      if(rho0[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstar0[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	}
      }
    }
  }
  else if(imode==3) {
    // channels for K0 pi- K0bar
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
      if(rho0[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstarc[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	}
      }
    }
  }
  else if(imode==4) {
    // channels for K- pi0 K0
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstar0[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
      if(Kstarc[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      if(rhoc[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rhoc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstar0[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	}
	if(Kstarc[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			    ires+2,iloc+1,ires+2,iloc+2));
	}
      }
    }
  }
  else if(imode==5) {  
    // channels for pi0 pi0 K-
    for(unsigned int ix=0;ix<3;++ix) {
      if(!Kstarc[ix]) continue;
      for(unsigned int ik=0;ik<2;++ik) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!Kstarc[iy]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==6) {
    // channels for K- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(rho0[ix])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,rho0[ix],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	if(Kstar0[ix])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstar0[ix],ires+1,iloc+2,
			    ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!Kstarc[ix]) continue;
	if(rho0[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,rho0[ix],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	}
	if(Kstar0[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstar0[ix],ires+1,iloc+2,
			    ires+2,iloc+1,ires+2,iloc+3));
	}
      }
    }
  }
  else if(imode==7) {
    // channels for pi- kbar0 pi0
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(rhoc[ix])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,rhoc[ix],ires+1,iloc+2,
			    ires+2,iloc+1,ires+2,iloc+3));
	if(Kstar0[ix])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstar0[ix],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	if(Kstarc[ix])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+3,
			    ires+2,iloc+1,ires+2,iloc+2));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!Kstarc[ix]) continue;
	if(Kstar0[iy])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	if(rhoc[iy])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,  rhoc[iy],ires+1,iloc+2,
			    ires+2,iloc+1,ires+2,iloc+3));
	if(Kstarc[iy])
	  mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			    ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  else if(imode==9||imode==10) {
    // channels for K_S0 pi- K_S0 and K_L0 pi- K_L0 
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstarc[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			    ires+2,iloc+1,ires+2,iloc+2));
	}
      }
    }
  }
  else if(imode==11) {
    // channels for K_S0 pi- K_L0
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      if(rho0[ix])
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstarc[iy]) {
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[ix],ires+1,iloc+1,
			    ires+2,iloc+2,ires+2,iloc+3));
	  mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[ix],ires+1,iloc+3,
			    ires+2,iloc+1,ires+2,iloc+2));
	}
      }
    }
  }
  if(_rhoparameters) {
    for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
      if(rhoc[ix]) mode->resetIntermediate(rhoc[ix],_rho1mass[ix],
					   _rho1width[ix]);
      if(rho0[ix]) mode->resetIntermediate(rho0[ix],_rho1mass[ix],
					   _rho1width[ix]);
    }
  }
  // K star parameters in the base class
  if(_kstarparameters) {
    for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
      if(Kstarc[ix]) mode->resetIntermediate(Kstarc[ix],_kstar1mass[ix],
					     _kstar1width[ix]);
      if(Kstar0[ix]) mode->resetIntermediate(Kstar0[ix],_kstar1mass[ix],
					     _kstar1width[ix]);
    }
  }
  return true;
}

void TwoKaonOnePionCurrent::dataBaseOutput(ofstream & os,
					   bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::TwoKaonOnePionCurrent " 
		<< name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rho1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
    else {
      os << "insert " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rho2wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":VectorRhoWeight " << ix 
	 << " " << _rho2wgts[ix] << "\n";
    }
    else {
      os << "insert " << name() << ":VectorRhoWeight " << ix 
	 << " " << _rho2wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar2wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":VectorKStarWeight " << ix 
	 << " " << _kstar2wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":VectorKStarWeight " << ix 
	 << " " << _kstar2wgts[ix] << "\n";
    }
  }
  os << "newdef " << name() << ":OmegaKStarWeight " << _omegaKstarwgt << "\n";
  os << "newdef " << name() << ":EpsOmega " << _epsomega << "\n";
  os << "newdef " << name() << ":Initializea1 " << _initializea1 << "\n";
  os << "newdef " << name() << ":RhoParameters " << _rhoparameters << "\n";
  os << "newdef " << name() << ":KstarParameters " << _kstarparameters << "\n";
  os << "newdef " << name() << ":a1Parameters " << _a1parameters << "\n";
  os << "newdef " << name() << ":K1Parameters " << _k1parameters << "\n";
  os << "newdef " << name() << ":OmegaParameters " << _omegaopt << "\n";
  os << "newdef " << name() << ":a1WidthOption " << _a1opt << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    os << "newdef " << name() << ":a1RunningWidth " << ix 
	   << " " << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    os << "newdef " << name() << ":a1RunningQ2 " << ix 
	   << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  os << "newdef " << name() << ":A1Width " << _a1width/GeV << "\n";
  os << "newdef " << name() << ":A1Mass " << _a1mass/GeV << "\n";
  os << "newdef " << name() << ":OmegaWidth " << _omegawidth/GeV << "\n";
  os << "newdef " << name() << ":OmegaMass " << _omegamass/GeV << "\n";
  os << "newdef " << name() << ":PhiWidth " << _phiwidth/GeV << "\n";
  os << "newdef " << name() << ":PhiMass " << _phimass/GeV << "\n";
  os << "newdef " << name() << ":FPi " << _fpi/MeV << "\n";
  for(unsigned int ix=0;ix<_k1mass.size();++ix) {
    if(ix<2) {
      os << "newdef " << name() << ":K1Masses " << ix 
	 << " " << _k1mass[ix]/GeV << "\n";
    }
    else {
      os << "insert " << name() << ":K1Masses " << ix 
	 << " " << _k1mass[ix]/GeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_k1width.size();++ix) {
    if(ix<2) {
      os << "newdef " << name() << ":K1Widths " << ix 
	 << " " << _k1width[ix]/GeV << "\n";
    }
    else {
      os << "insert " << name() << ":K1Widths " << ix 
	 << " " << _k1width[ix]/GeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialMasses " << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": RhoAxialMasses" << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho2mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoVectorMasses " << ix 
		<< " " << _rho2mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": RhoVectorMasses" << ix 
		<< " " << _rho2mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho2width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoVectorWidths " << ix 
		    << " " << _rho2width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":RhoVectorWidths " << ix 
		    << " " << _rho2width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialMasses " << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarAxialMasses" << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar2mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarVectorMasses " << ix 
		<< " " << _kstar2mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarVectorMasses" << ix 
		<< " " << _kstar2mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar2width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarVectorWidths " << ix 
		    << " " << _kstar2width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarVectorWidths " << ix 
		    << " " << _kstar2width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_k1wgta.size();++ix) {
    if(ix<2) os << "newdef " << name() << ":K1WeightKStarPi " << ix
		<< " " << _k1wgta[ix] << "\n";
    else     os << "insert " << name() << ":K1WeightKStarPi " << ix
		<< " " << _k1wgta[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_k1wgtb.size();++ix) {
    if(ix<2) os << "newdef " << name() << ":K1WeightRhoK " << ix
		<< " " << _k1wgtb[ix] << "\n";
    else     os << "insert " << name() << ":K1WeightRhoK " << ix
		<< " " << _k1wgtb[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(os,false,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}  

void TwoKaonOnePionCurrent::doinit() {
  WeakCurrent::doinit();
  // the particles we will use a lot
  tPDPtr a1(getParticleData(ParticleID::a_1minus)),
    pi0(getParticleData(ParticleID::pi0)),
    piplus(getParticleData(ParticleID::piplus)),
    piminus(getParticleData(ParticleID::piminus));
  tPDPtr k1[2]={getParticleData(ParticleID::K_1minus),
		getParticleData(ParticleID::Kstar_1minus)};
  // masses for the running widths
  _mpi=piplus->mass();
  _mK=getParticleData(ParticleID::K0)->mass();
  // the charged rho resonances
  tPDPtr rhoc[3]={getParticleData(-213),getParticleData(-100213),
 		  getParticleData(-30213)};
  // the charged K* resonances
  tPDPtr Kstarc[3]={getParticleData(-323),getParticleData(-100323),
 		    getParticleData(-30323)};
  if(!_a1parameters) {
    _a1mass=a1->mass();
    _a1width=a1->width();
  }
  // mass and width of the k_1
  if(!_k1parameters) {
    for(unsigned int ix=0;ix<2;++ix) {
      _k1mass[ix]  = k1[ix]->mass();
      _k1width[ix] = k1[ix]->width();
    }
  }
  // initialise the a_1 running width calculation
  inita1Width(-1);
  // rho parameters in the base classs
  tcPDPtr temp;
  unsigned int ix;
  if(_rhoparameters&&_rho1mass.size()<3) {
    ix = _rho1mass.size();
    _rho1mass.resize(3);
    _rho1width.resize(3);
    for(;ix<3;++ix) {
      if(rhoc[ix]) {
	_rho1mass [ix]=rhoc[ix]->mass();
	_rho1width[ix]=rhoc[ix]->width();
      }
    }
  }
  else if(!_rhoparameters) {
    _rho1mass.resize(3);_rho1width.resize(3);
    for(ix=0;ix<3;++ix) {
      if(rhoc[ix]) {
	_rho1mass[ix]=rhoc[ix]->mass();
	_rho1width[ix]=rhoc[ix]->width();
      }
    }
  }
  // K star parameters in the base class
  if(_kstarparameters&&_kstar1mass.size()<3) {
    ix = _kstar1mass.size();
    _kstar1mass.resize(3);_kstar1width.resize(3);
    for(;ix<3;++ix) {
      if(Kstarc[ix]) {
	_kstar1mass[ix]=Kstarc[ix]->mass();
	_kstar1width[ix]=Kstarc[ix]->width();
      }
    }
  }
  else if(!_kstarparameters) {
    _kstar1mass.resize(3);_kstar1width.resize(3);
    for(ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	_kstar1mass[ix]=Kstarc[ix]->mass();
	_kstar1width[ix]=Kstarc[ix]->width();
      }
    }
  }
  // rho parameters here
  if(_rhoparameters&&_rho2mass.size()<3) {
    ix = _rho2mass.size();
    _rho2mass.resize(3);_rho2width.resize(3);
    for(;ix<3;++ix) {
      if(rhoc[ix]) {
	_rho2mass[ix]=rhoc[ix]->mass();
	_rho2width[ix]=rhoc[ix]->width();
      }
    }
  }
  else if(!_rhoparameters) {
    _rho2mass.resize(3);_rho2width.resize(3);
    for(ix=0;ix<3;++ix) {
      if(rhoc[ix]) {
	_rho2mass[ix]=rhoc[ix]->mass();
	_rho2width[ix]=rhoc[ix]->width();
      }
    }
  }
  // Kstar parameters here
  if(_kstarparameters&&_kstar2width.size()<3) {
    ix = _kstar2mass.size();
    _kstar2mass.resize(3);_kstar2width.resize(3);
    for(;ix<3;++ix) {
      if(Kstarc[ix]) {
	_kstar2mass[ix]=Kstarc[ix]->mass();
	_kstar2width[ix]=Kstarc[ix]->width();
      }
    }
  }
  else if(!_kstarparameters) {
    _kstar2mass.resize(3);_kstar2width.resize(3);
    for(ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	_kstar2mass[ix]=Kstarc[ix]->mass();
	_kstar2width[ix]=Kstarc[ix]->width();
      }
    }
  }
  inita1Width(0);
}

TwoKaonOnePionCurrent::FormFactors
TwoKaonOnePionCurrent::calculateFormFactors(const int ichan,const int imode,
						 Energy2 q2,Energy2 s1,
						 Energy2 s2,Energy2 s3) const {
  useMe();
  Complex F1, F2, F5;
  // calculate the K- pi - K+ factor
  if(imode==2) {
    Complex a1fact(a1BreitWigner(q2)*sqrt(2.)/3.);
    if(ichan<0) {
      F1 = -a1fact*TKstar1(s1,-1);
      F2 =  a1fact*Trho1(s2,-1);
      F5 = Trho2(q2,-1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
    }
    else if(ichan%5==0) F1 = -a1fact*TKstar1(s1,    ichan/5);
    else if(ichan%5==1) F2 =  a1fact*Trho1(  s2,(ichan-1)/5);
    else if(ichan%5>=2) F5 = Trho2(q2,ichan/5)*TOmegaKStar(s1,s2,2*((ichan-2)%5))
      *sqrt(2.);
  }
  // calculate the K0 pi- K0bar
  else if(imode==3) {
    Complex a1fact(a1BreitWigner(q2)*sqrt(2.)/3.);
    if(ichan<0) {
      F1 =-a1fact*TKstar1(s1,-1);
      F2 = a1fact*Trho1  (s2,-1);
      F5 =-Trho2(q2,-1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
    }
    else if(ichan%5==0) F1 = -a1fact*TKstar1(s1,    ichan/5);
    else if(ichan%5==1) F2 =  a1fact*Trho1  (s2,(ichan-1)/5);
    else if(ichan%5>=2) F5 = -Trho2(q2,ichan/5)*TOmegaKStar(s1,s2,2*((ichan-2)%5))
      *sqrt(2.);
  }
  // calculate the K- pi0 k0
  else if(imode==4) {
    Complex a1fact(a1BreitWigner(q2)/3.);
    if(ichan<0) {
      F1 =  a1fact*( TKstar1(s1,-1)-TKstar1(s3,-1));
      F2 = -a1fact*(2.*Trho1(s2,-1)+TKstar1(s3,-1));
      F5 = Trho2(q2,-1)*(TKstar1(s3,-1)-TKstar1(s1,-1))/(1.+_omegaKstarwgt)/sqrt(2.);
    }
    else if(ichan%9==0) F1 =  a1fact*TKstar1(s1,ichan/9)/3.;
    else if(ichan%9==1) {
      F1 = +a1fact*TKstar1(s3,(ichan-1)/9)/3.;
      F2 = -a1fact*TKstar1(s3,(ichan-1)/9)/3.;
    }
    else if(ichan%9==2) F2 = -a1fact*2.*Trho1(s2,(ichan-2)/9)/3.;
    else if(ichan%9<6)  F5 =-Trho2(q2,ichan/9)*TKstar1(s1,(ichan-3)%9)
      /(1.+_omegaKstarwgt)/sqrt(2.);
    else                F5 = Trho2(q2,ichan/9)*TKstar1(s3,(ichan-6)%9)
      /(1.+_omegaKstarwgt)/sqrt(2.);
  }
  // calculate the pi0 pi0 K-
  else if(imode==5) {
    if(ichan<0) {
      Complex K1fact(TK1(q2,0,-1)/6.);
      F1 = K1fact*TKstar1(s1,-1);
      F2 =-K1fact*TKstar1(s2,-1);
      F5 =-0.25*TKstar2(q2,-1)*(TKstar1(s1,-1)-TKstar1(s2,-1));
    }
    else if(ichan%10==0) F1=  TK1(q2,0,0)/6.*TKstar1(s1,ichan/10);
    else if(ichan%10==1) F2= -TK1(q2,0,0)/6.*TKstar1(s2,ichan/10);
    else if(ichan%10==2) F1=  TK1(q2,0,1)/6.*TKstar1(s1,ichan/10);
    else if(ichan%10==3) F2= -TK1(q2,0,1)/6.*TKstar1(s2,ichan/10);
    else if(ichan%10<7 ) F5 =-sqrt(2.)/4*TKstar2(q2,ichan/10)*TKstar1(s1,(ichan-4)%10);
    else                 F5 = sqrt(2.)/4*TKstar2(q2,ichan/10)*TKstar1(s2,(ichan-7)%10);
  }
  // calculate the K- pi- pi+
  else if(imode==6) {
    double fact=sqrt(2.)/3.;
    if(ichan<0) {
      F1 = -fact*TK1(q2,1,-1)*Trho1(s1,-1);
      F2 =  fact*TK1(q2,0,-1)*TKstar1(s2,-1);
      F5 = -sqrt(0.5)*TKstar2(q2,-1)*(Trho1(s1,-1)+TKstar1(s2,-1));
    }
    else if(ichan%10==0) F1 = -fact*TK1(q2,1,0)*Trho1  (s1,ichan/10);
    else if(ichan%10==1) F2 =  fact*TK1(q2,0,0)*TKstar1(s2,ichan/10);
    else if(ichan%10==2) F1 = -fact*TK1(q2,1,1)*Trho1(  s1,ichan/10);
    else if(ichan%10==3) F2 =  fact*TK1(q2,0,1)*TKstar1(s2,ichan/10);
    else if(ichan%10<7)  F5 = -sqrt(0.5)*TKstar2(q2,ichan/10)*Trho1(  s1,(ichan-4)%10);
    else                 F5 = -sqrt(0.5)*TKstar2(q2,ichan/10)*TKstar1(s2,(ichan-7)%10);
  }
  // calculate the pi- K0bar pi0
  else if(imode==7) {
    if(ichan<0) {
      Complex K1facta(TK1(q2,0,-1)),K1factb(TK1(q2,1,-1));
      F1 = K1facta*(TKstar1(s1,-1)-TKstar1(s3,-1))/3.;
      F2 =-(2.*K1factb*Trho1(s2,-1)+K1facta*TKstar1(s3,-1))/3.;
      F5 = -0.5*TKstar2(q2,-1)*(2.*Trho1(s2,-1)+TKstar1(s1,-1)+TKstar1(s3,-1));
    }
    else if(ichan%15==0) F2 =-2.*TK1(q2,0,0)*Trho1  (s2,ichan/15)/3.;
    else if(ichan%15==1) F1 =    TK1(q2,1,0)*TKstar1(s1,ichan/15)/3.;
    else if(ichan%15==2) {
      F1 =-TK1(q2,1,0)*TKstar1(s3,ichan/15)/3.;
      F2 =-TK1(q2,1,0)*TKstar1(s3,ichan/15)/3.;
    }
    else if(ichan%15==3) F2 =-2.*TK1(q2,0,1)*Trho1  (s2,ichan/15)/3.;
    else if(ichan%15==4) F1 =    TK1(q2,1,1)*TKstar1(s1,ichan/15)/3.;
    else if(ichan%15==5) {
      F1 =-TK1(q2,1,1)*TKstar1(s3,ichan/15)/3.;
      F2 =-TK1(q2,1,1)*TKstar1(s3,ichan/15)/3.;
    }
    else if(ichan%15<9 ) F5 = -0.5*TKstar2(q2,ichan/15)*TKstar1(s1,(ichan- 6)%15);
    else if(ichan%15<12) F5 = -    TKstar2(q2,ichan/15)*Trho1  (s2,(ichan- 9)%15);
    else                 F5 = -0.5*TKstar2(q2,ichan/15)*TKstar1(s3,(ichan-12)%15);
  }
  // calculate the K_S0 pi- K_S0 or K_L0 pi- K_L0
  else if(imode==9||imode==10) {
    Complex a1fact(a1BreitWigner(q2)/6.);
    if(ichan<0) {
      F1 = a1fact*(TKstar1(s1,-1)+TKstar1(s3,-1));
      F2 = a1fact*TKstar1(s3,-1);
      F5 = 0.5*Trho2(q2,-1)*(TOmegaKStar(s1,s2,-1)-TOmegaKStar(s3,s2,-1));
    }
    else if(ichan%8==0) F1=a1fact*TKstar1(s1,ichan/8);
    else if(ichan%8==1) {
      F1 = a1fact*TKstar1(s3,ichan/8);
      F2 = a1fact*TKstar1(s3,ichan/8);
    }
    else if(ichan%8<5 ) F5 = -Trho2(q2,ichan/8)*TKstar1(s1,(ichan-2)%8)
      /(1.+_omegaKstarwgt)/2.;
    else                F5 =  Trho2(q2,ichan/8)*TKstar1(s3,(ichan-5)%8)
      /(1.+_omegaKstarwgt)/2.;
  }
  else if(imode==11) {
    Complex a1fact(a1BreitWigner(q2)/3./sqrt(2.));
    if(ichan<0) {
      F1 = -a1fact*(TKstar1(s1,-1)-TKstar1(s3,-1));
      F2 =  a1fact*(2.*Trho1(s2,-1)+TKstar1(s3,-1));
      F5 = -Trho2(q2,-1)*(TOmegaKStar(s1,s2,-1)+TOmegaKStar(s3,s2,-1))/sqrt(2.);
    }
    else if(ichan%9==0) F1 =-   a1fact*TKstar1(s1,ichan/9);
    else if(ichan%9==1) {
      F1 = a1fact*TKstar1(s3,ichan/9);
      F2 = a1fact*TKstar1(s3,ichan/9);
    }
    else if(ichan%9==2) F2 = 2.*a1fact*Trho1(  s2,ichan/9);
    else if(ichan%9<6 ) F5 = -sqrt(0.5)*TKstar2(q2,ichan/9)*
      TOmegaKStar(s1,s2,2*((ichan-3)%9))/sqrt(2.);
    else                F5 = -sqrt(0.5)*TKstar2(q2,ichan/9)*
      TOmegaKStar(s3,s2,2*((ichan-6)%9))/sqrt(2.);
  }
  return FormFactors(F1 / _fpi,
		     F2 / _fpi,
		     InvEnergy(),
		     InvEnergy(),
		     -F5 / sqr(Constants::twopi) / pow<3,1>(_fpi));
}

void TwoKaonOnePionCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  WeakCurrent::doinitrun();
}

void TwoKaonOnePionCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

double TwoKaonOnePionCurrent::
threeBodyMatrixElement(const int       , const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy    , 
		       const Energy    , const Energy    ) const {
  Energy2 mpi2(sqr(_mpi));
  Complex propb(Trho1(s1,-1)),propa(Trho1(s2,-1)); 
  // the matrix element
  Energy2 output(ZERO); 
  // first resonance
  output+= ((s1-4.*mpi2)+0.25*(s3-s2)*(s3-s2)/q2)*real(propb*conj(propb)); 
  // second resonance
  output+= ((s2-4.*mpi2)+0.25*(s3-s1)*(s3-s1)/q2)*real(propa*conj(propa)); 
  // the interference term 
  output+= (0.5*q2-s3-0.5*mpi2+0.25*(s3-s2)*(s3-s1)/q2)*real(propa*conj(propb)+
							     propb*conj(propa)); 
  return output / sqr(_rho1mass[0]);
}

Complex TwoKaonOnePionCurrent::Trho1(Energy2 q2,int ires) const {
  Complex output(0.);
  double norm(0.);
  for(unsigned int ix=0,N=_rho1wgts.size();ix<N;++ix) norm+=_rho1wgts[ix];
  if(ires<0) {
    for(unsigned int ix=0,N=_rho1wgts.size();ix<N;++ix) {
      output+=_rho1wgts[ix]*BWrho1(q2,ix);
    }
  }
  else {
    unsigned int temp(ires);
    if(temp<_rho1wgts.size()) output=_rho1wgts[temp]*BWrho1(q2,temp);
  }
  return output/norm;
}
  
Complex TwoKaonOnePionCurrent::Trho2(Energy2 q2,int ires) const {
  Complex output(0.);
  double norm(0.);
  for(unsigned int ix=0,N=_rho2wgts.size();ix<N;++ix) norm+=_rho2wgts[ix];
  if(ires<0) {
    for(unsigned int ix=0,N=_rho2wgts.size();ix<N;++ix) {
      output+=_rho2wgts[ix]*BWrho2(q2,ix);
    }
  }
  else {
    unsigned int temp(ires);
    if(temp<_rho2wgts.size()) output=_rho2wgts[temp]*BWrho2(q2,temp);
  }
  return output/norm;
}
  
Complex TwoKaonOnePionCurrent::TKstar1(Energy2 q2,int ires) const  {
  Complex output(0.);
  double norm(0.);
  for(unsigned int ix=0,N=_kstar1wgts.size();ix<N;++ix) norm+=_kstar1wgts[ix];
  if(ires<0) {
    for(unsigned int ix=0,N=_kstar1wgts.size();ix<N;++ix) {
      output+=_kstar1wgts[ix]*BWKstar1(q2,ix);
    }
  }
  else {
    unsigned int temp(ires);
    if(temp<_kstar1wgts.size()) output=_kstar1wgts[temp]*BWKstar1(q2,temp);
  }
  return output/norm;
}
  
Complex TwoKaonOnePionCurrent::TKstar2(Energy2 q2,int ires) const {
  Complex output(0.);
  double norm(0.);
  for(unsigned int ix=0,N=_kstar2wgts.size();ix<N;++ix) norm+=_kstar2wgts[ix];
  if(ires<0) {
    for(unsigned int ix=0,N=_kstar2wgts.size();ix<N;++ix) {
      output+=_kstar2wgts[ix]*BWKstar2(q2,ix);
    }
  }
  else {
    unsigned int temp(ires);
    if(temp<_kstar2wgts.size()) output=_kstar2wgts[temp]*BWKstar2(q2,temp);
  }
  return output/norm;
}

Complex TwoKaonOnePionCurrent::BWrho1(Energy2 q2, unsigned int ires) const {
  if(ires>=_rho1mass.size()) return 0.;
  Energy mass  = _rho1mass [ires];
  Energy width = _rho1width[ires];
  Energy q=sqrt(q2);
  Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mpi,_mpi);
  Energy pcm  = q<=2.*_mpi ? ZERO : Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  double ratio = Math::powi(pcm/pcm0, 3);
  Energy gam(width*mass*ratio/q);
  return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
} 

Complex TwoKaonOnePionCurrent::BWrho2(Energy2 q2, unsigned int ires) const {
  if(ires>=_rho2mass.size()) return 0.;
  Energy mass  = _rho2mass [ires];
  Energy width = _rho2width[ires];
  Energy q=sqrt(q2);
  Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mpi,_mpi);
  Energy pcm  = q<=2.*_mpi ? ZERO : Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  double ratio(pcm/pcm0);ratio*=ratio*ratio;
  Energy gam(width*mass*ratio/q);
  return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
}
  
Complex TwoKaonOnePionCurrent::BWKstar1(Energy2 q2, unsigned int ires) const {
  if(ires>=_kstar1mass.size()) return 0.;
  Energy mass  = _kstar1mass [ires];
  Energy width = _kstar1width[ires];
  Energy q=sqrt(q2);
  Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mK,_mpi);
  Energy pcm  = q<=_mpi+_mK ? ZERO : Kinematics::pstarTwoBodyDecay(q,_mK,_mpi);
  double ratio(pcm/pcm0);ratio*=ratio*ratio;
  Energy gam(width*mass*ratio/q);
  return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
}

Complex TwoKaonOnePionCurrent::BWKstar2(Energy2 q2, unsigned int ires) const  {
  if(ires>=_kstar2mass.size()) return 0.;
  Energy mass  = _kstar2mass [ires];
  Energy width = _kstar2width[ires];
  Energy q=sqrt(q2);
  Energy pcm0 = Kinematics::pstarTwoBodyDecay(mass,_mK,_mpi);
  Energy pcm  = q<=_mpi+_mK ? ZERO : Kinematics::pstarTwoBodyDecay(q,_mK,_mpi);
  double ratio(pcm/pcm0);ratio*=ratio*ratio;
  Energy gam(width*mass*ratio/q);
  return sqr(mass)/(sqr(mass)-q2-Complex(0.,1.)*mass*gam);
}

Complex TwoKaonOnePionCurrent::a1BreitWigner(Energy2 q2) const {
  Complex ii(0.,1.);
  Energy2 m2(_a1mass*_a1mass);
  Energy  q(sqrt(q2));
  return m2/(m2-q2-ii*q*a1Width(q2));
}
  
Complex TwoKaonOnePionCurrent::TK1(Energy2 q2,unsigned int iopt,int ires) const {
  Complex denom(0),num(0.);
  if(iopt==0) {
    for(unsigned int ix=0;ix<_k1wgta.size();++ix) denom+=_k1wgta[ix];
    if(ires==-1) {
      for(unsigned int ix=0;ix<_k1wgta.size();++ix) 
	num+=_k1wgta[ix]*K1BreitWigner(q2,ix);
    }
    else {	
      num+=_k1wgta[ires]*K1BreitWigner(q2,ires);
    }
  }
  else if(iopt==1) {
    for(unsigned int ix=0;ix<_k1wgtb.size();++ix) denom+=_k1wgtb[ix];
    if(ires==-1) {
      for(unsigned int ix=0;ix<_k1wgtb.size();++ix) 
	num+=_k1wgtb[ix]*K1BreitWigner(q2,ix);
    }
    else {	
      num+=_k1wgtb[ires]*K1BreitWigner(q2,ires);
    }
  }
  else {
    return 0.;
  }
  return num/denom;
}

Complex TwoKaonOnePionCurrent::K1BreitWigner(Energy2 q2,unsigned int ires) const {
  if(ires>=_k1mass.size()) return 0.;
  Energy2 m2=sqr(_k1mass[ires]),mg=_k1mass[ires]*_k1width[ires];
  return (-m2+Complex(0.,1.)*mg)/(q2-m2+Complex(0.,1.)*mg);
}

Energy TwoKaonOnePionCurrent::a1Width(Energy2 q2) const {
  if(!_a1opt) return _a1mass*_a1width*g(q2)/g(_a1mass*_a1mass)/sqrt(q2);
  else        return (*_a1runinter)(q2);
}
  
double TwoKaonOnePionCurrent::g(Energy2 q2) const {
  double output;
  if(q2<9.*_mpi*_mpi) {
    output=0.;
  }
  else if(q2<sqr(_rho1mass[0]+_mpi)) {
    double diff=(q2-9.*_mpi*_mpi)/GeV2;
      
    output=4.1*sqr(diff)*diff*(1.-3.3*diff+5.8*sqr(diff));
  }
  else {
    double ratio = q2/GeV2;
    output = ratio*(1.623+10.38/ratio-9.32/sqr(ratio)+0.65/(ratio*sqr(ratio)));
  }
  return output;
}
  
Complex TwoKaonOnePionCurrent::Tomega(Energy2 q2, int ires) const {
  double denom=(1.+_epsomega);
  Complex num(0.);
  if(ires<0) num=OmegaPhiBreitWigner(q2,0)+_epsomega*OmegaPhiBreitWigner(q2,1);
  else if(ires==0) num=OmegaPhiBreitWigner(q2,0);
  else             num=OmegaPhiBreitWigner(q2,1);
  return num/denom;
}

Complex TwoKaonOnePionCurrent::OmegaPhiBreitWigner(Energy2 q2, unsigned int ires) const {
  Energy2 m2,mg;
  if(ires==0) {
    m2=sqr(_omegamass);
    mg=_omegamass*_omegawidth;
  }
  else {
    m2=sqr(_phimass);
    mg=_phimass*_phiwidth;
  }
  return (-m2+Complex(0.,1.)*mg)/(q2-m2+Complex(0.,1.)*mg);
}

Complex TwoKaonOnePionCurrent::TOmegaKStar(Energy2 s1,Energy2 s2,int ires) const {
  Complex output;
  if(ires<0)         output = _omegaKstarwgt*TKstar1(s1,-1)+Tomega(s2,-1);
  else if(ires%2==0) output = _omegaKstarwgt*TKstar1(s1,ires/2);
  else if(ires%2==1) output = Tomega(s2,ires/2);
  return output/(1.+_omegaKstarwgt);
}


// the hadronic currents    
vector<LorentzPolarizationVectorE> 
TwoKaonOnePionCurrent::current(tcPDPtr resonance,
			      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & ,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  FormFactors F = calculateFormFactors(ichan,imode,q2,s1,s2,s3);
  //if(inpart.id()==ParticleID::tauplus){F.F5=conj(F.F5);}
  // the first three form-factors
  LorentzPolarizationVector vect;
  vect = (F.F2-F.F1)*momenta[2]
        +(F.F1-F.F3)*momenta[1]
        +(F.F3-F.F2)*momenta[0];
  // multiply by the transverse projection operator
  complex<InvEnergy> dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect += (F.F4-dot)*q;
  if(F.F5!=complex<InvEnergy3>()) 
    vect += Complex(0.,1.)*F.F5*Helicity::epsilon(momenta[0],
						  momenta[1],
						  momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()*vect);
}

bool TwoKaonOnePionCurrent::accept(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nk0(0),nk0bar(0),neta(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::eta)     ++neta;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  int imode(-1);
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) imode= 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) imode= 1;
  else if( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))                 imode= 2;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))              imode= 3;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))                 imode= 4;
  else if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   imode= 5;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )               imode= 6;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))             imode= 7;
  else if( (npip==1&&npi0==1&&neta==1) ||
	   (npim==1&&npi0==1&&neta==1))               imode= 8;
  else if( nks==2 && (npip==1||npim==1) )             imode= 9;
  else if( nkl==2 && (npip==1||npim==1) )             imode=10;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) )     imode=11;
  return imode==-1 ? false : acceptMode(imode);
}

unsigned int TwoKaonOnePionCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nk0(0),nk0bar(0),neta(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::eta)     ++neta;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  int imode(-1);
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) imode= 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) imode= 1;
  else if( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))                 imode= 2;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))              imode= 3;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))                 imode= 4;
  else if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   imode= 5;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )               imode= 6;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))             imode= 7;
  else if( (npip==1&&npi0==1&&neta==1) ||
	   (npim==1&&npi0==1&&neta==1))               imode= 8;
  else if( nks==2 && (npip==1||npim==1) )             imode= 9;
  else if( nkl==2 && (npip==1||npim==1) )             imode=10;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) )     imode=11;
  return imode;
}


tPDVector TwoKaonOnePionCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::piminus);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kplus);
  }
  else if(imode==3) {
    extpart[0]=getParticleData(ParticleID::K0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==4) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::K0);
  }
  else if(imode==5) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::Kminus);
  }
  else if(imode==6) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==7) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::Kbar0);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else if(imode==8) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::eta);
  }
  else if(imode==9) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_S0);
  }
  else if(imode==10) {
    extpart[0]=getParticleData(ParticleID::K_L0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  else if(imode==11) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}

