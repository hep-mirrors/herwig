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
  addDecayMode(2,-1);
  setInitialModes(7);
  // rho parameters
  // rho parameters for axial-vector pieces
  _rho1wgts  = {1.0,-0.145,0.};
  _rho1mass  = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rho1width = {0.145*GeV,0.510*GeV,0.120*GeV};
  // rho parameters for vector pieces
  _rho2wgts  = {1.0,-0.25,-0.038};
  _rho2mass  = {0.773*GeV,1.500*GeV,1.750*GeV};
  _rho2width = {0.145*GeV,0.220*GeV,0.120*GeV};
  // K* parameters
  // K* parameters for the axial-vector pieces
  _kstar1wgts  = {1.0,-0.135,0.};
  _kstar1mass  = {0.892*GeV,1.412*GeV,1.714*GeV};
  _kstar1width = {0.050*GeV,0.227*GeV,0.323*GeV};
  // a_1 parameters
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
  // parameters for the T_omega function
  _epsomega   = 0.05;
  _omegamass  = 0.782*GeV;
  _omegawidth = 0.00843*GeV;
  _phimass    = 1.020*GeV;
  _phiwidth   = 0.00443*GeV;
  _omegaKstarwgt=1./sqrt(2.);
  // the pion decay constant
  _fpi     = 130.7*MeV/sqrt(2.);
  _mpi     = ZERO;
  _mK      = ZERO;
  _maxmass = ZERO;
  _maxcalc = ZERO;
}


void TwoKaonOnePionCurrent::persistentOutput(PersistentOStream & os) const {
  os << _a1runinter
     << _rho1wgts << ounit(_rho1mass,GeV) << ounit(_rho1width,GeV) 
     << _rho2wgts << ounit(_rho2mass,GeV) << ounit(_rho2width,GeV)
     << _kstar1wgts << ounit(_kstar1mass,GeV) << ounit(_kstar1width,GeV)
     << ounit(_a1mass,GeV) << ounit(_a1width,GeV)
     << ounit(_a1runwidth,GeV) << ounit(_a1runq2,GeV2) << _epsomega 
     << ounit(_omegamass,GeV) << ounit(_omegawidth,GeV) 
     << ounit(_phimass,GeV) << ounit(_phiwidth,GeV) << _omegaKstarwgt 
     << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV) 
     << _initializea1 << _a1opt
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV);
}

void TwoKaonOnePionCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _a1runinter
     >> _rho1wgts >> iunit(_rho1mass,GeV) >> iunit(_rho1width,GeV) 
     >> _rho2wgts >> iunit(_rho2mass,GeV) >> iunit(_rho2width,GeV) 
     >> _kstar1wgts >> iunit(_kstar1mass,GeV) >> iunit(_kstar1width,GeV)
     >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV)
     >> iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) >> _epsomega 
     >> iunit(_omegamass,GeV) >> iunit(_omegawidth,GeV) 
     >> iunit(_phimass,GeV) >> iunit(_phiwidth,GeV) >> _omegaKstarwgt 
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV) 
     >> _initializea1 >> _a1opt
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

// complete the construction of the decay mode for integration
bool TwoKaonOnePionCurrent::createMode(int icharge, tcPDPtr resonance,
				       IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    switch(i3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return false;
      break;
    case IsoSpin::I3One:
      if( imode>1 || icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusOne:
      if( imode>1 || icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  // get the particles and check the mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1    = getParticleData(ParticleID::a_1minus);
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
    for(unsigned int ix=0;ix<3;++ix) {
      if(rhoc[ix]) rhoc[ix]=rhoc[ix]->CC();
      if(Kstar0[ix]) Kstar0[ix]=Kstar0[ix]->CC();
      if(Kstarc[ix]) Kstarc[ix]=Kstarc[ix]->CC();
    }
  }
  if(imode==0) {
    // channels for K- pi- K+
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
    }
  }
  else if(imode==1) {
    // channels for K0 pi- K0bar
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for K- pi0 K0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rhoc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  else if(imode==3||imode==4) {
    // channels for K_S0 pi- K_S0 and K_L0 pi- K_L0 
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  else if(imode==5) {
    // channels for K_S0 pi- K_L0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==a1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
	mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,rho0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=rhoc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,rhoc[ix],ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  // update the integration parameters
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rho1mass[ix],
			    _rho1width[ix]);
    mode->resetIntermediate(rho0[ix],_rho1mass[ix],
			    _rho1width[ix]);
  }
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
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
  os << "newdef " << name() << ":OmegaKStarWeight " << _omegaKstarwgt << "\n";
  os << "newdef " << name() << ":EpsOmega " << _epsomega << "\n";
  os << "newdef " << name() << ":Initializea1 " << _initializea1 << "\n";
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
  WeakCurrent::dataBaseOutput(os,false,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}  

void TwoKaonOnePionCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::K0)    ->mass();
  // initialise the a_1 running width calculation
  inita1Width(-1);
  inita1Width(0);
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
  
Complex TwoKaonOnePionCurrent::Tomega(Energy2 q2, int ires) const {
  double denom=(1.+_epsomega);
  Complex num(0.);
  if(ires<0) num=OmegaPhiBreitWigner(q2,0)+_epsomega*OmegaPhiBreitWigner(q2,1);
  else if(ires==0) num=OmegaPhiBreitWigner(q2,0);
  else             num=OmegaPhiBreitWigner(q2,1);
  return num/denom;
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
  // check the isospin
  if(Itotal!=IsoSpin::IUnknown && Itotal!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    switch(i3) {
    case IsoSpin::I3Zero:
      return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One: case IsoSpin::I3MinusOne:
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  // check the resonance
  int ires1=-1;
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      ires1=0; break;
    case 100:
      ires1=1; break;
    case  30:
      ires1=2; break;
    case  10:
      ires1=3; break;
    default:
      assert(false);
    }
  }
  useMe();
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
  // calculate the form factors
  useMe();
  Complex F1(0.), F2(0.), F5(0.);
  Complex a1fact = ires1<0 || ires1==3 ? a1BreitWigner(q2) : 0.;
  // calculate the K- pi - K+ factor
  if(imode==0) {
    a1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 = -a1fact*TKstar1(s1,-1);
      F2 =  a1fact*Trho1(s2,-1);
      if(ires1<0)
	F5 = Trho2(q2,   -1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else if(ires1<3)
	F5 = Trho2(q2,ires1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%5==0) F1 = -a1fact*TKstar1(s1,    ichan/5);
    else if(ichan%5==1) F2 =  a1fact*Trho1(  s2,(ichan-1)/5);
    else if(ichan%5>=2) F5 = Trho2(q2,ichan/5)*TOmegaKStar(s1,s2,2*((ichan-2)%5))
      *sqrt(2.);
  }
  // calculate the K0 pi- K0bar
  else if(imode==1) {
    a1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 =-a1fact*TKstar1(s1,-1);
      F2 = a1fact*Trho1  (s2,-1);
      if(ires1<0)
	F5 =-Trho2(q2,   -1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else if(ires1<3)
	F5 =-Trho2(q2,ires1)*TOmegaKStar(s1,s2,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%5==0) F1 = -a1fact*TKstar1(s1,    ichan/5);
    else if(ichan%5==1) F2 =  a1fact*Trho1  (s2,(ichan-1)/5);
    else if(ichan%5>=2) F5 = -Trho2(q2,ichan/5)*TOmegaKStar(s1,s2,2*((ichan-2)%5))
      *sqrt(2.);
  }
  // calculate the K- pi0 k0
  else if(imode==2) {
    a1fact /= 3.;
    if(ichan<0) {
      F1 =  a1fact*( TKstar1(s1,-1)-TKstar1(s3,-1));
      F2 = -a1fact*(2.*Trho1(s2,-1)+TKstar1(s3,-1));
      if(ires1<0)
	F5 = Trho2(q2,   -1)*(TKstar1(s3,-1)-TKstar1(s1,-1))/(1.+_omegaKstarwgt)/sqrt(2.);
      else if(ires1<3)
	F5 = Trho2(q2,ires1)*(TKstar1(s3,-1)-TKstar1(s1,-1))/(1.+_omegaKstarwgt)/sqrt(2.);
      else
	F5 = 0.;
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
  // calculate the K_S0 pi- K_S0 or K_L0 pi- K_L0
  else if(imode==3||imode==4) {
    a1fact /=6;
    if(ichan<0) {
      F1 = a1fact*(TKstar1(s1,-1)+TKstar1(s3,-1));
      F2 = a1fact*TKstar1(s3,-1);
      if(ires1<0)
	F5 = 0.5*Trho2(q2,   -1)*(TOmegaKStar(s1,s2,-1)-TOmegaKStar(s3,s2,-1));
      else if(ires1<3)
	F5 = 0.5*Trho2(q2,ires1)*(TOmegaKStar(s1,s2,-1)-TOmegaKStar(s3,s2,-1));
      else
	F5 = 0.;
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
  else if(imode==5) {
    a1fact *= 1./3./sqrt(2.);
    if(ichan<0) {
      F1 = -a1fact*(TKstar1(s1,-1)-TKstar1(s3,-1));
      F2 =  a1fact*(2.*Trho1(s2,-1)+TKstar1(s3,-1));
      if(ires1<0)
	F5 = -Trho2(q2,   -1)*(TOmegaKStar(s1,s2,-1)+TOmegaKStar(s3,s2,-1))/sqrt(2.);
      else if(ires1<3)
	F5 = -Trho2(q2,ires1)*(TOmegaKStar(s1,s2,-1)+TOmegaKStar(s3,s2,-1))/sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%9==0) F1 =-   a1fact*TKstar1(s1,ichan/9);
    else if(ichan%9==1) {
      F1 = a1fact*TKstar1(s3,ichan/9);
      F2 = a1fact*TKstar1(s3,ichan/9);
    }
    else if(ichan%9==2) F2 = 2.*a1fact*Trho1(  s2,ichan/9);
    else if(ichan%9<6 ) F5 = -sqrt(0.5)*Trho2(q2,ichan/9)*
      TOmegaKStar(s1,s2,2*((ichan-3)%9))/sqrt(2.);
    else                F5 = -sqrt(0.5)*Trho2(q2,ichan/9)*
      TOmegaKStar(s3,s2,2*((ichan-6)%9))/sqrt(2.);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect = (F2-F1)*momenta[2] + F1*momenta[1] - F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  if(F5!=0.) 
    vect -= Complex(0.,1.)*F5/sqr(Constants::twopi)/sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool TwoKaonOnePionCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if     ( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))             return true;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))          return true;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))             return true;
  else if( nks==2 && (npip==1||npim==1) )         return true;
  else if( nkl==2 && (npip==1||npim==1) )         return true;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) ) return true;
  return false;
}

unsigned int TwoKaonOnePionCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
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
  if     ( (nkp==1&&nkm==1&&npip==1) ||
	   (nkp==1&&nkm==1&&npim==1))                 return 0;
  else if( (nk0==1&&nk0bar==1&&npip==1) ||
	   (nk0==1&&nk0bar==1&&npim==1))              return 1;
  else if( (nkp==1&&nk0bar==1&&npi0==1) ||
	   (nkm==1&&npi0==1&&nk0==1))                 return 2;
  else if( nks==2 && (npip==1||npim==1) )             return 3;
  else if( nkl==2 && (npip==1||npim==1) )             return 4;
  else if( nks==1&&nkl==1 && (npip==1||npim==1) )     return 5;
  assert(false);
}


tPDVector TwoKaonOnePionCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kplus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::K0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::Kbar0);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::K0);
  }
  else if(imode==3) {
    extpart[0]=getParticleData(ParticleID::K_S0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_S0);
  }
  else if(imode==4) {
    extpart[0]=getParticleData(ParticleID::K_L0);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::K_L0);
  }
  else if(imode==5) {
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
