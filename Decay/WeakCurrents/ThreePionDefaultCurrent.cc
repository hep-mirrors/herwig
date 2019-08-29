// -*- C++ -*-
//
// ThreePionDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreePionDefaultCurrent class.
//

#include "ThreePionDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using namespace ThePEG;

DescribeClass<ThreePionDefaultCurrent,WeakCurrent>
describeHerwigThreePionDefaultCurrent("Herwig::ThreePionDefaultCurrent",
				       "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(ThreePionDefaultCurrent,Energy,Energy2)


ThreePionDefaultCurrent::ThreePionDefaultCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(4);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi=ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = {1.0,-0.145,0.};
  // local values of the a_1 parameters
  _a1mass  = 1.251*GeV;
  _a1width = 0.599*GeV;
  _a1opt   = true;
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  // initialization of the a_1 running width
  _initializea1=false;
  double a1q2in[200]={0,15788.6,31577.3,47365.9,63154.6,78943.2,94731.9,110521,
		       126309,142098,157886,173675,189464,205252,221041,236830,
		       252618,268407,284196,299984,315773,331562,347350,363139,
		       378927,394716,410505,426293,442082,457871,473659,489448,
		       505237,521025,536814,552603,568391,584180,599969,615757,
		       631546,647334,663123,678912,694700,710489,726278,742066,
		       757855,773644,789432,805221,821010,836798,852587,868375,
		       884164,899953,915741,931530,947319,963107,978896,994685,
		       1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06,1.07363e+06,
		       1.08942e+06,1.10521e+06,1.12099e+06,1.13678e+06,1.15257e+06,
		       1.16836e+06,1.18415e+06,1.19994e+06,1.21573e+06,1.23151e+06,
		       1.2473e+06,1.26309e+06,1.27888e+06,1.29467e+06,1.31046e+06,
		       1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06,1.3894e+06,
		       1.40519e+06,1.42098e+06,1.43677e+06,1.45256e+06,1.46834e+06
		       ,1.48413e+06,1.49992e+06,1.51571e+06,1.5315e+06,1.54729e+06,
		       1.56308e+06,1.57886e+06,1.59465e+06,1.61044e+06,1.62623e+06,
		       1.64202e+06,1.65781e+06,1.6736e+06,1.68939e+06,1.70517e+06,
		       1.72096e+06,1.73675e+06,1.75254e+06,1.76833e+06,1.78412e+06,
		       1.79991e+06,1.81569e+06,1.83148e+06,1.84727e+06,1.86306e+06,
		       1.87885e+06,1.89464e+06,1.91043e+06,1.92621e+06,1.942e+06,
		       1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06,2.02095e+06,
		       2.03674e+06,2.05252e+06,2.06831e+06,2.0841e+06,2.09989e+06,
		       2.11568e+06,2.13147e+06,2.14726e+06,2.16304e+06,2.17883e+06,
		       2.19462e+06,2.21041e+06,2.2262e+06,2.24199e+06,2.25778e+06,
		       2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06,2.33672e+06,
		       2.35251e+06,2.3683e+06,2.38409e+06,2.39987e+06,2.41566e+06,
		       2.43145e+06,2.44724e+06,2.46303e+06,2.47882e+06,2.49461e+06,
		       2.51039e+06,2.52618e+06,2.54197e+06,2.55776e+06,2.57355e+06,
		       2.58934e+06,2.60513e+06,2.62092e+06,2.6367e+06,2.65249e+06,
		       2.66828e+06,2.68407e+06,2.69986e+06,2.71565e+06,2.73144e+06,
		       2.74722e+06,2.76301e+06,2.7788e+06,2.79459e+06,2.81038e+06,
		       2.82617e+06,2.84196e+06,2.85774e+06,2.87353e+06,2.88932e+06,
		       2.90511e+06,2.9209e+06,2.93669e+06,2.95248e+06,2.96827e+06,
		       2.98405e+06,2.99984e+06,3.01563e+06,3.03142e+06,3.04721e+06,
		       3.063e+06,3.07879e+06,3.09457e+06,3.11036e+06,3.12615e+06,
		       3.14194e+06};
  double a1widthin[200]={0,0,0,0,0,0,0,0,0,0,0,0,0.00153933,0.0136382,0.0457614,
			 0.105567,0.199612,0.333825,0.513831,0.745192,1.0336,1.38501,
			 1.80581,2.30295,2.88403,3.5575,4.33278,5.22045,6.23243,
			 7.38223,8.68521,10.1589,11.8234,13.7018,15.8206,18.2107,
			 20.9078,23.9533,27.3954,31.2905,35.7038,40.7106,46.3984,
			 52.8654,60.2207,68.581,78.0637,88.7754,100.794,114.145,
			 128.783,144.574,161.299,178.683,196.426,214.248,231.908,
			 249.221,266.059,282.336,298.006,313.048,327.46,341.254,
			 354.448,367.066,379.133,390.677,401.726,412.304,422.439,
			   432.155,441.474,450.419,459.01,467.267,475.207,482.847,
			 490.203,497.29,504.121,510.71,517.068,523.207,529.138,
			 534.869,540.411,545.776,550.961,556.663,560.851,565.566,
			 570.137,574.569,578.869,583.041,587.091,591.023,594.843,
			 598.553,602.16,605.664,609.072,612.396,615.626,618.754,
			 621.796,624.766,627.656,630.47,633.21,635.878,638.5,
			 641.006,643.471,645.873,648.213,650.493,652.715,654.88,
			 656.99,659.047,661.052,663.007,664.963,666.771,668.6,
			 670.351,672.075,673.828,675.397,676.996,678.567,680.083,
			 681.589,683.023,684.457,685.825,687.18,688.499,689.789,
			 691.058,692.284,693.501,694.667,695.82,696.947,698.05,
			 699.129,700.186,701.221,702.234,703.226,704.198,705.158,
			 706.085,707.001,707.899,708.78,709.644,710.474,711.334,
			 712.145,712.943,713.727,714.505,715.266,716.015,716.751,
			 717.474,718.183,718.88,719.645,720.243,720.91,721.565,
			 722.211,722.851,723.473,724.094,724.697,725.296,725.886,
			 726.468,727.041,727.608,728.166,728.718,729.262,729.808,
			 730.337,730.856,731.374,731.883,732.386,732.884,733.373,
			 733.859,734.339,734.813};

  vector<double> tmp1(a1widthin,a1widthin+200);
  _a1runwidth.clear();
  std::transform(tmp1.begin(), tmp1.end(),
		 back_inserter(_a1runwidth),
		 [](double x){return x*MeV;});
  
  vector<double> tmp2(a1q2in,a1q2in+200);
  _a1runq2.clear();
  std::transform(tmp2.begin(), tmp2.end(),
		 back_inserter(_a1runq2),
		 [](double x){return x*MeV2;});

  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void ThreePionDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  // initialise the a_1 running width calculation
  inita1Width(-1);
}

void ThreePionDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts <<  ounit(_a1runwidth,GeV)<< ounit(_a1runq2,GeV2)
     << _initializea1 << ounit(_a1mass,GeV) << ounit(_a1width,GeV)
     << ounit(_fpi,GeV) << ounit(_mpi,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF123widths,GeV)
     << _a1opt << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV) << _a1runinter;
}

void ThreePionDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >>  iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) 
     >> _initializea1 >> iunit(_a1mass,GeV) >> iunit(_a1width,GeV)
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF123widths,GeV)
     >> _a1opt >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV) >> _a1runinter;
}

void ThreePionDefaultCurrent::Init() {
        
  static ClassDocumentation<ThreePionDefaultCurrent> documentation
    ("The ThreePionDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );

  
  static ParVector<ThreePionDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &ThreePionDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<ThreePionDefaultCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreePionDefaultCurrent::_initializea1, false, false, false);
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
  
  static Switch<ThreePionDefaultCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &ThreePionDefaultCurrent::_a1opt, true, false, false);
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

  static ParVector<ThreePionDefaultCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreePionDefaultCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<ThreePionDefaultCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreePionDefaultCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);
    
  static Parameter<ThreePionDefaultCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &ThreePionDefaultCurrent::_a1width, GeV, 0.599*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<ThreePionDefaultCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &ThreePionDefaultCurrent::_a1mass, GeV, 1.251*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static ParVector<ThreePionDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &ThreePionDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreePionDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &ThreePionDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ThreePionDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &ThreePionDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool ThreePionDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
					  FlavourInfo flavour,
					  unsigned int imode,PhaseSpaceModePtr mode,
					  unsigned int iloc,int ires,
					  PhaseSpaceChannel phase, Energy upp ) {
  // check the charge and resonance
  if(imode==2||imode==3) {
    if(icharge!=0) return false;
    if(resonance && resonance->id()!=ParticleID::a_10) return false;
  }
  else if(imode<2) {
    if(abs(icharge)!=3) return false;
    if(resonance && abs(resonance->id())!=ParticleID::a_1plus) return false;
  }
  else
    assert(false);
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IOne) return false;
  }
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
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
  // and other flavour
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return false;
  // get the particles and check the masses
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1;
  if(icharge==-3) {
    a1=getParticleData(ParticleID::a_1minus);
  }
  else if(icharge==3) {
    a1=getParticleData(ParticleID::a_1plus);
  }
  else {
    a1=getParticleData(ParticleID::a_10);
  }
  _maxmass=max(_maxmass,upp);
  // the rho0 resonances
  tPDPtr rho0[3]   = { getParticleData(113), getParticleData(100113), getParticleData(30113)};
  tPDPtr rhoc[3]   = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  if(icharge==3)
    for(unsigned int ix=0;ix<3;++ix) rhoc  [ix] =   rhoc[ix]->CC();
  // create the mode
  if(imode==0) {
    // channels for pi- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,rho0[ix],
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rho0[ix],
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  else if(imode==1) {
    // channels for pi0 pi0 pi-
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,rhoc[ix],
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rhoc[ix],
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  else if(imode>=2) {
    // channels  for pi+ pi- pi0
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+2,ires+1,rhoc[ix]->CC(),
			ires+2,iloc+1,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,a1,ires+1,iloc+1,ires+1,rhoc[ix],
			ires+2,iloc+2,ires+2,iloc+3));
    }
  }
  // reset the parameters for the resonances in the integration
  mode->resetIntermediate(a1,_a1mass,_a1width);
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],
			    _rhoF123widths[ix]);
    mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],
			    _rhoF123widths[ix]);
  }
  return true;
}

// initialisation of the a_1 width
// (iopt=-1 initialises, iopt=0 starts the interpolation)
void ThreePionDefaultCurrent::inita1Width(int iopt) {
  if(iopt==-1) {
    _maxcalc=_maxmass;
    if(!_initializea1||_maxmass==ZERO) return;
    // parameters for the table of values
    Energy2 step(sqr(_maxcalc)/199.);
    // integrator to perform the integral
    vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
    vector<int> intype;intype.push_back(2);intype.push_back(3);
    Energy mrho(getParticleData(ParticleID::rhoplus)->mass()),
      wrho(getParticleData(ParticleID::rhoplus)->width());
    vector<Energy> inmass(2,mrho),inwidth(2,wrho);
    vector<double> inpow(2,0.0);
    ThreeBodyAllOnCalculator<ThreePionDefaultCurrent> 
      widthgen(inweights,intype,inmass,inwidth,inpow,*this,0,_mpi,_mpi,_mpi);
    // normalisation constant to give physical width if on shell
    double a1const(_a1width/(widthgen.partialWidth(sqr(_a1mass))));
    // loop to give the values
    _a1runq2.clear(); _a1runwidth.clear();
    for(Energy2 moff2(ZERO); moff2<=sqr(_maxcalc); moff2+=step) {
      _a1runwidth.push_back(widthgen.partialWidth(moff2)*a1const);
      _a1runq2.push_back(moff2);
    }
  }
  // set up the interpolator
  else if(iopt==0) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
  }
}

void ThreePionDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ThreePionDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  output << "newdef " << name() << ":Initializea1 "    << _initializea1    << "\n";
  output << "newdef " << name() << ":a1WidthOption "   << _a1opt           << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    output << "newdef " << name() << ":a1RunningWidth " << ix 
	   << " " << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    output << "newdef " << name() << ":a1RunningQ2 " << ix 
	   << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  output << "newdef " << name() << ":A1Width " << _a1width/GeV << "\n";
  output << "newdef " << name() << ":A1Mass "  << _a1mass/GeV  << "\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void ThreePionDefaultCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  WeakCurrent::doinitrun();
}

void ThreePionDefaultCurrent::doupdate() {
  WeakCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

double ThreePionDefaultCurrent::
threeBodyMatrixElement(const int       , const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy    , 
		       const Energy    , const Energy    ) const {
  Energy2 mpi2(sqr(_mpi));
  Complex propb(BrhoF123(s1,-1)),propa(BrhoF123(s2,-1)); 
  // the matrix element
  Energy2 output(ZERO); 
  // first resonance
  output += ((s1-4.*mpi2) + 0.25*(s3-s2)*(s3-s2)/q2) * real(propb*conj(propb)); 
  // second resonance
  output += ((s2-4.*mpi2) + 0.25*(s3-s1)*(s3-s1)/q2) * real(propa*conj(propa)); 
  // the interference term 
  output += (0.5*q2-s3-0.5*mpi2+0.25*(s3-s2)*(s3-s1)/q2)*real(propa*conj(propb)+
							      propb*conj(propa)); 
  return output/sqr(_rhoF123masses[0]);
}


// the hadronic currents    
vector<LorentzPolarizationVectorE> 
ThreePionDefaultCurrent::current(tcPDPtr resonance,
			      FlavourInfo flavour,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & ,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IOne)
    return vector<LorentzPolarizationVectorE>();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
    case IsoSpin::I3Zero:
      if(imode<=1) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3One:
      if(imode>=2) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusOne:
      if(imode>=2) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return vector<LorentzPolarizationVectorE>();
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return vector<LorentzPolarizationVectorE>();
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return vector<LorentzPolarizationVectorE>();
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
  Complex F1(0.), F2(0.);
  Complex a1fact(2./3.);
  if(!resonance) a1fact *= a1BreitWigner(q2);
  if(ichan<0) {
    F1 = a1fact*BrhoF123(s1,-1);
    F2 =-a1fact*BrhoF123(s2,-1);
  }
  else if(ichan%2==0) F1 = a1fact*BrhoF123(s1,    ichan/2);
  else if(ichan%2==1) F2 =-a1fact*BrhoF123(s2,(ichan-1)/2);
  // the first three form-factors
  LorentzPolarizationVectorE vect = (F2-F1)*momenta[2] +F1*momenta[1] -F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool ThreePionDefaultCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),npi0(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) return true;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) return true;
  else if( npip==1 && npim==1 && npi0 ==1 )           return true;
  return false;
}

unsigned int ThreePionDefaultCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),npi0(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(      (npip==2&&npim==1) || (npim==2&&npip==1) ) return 0;
  else if( (npip==1&&npi0==2) || (npim==1&&npi0==2) ) return 1;
  else if( npip==1 && npim==1 && npi0 ==1 )           return 2;
  else assert(false);
}

tPDVector ThreePionDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
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
  else if(imode==2||imode==3) {
    extpart[0]=getParticleData(ParticleID::piplus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else
    assert(false);
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
