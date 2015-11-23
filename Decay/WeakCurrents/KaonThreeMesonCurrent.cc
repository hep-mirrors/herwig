// -*- C++ -*-
//
// KaonThreeMesonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KaonThreeMesonCurrent class.
//

#include "KaonThreeMesonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeClass<KaonThreeMesonCurrent,ThreeMesonCurrentBase>
describeHerwigKaonThreeMesonCurrent("Herwig::KaonThreeMesonCurrent",
				    "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(KaonThreeMesonCurrent,Energy,Energy2)

namespace {
  inline Energy  timesGeV (double x) { return x * GeV; }
  inline Energy2 timesGeV2(double x) { return x * GeV2; }
}
 
IBPtr KaonThreeMesonCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr KaonThreeMesonCurrent::fullclone() const {
  return new_ptr(*this);
}

KaonThreeMesonCurrent::KaonThreeMesonCurrent() {
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
		 timesGeV);
  
  vector<double> tmp2(a1q2in,a1q2in+200);
  _a1runq2.clear();
  std::transform(tmp2.begin(), tmp2.end(),
		 back_inserter(_a1runq2),
		 timesGeV2);


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


void KaonThreeMesonCurrent::persistentOutput(PersistentOStream & os) const {
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

void KaonThreeMesonCurrent::persistentInput(PersistentIStream & is, int) {
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


void KaonThreeMesonCurrent::Init() {

  static ClassDocumentation<KaonThreeMesonCurrent> documentation
    ("The KaonThreeMesonCurrent class implements the model of "
     "Z. Phys.  C 69 (1996) 243 [arXiv:hep-ph/9503474]"
     " for the weak current with three "
     "mesons, at least one of which is a kaon",
     "The KaonThreeMesonCurrent class implements the model of "
     "\\cite{Finkemeier:1995sr} for the weak current with three "
     "mesons, at least one of which is a kaon.",
     "\\bibitem{Finkemeier:1995sr}\n"
     "M.~Finkemeier and E.~Mirkes,\n"
     "Z.\\ Phys.\\  C {\\bf 69} (1996) 243 [arXiv:hep-ph/9503474].\n"
     " %%CITATION = ZEPYA,C69,243;%%\n"

);

  static Switch<KaonThreeMesonCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &KaonThreeMesonCurrent::_initializea1, false, false, false);
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

  static Parameter<KaonThreeMesonCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &KaonThreeMesonCurrent::_a1width, GeV, 0.599*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<KaonThreeMesonCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &KaonThreeMesonCurrent::_a1mass, GeV, 1.251*GeV, ZERO, 10.0*GeV,
     false, false, false);

  static Parameter<KaonThreeMesonCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &KaonThreeMesonCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static ParVector<KaonThreeMesonCurrent,Energy> interfaceRhoAxialMasses
    ("RhoAxialMasses",
     "The masses for the rho resonances if used local values",
     &KaonThreeMesonCurrent::_rho1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,Energy> interfaceRhoAxialWidths
    ("RhoAxialWidths",
     "The widths for the rho resonances if used local values",
     &KaonThreeMesonCurrent::_rho1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,Energy> interfaceRhoVectorMasses
    ("RhoVectorMasses",
     "The masses for the rho resonances if used local values",
     &KaonThreeMesonCurrent::_rho2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,Energy> interfaceRhoVectorWidths
    ("RhoVectorWidths",
     "The widths for the rho resonances if used local values",
     &KaonThreeMesonCurrent::_rho2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,Energy> interfaceKstarAxialMasses
    ("KstarAxialMasses",
     "The masses for the Kstar resonances if used local values",
     &KaonThreeMesonCurrent::_kstar1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,Energy> interfaceKstarAxialWidths
    ("KstarAxialWidths",
     "The widths for the Kstar resonances if used local values",
     &KaonThreeMesonCurrent::_kstar1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<KaonThreeMesonCurrent,Energy> interfaceKstarVectorMasses
    ("KstarVectorMasses",
     "The masses for the Kstar resonances if used local values",
     &KaonThreeMesonCurrent::_kstar2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,Energy> interfaceKstarVectorWidths
    ("KstarVectorWidths",
     "The widths for the Kstar resonances if used local values",
     &KaonThreeMesonCurrent::_kstar2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<KaonThreeMesonCurrent,double> interfaceAxialRhoWeight
    ("AxialRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &KaonThreeMesonCurrent::_rho1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,double> interfaceAxialKStarWeight
    ("AxialKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &KaonThreeMesonCurrent::_kstar1wgts,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<KaonThreeMesonCurrent,double> interfaceVectorRhoWeight
    ("VectorRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &KaonThreeMesonCurrent::_rho2wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<KaonThreeMesonCurrent,double> interfaceVectorKStarWeight
    ("VectorKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &KaonThreeMesonCurrent::_kstar2wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Switch<KaonThreeMesonCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values of the rho meson masses and widths",
     &KaonThreeMesonCurrent::_rhoparameters, true, false, false);
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
  
  static Switch<KaonThreeMesonCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values of the rho meson masses and widths",
     &KaonThreeMesonCurrent::_kstarparameters, true, false, false);
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
  
  static Switch<KaonThreeMesonCurrent,bool> interfacea1Parameters
    ("a1Parameters",
     "Use local values of the rho meson masses and widths",
     &KaonThreeMesonCurrent::_a1parameters, true, false, false);
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
  
  static Switch<KaonThreeMesonCurrent,bool> interfaceK1Parameters
    ("K1Parameters",
     "Use local values of the rho meson masses and widths",
     &KaonThreeMesonCurrent::_k1parameters, true, false, false);
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
  
  static Switch<KaonThreeMesonCurrent,bool> interfacea1WidthOption
    ("a1WidthOption",
     "Option for the treatment of the a1 width",
     &KaonThreeMesonCurrent::_a1opt, true, false, false);
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

  static ParVector<KaonThreeMesonCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &KaonThreeMesonCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<KaonThreeMesonCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &KaonThreeMesonCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

  static ParVector<KaonThreeMesonCurrent,Energy> interfaceK1Masses
    ("K1Masses",
     "Masses of the K_1 mesons",
     &KaonThreeMesonCurrent::_k1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KaonThreeMesonCurrent,Energy> interfaceK1Widths
    ("K1Widths",
     "Widths of the K_1 mesons",
     &KaonThreeMesonCurrent::_k1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KaonThreeMesonCurrent,double> interfaceK1WeightKStarPi
    ("K1WeightKStarPi",
     "The relative weights for the K_1 resonances in the K* pi final-state",
     &KaonThreeMesonCurrent::_k1wgta, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);

  static ParVector<KaonThreeMesonCurrent,double> interfaceK1WeightRhoK
    ("K1WeightRhoK",
     "The relative weights for the K_1 resonances in the rho K final-state",
     &KaonThreeMesonCurrent::_k1wgtb, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);

  static Parameter<KaonThreeMesonCurrent,double> interfaceEpsOmega
    ("EpsOmega",
     "The omega-phi mixing ",
     &KaonThreeMesonCurrent::_epsomega, 0.05, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<KaonThreeMesonCurrent,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &KaonThreeMesonCurrent::_omegamass, GeV, 0.782*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<KaonThreeMesonCurrent,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &KaonThreeMesonCurrent::_omegawidth, GeV, 0.00843*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<KaonThreeMesonCurrent,Energy> interfacePhiMass
    ("PhiMass",
     "The mass of the phi meson",
     &KaonThreeMesonCurrent::_phimass, GeV, 1.020*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<KaonThreeMesonCurrent,Energy> interfacePhiWidth
    ("PhiWidth",
     "The width of the phi meson",
     &KaonThreeMesonCurrent::_phiwidth, GeV, 0.00443*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<KaonThreeMesonCurrent,double> interfaceOmegaKStarWeight
    ("OmegaKStarWeight",
     "The relative weight of the omega-phi and K* terms",
     &KaonThreeMesonCurrent::_omegaKstarwgt, 1./sqrt(2.), 0.0, 100.0,
     false, false, Interface::limited);
  
  static Switch<KaonThreeMesonCurrent,bool> interfaceOmegaParameters
    ("OmegaParameters",
     "Use local values of the omega/phi meson masses and widths",
     &KaonThreeMesonCurrent::_omegaopt, true, false, false);
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

void KaonThreeMesonCurrent::inita1Width(int iopt) {
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
    ThreeBodyAllOnCalculator<KaonThreeMesonCurrent> 
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
bool KaonThreeMesonCurrent::acceptMode(int imode) const { 
  return imode>=2&&imode<=11&&imode!=8;
}



// complete the construction of the decay mode for integration
bool KaonThreeMesonCurrent::createMode(int icharge, unsigned int imode,
					  DecayPhaseSpaceModePtr mode,
					  unsigned int iloc,unsigned int ires,
					  DecayPhaseSpaceChannelPtr phase,Energy upp) {
  int iq(0),ia(0);
  if(!acceptMode(imode)) return false;
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1,k1[2];
  if(icharge==-3) {
    a1    = getParticleData(ParticleID::a_1minus);
    k1[0] = getParticleData(ParticleID::K_1minus);
    k1[1] = getParticleData(ParticleID::Kstar_1minus);
  }
  else if(icharge==3) {
    a1    = getParticleData(ParticleID::a_1plus);
    k1[0] = getParticleData(ParticleID::K_1plus);
    k1[1] = getParticleData(ParticleID::Kstar_1plus);
  }
  else {
    return false;
  }
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
    for(unsigned int ix=0;ix<3;++ix) {
      if(rhoc[ix]) rhoc[ix]=rhoc[ix]->CC();
      if(Kstar0[ix]) Kstar0[ix]=Kstar0[ix]->CC();
      if(Kstarc[ix]) Kstarc[ix]=Kstarc[ix]->CC();
    }
  }
  DecayPhaseSpaceChannelPtr newchannel;
  if(imode==2) {
    // channels for K- pi- K+
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstar0[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
      }
      if(rho0[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	newchannel->addIntermediate(rho0[ix],0,0.0,iloc,iloc+2);
	mode->addChannel(newchannel);
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstar0[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
      }
    }
  }
  else if(imode==3) {
    // channels for K0 pi- K0bar
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
      }
      if(rho0[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	newchannel->addIntermediate(rho0[ix],0,0.0, iloc,iloc+2);
	mode->addChannel(newchannel);
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstarc[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
      }
    }
  }
  else if(imode==4) {
    // channels for K- pi0 K0
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstar0[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
      }
      if(Kstarc[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc+2);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc  ,iloc+1);
	mode->addChannel(newchannel);
      }
      if(rhoc[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	mode->addChannel(newchannel);
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstar0[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(Kstarc[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc+2);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc  ,iloc+1);
	  mode->addChannel(newchannel);
	}
      }
    }
  }
  else if(imode==5) {  
    // channels for pi0 pi0 K-
    for(unsigned int ix=0;ix<3;++ix) {
      if(!Kstarc[ix]) continue;
      for(unsigned int ik=0;ik<2;++ik) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(k1[ik]    ,0,0.0,-ires-1,iloc);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(k1[ik]    ,0,0.0,-ires-1,iloc+1);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc,iloc+2);
	mode->addChannel(newchannel);
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!Kstarc[iy]) continue;
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc  );
	newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc  ,iloc+2);
	mode->addChannel(newchannel);
      }
    }
  }
  else if(imode==6) {
    // channels for K- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(rho0[ix]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1[ik]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(Kstar0[ix]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1[ik]    ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!Kstarc[ix]) continue;
	if(rho0[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[iy]  ,0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(Kstar0[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
      }
    }
  }
  else if(imode==7) {
    // channels for pi- kbar0 pi0
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(rhoc[ix]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1[ik]  ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(Kstar0[ix]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1[ik]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstar0[ix],0,0.0,iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(Kstarc[ix]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1[ik]  ,0,0.0,-ires-1,iloc+2);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc,iloc+1);
	  mode->addChannel(newchannel);
	}
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!Kstarc[ix]) continue;
	if(Kstar0[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(rhoc[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[iy]  ,0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
	if(Kstarc[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+2);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc  ,iloc+1);
	  mode->addChannel(newchannel);
	}
      }
    }
  }
  else if(imode==9||imode==10) {
    // channels for K_S0 pi- K_S0 and K_L0 pi- K_L0 
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc+2);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc  ,iloc+1);
	mode->addChannel(newchannel);
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstarc[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc+2);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc  ,iloc+1);
	  mode->addChannel(newchannel);
	}
      }
    }
  }
  else if(imode==11) {
    // channels for K_S0 pi- K_L0
    for(unsigned int ix=0;ix<3;++ix) {
      if(Kstarc[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	mode->addChannel(newchannel);
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc+2);
	newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc  ,iloc+1);
	mode->addChannel(newchannel);
      }
      if(rho0[ix]) {
	newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	newchannel->addIntermediate(rho0[ix],0,0.0, iloc,iloc+2);
	mode->addChannel(newchannel);
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(!rhoc[ix]) continue;
	if(Kstarc[iy]) {
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc+2);
	  newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc  ,iloc+1);
	  mode->addChannel(newchannel);
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

void KaonThreeMesonCurrent::dataBaseOutput(ofstream & os,
					   bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::KaonThreeMesonCurrent " 
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
  ThreeMesonCurrentBase::dataBaseOutput(os,false,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}  

void KaonThreeMesonCurrent::doinit() {
  ThreeMesonCurrentBase::doinit();
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

KaonThreeMesonCurrent::FormFactors
KaonThreeMesonCurrent::calculateFormFactors(const int ichan,const int imode,
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

void KaonThreeMesonCurrent::doinitrun() {
  // set up the running a_1 width
  inita1Width(0);
  ThreeMesonCurrentBase::doinitrun();
}

void KaonThreeMesonCurrent::doupdate() {
  ThreeMesonCurrentBase::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1Width(-1);
}

double KaonThreeMesonCurrent::
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

Complex KaonThreeMesonCurrent::Trho1(Energy2 q2,int ires) const {
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
  
Complex KaonThreeMesonCurrent::Trho2(Energy2 q2,int ires) const {
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
  
Complex KaonThreeMesonCurrent::TKstar1(Energy2 q2,int ires) const  {
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
  
Complex KaonThreeMesonCurrent::TKstar2(Energy2 q2,int ires) const {
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

Complex KaonThreeMesonCurrent::BWrho1(Energy2 q2, unsigned int ires) const {
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

Complex KaonThreeMesonCurrent::BWrho2(Energy2 q2, unsigned int ires) const {
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
  
Complex KaonThreeMesonCurrent::BWKstar1(Energy2 q2, unsigned int ires) const {
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

Complex KaonThreeMesonCurrent::BWKstar2(Energy2 q2, unsigned int ires) const  {
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

Complex KaonThreeMesonCurrent::a1BreitWigner(Energy2 q2) const {
  Complex ii(0.,1.);
  Energy2 m2(_a1mass*_a1mass);
  Energy  q(sqrt(q2));
  return m2/(m2-q2-ii*q*a1Width(q2));
}
  
Complex KaonThreeMesonCurrent::TK1(Energy2 q2,unsigned int iopt,int ires) const {
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

Complex KaonThreeMesonCurrent::K1BreitWigner(Energy2 q2,unsigned int ires) const {
  if(ires>=_k1mass.size()) return 0.;
  Energy2 m2=sqr(_k1mass[ires]),mg=_k1mass[ires]*_k1width[ires];
  return (-m2+Complex(0.,1.)*mg)/(q2-m2+Complex(0.,1.)*mg);
}

Energy KaonThreeMesonCurrent::a1Width(Energy2 q2) const {
  if(!_a1opt) return _a1mass*_a1width*g(q2)/g(_a1mass*_a1mass)/sqrt(q2);
  else        return (*_a1runinter)(q2);
}
  
double KaonThreeMesonCurrent::g(Energy2 q2) const {
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
  
Complex KaonThreeMesonCurrent::Tomega(Energy2 q2, int ires) const {
  double denom=(1.+_epsomega);
  Complex num(0.);
  if(ires<0) num=OmegaPhiBreitWigner(q2,0)+_epsomega*OmegaPhiBreitWigner(q2,1);
  else if(ires==0) num=OmegaPhiBreitWigner(q2,0);
  else             num=OmegaPhiBreitWigner(q2,1);
  return num/denom;
}

Complex KaonThreeMesonCurrent::OmegaPhiBreitWigner(Energy2 q2, unsigned int ires) const {
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

Complex KaonThreeMesonCurrent::TOmegaKStar(Energy2 s1,Energy2 s2,int ires) const {
  Complex output;
  if(ires<0)         output = _omegaKstarwgt*TKstar1(s1,-1)+Tomega(s2,-1);
  else if(ires%2==0) output = _omegaKstarwgt*TKstar1(s1,ires/2);
  else if(ires%2==1) output = Tomega(s2,ires/2);
  return output/(1.+_omegaKstarwgt);
}


