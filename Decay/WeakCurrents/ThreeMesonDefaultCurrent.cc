// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeMesonDefaultCurrent class.
//

#include "ThreeMesonDefaultCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeMesonDefaultCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

ThreeMesonDefaultCurrent::ThreeMesonDefaultCurrent() 
{
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi=0.;_mK=0.;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts.push_back(1.0);_rhoF123wgts.push_back(-0.145);
  _rhoF123wgts.push_back(0.);
  _rhoF5wgts.push_back(-26.);_rhoF5wgts.push_back(6.5);
  _rhoF5wgts.push_back(1.);
  // the Kstar weights
  _KstarF123wgts.push_back(1.);_KstarF123wgts.push_back(0.);
  _KstarF123wgts.push_back(0.);
  _KstarF5wgts.push_back(1.);_KstarF5wgts.push_back(0.);
  _KstarF5wgts.push_back(0.);
  // relative rho/Kstar weights
  _rhoKstarwgt=-0.2;
  // local values of the a_1 parameters
  _a1parameters=true;_a1mass=1.251*GeV;_a1width=0.599*GeV;
  // local values of the K_1 parameters
  _K1parameters=true;_K1mass=1.402*GeV,_K1width=0.174*GeV;
  // local values of the rho parameters
  _rhoparameters=true;
  _rhoF123masses.push_back(0.773*GeV);_rhoF123masses.push_back(1.370*GeV);
  _rhoF123masses.push_back(1.750*GeV);
  _rhoF123widths.push_back(0.145*GeV);_rhoF123widths.push_back(0.510*GeV);
  _rhoF123widths.push_back(0.120*GeV);
  _rhoF5masses.push_back(0.773*GeV);_rhoF5masses.push_back(1.500*GeV);
  _rhoF5masses.push_back(1.750*GeV);
  _rhoF5widths.push_back(0.145*GeV);_rhoF5widths.push_back(0.220*GeV);
  _rhoF5widths.push_back(0.120*GeV);
  // local values for the Kstar parameters
  _Kstarparameters=true;
  _KstarF123masses.push_back(0.8921*GeV);
  _KstarF123widths.push_back(0.0513*GeV);
  _KstarF5masses.push_back(0.892*GeV);
  _KstarF5widths.push_back(0.0513*GeV);
  // initialization of the a_1 running width
  _initializea1=false;
  Energy2 a1q2in[200]={0,15788.6,31577.3,47365.9,63154.6,78943.2,94731.9,110521,
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
  Energy a1widthin[200]={0,0,0,0,0,0,0,0,0,0,0,0,0.00153933,0.0136382,0.0457614,
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
  _a1runwidth=vector<Energy>(a1widthin,a1widthin+200);
  _a1runq2=vector<Energy2>(a1q2in,a1q2in+200);
}

void ThreeMesonDefaultCurrent::doinit() throw(InitException) {
  ThreeMesonCurrentBase::doinit();
  // the particles we will use a lot
  tPDPtr a1(getParticleData(ParticleID::a_1minus)),
    k1(getParticleData(ParticleID::K_1minus)),pi0(getParticleData(ParticleID::pi0)),
    piplus(getParticleData(ParticleID::piplus)),
    piminus(getParticleData(ParticleID::piminus));
  // masses for the running widths
  _mpi=piplus->mass();
  _mK=getParticleData(ParticleID::Kminus)->mass();
  // the charged rho resonances
  tPDPtr rhoc[3]={getParticleData(-213),getParticleData(-100213),
		  getParticleData(-30213)};
  // the charged K* resonances
  tPDPtr Kstarc[3]={getParticleData(-323),getParticleData(-100323),
		    getParticleData(-30323)};
  if(!_a1parameters){_a1mass=a1->mass();_a1width=a1->width();}
  // mass and width of the k_1
  if(!_K1parameters){_K1mass=k1->mass();_K1width=k1->width();}
  // initialise the a_1 running width calculation
  if(_initializea1){inita1width(-1);}
  // rho parameters in the base classs
  tcPDPtr temp;
  unsigned int ix;
  if(_rhoparameters&&_rhoF123masses.size()<3)
    {
      ix = _rhoF123masses.size();
      _rhoF123masses.resize(3);_rhoF123widths.resize(3);
      for(;ix<3;++ix)
	{
	  _rhoF123masses[ix]=rhoc[ix]->mass();
	  _rhoF123widths[ix]=rhoc[ix]->width();
	}
    }
  else if(!_rhoparameters)
    {
      _rhoF123masses.resize(3);_rhoF123widths.resize(3);
      for(ix=0;ix<3;++ix)
	{
	  _rhoF123masses[ix]=rhoc[ix]->mass();
	  _rhoF123widths[ix]=rhoc[ix]->width();
	}
    }
  // K star parameters in the base class
  if(_Kstarparameters&&_KstarF123masses.size()<3)
    {
      ix = _KstarF123masses.size();
      _KstarF123masses.resize(3);_KstarF123widths.resize(3);
      for(;ix<3;++ix)
	{
	  _KstarF123masses[ix]=Kstarc[ix]->mass();
	  _KstarF123widths[ix]=Kstarc[ix]->width();
	}
    }
  else if(!_Kstarparameters)
    {
      _KstarF123masses.resize(3);_KstarF123widths.resize(3);
      for(ix=0;ix<3;++ix)
	{
	  _KstarF123masses[ix]=Kstarc[ix]->mass();
	  _KstarF123widths[ix]=Kstarc[ix]->width();
	}
    }
  // rho parameters here
  if(_rhoparameters&&_rhoF5masses.size()<3)
    {
      ix = _rhoF5masses.size();
      _rhoF5masses.resize(3);_rhoF5widths.resize(3);
      for(;ix<3;++ix)
	{
	  _rhoF5masses[ix]=rhoc[ix]->mass();
	  _rhoF5widths[ix]=rhoc[ix]->width();
	}
    }
  else if(!_rhoparameters)
    {
      _rhoF5masses.resize(3);_rhoF5widths.resize(3);
      for(ix=0;ix<3;++ix)
	{
	  _rhoF5masses[ix]=rhoc[ix]->mass();
	  _rhoF5widths[ix]=rhoc[ix]->width();
	}
    }
  // Kstar parameters here
  if(_Kstarparameters&&_KstarF5widths.size()<3)
    {
      ix = _KstarF5masses.size();
      _KstarF5masses.resize(3);_KstarF5widths.resize(3);
      for(;ix<3;++ix)
	{
	  _KstarF5masses[ix]=Kstarc[ix]->mass();
	  _KstarF5widths[ix]=Kstarc[ix]->width();
	}
    }
  else if(!_Kstarparameters)
    {
      _KstarF5masses.resize(3);_KstarF5widths.resize(3);
      for(ix=0;ix<3;++ix)
	{
	  _KstarF5masses[ix]=Kstarc[ix]->mass();
	  _KstarF5widths[ix]=Kstarc[ix]->width();
	}
    }
}

void ThreeMesonDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _KstarF123wgts << _rhoF5wgts << _KstarF5wgts
     << _rhoKstarwgt <<  _a1runwidth << _a1runq2 <<  _initializea1
     << _a1mass << _a1width << _K1mass << _K1width << _fpi << _mpi << _mK
     <<_rhoparameters << _rhoF123masses << _rhoF5masses << _rhoF123widths 
     << _rhoF5widths << _Kstarparameters << _KstarF123masses <<_KstarF5masses
     << _KstarF123widths << _KstarF5widths << _a1parameters << _K1parameters;
}

void ThreeMesonDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _KstarF123wgts >> _rhoF5wgts >> _KstarF5wgts
     >> _rhoKstarwgt >>  _a1runwidth >> _a1runq2 >>  _initializea1
     >> _a1mass >> _a1width >> _K1mass >> _K1width >> _fpi >> _mpi >> _mK
     >>_rhoparameters >> _rhoF123masses >> _rhoF5masses >> _rhoF123widths 
     >> _rhoF5widths >> _Kstarparameters >> _KstarF123masses >>_KstarF5masses
     >> _KstarF123widths >> _KstarF5widths >> _a1parameters >> _K1parameters;
}

ClassDescription<ThreeMesonDefaultCurrent> ThreeMesonDefaultCurrent::initThreeMesonDefaultCurrent;
// Definition of the static class description member.

void ThreeMesonDefaultCurrent::Init() {
        
  static ClassDocumentation<ThreeMesonDefaultCurrent> documentation
    ("The ThreeMesonDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.");
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF123KstarWgt
    ("F123KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_KstarF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF5RhoWgt
    ("F5RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_rhoF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,double> interfaceF5KstarWgt
    ("F5KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &ThreeMesonDefaultCurrent::_KstarF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Parameter<ThreeMesonDefaultCurrent,double> interfaceRhoKstarWgt
    ("RhoKstarWgt",
     "The relative weights of the rho and K* in the F5 form factor",
     &ThreeMesonDefaultCurrent::_rhoKstarwgt, -0.2, -10., 10.,
     false, false, false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &ThreeMesonDefaultCurrent::_initializea1, false, false, false);
  static SwitchOption interfaceInitializea1Initialization
    (interfaceInitializea1,
     "Initialization",
     "Initialize the calculation",
     true);
  static SwitchOption interfaceInitializea1NoInitialization
    (interfaceInitializea1,
     "NoInitialization",
     "Use the default values",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceRhoParameters
    ("RhoParameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_rhoparameters, true, false, false);
  static SwitchOption interfaceRhoParameterstrue
    (interfaceRhoParameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceRhoParametersParticleData
    (interfaceRhoParameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_Kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
       "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfacea1Parameters
    ("a1Parameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_a1parameters, true, false, false);
  static SwitchOption interfacea1Parameterstrue
    (interfacea1Parameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfacea1ParametersParticleData
    (interfacea1Parameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static Switch<ThreeMesonDefaultCurrent,bool> interfaceK1Parameters
    ("K1Parameters",
     "Use local values of the rho meson masses and widths",
     &ThreeMesonDefaultCurrent::_K1parameters, true, false, false);
  static SwitchOption interfaceK1Parameterstrue
    (interfaceK1Parameters,
     "Local",
     "Use local values of the parameters",
     true);
  static SwitchOption interfaceK1ParametersParticleData
    (interfaceK1Parameters,
     "ParticleData",
     "Use the masses and wdiths from the particle data objects",
     false);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &ThreeMesonDefaultCurrent::_a1runwidth, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<ThreeMesonDefaultCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &ThreeMesonDefaultCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, 0.0*GeV2, 10.0*GeV2,
     false, false, true);
    
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceA1Width
    ("A1Width",
     "The a_1 width if using local values.",
     &ThreeMesonDefaultCurrent::_a1width, GeV, 0.599*GeV, 0*GeV, 10.0*GeV,
     false, false, false);
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceA1Mass
    ("A1Mass",
     "The a_1 mass if using local values.",
     &ThreeMesonDefaultCurrent::_a1mass, GeV, 1.251*GeV, 0*GeV, 10.0*GeV,
     false, false, false);
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceK1Width
    ("K1Width",
     "The K_1 width if using local values.",
     &ThreeMesonDefaultCurrent::_K1width, GeV, 0.174*GeV, 0*GeV, 10.0*GeV,
     false, false, false);
  
  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceK1Mass
    ("K1Mass",
     "The K_1 mass if using local values.",
     &ThreeMesonDefaultCurrent::_K1mass, GeV, 1.402*GeV, 0*GeV, 10.0*GeV,
     false, false, false);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfacerhoF5masses
    ("rhoF5masses",
     "The masses for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF5masses, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfacerhoF5widths
    ("rhoF5widths",
     "The widths for the rho resonances if used local values",
     &ThreeMesonDefaultCurrent::_rhoF5widths, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<ThreeMesonDefaultCurrent,Energy> interfaceKstarF123masses
    ("KstarF123masses",
     "The masses for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF123masses, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfaceKstarF123widths
    ("KstarF123widths",
     "The widths for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF123widths, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfaceKstarF5masses
    ("KstarF5masses",
     "The masses for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF5masses, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);
  
  static ParVector<ThreeMesonDefaultCurrent,Energy> interfaceKstarF5widths
    ("KstarF5widths",
     "The widths for the Kstar resonances if used local values",
     &ThreeMesonDefaultCurrent::_KstarF5widths, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ThreeMesonDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &ThreeMesonDefaultCurrent::_fpi, MeV, 92.4*MeV, 0.0*MeV, 200.0*MeV,
     false, false, true);
}
  
// modes handled by this class
bool ThreeMesonDefaultCurrent::acceptMode(int imode) const{return imode>=0&&imode<=8;}

// calculate the form-factors
void ThreeMesonDefaultCurrent::
calculateFormFactors(const int ichan, const int imode,
		     Energy2 q2, Energy2 s1, Energy2 s2, Energy2 s3,
		     Complex & F1, Complex & F2, Complex & F3, Complex & F4,
		     Complex & F5) const
{
  F1=0.;F2=0.;F3=0.;F4=0.;F5=0.;
  // calculate the pi- pi- pi+ factor
  if(imode==0)
    {
      Complex a1fact(a1BreitWigner(q2)*2./3.);
      if(ichan<0){F1= a1fact*BrhoF123(s1,-1);F2 =-a1fact*BrhoF123(s2,-1);}
      else if(ichan%2==0){F1 = a1fact*BrhoF123(s1,ichan/2);}
      else if(ichan%2==1){F2 =-a1fact*BrhoF123(s2,(ichan-1)/2);}
    }
  // calculate the pi0 pi0 pi- factor
  else if(imode==1)
    {
      Complex a1fact(a1BreitWigner(q2)*2./3.);
      if(ichan<0){F1 = a1fact*BrhoF123(s1,-1);F2 =-a1fact*BrhoF123(s2,-1);}
      else if(ichan%2==0){F1 = a1fact*BrhoF123(s1,ichan/2);}
      else if(ichan%2==1){F2 =-a1fact*BrhoF123(s2,(ichan-1)/2);}
    }
  // calculate the K- pi - K+ factor
  else if(imode==2)
    {
      Complex a1fact(a1BreitWigner(q2)*sqrt(2.)/3.);
      if(ichan<0)
	{
	  F1 =-a1fact*BKstarF123(s1,-1); F2 = a1fact*BrhoF123(s2,-1);
	  F5 = BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
	}
      else if(ichan%8==0){F1 =-a1fact*BKstarF123(s1,ichan/8);}
      else if(ichan%8==1){F2 = a1fact*BrhoF123(s2,(ichan-1)/8);}
      else if(ichan%8>=2){F5 = BrhoF5(q2,ichan/8)*FKrho(s1,s2,(ichan-2)%8)*sqrt(2.);}
    }
  // calculate the K0 pi- K0bar
  else if(imode==3)
    {
      Complex a1fact(a1BreitWigner(q2)*sqrt(2.)/3.);
      if(ichan<0)
	{
	  F1 =-a1fact*BKstarF123(s1,-1);F2 = a1fact*BrhoF123(s2,-1);
	  F5 =-BrhoF5(q2,-1)*FKrho(s1,s2,-1)*sqrt(2.);
	}
      else if(ichan%8==0){F1 = -a1fact*BKstarF123(s1,ichan/8);}
      else if(ichan%8==1){F2 = a1fact*BrhoF123(s2,(ichan-1)/8);}
      else if(ichan%8>=2){F5 = -BrhoF5(q2,ichan/8)*FKrho(s1,s2,(ichan-2)%8)*sqrt(2.);}
    }
  // calculate the K- pi0 k0
  else if(imode==4)
    {
      Complex a1fact(a1BreitWigner(q2));
      if(ichan<0){F2 =-a1fact*BrhoF123(s2,-1);}
      else{F2 =-a1fact*BrhoF123(s2,ichan);}
    }
  // calculate the pi0 pi0 K-
  else if(imode==5)
    {
      Complex K1fact(K1BreitWigner(q2)/6.);
      if(ichan<0){F1 = K1fact*BKstarF123(s1,-1);F2 =-K1fact*BKstarF123(s2,-1);}
      else if(ichan%2==0){F1 = K1fact*BKstarF123(s1,ichan/2);}
      else{F2 =-K1fact*BKstarF123(s2,(ichan-1)/2);}
    }
  // calculate the K- pi- pi+
  else if(imode==6)
    {
      Complex K1fact(K1BreitWigner(q2)*sqrt(2.)/3.);
      if(ichan<0)
	{
	  F1 =-K1fact*BrhoF123(s1,-1);F2 = K1fact*BKstarF123(s2,-1);
	  F5 =-BKstarF123(q2,-1)*FKrho(s2,s1,-1)*sqrt(2.);
	}
      else if(ichan%8==0){F1 =-K1fact*BrhoF123(s1,ichan/8);}
      else if(ichan%8==1){F2 = K1fact*BKstarF123(s2,(ichan-1)/8);}
      else{F5 = -BKstarF123(q2,ichan/8)*FKrho(s2,s1,(ichan-2)%8)*sqrt(2.);}
    }
  // calculate the pi- K0bar pi0
  else if(imode==7)
    {
      Complex K1fact(K1BreitWigner(q2));
      if(ichan<0){F2 =-K1fact*BrhoF123(s2,-1);F5 =-2.*BKstarF123(q2,-1)*FKrho(s1,s2,-1);}
      else if(ichan%7==0){F2 =-K1fact*BrhoF123(s2,ichan/7);}
      else {F5 =-2.*BKstarF123(q2,ichan/7)*FKrho(s1,s2,(ichan-1)%7);}
    }
  // calculate the pi- pi0 eta
  else if(imode==8)
    {
      if(ichan<0){F5 = BrhoF5(q2,-1)*BrhoF123(s3,-1)*sqrt(2./3.);}
      else{F5 = BrhoF5(q2,ichan/3)*BrhoF123(s3,ichan%3)*sqrt(2./3.);}
    }
  // multiply by the prefactors
  F1/=_fpi;F2/=_fpi;F3/=_fpi;F4/=_fpi;F5/=_fpi;
  F5 =-F5*Complex(0.,1.)/4./pi/pi/_fpi/_fpi;
}

// complete the construction of the decay mode for integration
bool ThreeMesonDefaultCurrent::createMode(int icharge, unsigned int imode,
					  DecayPhaseSpaceModePtr mode,
					  unsigned int iloc,unsigned int ires,
					  DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  int iq(0),ia(0);
  bool kineallowed(true);
  if(!acceptMode(imode)){return false;}
  PDVector extpart(particles(1,imode,iq,ia));
  Energy min(0.);
  for(unsigned int ix=0;ix<extpart.size();++ix){min+=extpart[ix]->massMin();}
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return kineallowed;}
  // the particles we will use a lot
  tPDPtr a1,k1;
  if(icharge==-3)
    {
      a1=getParticleData(ParticleID::a_1minus);
      k1=getParticleData(ParticleID::K_1minus);
    }
  else if(icharge==3)
    {
      a1=getParticleData(ParticleID::a_1plus);
      k1=getParticleData(ParticleID::K_1plus);
    }
  else
    {return false;}
  // the rho0 resonances
  tPDPtr rho0[3]={getParticleData(113),getParticleData(100113),getParticleData(30113)};
  tPDPtr rhoc[3],Kstar0[3],Kstarc[3];
  if(icharge==-3)
    {
      // the charged rho resonances
      rhoc[0] = getParticleData(-213);
      rhoc[1] = getParticleData(-100213);
      rhoc[2] = getParticleData(-30213);
      // the K*0 resonances
      Kstar0[0] = getParticleData(313);
      Kstar0[1] = getParticleData(100313);
      Kstar0[2] = getParticleData(30313);
      // the charged K* resonances
      Kstarc[0] = getParticleData(-323);
      Kstarc[1] = getParticleData(-100323);
      Kstarc[2] = getParticleData(-30323);
    }
  else
    {
      // the charged rho resonances
      rhoc[0] = getParticleData(213);
      rhoc[1] = getParticleData(100213);
      rhoc[2] = getParticleData(30213);
      // the K*0 resonances
      Kstar0[0] = getParticleData(-313);
      Kstar0[1] = getParticleData(-100313);
      Kstar0[2] = getParticleData(-30313);
      // the charged K* resonances
      Kstarc[0] = getParticleData(323);
      Kstarc[1] = getParticleData(100323);
      Kstarc[2] = getParticleData(30323);
    }
  DecayPhaseSpaceChannelPtr newchannel;
  if(imode==0)
    {
      // channels for pi- pi- pi+
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==1)
    {
      // channels for pi0 pi0 pi-
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==2)
    {
      // channels for K- pi- K+
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0,iloc,iloc+2);
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(rho0[iy],0,0.0,iloc,iloc+2);
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==3)
    {
      // channels for K0 pi- K0bar
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1        ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix]  ,0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(Kstarc[iy],0,0.0, iloc+1,iloc+2);
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(rho0[iy],0,0.0, iloc,iloc+2);
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==4)
    {
      // channels for K- pi0 K0
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(a1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==5)
    {  
      // channels for pi0 pi0 K-
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1        ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1        ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(Kstarc[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	}
    }
  else if(imode==6)
    {
      // channels for K- pi- pi+
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1      ,0,0.0,-ires-1,iloc);
	  newchannel->addIntermediate(rho0[ix],0,0.0, iloc+1,iloc+2);
	  mode->addChannel(newchannel);
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1        ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(Kstar0[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(rho0[iy]  ,0,0.0, iloc+1,iloc+2);
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc,iloc+2);
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==7)
    {
      // channels for pi- kbar0 pi0
      for(unsigned int ix=0;ix<3;++ix)
	{
	  newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	  newchannel->addIntermediate(k1      ,0,0.0,-ires-1,iloc+1);
	  newchannel->addIntermediate(rhoc[ix],0,0.0, iloc,iloc+2);
	  mode->addChannel(newchannel);
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(Kstar0[iy],0,0.0, iloc+1,iloc+2);
	      mode->addChannel(newchannel);
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(Kstarc[ix],0,0.0,-ires-1,iloc+1);
	      newchannel->addIntermediate(rhoc[iy]  ,0,0.0, iloc,iloc+2);
	      mode->addChannel(newchannel);
	    }
	}
    }
  else if(imode==8)
    {
      // channels for pi- pi0 eta
      for(unsigned int ix=0;ix<3;++ix)
	{
	  for(unsigned int iy=0;iy<3;++iy)
	    {
	      newchannel= new_ptr(DecayPhaseSpaceChannel(*phase));
	      newchannel->addIntermediate(rhoc[ix],0,0.0,-ires-1,iloc);
	      newchannel->addIntermediate(rho0[iy],0,0.0, iloc+1,iloc+2);
	      mode->addChannel(newchannel);
	    }
	}
    }
  if(_rhoparameters)
    {
      for(unsigned int ix=0;ix<_rhoF123masses.size();++ix)
	{
	  mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
	  mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
	}
    }
  // K star parameters in the base class
  if(_Kstarparameters)
    {
      for(unsigned int ix=0;ix<_KstarF123masses.size();++ix)
	{
	  mode->resetIntermediate(Kstarc[ix],_KstarF123masses[ix],_KstarF123widths[ix]);
	  mode->resetIntermediate(Kstar0[ix],_KstarF123masses[ix],_KstarF123widths[ix]);
	}
    }
  return kineallowed;
}


PDVector ThreeMesonDefaultCurrent::particles(int icharge, unsigned int imode,int iq,
					       int ia)
{
  PDVector extpart(3);
  if(imode==0)
    {
      extpart[0]=getParticleData(ParticleID::piminus);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::piplus);
    }
  else if(imode==1)
    {
      extpart[0]=getParticleData(ParticleID::pi0);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::piminus);
    }
  else if(imode==2)
    {
      extpart[0]=getParticleData(ParticleID::Kminus);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::Kplus);
    }
  else if(imode==3)
    {
      extpart[0]=getParticleData(ParticleID::K0);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::Kbar0);
    }
  else if(imode==4)
    {
      extpart[0]=getParticleData(ParticleID::Kminus);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::K0);
    }
  else if(imode==5)
    {
      extpart[0]=getParticleData(ParticleID::pi0);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::Kminus);
    }
  else if(imode==6)
    {
      extpart[0]=getParticleData(ParticleID::Kminus);
      extpart[1]=getParticleData(ParticleID::piminus);
      extpart[2]=getParticleData(ParticleID::piplus);
    }
  else if(imode==7)
    {
      extpart[0]=getParticleData(ParticleID::piminus);
      extpart[1]=getParticleData(ParticleID::Kbar0);
      extpart[2]=getParticleData(ParticleID::pi0);
    }
  else if(imode==8)
    {
      extpart[0]=getParticleData(ParticleID::piminus);
      extpart[1]=getParticleData(ParticleID::pi0);
      extpart[2]=getParticleData(ParticleID::eta);
    }
  // conjugate the particles if needed
  if(icharge==3)
    {for(unsigned int ix=0;ix<3;++ix)
	{if(extpart[0]->CC()){extpart[0]=extpart[0]->CC();}}}
  // return the answer
  return extpart;
}

void ThreeMesonDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::ThreeMesonDefaultCurrent " << fullName() << " \n";}
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":F123RhoWeight " << ix 
		      << " " << _rhoF123wgts[ix] << "\n";}
      else{output << "insert " << fullName() << ":F123RhoWeight " << ix 
		  << " " << _rhoF123wgts[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_KstarF123wgts.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":F123KstarWeight " << ix 
		      << " " << _KstarF123wgts[ix] << "\n";}
      else{output << "insert " << fullName() << ":F123KstarWeight " << ix 
		  << " " << _KstarF123wgts[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_rhoF5wgts.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":F5RhoWeight " << ix 
		      << " " << _rhoF5wgts[ix] << "\n";}
      else{output << "insert " << fullName() << ":F5RhoWeight " << ix 
		      << " " << _rhoF5wgts[ix] << "\n";}
    }
  for(unsigned int ix=0;ix<_KstarF5wgts.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":F5KstarWeight " << ix 
		      << " " << _KstarF5wgts[ix] << "\n";}
      else{output << "insert " << fullName() << ":F5KstarWeight " << ix 
		      << " " << _KstarF5wgts[ix] << "\n";}
    }
  output << "set " << fullName() << ":RhoKstarWgt " << _rhoKstarwgt << "\n";
  output << "set " << fullName() << ":Initializea1 " << _initializea1 << "\n";
  output << "set " << fullName() << ":RhoParameters " << _rhoparameters << "\n";
  output << "set " << fullName() << ":KstarParameters " << _Kstarparameters << "\n";
  output << "set " << fullName() << ":a1Parameters " << _a1parameters << "\n";
  output << "set " << fullName() << ":K1Parameters " << _K1parameters << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix)
    {output << "set " << fullName() << ":a1RunningWidth " << ix 
	    << " " << _a1runwidth[ix]/GeV << "\n";}
  for(unsigned int ix=0;ix<_a1runq2.size();++ix)
    {output << "set " << fullName() << ":a1RunningQ2 " << ix 
	    << " " << _a1runq2[ix]/GeV2 << "\n";}
  output << "set " << fullName() << ":A1Width " << _a1width/GeV << "\n";
  output << "set " << fullName() << ":A1Mass " << _a1mass/GeV << "\n";
  output << "set " << fullName() << ":K1Width " << _K1width/GeV << "\n";
  output << "set " << fullName() << ":K1Mass " << _K1mass/GeV << "\n";
  output << "set " << fullName() << ":FPi " << _fpi/MeV << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":rhoF123masses " << ix 
		      << " " << _rhoF123masses[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":rhoF123masses " << ix 
		      << " " << _rhoF123masses[ix]/GeV << "\n";}
    }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":rhoF123widths " << ix 
		      << " " << _rhoF123widths[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":rhoF123widths " << ix 
		  << " " << _rhoF123widths[ix]/GeV << "\n";}
    }
  for(unsigned int ix=0;ix<_rhoF5masses.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":rhoF5masses " << ix 
		      << " " << _rhoF5masses[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":rhoF5masses " << ix 
		      << " " << _rhoF5masses[ix]/GeV << "\n";}
    }
  for(unsigned int ix=0;ix<_rhoF5widths.size();++ix)
    {
      if(ix<3){output << "set " << fullName() << ":rhoF5widths " << ix 
		      << " " << _rhoF5widths[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":rhoF5widths " << ix 
		      << " " << _rhoF5widths[ix]/GeV << "\n";}
    }
  for(unsigned int ix=0;ix<_KstarF123masses.size();++ix)
    {
      if(ix<1){output << "set " << fullName() << ":KstarF123masses " << ix 
		      << " " << _KstarF123masses[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":KstarF123masses " << ix 
		      << " " << _KstarF123masses[ix]/GeV << "\n";}
    }
  for(unsigned int ix=0;ix<_KstarF123widths.size();++ix)
    {
      if(ix<1){output << "set " << fullName() << ":KstarF123widths " << ix 
		      << " " << _KstarF123widths[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":KstarF123widths " << ix 
		      << " " << _KstarF123widths[ix]/GeV << "\n";}
	}
  for(unsigned int ix=0;ix<_KstarF5masses.size();++ix)
    {
      if(ix<1){output << "set " << fullName() << ":KstarF5masses " << ix 
		      << " " << _KstarF5masses[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":KstarF5masses " << ix 
		      << " " << _KstarF5masses[ix]/GeV << "\n";}
    }
  for(unsigned int ix=0;ix<_KstarF5widths.size();++ix)
    {
      if(ix<1){output << "set " << fullName() << ":KstarF5widths " << ix 
		      << " " << _KstarF5widths[ix]/GeV << "\n";}
      else{output << "insert " << fullName() << ":KstarF5widths " << ix 
		      << " " << _KstarF5widths[ix]/GeV << "\n";}
    }
  ThreeMesonCurrentBase::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
  
}

// the functions for the integrands of the a_1 width
namespace Herwig {
using namespace Genfun;

FUNCTION_OBJECT_IMP(Defaulta1MatrixElement)

Defaulta1MatrixElement::Defaulta1MatrixElement(Ptr<Herwig::ThreeMesonDefaultCurrent>::pointer in)
  {_decayer=in;}

  
Defaulta1MatrixElement::Defaulta1MatrixElement(const Defaulta1MatrixElement & right)  {  }
  
unsigned int Defaulta1MatrixElement::dimensionality() const {return 7;}

double Defaulta1MatrixElement::operator ()(const Argument & a) const 
{return _decayer->a1MatrixElement(a[0],a[1],a[2],a[3],a[4],a[5],a[6]);}

}
