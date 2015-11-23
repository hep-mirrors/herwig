// -*- C++ -*-
//
// FourPionNovosibirskCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourPionNovosibirskCurrent class.
//

#include "FourPionNovosibirskCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include <functional>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

DescribeClass<FourPionNovosibirskCurrent,WeakDecayCurrent>
describeHerwigFourPionNovosibirskCurrent("Herwig::FourPionNovosibirskCurrent",
					 "HwWeakCurrents.so");

HERWIG_INTERPOLATOR_CLASSDESC(FourPionNovosibirskCurrent1,Energy,Energy2)

HERWIG_INTERPOLATOR_CLASSDESC(FourPionNovosibirskCurrent2,double,Energy)

HERWIG_INTERPOLATOR_CLASSDESC(FourPionNovosibirskCurrent3,double,Energy2)

namespace {
  inline Energy  timesGeV (double x) { return x * GeV; }
  inline Energy2 timesGeV2(double x) { return x * GeV2; }
}
 
IBPtr FourPionNovosibirskCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr FourPionNovosibirskCurrent::fullclone() const {
  return new_ptr(*this);
}

void FourPionNovosibirskCurrent::doupdate() {
  WeakDecayCurrent::doupdate();
  // update running width if needed
  if ( !touched() ) return;
  if(_maxmass!=_maxcalc) inita1width(-1);
}

FourPionNovosibirskCurrent::FourPionNovosibirskCurrent() : _mpic(), _mpi0(),
							   _mpic2(), _mpi02(), _prho()
{
  // set the number of modes
  addDecayMode(2,-1);
  addDecayMode(2,-1);
  setInitialModes(2);
  // masses of the particles used in the current
  _rhomass    = 0.7761*GeV;
  _a1mass     = 1.2300*GeV;
  _omegamass  = 0.7820*GeV;
  _sigmamass  = 0.8000*GeV;
  // widths of the particles used in the current
  _rhowidth   = 0.14450*GeV;
  _a1width    = 0.45000*GeV;
  _omegawidth = 0.00841*GeV;
  _sigmawidth = 0.80000*GeV;
  // parameters for the resonance used in the integration
  _intmass = 1.4*GeV;
  _intwidth = 0.5*GeV;
  // relative coupling of the sigma and rho in the a_1 decay
  _zmag   = 1.3998721;
  _zphase = 0.43585036;
  _zsigma=0.;
  // parameter for f_a1
  _lambda2 = 1.2*GeV2;
  _onedlam2 = 1./_lambda2;
  _a1massolam2 = _a1mass*_a1mass*_onedlam2;
  _hm2=ZERO; 
  _rhoD=ZERO;
  _dhdq2m2=0.;
  // use local values of the parameters
  _localparameters=true;
  // magic numbers from TAUOLA (N.B. conversion from GeV to MeV)
  _athreec = 76.565/GeV;
  _bthreec = 0.71709;
  _cthreec = 0.27505;
  _aomega  = 886.84/GeV;
  _bomega  = 0.70983;
  _comega  = 0.26689;
  _aonec   = 96.867/GeV;
  _bonec   = 0.70907;
  _conec   = 0.26413;
  //parameters for the running omega width
  _omegaparam.resize(10);
  _omegaparam[0] = 17.560;
  _omegaparam[1] = 141.110;
  _omegaparam[2] = 894.884;
  _omegaparam[3] = 4977.35;
  _omegaparam[4] = 7610.66;
  _omegaparam[5] =-42524.4;
  _omegaparam[6] =-1333.26;
  _omegaparam[7] = 4860.19;
  _omegaparam[8] =-6000.81;
  _omegaparam[9] = 2504.97;
  // values of the g omega function from hep-ph/0201149
  double Fomegainit[98]={ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 2.2867811,
			  2.9710648, 2.9344304, 2.6913538, 2.5471206, 2.3557470,
			  2.2448280, 2.1074708, 2.0504866, 1.9270257, 1.8669430,
			  1.7907301, 1.7184515, 1.6535717, 1.6039416, 1.5535343,
			  1.5065620, 1.4608675, 1.4215596, 1.3849826, 1.3480113,
			  1.3147917, 1.2793381, 1.2487282, 1.2184237, 1.1952927,
			  1.1683835, 1.1458827, 1.1145806, 1.0935820, 1.0608720,
			  1.0390474, 1.0164336, 0.9908721, 0.9585276, 0.9307971,
			  0.9017274, 0.8731154, 0.8452763, 0.8145532, 0.7817339,
			  0.7493086, 0.7199919, 0.6887290, 0.6568120, 0.6255773,
			  0.5944664, 0.5661956, 0.5391204, 0.5102391, 0.4786543,
			  0.4546428, 0.4316779, 0.4063754, 0.3769831, 0.3561141,
			  0.3333555, 0.3139160, 0.2949214, 0.2814728, 0.2602444,
			  0.2349602, 0.2269845, 0.2192318, 0.2286938, 0.2839763,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000};
  // values of the three charged pion G function from hep-ph/0201149
  double Fthreeinit[98]={ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  13.1664906,10.7234087, 8.8219614,10.7989664, 9.1883001,
			  7.8526378, 7.7481031, 8.2633696, 5.5042820, 4.9029269,
			  4.4794345, 3.9654009, 4.5254011, 3.6509495, 3.5005512,
			  3.2274280, 3.1808922, 2.9925177, 2.6886659, 2.5195024,
			  2.4678771, 2.3540580, 2.2123868, 2.1103525, 2.0106986,
			  1.8792295, 1.8250662, 1.7068460, 1.6442842, 1.5503920,
			  1.4814349, 1.4225838, 1.3627135, 1.3205355, 1.2784383,
			  1.2387408, 1.1975995, 1.1633024, 1.1318133, 1.1114354,
			  1.0951439, 1.0691465, 1.0602311, 1.0392803, 1.0220672,
			  1.0154786, 1.0010130, 0.9908018, 0.9710845, 0.9602382,
			  0.9488459, 0.9316537, 0.9118049, 0.8920435, 0.8719332,
			  0.8520256, 0.8280582, 0.8064085, 0.7767881, 0.7570597,
			  0.7382626, 0.7100251, 0.6846500, 0.6666913, 0.6372250,
			  0.6162248, 0.6007728, 0.5799103, 0.5674670, 0.5446148,
			  0.5352115, 0.5128809, 0.4932536, 0.5310397, 0.8566489,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000};
  double   Foneinit[98]={ 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  1.4819183, 1.7086354, 1.6958492, 1.6172935, 1.6301320,
			  1.5719706, 1.5459771, 1.5377471, 1.5008864, 1.4924121,
			  1.4720882, 1.4371741, 1.3990080, 1.3879193, 1.4030601,
			  1.3768673, 1.3493533, 1.3547127, 1.3275831, 1.3167892,
			  1.3035913, 1.2968298, 1.2801558, 1.2650299, 1.2557997,
			  1.2325822, 1.2210644, 1.1935984, 1.1746194, 1.1510350,
			  1.1358515, 1.1205584, 1.1010553, 1.0903869, 1.0731295,
			  1.0578678, 1.0438409, 1.0377911, 1.0253277, 1.0103551,
			  1.0042409, 0.9937978, 0.9858117, 0.9770346, 0.9724492,
			  0.9656686, 0.9606671, 0.9525813, 0.9488522, 0.9417335,
			  0.9399430, 0.9323438, 0.9281269, 0.9244171, 0.9237418,
			  0.9174354, 0.9177181, 0.9120840, 0.9047825, 0.9065579,
			  0.9034142, 0.8992961, 0.9011586, 0.9036470, 0.8954964,
			  0.8898208, 0.8911991, 0.8854824, 0.8888282, 0.8868449,
			  0.9004632, 0.8981572, 0.9096183, 0.9046990, 1.7454215,
			  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000,
			  0.0000000, 0.0000000, 0.0000000};
  // eninit in GeV
  double     eninit[98]={ 0.6000000, 0.6131313, 0.6262626, 0.6393939, 0.6525252,
			  0.6656566, 0.6787879, 0.6919192, 0.7050505, 0.7181818,
			  0.7313131, 0.7444444, 0.7575758, 0.7707071, 0.7838384,
			  0.7969697, 0.8101010, 0.8232324, 0.8363636, 0.8494949,
			  0.8626263, 0.8757576, 0.8888889, 0.9020202, 0.9151515,
			  0.9282829, 0.9414141, 0.9545454, 0.9676768, 0.9808081,
			  0.9939394, 1.0070707, 1.0202020, 1.0333333, 1.0464647,
			  1.0595959, 1.0727273, 1.0858586, 1.0989898, 1.1121212,
			  1.1252525, 1.1383839, 1.1515151, 1.1646465, 1.1777778,
			  1.1909091, 1.2040404, 1.2171717, 1.2303030, 1.2434343,
			  1.2565657, 1.2696970, 1.2828283, 1.2959596, 1.3090909,
			  1.3222222, 1.3353535, 1.3484849, 1.3616161, 1.3747475,
			  1.3878788, 1.4010102, 1.4141414, 1.4272727, 1.4404041,
			  1.4535353, 1.4666667, 1.4797980, 1.4929293, 1.5060606,
			  1.5191919, 1.5323232, 1.5454545, 1.5585859, 1.5717171,
			  1.5848485, 1.5979798, 1.6111112, 1.6242424, 1.6373737,
			  1.6505051, 1.6636363, 1.6767677, 1.6898990, 1.7030303,
			  1.7161616, 1.7292930, 1.7424242, 1.7555555, 1.7686869,
			  1.7818182, 1.7949495, 1.8080808, 1.8212122, 1.8343434,
			  1.8474747, 1.8606061, 1.8737373};          
  // ensigma in GeV2
  double   ensigma[100]={ 0.2916000, 0.3206586, 0.3497172, 0.3787757, 0.4078344,
			  0.4368929, 0.4659515, 0.4950101, 0.5240687, 0.5531273,
			  0.5821859, 0.6112444, 0.6403030, 0.6693616, 0.6984202,
			  0.7274788, 0.7565374, 0.7855960, 0.8146545, 0.8437131,
			  0.8727717, 0.9018303, 0.9308889, 0.9599475, 0.9890060,
			  1.0180646, 1.0471232, 1.0761818, 1.1052403, 1.1342990,
			  1.1633576, 1.1924162, 1.2214748, 1.2505333, 1.2795919,
			  1.3086505, 1.3377091, 1.3667676, 1.3958262, 1.4248848,
			  1.4539435, 1.4830021, 1.5120606, 1.5411192, 1.5701778,
			  1.5992364, 1.6282949, 1.6573535, 1.6864121, 1.7154707,
			  1.7445292, 1.7735878, 1.8026465, 1.8317051, 1.8607637,
			  1.8898222, 1.9188808, 1.9479394, 1.9769980, 2.0060565,
			  2.0351152, 2.0641737, 2.0932324, 2.1222908, 2.1513495,
			  2.1804080, 2.2094667, 2.2385252, 2.2675838, 2.2966425,
			  2.3257010, 2.3547597, 2.3838181, 2.4128768, 2.4419353,
			  2.4709940, 2.5000525, 2.5291111, 2.5581696, 2.5872283,
			  2.6162868, 2.6453454, 2.6744041, 2.7034626, 2.7325213,
			  2.7615798, 2.7906384, 2.8196969, 2.8487556, 2.8778141,
			  2.9068727, 2.9359312, 2.9649899, 2.9940486, 3.0231071,
			  3.0521657, 3.0812242, 3.1102829, 3.1393414, 3.1684000};
  double    Fsigma[100]={ 2.0261996, 2.2349865, 2.4839740, 2.7840748, 3.1488798,
			  3.5936222, 4.1301847, 4.7517977, 5.3984838, 5.9147439,
			  6.0864558, 5.8283591, 5.2841811, 4.6615186, 4.0839195,
			  3.5914702, 3.1841860, 2.8494759, 2.5732665, 2.3434010,
			  2.1502059, 1.9862038, 1.8456544, 1.7241427, 1.6182493,
			  1.5253036, 1.4432002, 1.3702650, 1.3051554, 1.2467849,
			  1.1942677, 1.1468738, 1.1039963, 1.0651271, 1.0298390,
			  0.9977714, 0.9686196, 0.9421255, 0.9180685, 0.8962603,
			  0.8765374, 0.8587573, 0.8427927, 0.8285285, 0.8158574,
			  0.8046767, 0.7948853, 0.7863811, 0.7790571, 0.7728010,
			  0.7674922, 0.7630011, 0.7591889, 0.7559078, 0.7530031,
			  0.7503147, 0.7476809, 0.7449428, 0.7419487, 0.7385587,
			  0.7346500, 0.7301207, 0.7248930, 0.7189151, 0.7121620,
			  0.7046344, 0.6963565, 0.6873729, 0.6777444, 0.6675445,
			  0.6568548, 0.6457604, 0.6343476, 0.6227004, 0.6108983,
			  0.5990148, 0.5871165, 0.5752623, 0.5635037, 0.5518846,
			  0.5404415, 0.5292045, 0.5181981, 0.5074410, 0.4969472,
			  0.4867267, 0.4767860, 0.4671288, 0.4577557, 0.4486661,
			  0.4398569, 0.4313242, 0.4230627, 0.4150662, 0.4073282,
			  0.3998415, 0.3925985, 0.3855914, 0.3788125, 0.3722538};
  // set up the interpolators
  _Fomega  = make_InterpolatorPtr( 98,Fomegainit,1.0,eninit, GeV, 1);
  _Fthreec = make_InterpolatorPtr( 98,Fthreeinit,1.0,eninit, GeV, 1);
  _Fonec   = make_InterpolatorPtr( 98,Foneinit  ,1.0,eninit, GeV, 1);
  _Fsigma  = make_InterpolatorPtr(100,Fsigma    ,1.0,ensigma,GeV2,1);
  // initialise the calculation of the a_1 width
  _initializea1=false;
  // in GeV2
  double  a1q2in[200]={0,15788.6,31577.3,47365.9,63154.6,78943.2,94731.9,110521,
		       126309,142098,157886,173675,189464,205252,221041,236830,252618,
		       268407,284196,299984,315773,331562,347350,363139,378927,394716,
		       410505,426293,442082,457871,473659,489448,505237,521025,536814,
		       552603,568391,584180,599969,615757,631546,647334,663123,678912,
		       694700,710489,726278,742066,757855,773644,789432,805221,821010,
		       836798,852587,868375,884164,899953,915741,931530,947319,963107,
		       978896,994685,1.01047e+06,1.02626e+06,1.04205e+06,1.05784e+06,
		       1.07363e+06,1.08942e+06,1.10521e+06,1.12099e+06,1.13678e+06,
		       1.15257e+06,1.16836e+06,1.18415e+06,1.19994e+06,1.21573e+06,
		       1.23151e+06,1.2473e+06,1.26309e+06,1.27888e+06,1.29467e+06,
		       1.31046e+06,1.32625e+06,1.34203e+06,1.35782e+06,1.37361e+06,
		       1.3894e+06,1.40519e+06,1.42098e+06,1.43677e+06,1.45256e+06,
		       1.46834e+06,1.48413e+06,1.49992e+06,1.51571e+06,1.5315e+06,
		       1.54729e+06,1.56308e+06,1.57886e+06,1.59465e+06,1.61044e+06,
		       1.62623e+06,1.64202e+06,1.65781e+06,1.6736e+06,1.68939e+06,
		       1.70517e+06,1.72096e+06,1.73675e+06,1.75254e+06,1.76833e+06,
		       1.78412e+06,1.79991e+06,1.81569e+06,1.83148e+06,1.84727e+06,
		       1.86306e+06,1.87885e+06,1.89464e+06,1.91043e+06,1.92621e+06,
		       1.942e+06,1.95779e+06,1.97358e+06,1.98937e+06,2.00516e+06,
		       2.02095e+06,2.03674e+06,2.05252e+06,2.06831e+06,2.0841e+06,
		       2.09989e+06,2.11568e+06,2.13147e+06,2.14726e+06,2.16304e+06,
		       2.17883e+06,2.19462e+06,2.21041e+06,2.2262e+06,2.24199e+06,
		       2.25778e+06,2.27356e+06,2.28935e+06,2.30514e+06,2.32093e+06,
		       2.33672e+06,2.35251e+06,2.3683e+06,2.38409e+06,2.39987e+06,
		       2.41566e+06,2.43145e+06,2.44724e+06,2.46303e+06,2.47882e+06,
		       2.49461e+06,2.51039e+06,2.52618e+06,2.54197e+06,2.55776e+06,
		       2.57355e+06,2.58934e+06,2.60513e+06,2.62092e+06,2.6367e+06,
		       2.65249e+06,2.66828e+06,2.68407e+06,2.69986e+06,2.71565e+06,
		       2.73144e+06,2.74722e+06,2.76301e+06,2.7788e+06,2.79459e+06,
		       2.81038e+06,2.82617e+06,2.84196e+06,2.85774e+06,2.87353e+06,
		       2.88932e+06,2.90511e+06,2.9209e+06,2.93669e+06,2.95248e+06,
		       2.96827e+06,2.98405e+06,2.99984e+06,3.01563e+06,3.03142e+06,
		       3.04721e+06,3.063e+06,3.07879e+06,3.09457e+06,3.11036e+06,
		       3.12615e+06,3.14194e+06};
  // in GeV
  double a1widthin[200]={0,0,0,0,0,0,0,0,
			 0,0,0,0,0.000634625,0.00686721,0.026178,0.066329,
			 0.134996,0.239698,0.387813,0.586641,0.843471,1.16567,
			 1.56076,2.03654,2.60115,3.26324,4.03202,4.91749,
			 5.93053,7.08313,8.38858,9.86176,11.5194,13.3805,
			 15.4667,17.8029,20.4175,23.3438,26.6202,30.2917,
			 34.4108,39.0384,44.2457,50.1143,56.7369,64.2147,
			 72.6566,82.1666,92.8329,104.708,117.786,131.981,
			 147.124,162.974,179.244,195.64,211.904,227.818,
			 243.223,257.991,272.06,285.392,297.971,309.8,
			 320.894,331.278,340.979,350.03,358.463,366.31,
			 373.608,380.386,386.677,392.511,397.945,402.935,
			 407.563,411.841,415.79,419.433,422.766,425.853,
			 428.695,431.302,433.715,435.883,437.887,439.716,
			 441.426,442.947,444.326,445.575,446.65,447.666,
			 448.578,449.395,450.123,450.821,451.343,451.847,
			 452.283,452.859,452.987,453.266,453.496,453.686,
			 453.839,453.958,454.05,454.113,454.149,454.16,
			 454.154,454.13,454.091,454.037,453.966,453.9,
			 453.814,453.724,453.628,453.528,453.417,453.314,
			 453.206,453.1,452.995,452.891,452.79,452.697,
			 452.598,452.509,452.423,452.343,452.269,452.201,
			 452.141,452.086,452.039,452.004,451.966,451.941,
			 451.926,451.888,451.919,451.928,451.945,451.971,
			 452.006,452.05,452.102,452.163,452.234,452.421,
			 452.401,452.498,452.605,452.718,452.84,452.971,
			 453.111,453.261,453.417,453.583,453.756,453.937,
			 454.126,454.324,454.529,455.023,454.964,455.719,
			 455.428,455.671,455.921,456.179,456.444,456.695,
			 456.996,457.276,457.57,457.867,458.171,458.478,
			 458.793,459.115,459.442,459.777,460.115,460.458,
			 460.809,461.161,461.52,461.884,462.253,462.626,
			 463.004,463.832,463.778,464.166};
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

  _maxmass=ZERO;
  _maxcalc=ZERO;
}

void FourPionNovosibirskCurrent::doinit() {
  WeakDecayCurrent::doinit();
  // pion masses
  _mpic=getParticleData(ParticleID::piplus)->mass();
  _mpic2=sqr(_mpic);
  _mpi0=getParticleData(ParticleID::pi0)->mass();
  _mpi02=sqr(_mpi0);
  if(!_localparameters) {
    _rhomass    = getParticleData(ParticleID::rhominus)->mass();
    _rhowidth   = getParticleData(ParticleID::rhominus)->width();
    _omegamass  = getParticleData(ParticleID::omega)->mass();
    _omegawidth = getParticleData(ParticleID::omega)->width();
    _sigmamass  = getParticleData(9000221)->mass();
    _sigmawidth = getParticleData(9000221)->width();
    _a1mass    = getParticleData(ParticleID::a_1minus)->mass();
    _a1width   = getParticleData(ParticleID::a_1minus)->width();
  }
  // calculate the constants for the a_1 form factor
  _onedlam2 = 1./_lambda2;
  _a1massolam2 = _a1mass*_a1mass*_onedlam2;
  // parameter for the sigma breit-wigner
  _psigma.push_back(Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi0,_mpi0));
  _psigma.push_back(Kinematics::pstarTwoBodyDecay(_sigmamass,_mpic,_mpic));
  // parameters for the rho breit wigner
  _prho=Kinematics::pstarTwoBodyDecay(_rhomass,_mpic,_mpic);
  _hm2 = hFunction(_rhomass);
  _dhdq2m2=dhdq2Parameter();
  _rhoD=DParameter();
  // convert the magnitude and phase of z into a phase
  _zsigma = _zmag*(cos(_zphase)+Complex(0.,1.)*sin(_zphase));
  // initialize the a_1 width
  inita1width(-1);
}

void FourPionNovosibirskCurrent::doinitrun() {
  // set up the running a_1 width
  inita1width(0);
  WeakDecayCurrent::doinitrun();
}

void FourPionNovosibirskCurrent::persistentOutput(PersistentOStream & os) const {
  os << _a1runinter << _Fomega << _Fthreec << _Fonec << _Fsigma
     << ounit(_rhomass,GeV) << ounit(_a1mass,GeV) << ounit(_omegamass,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_rhowidth,GeV) << ounit(_a1width,GeV)
     << ounit(_omegawidth,GeV) << ounit(_sigmawidth,GeV) 
     << _zsigma << ounit(_lambda2,GeV2)
     << _initializea1 << _localparameters 
     << ounit(_a1runwidth,GeV) << ounit(_a1runq2,GeV2) << ounit(_onedlam2,1/GeV2) 
     << _a1massolam2 << ounit(_psigma,GeV) << ounit(_mpic,GeV) << ounit(_mpi0,GeV)
     << ounit(_aomega,1/GeV) << ounit(_athreec,1/GeV) << ounit(_aonec,1/GeV) 
     << _bomega << _bthreec << _bonec 
     << _comega << _cthreec <<_conec << _omegaparam 
     << ounit(_intwidth,GeV) << ounit(_intmass,GeV)
     << ounit(_mpic2,GeV2) << ounit(_mpi02,GeV2) << ounit(_hm2,GeV2) << _dhdq2m2 
     << ounit(_prho,GeV) << ounit(_rhoD,GeV2) << _zmag << _zphase
     << ounit(_maxmass,GeV) << ounit(_maxcalc,GeV);
}

void FourPionNovosibirskCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _a1runinter >> _Fomega >> _Fthreec >> _Fonec >> _Fsigma
     >> iunit(_rhomass,GeV) >> iunit(_a1mass,GeV) >> iunit(_omegamass,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_rhowidth,GeV) >> iunit(_a1width,GeV)
     >> iunit(_omegawidth,GeV) >> iunit(_sigmawidth,GeV) 
     >> _zsigma >> iunit(_lambda2,GeV2)
     >> _initializea1 >> _localparameters 
     >> iunit(_a1runwidth,GeV) >> iunit(_a1runq2,GeV2) >> iunit(_onedlam2,1/GeV2)
     >> _a1massolam2 >> iunit(_psigma,GeV) >> iunit(_mpic,GeV) >> iunit(_mpi0,GeV)
     >> iunit(_aomega,1/GeV) >> iunit(_athreec,1/GeV) >> iunit(_aonec,1/GeV) 
     >> _bomega >> _bthreec >> _bonec 
     >> _comega >> _cthreec >>_conec >> _omegaparam 
     >> iunit(_intwidth,GeV) >> iunit(_intmass,GeV)
     >> iunit(_mpic2,GeV2) >> iunit(_mpi02,GeV2)>> iunit(_hm2,GeV2) >> _dhdq2m2
     >> iunit(_prho,GeV) >> iunit(_rhoD,GeV2) >> _zmag >> _zphase
     >> iunit(_maxmass,GeV) >> iunit(_maxcalc,GeV);
}

// Definition of the static class description member.

void FourPionNovosibirskCurrent::Init() {

  static ClassDocumentation<FourPionNovosibirskCurrent> documentation
    ("The FourPionNovosibirskCurrent class performs the decay"
     " of the tau to four pions using currents based on the the"
     " Novosibirsk e+e- data",
     "The decay of the tau to four pions uses currents based on \\cite{Bondar:2002mw}.",
     "%\\cite{Bondar:2002mw}\n"
     "\\bibitem{Bondar:2002mw}\n"
     "  A.~E.~Bondar, S.~I.~Eidelman, A.~I.~Milstein, T.~Pierzchala, N.~I.~Root, Z.~Was and M.~Worek,\n"
     "   ``Novosibirsk hadronic currents for tau --> 4pi channels of tau decay\n"
     "  %library TAUOLA,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 146}, 139 (2002)\n"
     "  [arXiv:hep-ph/0201149].\n"
     "  %%CITATION = CPHCB,146,139;%%\n"
     );

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacerhoMass
    ("rhoMass",
     "The local value of the rho mass",
     &FourPionNovosibirskCurrent::_rhomass, GeV,0.7761*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacea1mass
    ("a1Mass",
     "The local value of the square of the a_1 mass",
     &FourPionNovosibirskCurrent::_a1mass, GeV, 1.2300*GeV, 0.5*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceSigmaMass
    ("sigmaMass",
     "The local value of the sigma mass",
     &FourPionNovosibirskCurrent::_sigmamass, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceOmegaMass
    ("omegaMass",
     "The local value of the omega mass",
     &FourPionNovosibirskCurrent::_omegamass, GeV, 0.7820*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacerhoWidth
    ("rhoWidth",
     "The local value of the rho width",
     &FourPionNovosibirskCurrent::_rhowidth, GeV,0.1445*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfacea1width
    ("a1Width",
     "The local value of the square of the a_1 width",
     &FourPionNovosibirskCurrent::_a1width, GeV, 0.45*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceSigmaWidth
    ("sigmaWidth",
     "The local value of the sigma width",
     &FourPionNovosibirskCurrent::_sigmawidth, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceOmegaWidth
    ("omegaWidth",
     "The local value of the omega width",
     &FourPionNovosibirskCurrent::_omegawidth, GeV, 0.00841*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceIntegrationMass
    ("IntegrationMass",
     "Mass of the pseudoresonance used to improve integration effciency",
     &FourPionNovosibirskCurrent::_intmass, GeV, 1.4*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy> interfaceIntegrationWidth
    ("IntegrationWidth",
     "Width of the pseudoresonance used to improve integration effciency",
     &FourPionNovosibirskCurrent::_intwidth, GeV, 0.5*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,double> interfaceSigmaMagnitude
    ("SigmaMagnitude",
     "magnitude of the relative sigma coupling",
     &FourPionNovosibirskCurrent::_zmag, 1.3998721, 0.0, 10.0e20,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,double> interfaceSigmaPhase
    ("SigmaPhase",
     "phase of the relative sigma coupling",
     &FourPionNovosibirskCurrent::_zphase, 0.43585036, 0.0, Constants::twopi,
     false, false, true);

  static Parameter<FourPionNovosibirskCurrent,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of the mass scale squared to use in the form-factor",
     &FourPionNovosibirskCurrent::_lambda2, GeV2, 1.2*GeV2, 0.0001*GeV2, 10.0*GeV2,
     false, false, true);

  static Switch<FourPionNovosibirskCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &FourPionNovosibirskCurrent::_localparameters, true, false, false);
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

  static Switch<FourPionNovosibirskCurrent,bool> interfaceInitializea1
    ("Initializea1",
     "Initialise the calculation of the a_1 running width",
     &FourPionNovosibirskCurrent::_initializea1, false, false, false);
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
  
  static ParVector<FourPionNovosibirskCurrent,Energy> interfacea1RunningWidth
    ("a1RunningWidth",
     "The values of the a_1 width for interpolation to giving the running width.",
     &FourPionNovosibirskCurrent::_a1runwidth, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<FourPionNovosibirskCurrent,Energy2> interfacea1RunningQ2
    ("a1RunningQ2",
     "The values of the q^2 for interpolation to giving the running width.",
     &FourPionNovosibirskCurrent::_a1runq2, GeV2, -1, 1.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

}

// initialisation of the a_1 running width 
void FourPionNovosibirskCurrent::inita1width(int iopt) {
  // set up the interpolator
  if(iopt==0||!_initializea1) {
    _a1runinter = make_InterpolatorPtr(_a1runwidth,_a1runq2,3);
    return;
  }
  _maxcalc=_maxmass;
  if(_maxmass==ZERO) return;
  // parameters for the table of values
  Energy2 step(sqr(_maxmass)/200.);
  // function to be integrated to give the matrix element
  // integrator to perform the integral
  // weights for the integration channels
  vector<double> inweights;
  inweights.push_back(0.3);inweights.push_back(0.3);inweights.push_back(0.3);
  vector<double> inpower(3, 0.0);
  // types of integration channels
  vector<int> intype;
  intype.push_back(2);intype.push_back(3);intype.push_back(1);
  // masses for the integration channels
  vector<Energy> inmass(2,_rhomass);inmass.push_back(_sigmamass);
  // widths for the integration channels
  vector<Energy> inwidth(2,_rhowidth);inwidth.push_back(_sigmawidth);
  ThreeBodyAllOnCalculator<FourPionNovosibirskCurrent> 
    widthgen1(inweights,intype,inmass,inwidth,inpower,*this,0,_mpi0,_mpic,_mpic); 
  ThreeBodyAllOnCalculator<FourPionNovosibirskCurrent>
    widthgen2(inweights,intype,inmass,inwidth,inpower,*this,1,_mpi0,_mpi0,_mpi0); 
  // normalisation constant to give physical width if on shell
  double a1const(_a1width/(widthgen1.partialWidth(sqr(_a1mass))+
			   widthgen2.partialWidth(sqr(_a1mass))));
  // loop to give the values
  Energy2 moff2(ZERO);
  _a1runwidth.clear();_a1runq2.clear();
  for(;moff2<=sqr(_maxmass);moff2+=step) {
    Energy total = a1const*(widthgen1.partialWidth(moff2)+widthgen2.partialWidth(moff2));
    _a1runwidth.push_back(total);
    _a1runq2.push_back(moff2);
  }
}

// complete the construction of the decay mode for integration
bool FourPionNovosibirskCurrent::createMode(int icharge, unsigned int imode,
					    DecayPhaseSpaceModePtr mode,
					    unsigned int iloc,unsigned int ires,
					    DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  // check the charge
  if(abs(icharge)!=3) return false;
  // check that the modes are kinematical allowed
  Energy min(ZERO);
  if(imode==0) {
    min=   getParticleData(ParticleID::piplus)->mass()
      +3.*getParticleData(ParticleID::pi0)->mass();
  }
  else {
    min=3.*getParticleData(ParticleID::piplus)->mass()
      +getParticleData(ParticleID::pi0)->mass();
  }
  if(min>upp) return false;
  _maxmass=max(upp,_maxmass);
  // intermediates for the channels
  tPDPtr omega(getParticleData(ParticleID::omega)),rhop,rhom,
    rho0(getParticleData(ParticleID::rho0)),a1m,a10(getParticleData(ParticleID::a_10)),
    sigma(getParticleData(9000221)),rhot;
  if(icharge==3) {
    rhop = getParticleData(ParticleID::rhominus);
    rhom = getParticleData(ParticleID::rhoplus);
    a1m  = getParticleData(ParticleID::a_1plus);
    rhot = getParticleData(24);
  }
  else {
    rhop = getParticleData(ParticleID::rhoplus);
    rhom = getParticleData(ParticleID::rhominus);
    a1m  = getParticleData(ParticleID::a_1minus);
    rhot = getParticleData(-24);
  }
  DecayPhaseSpaceChannelPtr newchannel;
  // the omega channels for the three charged pion mode
  // first  channel two channels with rho0
  if(imode==1) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
    newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rho0    ,0,0.0, iloc+1,iloc+2);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rho0    ,0,0.0, iloc,iloc+2);
    mode->addChannel(newchannel);
    // second two channels with rho -
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
    newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc+3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc,iloc+3);
    mode->addChannel(newchannel);
    // third two channels with rho +
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
    newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(omega   ,0,0.0,-ires-2,iloc);
    newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
    mode->addChannel(newchannel);
    //  a_1 channels with rhos
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
    newchannel->addIntermediate(rho0    ,0,0.0, iloc+1,iloc+2);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rho0    ,0,0.0, iloc,iloc+2);
    mode->addChannel(newchannel);
    // neutral a_1 channels with rhos
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
    newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc+3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc,iloc+3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
    newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc);
    newchannel->addIntermediate(rhop    ,0,0.0, iloc+2,iloc+3);
    mode->addChannel(newchannel);
    //  a_1 channels with sigmas
    if(sigma) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc,iloc+2);
      mode->addChannel(newchannel);
      // neutral a_1 channels with sigma
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc,iloc+2);
      mode->addChannel(newchannel);
    }
  }
  else {
    // channels with an a1- and a rho -
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+3,iloc);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+2,iloc);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+3,iloc);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+2);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+3);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+1);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+2,iloc);
    mode->addChannel(newchannel);
    newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
    newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
    newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc+2);
    newchannel->addIntermediate(rhom    ,0,0.0, iloc+1,iloc);
    mode->addChannel(newchannel);
    // channels with a sigma and a10
    if(sigma ) {
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+1);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+2,iloc+3);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+2);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+3);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc);
      newchannel->addIntermediate(a10     ,0,0.0,-ires-2,iloc+3);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      mode->addChannel(newchannel);
      // channels with a1- and sigma
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+1);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+2,iloc+3);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+2);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+3);
      mode->addChannel(newchannel);
      newchannel = new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(rhot    ,0,0.0,-ires-1,iloc+3);
      newchannel->addIntermediate(a1m     ,0,0.0,-ires-2,iloc);
      newchannel->addIntermediate(sigma   ,0,0.0, iloc+1,iloc+2);
      mode->addChannel(newchannel);
    }
  }
  // reset the parameters of the dummy resonance used for integration
  mode->resetIntermediate(rhot,_intmass,_intwidth);
  // reset the parameters of the resonances if using local values
  if(_localparameters) {
    mode->resetIntermediate(rhom,_rhomass,_rhowidth);
    mode->resetIntermediate(rhop,_rhomass,_rhowidth);
    mode->resetIntermediate(rho0,_rhomass,_rhowidth);
    mode->resetIntermediate(omega,_omegamass,_omegawidth);
    if(sigma) mode->resetIntermediate(sigma,_sigmamass,_sigmawidth);
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector FourPionNovosibirskCurrent::particles(int icharge, unsigned int imode,int,int) {
  if(abs(icharge)!=3) return tPDVector();
  tPDVector output(4);
  if(imode==1) {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::piplus);
    output[2]=getParticleData(ParticleID::piminus);
  }
  else {
    output[0]=getParticleData(ParticleID::piplus);
    output[1]=getParticleData(ParticleID::pi0);
    output[2]=getParticleData(ParticleID::pi0);
  }
  output[3]=getParticleData(ParticleID::pi0);
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

 
// the hadronic currents    
vector<LorentzPolarizationVectorE> 
FourPionNovosibirskCurrent::current(const int imode, const int ichan,
				    Energy & scale,const ParticleVector & decay,
				    DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<4;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  LorentzVector<complex<InvEnergy> > output;
  double fact(1.);
  // the momenta of the particles
  Lorentz5Momentum q1(decay[0]->momentum()),q2(decay[2]->momentum()),
    q3(decay[1]->momentum()),q4(decay[3]->momentum());
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  scale = Q.mass();
  // decide which decay mode
  // three charged pions
  if(imode==1) {
    // momenta of the particles
    LorentzVector<complex<Energy5> > veca1rho,vecomega,veca1sig;
    if(ichan<0) {
      // a_1 rho current
      veca1rho = 
	t1(q1,q2,q3,q4)+t1(q3,q2,q1,q4)+t1(q1,q3,q2,q4)
	+t1(q3,q1,q2,q4)+t1(q4,q3,q1,q2)+t1(q4,q1,q3,q2);
      // a_1 sigma current
      veca1sig = 
	t2(q4,q3,q1,q2,1)+t2(q4,q1,q3,q2,1)
	-t2(q1,q4,q3,q2,1)-t2(q3,q4,q1,q2,1);
      // omega current
      vecomega = 
	t3(q1,q2,q3,q4)+t3(q3,q2,q1,q4)-t3(q1,q3,q2,q4)
	-t3(q3,q1,q2,q4)-t3(q1,q4,q3,q2)-t3(q3,q4,q1,q2);
    }
    else if(ichan== 0) vecomega = t3(q1,q4,q3,q2);
    else if(ichan== 1) vecomega = t3(q3,q4,q1,q2);
    else if(ichan== 2) vecomega = t3(q1,q2,q3,q4);
    else if(ichan== 3) vecomega = t3(q3,q2,q1,q4);
    else if(ichan== 4) vecomega = t3(q1,q3,q2,q4);
    else if(ichan== 5) vecomega = t3(q3,q1,q2,q4);
    else if(ichan== 6) veca1rho = t1(q4,q1,q3,q2);
    else if(ichan== 7) veca1rho = t1(q4,q3,q1,q2);
    else if(ichan== 8) veca1rho = t1(q1,q2,q3,q4);
    else if(ichan== 9) veca1rho = t1(q3,q2,q1,q4);
    else if(ichan==10) veca1rho = t1(q1,q3,q2,q4);
    else if(ichan==11) veca1rho = t1(q3,q1,q2,q4);
    else if(ichan==12) veca1sig = t2(q4,q1,q3,q2,1);
    else if(ichan==13) veca1sig = t2(q4,q3,q1,q2,1);
    else if(ichan==14) veca1sig = t2(q1,q4,q3,q2,1);
    else if(ichan==15) veca1sig = t2(q3,q4,q1,q2,1);
    // final manipulations
    veca1rho += veca1sig;
    LorentzVector<complex<InvEnergy> > 
      veca1rho1 = veca1rho * gFunction(Q.mass2(),1);
    LorentzVector<complex<InvEnergy> > 
      vecomega1 = vecomega * gFunction(Q.mass2(),2);
    output = vecomega1 + veca1rho1;
    // this is 1/sqrt(2) for identical particles
    fact *= 1./sqrt(2.);
  }
  else if(imode==0) {
    // momenta of the particles
    LorentzVector<complex<Energy5> > veca1rho,veca1sig;
    if(ichan<0) {
      // a_1 rho current
      veca1rho= t1(q2,q3,q1,q4)+t1(q2,q4,q1,q3)+t1(q3,q2,q1,q4)
	+t1(q3,q4,q1,q2)+t1(q4,q2,q1,q3)+t1(q4,q3,q1,q2);
      // a_1 sigma current
      veca1sig=
 	t2(q2,q1,q3,q4,0)+t2(q3,q1,q2,q4,0)+t2(q4,q1,q3,q2,0)
	-t2(q1,q2,q3,q4,0)-t2(q1,q3,q2,q4,0)-t2(q1,q4,q3,q2,0);
    }
    else if(ichan== 0) veca1rho = t1(q2,q3,q1,q4);
    else if(ichan== 1) veca1rho = t1(q2,q4,q1,q3);
    else if(ichan== 2) veca1rho = t1(q3,q2,q1,q4);
    else if(ichan== 3) veca1rho = t1(q3,q4,q1,q2);
    else if(ichan== 4) veca1rho = t1(q4,q2,q1,q3);
    else if(ichan== 5) veca1rho = t1(q4,q3,q1,q2);
    else if(ichan== 6) veca1sig = t2(q2,q1,q3,q4,0);
    else if(ichan== 7) veca1sig = t2(q3,q1,q2,q4,0);
    else if(ichan== 8) veca1sig = t2(q4,q1,q3,q2,0);
    else if(ichan== 9) veca1sig = t2(q1,q2,q3,q4,0);
    else if(ichan==10) veca1sig = t2(q1,q3,q2,q4,0);
    else if(ichan==11) veca1sig = t2(q1,q4,q3,q2,0);
    // add them up 
    output = (veca1rho + veca1sig) * gFunction(Q.mass2(),0);
    // this is sqrt(1/3!) for identical particles
    fact *= 1./sqrt(6.);
  }     
  else {
    throw DecayIntegratorError() << "Unknown decay mode in the " 
				 << "FourPionNovosibirskCurrent::"
				 << "hadronCurrent()" << Exception::abortnow;
  }  
  vector<LorentzPolarizationVectorE> temp(1, output * fact * Q.mass2()); 
  return temp;
}

bool FourPionNovosibirskCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check four products
  if(id.size()!=4){return false;}
  int npiminus=0,npiplus=0,npi0=0;
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID:: piplus)      ++npiplus;
    else if(id[ix]==ParticleID::piminus) ++npiminus;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
  }
  if(npiminus==2&&npiplus==1&&npi0==1)      allowed=true;
  else if(npiminus==1&&npi0==3)             allowed=true;
  else if(npiplus==2&&npiminus==1&&npi0==1) allowed=true;
  else if(npiplus==1&&npi0==3)              allowed=true;
  return allowed;
}

// the decay mode
unsigned int FourPionNovosibirskCurrent::decayMode(vector<int> idout) {
  unsigned int npi(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::piplus) ++npi;
  }
  if(npi==3) return 1;
  return 0;
}


// output the information for the database
void FourPionNovosibirskCurrent::dataBaseOutput(ofstream & output,bool header,
						bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::FourPionNovosibirskCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  output << "newdef " << name() << ":rhoMass "    << _rhomass/GeV << "\n";
  output << "newdef " << name() << ":a1Mass  "    << _a1mass/GeV  << "\n";
  output << "newdef " << name() << ":sigmaMass  " << _sigmamass/GeV  << "\n";
  output << "newdef " << name() << ":omegaMass  " << _omegamass/GeV  << "\n";
  output << "newdef " << name() << ":rhoWidth "    << _rhowidth/GeV << "\n";
  output << "newdef " << name() << ":a1Width  "    << _a1width/GeV  << "\n";
  output << "newdef " << name() << ":sigmaWidth  " << _sigmawidth/GeV  << "\n";
  output << "newdef " << name() << ":omegaWidth  " << _omegawidth/GeV  << "\n";
  output << "newdef " << name() << ":IntegrationMass "  << _intmass/GeV  << "\n";
  output << "newdef " << name() << ":IntegrationWidth " << _intwidth/GeV  << "\n";
  output << "newdef " << name() << ":SigmaMagnitude "  <<  _zmag << "\n";
  output << "newdef " << name() << ":SigmaPhase " << _zphase  << "\n";
  output << "newdef " << name() << ":Lambda2 "  <<  _lambda2/GeV2 << "\n";
  output << "newdef " << name() << ":LocalParameters " <<  _localparameters << "\n";
  output << "newdef " << name() << ":Initializea1 " <<  _initializea1 << "\n";
  for(unsigned int ix=0;ix<_a1runwidth.size();++ix) {
    if(ix<200) output << "newdef ";
    else       output << "insert ";
    output << name() << ":a1RunningWidth " << ix << " " 
	   << _a1runwidth[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_a1runq2.size();++ix) {
    if(ix<200) output << "newdef ";
    else       output << "insert ";
    output << name() << ":a1RunningQ2 " << ix << " " << _a1runq2[ix]/GeV2 << "\n";
  }
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

double FourPionNovosibirskCurrent::
threeBodyMatrixElement(const int iopt, const Energy2 q2,
		       const Energy2 s3, const Energy2 s2, 
		       const Energy2 s1, const Energy,
		       const Energy, const Energy) const {
  unsigned int ix;
  // construct the momenta of the decay products
  Energy p1[5],p2[5],p3[5];
  Energy2 p1sq, p2sq, p3sq;
  Energy q(sqrt(q2));
  if(iopt==0) {
    p1[0] = 0.5*(q2+_mpi02-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-_mpi02);
    p2[0] = 0.5*(q2+_mpic2-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-_mpic2);
    p3[0] = 0.5*(q2+_mpic2-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-_mpic2);
  }
  else {
    p1[0] = 0.5*(q2+_mpi02-s1)/q; p1sq=p1[0]*p1[0]; p1[4]=sqrt(p1sq-_mpi02);
    p2[0] = 0.5*(q2+_mpi02-s2)/q; p2sq=p2[0]*p2[0]; p2[4]=sqrt(p2sq-_mpi02);
    p3[0] = 0.5*(q2+_mpi02-s3)/q; p3sq=p3[0]*p3[0]; p3[4]=sqrt(p3sq-_mpi02);
  }
  // take momentum of 1 parallel to z axis
  p1[1]=ZERO;p1[2]=ZERO;p1[3]=p1[4];
  // construct 2 
  double cos2(0.5*(sqr(p1[4])+sqr(p2[4])-sqr(p3[4]))/p1[4]/p2[4]);
  p2[1] = p2[4]*sqrt(1.-sqr(cos2)); p2[2]=ZERO; p2[3]=-p2[4]*cos2;
  // construct 3
  double cos3(0.5*(sqr(p1[4])-sqr(p2[4])+sqr(p3[4]))/p1[4]/p3[4]);
  p3[1] =-p3[4]*sqrt(1.-sqr(cos3)); p3[2]=ZERO; p3[3]=-p3[4]*cos3;
  // pi+pi-pi0 term
  complex<Energy4> output(0.*sqr(MeV2));
  if(iopt==0) {
    // values for the different Breit-Wigner terms
    Complex rho1(2.365*rhoBreitWigner(s2)),
      rho2(2.365*rhoBreitWigner(s3)),
      sig1(sigmaBreitWigner(s1,1));
    // compute the vector
    complex<Energy2> term;
    for(ix=1;ix<4;++ix) { 
      term = (p1[0]*p2[ix]-p2[0]*p1[ix])*rho2+(p1[0]*p3[ix]-p3[0]*p1[ix])*rho1
	+_zsigma*q*p1[ix]*sig1;
      output+=term*conj(term);
    }
  }
  // pi0pi0pi0 term
  else if(iopt==1) {
    // values for the different Breit-Wigner terms
    Complex sig1(sigmaBreitWigner(s1,0)),
      sig2(sigmaBreitWigner(s2,0)),
      sig3(sigmaBreitWigner(s3,0));
    // compute the vector
    complex<Energy2> term;
    for(ix=1;ix<4;++ix) {
      term = _zsigma * q * (p1[ix]*sig1 + p2[ix]*sig2 + p3[ix]*sig3);
      output += term*conj(term);
    }
    output/=6.;
  }
  output *= a1FormFactor(q2);
  return output.real() / pow<4,1>(_rhomass);
}

Complex FourPionNovosibirskCurrent::sigmaBreitWigner(Energy2 q2,
						     unsigned int iopt) const {
  Energy q(sqrt(q2));
  Energy pcm = iopt==0 ? 
    Kinematics::pstarTwoBodyDecay(q,_mpi0,_mpi0) :
    Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic);
  if(pcm<ZERO) pcm=ZERO;
  Energy  width(_sigmawidth*pcm/_psigma[iopt]);
  Energy2 msigma2 = sqr(_sigmamass);
  return msigma2/(q2-msigma2+Complex(0.,1.)*msigma2*width/q);
}

Complex FourPionNovosibirskCurrent::a1BreitWigner(Energy2 q2) const {
  Complex ii(0.,1.);
  Energy2 m2 = sqr(_a1mass);
  Energy q = sqrt(q2);
  return (m2/complex<Energy2>(q2 - m2 + ii*q*a1width(q2)));
}

Complex FourPionNovosibirskCurrent::omegaBreitWigner(Energy2 q2) const {
  Energy q(sqrt(q2));
  // calcluate the running width
  double diff((q-_omegamass)/GeV),temp(diff);
  double gomega(1.);
  Complex ii(0.,1.);
  if(q<=1.*GeV) {
    for(unsigned int ix=0;ix<6;++ix) {
      gomega +=temp*_omegaparam[ix];
      temp*=diff;
    }
  }
  else {
    gomega=_omegaparam[6]+q/GeV*(_omegaparam[7]+q/GeV*_omegaparam[8]
				 +q2/GeV2*_omegaparam[9]);
  }
  if(gomega<0.){gomega=0.;}
  Energy2 numer=_omegamass*_omegamass;
  complex<Energy2> denom=q2-_omegamass*_omegamass+ii*_omegamass*_omegawidth*gomega;
  return numer/denom;
}

Complex FourPionNovosibirskCurrent::rhoBreitWigner(Energy2 q2) const {
  Energy q(sqrt(q2));
  Energy2 grhom(8.*_prho*_prho*_prho/_rhomass);
  complex<Energy2> denom;
  Complex ii(0.,1.);
  if(q2<4.*_mpic2) {
    denom=q2-_rhomass*_rhomass
      -_rhowidth*_rhomass*(hFunction(q)-_hm2-(q2-_rhomass*_rhomass)*_dhdq2m2)/grhom;
  }
  else {
    Energy pcm(2.*Kinematics::pstarTwoBodyDecay(q,_mpic,_mpic));
    Energy2 grho(pcm*pcm*pcm/q);
    denom=q2-_rhomass*_rhomass
      -_rhowidth*_rhomass*(hFunction(q)-_hm2-(q2-_rhomass*_rhomass)*_dhdq2m2)/grhom
      +ii*_rhomass*_rhowidth*grho/grhom;
  }
  return _rhoD/denom;
}

LorentzVector<complex<Energy5> > 
FourPionNovosibirskCurrent::t1(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
			       Lorentz5Momentum & q3,Lorentz5Momentum & q4) const {
  // momentum of the whole sysytem
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  // compute the virtuality of the a_1
  Lorentz5Momentum a1(q2+q3+q4);a1.rescaleMass();
  // compute the virtuality of the  rho
  Lorentz5Momentum rho(q3+q4);rho.rescaleMass();
  // compute the prefactor    
  Complex pre(-a1FormFactor(a1.mass2())*a1BreitWigner(a1.mass2())*
	      rhoBreitWigner(rho.mass2()));
  // dot products we need
  Energy2 QdQmq1(Q*a1);
  complex<Energy4> consta(QdQmq1*(a1*q3)), constb(QdQmq1*(a1*q4)),
    constc(((Q*q4)*(q1*q3)-(Q*q3)*(q1*q4)));
  // compute the current
  return pre*(consta*q4-constb*q3+constc*a1);
}

LorentzVector<complex<Energy5> > 
FourPionNovosibirskCurrent::t2(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
			       Lorentz5Momentum & q3,Lorentz5Momentum & q4,
			       unsigned int iopt) const {
  // momentum of the whole system
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  // compute the virtuality of the a_1
  Lorentz5Momentum a1(q2+q3+q4);a1.rescaleMass();
  // compute the virtuality of the  sigma
  Lorentz5Momentum sigma(q3+q4);sigma.rescaleMass();
  // compute the prefactor
  Complex pre(_zsigma*a1FormFactor(a1.mass2())
	      *a1BreitWigner(a1.mass2())*
	      sigmaBreitWigner(sigma.mass2(),iopt));
  // dot products we need
  complex<Energy4> consta((Q*a1)*a1.mass2()),constb((Q*q2)*a1.mass2());
  // compute the current
  return pre*(consta*q2-constb*a1);
}

LorentzVector<complex<Energy5> > 
FourPionNovosibirskCurrent::t3(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
			       Lorentz5Momentum & q3,Lorentz5Momentum & q4) const {
  // momentum of the whole sysytem
  Lorentz5Momentum Q(q1+q2+q3+q4);Q.rescaleMass();
  // compute the virtuality of the omega
  Lorentz5Momentum omega(q2+q3+q4);omega.rescaleMass();
  // compute the virtuality of the  rho
  Lorentz5Momentum rho(q3+q4);rho.rescaleMass();
  // compute the prefactor
  Complex pre(omegaBreitWigner(omega.mass2())*rhoBreitWigner(rho.mass2()));
  // dot products we need
  complex<Energy4> consta((Q*q3)*(q1*q4)-(Q*q4)*(q1*q3)),
    constb(-(Q*q2)*(q1*q4)+(q1*q2)*(Q*q4)),
    constc((Q*q2)*(q1*q3)-(q1*q2)*(Q*q3));
  // compute the current
  return pre*(consta*q2+constb*q3+constc*q4);
}

InvEnergy6 FourPionNovosibirskCurrent::gFunction(Energy2 q2, int ichan) const {
  Energy q(sqrt(q2));
  InvEnergy4 invmrho4 = 1/sqr(sqr(_rhomass));
  // the one charged pion G function
  if(ichan==0) {
    return (*_Fonec)(q) * _aonec * (*_Fsigma)(q2) * sqrt(_bonec*q/GeV-_conec) *
      invmrho4/q;
  }
  // the three charged pion G function
  else if(ichan==1) {
    return (*_Fthreec)(q)*_athreec*sqrt(_bthreec*q/GeV-_cthreec)*invmrho4/q; 
  }
  // the omega G function
  else if(ichan==2) {
    return(*_Fomega)(q)*_aomega*sqrt(_bomega*q/GeV-_comega)*invmrho4/q;
  }
  assert(false);
  return InvEnergy6();
}

Energy2 FourPionNovosibirskCurrent::DParameter() const {
  Energy2 grhom(8.*_prho*_prho*_prho/_rhomass);
  return _rhomass*_rhomass+_rhowidth*_rhomass*
    (hFunction(ZERO)-_hm2+_rhomass*_rhomass*_dhdq2m2)/grhom;
}

double FourPionNovosibirskCurrent::dhdq2Parameter() const {
  Energy2 mrho2(_rhomass*_rhomass);
  double root(sqrt(1.-4.*_mpic2/mrho2));
  return root/Constants::pi*(root+(1.+2*_mpic2/mrho2)*log((1+root)/(1-root)));
}

Energy2 FourPionNovosibirskCurrent::hFunction(const Energy q) const {
  using Constants::pi;
  static const Energy2 eps(0.01*MeV2);

  Energy2 q2(q*q), output;
  if (q2 > 4*_mpic2) {
    double root = sqrt(1.-4.*_mpic2/q2);
    output = root*log((1.+root)/(1.-root))*(q2-4*_mpic2)/pi;
  }
  else if (q2 > eps) output = ZERO;
  else               output = -8.*_mpic2/pi;
  return output;
}

