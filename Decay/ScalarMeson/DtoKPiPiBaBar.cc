// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiBaBar class.
//

#include "DtoKPiPiBaBar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

DtoKPiPiBaBar::DtoKPiPiBaBar() {
  // choice of the model
  _imodel=0;
  // Parameters for the K-matrix based fit
  // g_pipi
  _gpipi.push_back( 0.229*GeV);
  _gpipi.push_back( 0.941*GeV);
  _gpipi.push_back( 0.369*GeV);
  _gpipi.push_back( 0.337*GeV);
  _gpipi.push_back( 0.182*GeV);
  // g_KK
  _gKK.push_back(-0.554*GeV);
  _gKK.push_back( 0.551*GeV);
  _gKK.push_back( 0.239*GeV);
  _gKK.push_back( 0.409*GeV);
  _gKK.push_back(-0.176*GeV);
  // g_4pi
  _g4pi.push_back( 0.000*GeV);
  _g4pi.push_back( 0.000*GeV);
  _g4pi.push_back( 0.556*GeV);
  _g4pi.push_back( 0.857*GeV);
  _g4pi.push_back(-0.797*GeV);
  // g_etaeta
  _getaeta.push_back(-0.399*GeV);
  _getaeta.push_back( 0.391*GeV);
  _getaeta.push_back( 0.183*GeV);
  _getaeta.push_back( 0.199*GeV);
  _getaeta.push_back(-0.004*GeV);
  // g_etaetap
  _getaetap.push_back(-0.346*GeV);
  _getaetap.push_back( 0.315*GeV);
  _getaetap.push_back( 0.187*GeV);
  _getaetap.push_back(-0.010*GeV);
  _getaetap.push_back( 0.224*GeV);
  // m_alpha
  _malpha.push_back(0.651*GeV);
  _malpha.push_back(1.204*GeV);
  _malpha.push_back(1.558*GeV);
  _malpha.push_back(1.210*GeV);
  _malpha.push_back(1.822*GeV);
  // f^scatt
  _fscatt.push_back( 0.234);
  _fscatt.push_back( 0.150);
  _fscatt.push_back(-0.206);
  _fscatt.push_back( 0.328);
  _fscatt.push_back( 0.354);
  // s^scatt_0
  _s0scatt=-3.926*GeV2;
  // s_a0
  _sA0 = -0.15*GeV2;
  // s_A
  _sA  = 1.*GeV2;
  // the beta production coupling
  _betare.push_back( -3.78*GeV);_betaim.push_back(1.23*GeV);
  _betare.push_back(  9.55*GeV);_betaim.push_back(3.43*GeV);
  _betare.push_back(  0.00*GeV);_betaim.push_back(0.00*GeV);
  _betare.push_back( 12.97*GeV);_betaim.push_back(1.27*GeV);
  _betare.push_back(  0.00*GeV);_betaim.push_back(0.00*GeV);
  // f^prod
  _fprodre = -10.22; _fprodim = -6.35;
  // zero the masses
  _mpi   = 0.;
  _mK    = 0.;
  _meta  = 0.;
  _metap = 0.;
  // rho_3 interpolation tables
  _initrho3 = false;
  double rscale[101]={0.311677,0.31856 ,0.325443,0.332326,0.339209,
		      0.346093,0.352976,0.359859,0.366742,0.373626,
		      0.380509,0.387392,0.394275,0.401159,0.408042,
		      0.414925,0.421808,0.428692,0.435575,0.442458,
		      0.449341,0.456224,0.463108,0.469991,0.476874,
		      0.483757,0.490641,0.497524,0.504407,0.51129 ,
		      0.518174,0.525057,0.53194 ,0.538823,0.545707,
		      0.55259 ,0.559473,0.566356,0.573239,0.580123,
		      0.587006,0.593889,0.600772,0.607656,0.614539,
		      0.621422,0.628305,0.635189,0.642072,0.648955,
		      0.655838,0.662722,0.669605,0.676488,0.683371,
		      0.690254,0.697138,0.704021,0.710904,0.717787,
		      0.724671,0.731554,0.738437,0.74532 ,0.752204,
		      0.759087,0.76597 ,0.772853,0.779736,0.78662 ,
		      0.793503,0.800386,0.807269,0.814153,0.821036,
		      0.827919,0.834802,0.841686,0.848569,0.855452,
		      0.862335,0.869219,0.876102,0.882985,0.889868,
		      0.896751,0.903635,0.910518,0.917401,0.924284,
		      0.931168,0.938051,0.944934,0.951817,0.958701,
		      0.965584,0.972467,0.97935 ,0.986234,0.993117,
		      1};
  double rvalue[101]={0          ,8.53899e-12,3.80878e-10,3.49296e-09,1.67667e-08,
		      5.64581e-08,1.51939e-07,3.50344e-07,7.21493e-07,1.36305e-06,
		      2.40589e-06,4.0196e-06 ,6.41823e-06,9.8661e-06 ,1.46838e-05,
		      2.12546e-05,3.00302e-05,4.15379e-05,5.6387e-05 ,7.52756e-05,
		      9.89975e-05,0.00012845 ,0.000164641,0.000208697,0.00026187 ,
		      0.000325548,0.00040126 ,0.00049069 ,0.00059568 ,0.000718245,
		      0.000860577,0.00102506 ,0.00121429 ,0.00143104 ,0.00167836 ,
		      0.00195948 ,0.00227793 ,0.00263745 ,0.00304208 ,0.00349616 ,
		      0.00400431 ,0.00457147 ,0.00520293 ,0.0059043  ,0.00668158 ,
		      0.00754116 ,0.00848981 ,0.00953475 ,0.0106836  ,0.0119445  ,
		      0.0133261  ,0.0148374  ,0.0164882  ,0.0182885  ,0.0202493  ,
		      0.022382   ,0.0246986  ,0.0272118  ,0.0299351  ,0.0328828  ,
		      0.0360698  ,0.0395118  ,0.0432257  ,0.0472288  ,0.0515398  ,
		      0.0561781  ,0.0611643  ,0.0665199  ,0.0722678  ,0.0784319  ,
		      0.0850375  ,0.092111   ,0.0996805  ,0.107775   ,0.116426   ,
		      0.125666   ,0.135529   ,0.146051   ,0.15727    ,0.169226   ,
		      0.18196    ,0.195517   ,0.209943   ,0.225286   ,0.241599   ,
		      0.258934   ,0.277349   ,0.296904   ,0.31766    ,0.339684   ,
		      0.363046   ,0.387818   ,0.414078   ,0.441905   ,0.471385   ,
		      0.502608   ,0.535667   ,0.570662   ,0.607696   ,0.646878   ,
		      0.688323};
  _rho3scale=vector<double>(rscale,rscale+101);
  _rho3     =vector<double>(rvalue,rvalue+101);
  // The amplitudes and phases for the K-matrix fit
  _k892mre   =-1.159     ; _k892mim   = 1.361     ;
  _k1430mre0 = 2.482*GeV2; _k1430mim0 =-0.653*GeV2;
  _k1430mre2 = 0.852/GeV2; _k1430mim2 =-0.729/GeV2;
  _k1410mre  =-0.402     ; _k1410mim  = 0.050     ;
  _k1680mre  =-1.000     ; _k1680mim  = 1.690     ;
  _k892pre   = 0.133     ; _k892pim   =-0.132     ;
  _k1430pre0 = 0.375*GeV2; _k1430pim0 =-0.143*GeV2;
  _k1430pre2 = 0.088/GeV2; _k1430pim2 =-0.057/GeV2;
  _rho770re  = 1.000     ; _rho770im  = 0.000     ;
  _omegare   =-0.0182    ; _omegaim   =-0.0367    ;
  _f2re      = 0.787/GeV2; _f2im      =-0.397/GeV2;
  _rho1450re = 0.405     ; _rho1450im =-0.458     ;
  // amplitudes and phase for the normal fit
  _k892mamp     = 1.781     ; _k892mphase   =  131.0;
  _k1430mamp0   = 2.45*GeV2 ; _k1430mphase0 = -  8.3;
  _k1430mamp2   = 1.05/GeV2 ; _k1430mphase2 = - 54.3;
  _k1410mamp    = 0.52      ; _k1410mphase  =  154. ;
  _k1680mamp    = 0.89      ; _k1680mphase  = -139. ;
  _k892pamp     = 0.180     ; _k892pphase   =  -44.1;
  _k1430pamp0   = 0.37*GeV2 ; _k1430pphase0 =   18. ;
  _k1430pamp2   = 0.075/GeV2; _k1430pphase2 = -104. ;
  _rho770amp    = 1.0       ; _rho770phase  =    0. ;
  _omegaamp     = 0.0391    ; _omegaphase   =  115.3;
  _f2amp        = 0.922/GeV2; _f2phase      = - 21.3;
  _rho1450amp   = 0.52      ; _rho1450phase =   38. ;
  _f980amp      = 0.482*GeV2; _f980phase    = -141.8;
  _f1370amp     = 2.25*GeV2 ; _f1370phase   =  113.2;
  _sigmaamp     = 1.36*GeV2 ; _sigmaphase   = -177.9;
  _sigmapamp    = 0.34*GeV2 ; _sigmapphase  =  153. ;
  _nonamp       = 3.53      ; _nonphase     =  128. ;
  // radial sizes
  _rD0  = 5.0/GeV;
  _rres = 1.5/GeV;
  // masses and widths of the resonaces
  _mK892    =  891.66*MeV; _wK892    =  50.80*MeV;      
  _mK14300  = 1412.00*MeV; _wK14300  = 294.00*MeV;    
  _mK14302  = 1425.60*MeV; _wK14302  =  98.50*MeV;    
  _mK1410   = 1414.00*MeV; _wK1410   = 232.00*MeV;     
  _mK1680   = 1717.00*MeV; _wK1680   = 322.00*MeV;     
  _mrho770  =  775.80*MeV; _wrho770  = 150.30*MeV;    
  _momega   =  782.59*MeV; _womega   =   8.49*MeV;     
  _mf2      = 1275.40*MeV; _wf2      = 185.10*MeV;	      
  _mrho1450 = 1465.00*MeV; _wrho1450 = 400.00*MeV;   
  _mf980    =  977.00*MeV; _wf980    =  44.00*MeV;      
  _mf1370   = 1434.00*MeV; _wf1370   = 173.00*MeV;     
  _msigma   =  484.00*MeV; _wsigma   = 383.00*MeV;     
  _msigmap  = 1014.00*MeV; _wsigmap  =  88.00*MeV;
  // zero all the amplitudes
  _aKm892   = 0.;
  _aKm14300 = 0.;
  _aKm14302 = 0.;
  _aKm1410  = 0.;
  _aKm1680  = 0.;
  _aKp892   = 0.;
  _aKp14300 = 0.;
  _aKp14302 = 0.;
  _arho770  = 0.;
  _aomega   = 0.;
  _af2      = 0.;
  _arho1450 = 0.;
  _af980    = 0.;
  _af1370   = 0.;
  _asigma   = 0.;
  _asigmap  = 0.;
  _aNR      = 0.;
  // generate the intermediates
  generateIntermediates(true);
}

void DtoKPiPiBaBar::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // calculate the amplitudes for the different terms
  if(_imodel==0) {
    double fact = pi/180.;
    _aKm892   = _k892mamp  *Complex(cos(_k892mphase  *fact),sin(_k892mphase  *fact));
    _aKm14300 = _k1430mamp0*Complex(cos(_k1430mphase0*fact),sin(_k1430mphase0*fact));
    _aKm14302 = _k1430mamp2*Complex(cos(_k1430mphase2*fact),sin(_k1430mphase2*fact));
    _aKm1410  = _k1410mamp *Complex(cos(_k1410mphase *fact),sin(_k1410mphase *fact));
    _aKm1680  = _k1680mamp *Complex(cos(_k1680mphase *fact),sin(_k1680mphase *fact));
    _aKp892   = _k892pamp  *Complex(cos(_k892pphase  *fact),sin(_k892pphase  *fact));
    _aKp14300 = _k1430pamp0*Complex(cos(_k1430pphase0*fact),sin(_k1430pphase0*fact));
    _aKp14302 = _k1430pamp2*Complex(cos(_k1430pphase2*fact),sin(_k1430pphase2*fact));
    _arho770  = _rho770amp *Complex(cos(_rho770phase *fact),sin(_rho770phase *fact));
    _aomega   = _omegaamp  *Complex(cos(_omegaphase  *fact),sin(_omegaphase  *fact));
    _af2      = _f2amp     *Complex(cos(_f2phase     *fact),sin(_f2phase     *fact));
    _arho1450 = _rho1450amp*Complex(cos(_rho1450phase*fact),sin(_rho1450phase*fact));
    _af980    = _f980amp   *Complex(cos(_f980phase   *fact),sin(_f980phase   *fact));
    _af1370   = _f1370amp  *Complex(cos(_f1370phase  *fact),sin(_f1370phase  *fact));
    _asigma   = _sigmaamp  *Complex(cos(_sigmaphase  *fact),sin(_sigmaphase  *fact));
    _asigmap  = _sigmapamp *Complex(cos(_sigmapphase *fact),sin(_sigmapphase *fact));
    _aNR      = _nonamp    *Complex(cos(_nonphase    *fact),sin(_nonphase    *fact));
  }
  else {
    _aKm892   = Complex(_k892mre  ,_k892mim  );
    _aKm14300 = Complex(_k1430mre0,_k1430mim0);
    _aKm14302 = Complex(_k1430mre2,_k1430mim2);
    _aKm1410  = Complex(_k1410mre ,_k1410mim );
    _aKm1680  = Complex(_k1680mre ,_k1680mim );
    _aKp892   = Complex(_k892pre  ,_k892pim  );
    _aKp14300 = Complex(_k1430pre0,_k1430pim0);
    _aKp14302 = Complex(_k1430pre2,_k1430pim2);
    _arho770  = Complex(_rho770re ,_rho770im );
    _aomega   = Complex(_omegare  ,_omegaim  );
    _af2      = Complex(_f2re     ,_f2re     );
    _arho1450 = Complex(_rho1450re,_rho1450im);
    _af980    = 0.;
    _af1370   = 0.;
    _asigma   = 0.;
    _asigmap  = 0.;
    _aNR      = 0.;
  }
  if(_malpha  .size()!=5||_gKK     .size()!=5||_betaim.size()!=5||
     _g4pi    .size()!=5||_getaeta .size()!=5||_betare.size()!=5||
     _getaetap.size()!=5||_gpipi   .size()!=5) 
    throw InitException() << "Inconsitent sizes for vectors in DtoKPiPiBaBar::doinit()"
			  << Exception::runerror; 
  // K-matrix parameters
  _gialpha.push_back(_gpipi   );
  _gialpha.push_back(_gKK     );
  _gialpha.push_back(_g4pi    );
  _gialpha.push_back(_getaeta );
  _gialpha.push_back(_getaetap);
  _fprod = Complex(_fprodre,_fprodim);
  _beta.resize(5);
  for(unsigned int ix=0;ix<5;++ix) {
    _beta[ix]=Complex(_betare[ix],_betaim[ix]);
  }
  // masses for the K-matrix
  _mpi   = getParticleData(ParticleID::piplus  )->mass();
  _mK    = getParticleData(ParticleID::Kplus   )->mass();
  _meta  = getParticleData(ParticleID::eta     )->mass();
  _metap = getParticleData(ParticleID::etaprime)->mass();
  // particle data objects for the intermediates
  tPDPtr k892m   = getParticleData(ParticleID::Kstarminus);
  tPDPtr k892p   = getParticleData(ParticleID::Kstarplus );
  tPDPtr k1430m0 = getParticleData(ParticleID::Kstar_0minus);
  tPDPtr k1430p0 = getParticleData(ParticleID::Kstar_0plus);
  tPDPtr k1430m2 = getParticleData(ParticleID::Kstar_2plus);
  tPDPtr k1430p2 = getParticleData(ParticleID::Kstar_2minus);
  tPDPtr k1410m  = getParticleData(-100323);
  tPDPtr k1680m  = getParticleData(-30323);
  tPDPtr rho770  = getParticleData(ParticleID::rho0);
  tPDPtr omega   = getParticleData(ParticleID::omega);
  tPDPtr f980    = getParticleData(9010221);
  tPDPtr f1370   = getParticleData(10221);
  tPDPtr f2      = getParticleData(ParticleID::f_2);
  tPDPtr rho1450 = getParticleData(100113);
  tPDPtr sigma   = getParticleData(9000221);
  // construct the channels for the decay
  DecayPhaseSpaceChannelPtr newchannel;
  PDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kbar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  unsigned int ix=0;
  // D0 -> K*- pi+
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K(1430)*_0- pi+
  if(k1430m0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430m0,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K(1430)*_2- pi+
  if(k1430m2) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430m2,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K(1410)*- pi+
  if(k1410m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1410m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K(1680)*- pi+
  if(k1680m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1680m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K*+ pi-
  if(k892p) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k892p,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K(1430)*_0+ pi-
  if(k1430p0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1430p0,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> K(1430)*_2+ pi-
  if(k1430p2) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1430p2,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> rho Kbar0
  if(rho770) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho770,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> omega Kbar0
  if(omega) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(omega,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> f0(980) Kbar0
  if(f980) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f980,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> f0(1370) Kbar0
  if(f1370) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f1370,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> f2(1270) Kbar0
  if(f2) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f2,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> rho Kbar0
  if(rho1450) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho1450,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // D0 -> sigma Kbar0
  if(sigma) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(sigma,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(_weights.empty()) {
    _weights.resize(ix,1./double(ix));
    _maxwgt=1.;
  }
  addMode(mode,_maxwgt,_weights);
  // calculate the rho 3 function if needed
  if(_initrho3) {
    _rho3.resize(0);
    _rho3scale.resize(0);
    DtoKPiPiBaBarOuterIntegrand integrand(this,_mpi);
    GaussianIntegrator integrator;
    Energy2 step = (GeV2-16.*sqr(_mpi))/100.;
    integrand.s(GeV2);
    double low = sqr(2.*_mpi/GeV);
    double upp = sqr(GeV-2.*_mpi)/GeV2;
    double norm = (GeV2-16.*sqr(_mpi))/GeV2/integrator.value(integrand,low,upp);
    for(Energy2 s=16.*sqr(_mpi);s<GeV2;s+=step) {
      integrand.s(s);
      low = sqr(2.*_mpi/GeV);
      upp = sqr(sqrt(s)-2.*_mpi)/GeV2;
      integrand.s(s);
      double output = norm*integrator.value(integrand,low,upp);
      _rho3     .push_back(output);
      _rho3scale.push_back(s/GeV2);
    }
  }
  _rho3inter = new_ptr(Interpolator(_rho3,_rho3scale,3));


  for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
    Complex test=F1(s);
    cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
  }
  cerr << "join\n";
   Energy mD=getParticleData(ParticleID::D0    )->mass();
   Energy mA=getParticleData(ParticleID::K0    )->mass();
   Energy mB=getParticleData(ParticleID::piplus)->mass();
  for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
    Complex test=
      _af980  /GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_mf980  ,_wf980  )+
      _af1370 /GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_mf1370 ,_wf1370 )+
      _asigma /GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_msigma ,_wsigma )+
      _asigmap/GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_msigmap,_wsigmap)+
      _aNR;
    cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
  }
  cerr << "join red\n";
//   for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
//     Complex test=
//       _af980  /GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_mf980  ,_wf980  );
//     cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
//   }
//   cerr << "join blue\n";
//   for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
//     Complex test=
//       _af1370 /GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_mf1370 ,_wf1370 );
//     cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
//   }
//   cerr << "join green\n";
//    for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
//      Complex test=
//        _asigma /GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_msigma ,_wsigma );
//      cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
//    }
//    cerr << "join magenta\n";
//   for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
//     Complex test=
//       _asigmap/GeV2*amplitude(0,false,mD,mB,mB,mA,sqrt(s),0.,0.,_msigmap,_wsigmap);
//     cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
//   }
//   cerr << "join cyan\n";
//   for(Energy2 s=4.*sqr(_mpi);s<2.*GeV2;s+=0.01*GeV2) {
//     Complex test=_aNR;
//     cerr << s/GeV2 << " " << real(test*conj(test)) << "\n";
//   }
//   cerr << "join yellow\n";
}

void DtoKPiPiBaBar::persistentOutput(PersistentOStream & os) const {
  os << _imodel << _gpipi << _gKK << _g4pi << _getaeta << _getaetap
     << _malpha << _fscatt << _s0scatt << _sA0 << _sA << _rho0 
     << _betare << _betaim << _beta << _fprodre << _fprodim << _fprod
     << _k892mre << _k892mim << _k1430mre0 << _k1430mim0 << _k1430mre2
     << _k1430mim2 << _k1410mre << _k1410mim << _k1680mre << _k1680mim
     << _k892pre << _k892pim << _k1430pre0 << _k1430pim0 << _k1430pre2
     << _k1430pim2 << _rho770re << _rho770im << _omegare << _omegaim
     << _f2re << _f2im << _rho1450re << _rho1450im << _k892mamp 
     << _k892mphase << _k1430mamp0 << _k1430mphase0 << _k1430mamp2 
     << _k1430mphase2 << _k1410mamp << _k1410mphase << _k1680mamp 
     << _k1680mphase << _k892pamp << _k892pphase << _k1430pamp0 
     << _k1430pphase0 << _k1430pamp2 << _k1430pphase2 << _rho770amp 
     << _rho770phase << _omegaamp << _omegaphase << _f2amp << _f2phase 
     << _rho1450amp << _rho1450phase << _f980amp << _f980phase << _f1370amp 
     << _f1370phase << _sigmaamp << _sigmaphase << _sigmapamp << _sigmapphase 
     << _nonamp  << _nonphase << _rD0 << _rres << _aKm892 << _aKm14300 
     << _aKm14302 << _aKm1410 << _aKm1680 << _aKp892 <<  _aKp14300 << _aKp14302 
     << _arho770 << _aomega << _af2 << _arho1450 << _af980 << _af1370 
     << _asigma << _asigmap << _aNR <<  _mK892 << _mK14300 << _mK14302 
     << _mK1410 << _mK1680 << _mrho770 << _momega << _mf2 << _mrho1450 
     << _mf980 << _mf1370 << _msigma << _msigmap << _wK892 << _wK14300 
     << _wK14302 << _wK1410 << _wK1680 << _wrho770 << _womega << _wf2 
     << _wrho1450 << _wf980 << _wf1370 << _wsigma << _wsigmap << _maxwgt 
     << _weights << _gialpha << _mpi << _mK << _meta << _metap << _rho3scale 
     << _rho3 << _initrho3 << _rho3inter;
}

void DtoKPiPiBaBar::persistentInput(PersistentIStream & is, int) {
  is >> _imodel >> _gpipi >> _gKK >> _g4pi >> _getaeta >> _getaetap
     >> _malpha >> _fscatt >> _s0scatt >> _sA0 >> _sA >> _rho0 
     >> _betare >> _betaim >> _beta >> _fprodre >> _fprodim >> _fprod
     >> _k892mre >> _k892mim >> _k1430mre0 >> _k1430mim0 >> _k1430mre2
     >> _k1430mim2 >> _k1410mre >> _k1410mim >> _k1680mre >> _k1680mim
     >> _k892pre >> _k892pim >> _k1430pre0 >> _k1430pim0 >> _k1430pre2
     >> _k1430pim2 >> _rho770re >> _rho770im >> _omegare >> _omegaim
     >> _f2re >> _f2im >> _rho1450re >> _rho1450im >> _k892mamp 
     >> _k892mphase >> _k1430mamp0 >> _k1430mphase0 >> _k1430mamp2 
     >> _k1430mphase2 >> _k1410mamp >> _k1410mphase >> _k1680mamp 
     >> _k1680mphase >> _k892pamp >> _k892pphase >> _k1430pamp0 
     >> _k1430pphase0 >> _k1430pamp2 >> _k1430pphase2 >> _rho770amp 
     >> _rho770phase >> _omegaamp >> _omegaphase >> _f2amp >> _f2phase 
     >> _rho1450amp >> _rho1450phase >> _f980amp >> _f980phase >> _f1370amp 
     >> _f1370phase >> _sigmaamp >> _sigmaphase >> _sigmapamp >> _sigmapphase 
     >> _nonamp  >> _nonphase >> _rD0 >> _rres >> _aKm892 >> _aKm14300 
     >> _aKm14302 >> _aKm1410 >> _aKm1680 >> _aKp892 >>  _aKp14300 >> _aKp14302 
     >> _arho770 >> _aomega >> _af2 >> _arho1450 >> _af980 >> _af1370 
     >> _asigma >> _asigmap >> _aNR >>  _mK892 >> _mK14300 >> _mK14302 
     >> _mK1410 >> _mK1680 >> _mrho770 >> _momega >> _mf2 >> _mrho1450 
     >> _mf980 >> _mf1370 >> _msigma >> _msigmap >> _wK892 >> _wK14300 
     >> _wK14302 >> _wK1410 >> _wK1680 >> _wrho770 >> _womega >> _wf2 
     >> _wrho1450 >> _wf980 >> _wf1370 >> _wsigma >> _wsigmap >> _maxwgt 
     >> _weights >> _gialpha >> _mpi >> _mK >> _meta >> _metap >> _rho3scale 
     >> _rho3 >> _initrho3 >> _rho3inter;
}

ClassDescription<DtoKPiPiBaBar> DtoKPiPiBaBar::initDtoKPiPiBaBar;
// Definition of the static class description member.

void DtoKPiPiBaBar::Init() {

  static ClassDocumentation<DtoKPiPiBaBar> documentation
    ("There is no documentation for the DtoKPiPiBaBar class");

  static Switch<DtoKPiPiBaBar,unsigned int> interfaceModel
    ("Model",
     "Which model of the scalar pipi resonances to use",
     &DtoKPiPiBaBar::_imodel, 1, false, false);
  static SwitchOption interfaceModelBreitWigner
    (interfaceModel,
     "BreitWigner",
     "Use the model with Breit-Wigners",
     0);
  static SwitchOption interfaceModelKMatrix
    (interfaceModel,
     "KMatrix",
     "Use the K-matrix model",
     1);

  static ParVector<DtoKPiPiBaBar,double> interfacegPiPi
    ("gPiPi",
     "The K-matrix coupling for the pi-pi channel",
     &DtoKPiPiBaBar::_gpipi, -1, 0., -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiBaBar,double> interfacegKK
    ("gKK",
     "The K-matrix coupling for the KK channel",
     &DtoKPiPiBaBar::_gKK, -1, 0.0, -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiBaBar,double> interfaceg4pi
    ("g4pi",
     "The K-matrix coupling for the 4pi channel",
     &DtoKPiPiBaBar::_g4pi, -1, 0.0, -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiBaBar,double> interfacegetaeta
    ("getaeta",
     "The K-matrix coupling for the etaeta channel",
     &DtoKPiPiBaBar::_getaeta, -1, 0.0, -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiBaBar,double> interfacegetaetap
    ("getaetap",
     "The K-matrix coupling for the eta etaprime channel",
     &DtoKPiPiBaBar::_getaeta, -1, 0.0, -100., 100.,
     false, false, Interface::limited);

  static ParVector<DtoKPiPiBaBar,Energy> interfaceKMatrixMasses
    ("KMatrixMasses",
     "The masses of the resonances in the K-matrix fix",
     &DtoKPiPiBaBar::_malpha, MeV, -1, 1.0*GeV, 0.0*MeV, 10.0*MeV,
     false, false, Interface::limited);

//   /**
//    *  The \f$f^{\rm scatt}_{1\alpha}\f$ parameters for the K-matrix
//    */
//   vector<double> _fscatt;

//   /**
//    *  \f$s^{\rm scatt}_0\f$
//    */
//   Energy2 _s0scatt;

//   /**
//    *  \f$s_{A_0}\f$
//    */
//   Energy2 _sA0;

//   /**
//    * \f$s_A\f$
//    */
//   Energy2 _sA;

//   /**
//    *  The \f$\rho_0\f$ parameter for the multi-body function
//    */
//   Energy2 _rho0;
//   //@}

//   /**
//    *  The \f$\beta\f$ production coupling
//    */
//   //@{
//   /**
//    *  The real part
//    */
//   vector<double> _betare;
  
//   /**
//    *  The imaginary part
//    */
//   vector<double> _betaim;
  
//   /**
//    *  The full coupling
//    */
//   vector<Complex> _beta;
//   //@}

//   /**
//    *  The \f$f^{\rm prod}\f$ parameters
//    */
//   //@{
//   /**
//    *  The real part
//    */
//   vector<double> _fprodre;

//   /**
//    *  The imaginary part
//    */
//   vector<double> _fprodim;

//   /**
//    *  The full coupling
//    */
//   vector<double> _fprod;
//   //@}

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar892minusReal
    ("KStar892minusReal",
     "The real part of the amplitude for the K*(892)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k892mre, -1.159, -10., 10.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar892minusImag
    ("KStar892minusImag",
     "The imaginary part of the amplitude for the K*(892)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k892mim, 1.361, -10., 10.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy2> interfaceK0Star1430minusReal
    ("K0Star1430minusReal",
     "The real part of the amplitude for the K*0(1430)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430mre0, GeV2, 2.482*GeV2, -10.*GeV2, 10.*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy2> interfaceK0Star1430minusImag
    ("K0Star1430minusImag",
     "The imaginary part of the amplitude for the K*0(1430)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430mim0, GeV2, -0.653*GeV2, -10.*GeV2, 10.*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,InvEnergy2> interfaceK2Star1430minusReal
    ("K2Star1430minusReal",
     "The real part of the amplitude for the K*2(1430)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430mre2, 1./GeV2, 0.852/GeV2, -10./GeV2, 10./GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,InvEnergy2> interfaceK2Star1430minusImag
    ("K2Star1430minusImag",
     "The imaginary part of the amplitude for the K2*(1430)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430mim2, 1./GeV2, -0.729/GeV2, -10./GeV2, 10./GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar1410minusReal
    ("KStar1410minusReal",
     "The real part of the amplitude for the K*(1410)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1410mre, -0.402, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar1410minusImag
    ("KStar1410minusImag",
     "The imaginary part of the amplitude for the K*(1410)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1410mim, 0.050, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar1680minusReal
    ("KStar1680minusReal",
     "The real part of the amplitude for the K*(1680)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1680mre, -1.000, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar1680minusImag
    ("KStar1680minusImag",
     "The imaginary part of the amplitude for the K*(1680)- in the K-matrix fit",
     &DtoKPiPiBaBar::_k1680mim, 1.690, -10., 10.0,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar892plusReal
    ("KStar892plusReal",
     "The real part of the amplitude for the K*(892)+ in the K-matrix fit",
     &DtoKPiPiBaBar::_k892pre,  0.133, -10., 10.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,double> interfaceKStar892plusImag
    ("KStar892plusImag",
     "The imaginary part of the amplitude for the K*(892)+ in the K-matrix fit",
     &DtoKPiPiBaBar::_k892pim, -0.132, -10., 10.,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy2> interfaceK0Star1430plusReal
    ("K0Star1430plusReal",
     "The real part of the amplitude for the K*0(1430)+ in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430pre0, GeV2, 0.375*GeV2, -10.*GeV2, 10.*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy2> interfaceK0Star1430plusImag
    ("K0Star1430plusImag",
     "The imaginary part of the amplitude for the K*0(1430)+ in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430pim0, GeV2, -0.143*GeV2, -10.*GeV2, 10.*GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,InvEnergy2> interfaceK2Star1430plusReal
    ("K2Star1430plusReal",
     "The real part of the amplitude for the K*2(1430)+ in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430pre2, 1./GeV2, 0.088/GeV2, -10./GeV2, 10./GeV2,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,InvEnergy2> interfaceK2Star1430plusImag
    ("K2Star1430plusImag",
     "The imaginary part of the amplitude for the K2*(1430)+ in the K-matrix fit",
     &DtoKPiPiBaBar::_k1430pim2, 1./GeV2, -0.057/GeV2, -10./GeV2, 10./GeV2,
     false, false, Interface::limited);

//   /**
//    *  Real part of the amplitude for \f$\rho(770)\f$
//    */
//   double _rho770re;

//   /**
//    *  Imaginary part of the amplitude for \f$\rho(770)\f$
//    */
//   double _rho770im;

//   /**
//    *  Real part of the amplitude for \f$\omega(782)\f$
//    */
//   double _omegare;

//   /**
//    *  Imaginary part of the amplitude for \f$\omega(782)\f$
//    */
//   double _omegaim;

//   /**
//    *  Real part of the amplitude for \f$f_2(1270)\f$
//    */
//   double _f2re;

//   /**
//    *  Imaginary part of the amplitude for \f$f_2(1270)\f$
//    */
//   double _f2im;

//   /**
//    *  Real part of the amplitude for \f$\rho(1450)\f$
//    */
//   double _rho1450re;

//   /**
//    *  Imaginary part of the amplitude for \f$\rho(1450)\f$
//    */
//   double _rho1450im;
//   //@}

//   /**
//    *  The amplitudes and phases for the normal fit
//    */
//   //@{
//   /**
//    *  Amplitude for \f$K^*(892)^-\f$
//    */
//   double _k892mamp;

//   /**
//    *  Phase for \f$K^*(892)^-\f$
//    */
//   double _k892mphase;

//   /**
//    *  Amplitude for \f$K^*_0(1430)^-\f$
//    */
//   double _k1430mamp0;

//   /**
//    *  Phase for \f$K^*_0(1430)^-\f$
//    */
//   double _k1430mphase0;

//   /**
//    *  Amplitude for \f$K^*_2(1430)^-\f$
//    */
//   double _k1430mamp2;

//   /**
//    *  Phase for \f$K^*_2(1430)^-\f$
//    */
//   double _k1430mphase2;

//   /**
//    *  Amplitude for \f$K^*(1410)^-\f$
//    */
//   double _k1410mamp;

//   /**
//    *  Phase for \f$K^*(1410)^-\f$
//    */
//   double _k1410mphase;

//   /**
//    *  Amplitude for \f$K^*(1680)^-\f$
//    */
//   double _k1680mamp;

//   /**
//    *  Phase for \f$K^*(1680)^-\f$
//    */
//   double _k1680mphase;

//   /**
//    *  Amplitude for \f$K^*(892)^+\f$
//    */
//   double _k892pamp;

//   /**
//    *  Phase for \f$K^*(892)^+\f$
//    */
//   double _k892pphase;

//   /**
//    *  Amplitude for \f$K^*_0(1430)^+\f$
//    */
//   double _k1430pamp0;

//   /**
//    *  Phase for \f$K^*_0(1430)^+\f$
//    */
//   double _k1430pphase0;

//   /**
//    *  Amplitude for \f$K^*_2(1430)^+\f$
//    */
//   double _k1430pamp2;

//   /**
//    *  Phase for \f$K^*_2(1430)^+\f$
//    */
//   double _k1430pphase2;

//   /**
//    *  Amplitude for \f$\rho(770)\f$
//    */
//   double _rho770amp;

//   /**
//    *  Phase for \f$\rho(770)\f$
//    */
//   double _rho770phase;

//   /**
//    *  Amplitude for \f$\omega(782)\f$
//    */
//   double _omegaamp;

//   /**
//    *  Phase for \f$\omega(782)\f$
//    */
//   double _omegaphase;

//   /**
//    *  Amplitude for \f$f_2(1270)\f$
//    */
//   double _f2amp;

//   /**
//    *  Phase for \f$f_2(1270)\f$
//    */
//   double _f2phase;

//   /**
//    *  Amplitude for \f$\rho(1450)\f$
//    */
//   double _rho1450amp;

//   /**
//    *  Phase for \f$\rho(1450)\f$
//    */
//   double _rho1450phase;

//   /**
//    *  Amplitude for \f$f_0(980)\f$
//    */
//   double _f980amp;

//   /**
//    *  Phase for \f$f_0(980)\f$
//    */
//   double _f980phase;

//   /**
//    *  Amplitude for \f$f_0(1370)\f$
//    */
//   double _f1370amp;

//   /**
//    *  Phase for \f$f_0(1370)\f$
//    */
//   double _f1370phase;

//   /**
//    *  Amplitude for \f$\sigma\f$
//    */
//   double _sigmaamp;

//   /**
//    *  Phase for \f$\sigma\f$
//    */
//   double _sigmaphase;

//   /**
//    *  Amplitude for \f$\sigma'\f$
//    */
//   double _sigmapamp;

//   /**
//    *  Phase for \f$\sigma'\f$
//    */
//   double _sigmapphase;

//   /**
//    *  Amplitude for the non-resonant component
//    */
//   double _nonamp;

//   /**
//    *  Phase for the non-resonant component
//    */
//   double _nonphase;
//   //@}


  static Parameter<DtoKPiPiBaBar,InvEnergy> interfaceRD0
    ("RD0",
     "The size of the D0 meson for the Blatt-Weisskopf form factor",
     &DtoKPiPiBaBar::_rD0, 1./GeV, 5./GeV, 0.0/GeV, 10.0/GeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,InvEnergy> interfaceRresonance
    ("Rresonance",
     "The size of the intermediate resonances for the Blatt-Weisskopf form factor",
     &DtoKPiPiBaBar::_rres, 1./GeV, 1.5/GeV, 0./GeV, 10./GeV,
     false, false, Interface::limited);
   
  static Parameter<DtoKPiPiBaBar,Energy> interfaceKStar892Mass
    ("KStar892Mass",
     "The mass of the K*(892)",
     &DtoKPiPiBaBar::_mK892, MeV, 891.66*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceK0Star1430Mass
    ("K0Star1430Mass",
     "The mass of the K*0(1430)",
     &DtoKPiPiBaBar::_mK14300, MeV, 1412.*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceK2Star1430Mass
    ("K2Star1430Mass",
     "The mass of the K*2(1430)",
     &DtoKPiPiBaBar::_mK14302, MeV, 1425.6*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);  

  static Parameter<DtoKPiPiBaBar,Energy> interfaceKStar1410Mass
    ("KStar1410Mass",
     "The mass of the K*(1410)",
     &DtoKPiPiBaBar::_mK1410, MeV, 1414.*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceKStar1680Mass
    ("KStar1680Mass",
     "The mass of the K*(1680)",
     &DtoKPiPiBaBar::_mK1680, MeV, 1717.*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceRho770Mass
    ("Rho770Mass",
     "The mass of the rho(770)",
     &DtoKPiPiBaBar::_mrho770, MeV, 775.8*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega(782)",
     &DtoKPiPiBaBar::_momega, MeV, 782.59*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacef2Mass
    ("f2Mass",
     "The mass of the f_2(1270)",
     &DtoKPiPiBaBar::_momega, MeV, 1275.40*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceRho1450Mass
    ("Rho1450Mass",
     "The mass of the rho(1450)",
     &DtoKPiPiBaBar::_mrho1450, MeV, 1465*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacef980Mass
    ("f980Mass",
     "The mass of the f_0(980)",
     &DtoKPiPiBaBar::_mf980, MeV, 977.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacef1370Mass
    ("f1370Mass",
     "The mass of the f_0(1370)",
     &DtoKPiPiBaBar::_mf1370, MeV, 1434.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma",
     &DtoKPiPiBaBar::_msigma, MeV, 484.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacesigmaprimeMass
    ("sigmaprimeMass",
     "The mass of the sigmaprime",
     &DtoKPiPiBaBar::_msigmap, MeV, 1014.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);
   
  static Parameter<DtoKPiPiBaBar,Energy> interfaceKStar892Width
    ("KStar892Width",
     "The width of the K*(892)",
     &DtoKPiPiBaBar::_wK892, MeV, 50.80*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceK0Star1430Width
    ("K0Star1430Width",
     "The width of the K*0(1430)",
     &DtoKPiPiBaBar::_wK14300, MeV,  294.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceK2Star1430Width
    ("K2Star1430Width",
     "The width of the K*2(1430)",
     &DtoKPiPiBaBar::_wK14302, MeV,  98.50*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);  

  static Parameter<DtoKPiPiBaBar,Energy> interfaceKStar1410Width
    ("KStar1410Width",
     "The width of the K*(1410)",
     &DtoKPiPiBaBar::_wK1410, MeV, 232.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceKStar1680Width
    ("KStar1680Width",
     "The width of the K*(1680)",
     &DtoKPiPiBaBar::_wK1680, MeV,  322.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceRho770Width
    ("Rho770Width",
     "The width of the rho(770)",
     &DtoKPiPiBaBar::_wrho770, MeV,  150.30*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega(782)",
     &DtoKPiPiBaBar::_womega, MeV,   8.49*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacef2Width
    ("f2Width",
     "The width of the f_2(1270)",
     &DtoKPiPiBaBar::_womega, MeV, 185.10*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfaceRho1450Width
    ("Rho1450Width",
     "The width of the rho(1450)",
     &DtoKPiPiBaBar::_wrho1450, MeV, 400.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacef980Width
    ("f980Width",
     "The width of the f_0(980)",
     &DtoKPiPiBaBar::_wf980, MeV,  44.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacef1370Width
    ("f1370Width",
     "The width of the f_0(1370)",
     &DtoKPiPiBaBar::_wf1370, MeV,  173.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma",
     &DtoKPiPiBaBar::_wsigma, MeV, 383.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<DtoKPiPiBaBar,Energy> interfacesigmaprimeWidth
    ("sigmaprimeWidth",
     "The width of the sigmaprime",
     &DtoKPiPiBaBar::_wsigmap, MeV,  88.00*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

//   /**
//    *  Integration parameters
//    */
//   //@{
//   /**
//    *  The maximum weight
//    */
//   double _maxwgt;

//   /**
//    *  The weights for the different channels
//    */
//   vector<double> _weights;
//   //@}
}

int DtoKPiPiBaBar::modeNumber(bool & cc,const DecayMode & dm) const {
  int id0(dm.parent()->id());
  // incoming particle must be D0
  if(abs(id0)!=ParticleID::D0) return -1;
  cc = id0==ParticleID::Dbar0;
  // must be three decay products
  if(dm.products().size()!=3) return -1;
  ParticleMSet::const_iterator pit = dm.products().begin();
  unsigned int npip(0),npim(0),nk(0);
  for( ;pit!=dm.products().end();++pit) {
    id0=(**pit).id();
    if(id0==ParticleID::piplus)       ++npip;
    else if(id0==ParticleID::piminus) ++npim;
    else if(abs(id0)==ParticleID::K0) ++nk;
    else if(id0==ParticleID::K_L0)    ++nk;
    else if(id0==ParticleID::K_S0)    ++nk;
  }
  if(npim==1&&npip==1&&nk==1) return  0;
  else                        return -1;
}

double DtoKPiPiBaBar::me2(bool vertex, const int ichan,
			    const Particle & inpart,
			    const ParticleVector & decay) const {
  // wavefunnction for the decaying particle
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // wavefunctions for the outgoing particles
  for(unsigned int ix=0;ix<3;++ix) {
    PPtr mytemp = decay[ix]; 
    ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }
  // compute the invariant masses needed to calulate the amplitudes
  Energy mA  = decay[0]->mass();
  Energy mB  = decay[1]->mass();
  Energy mC  = decay[2]->mass();
  Energy mD  = inpart.mass();
  Energy mAB = (decay[0]->momentum()+decay[1]->momentum()).m();
  Energy mAC = (decay[0]->momentum()+decay[2]->momentum()).m();
  Energy mBC = (decay[1]->momentum()+decay[2]->momentum()).m();
  // compute the amplitudes for the resonaces present in both models
  Complex amp(0);
  if(ichan<0) {
    amp =0.
      +_aKm892       *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK892   ,_wK892   )
      +_aKm14300/GeV2*amplitude(0,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300 ,_wK14300 )
      +_aKm14302*GeV2*amplitude(2,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14302 ,_wK14302 )
      +_aKm1410      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1410  ,_wK1410  )
      +_aKm1680      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680  ,_wK1680  )
      +_aKp892       *amplitude(1,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK892   ,_wK892   )
      +_aKp14300/GeV2*amplitude(0,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14300 ,_wK14300 )
      +_aKp14302*GeV2*amplitude(2,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14302 ,_wK14302 )
      +_arho770      *amplitude(1,true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrho770 ,_wrho770 )
      +_aomega       *amplitude(1,false,mD,mB,mC,mA,mBC,mAB,mAC,_momega  ,_womega  )
      +_af2     *GeV2*amplitude(2,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf2     ,_wf2     )
      +_arho1450     *amplitude(1,true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrho1450,_wrho1450)
      ;
    if(_imodel==0) {
      amp+= 
      _af980  /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf980  ,_wf980  )+
      _af1370 /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf1370 ,_wf1370 )+
      _asigma /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_msigma ,_wsigma )+
      _asigmap/GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_msigmap,_wsigmap)+
      _aNR
	;
    }
    else {
      amp+=F1(sqr(mBC));
    }
  }
  else if(ichan==0 ) {
    amp = _aKm892       *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK892   ,_wK892  );
  }
  else if(ichan==1 ) {
    amp = _aKm14300/GeV2*amplitude(0,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300 ,_wK14300);
  }
  else if(ichan==2 ) {
    amp = _aKm14302*GeV2*amplitude(2,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14302 ,_wK14302);
  }
  else if(ichan==3 ) {
    amp = _aKm1410      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1410  ,_wK1410 );
  }
  else if(ichan==4 ) {
    amp = _aKm1680      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680  ,_wK1680 );
  }
  else if(ichan==5 ) {
    amp = _aKp892       *amplitude(1,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK892   ,_wK892  );
  }
  else if(ichan==6 ) {
    amp = _aKp14300/GeV2*amplitude(0,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14300 ,_wK14300);
  }
  else if(ichan==7 ) {
    amp = _aKp14302*GeV2*amplitude(2,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14302 ,_wK14302);
  }
  else if(ichan==8 ) {
    amp = _af980  /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf980  ,_wf980  );
  }
  else if(ichan==9 ) {
    amp = _af1370 /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf1370 ,_wf1370 );
  }
  else if(ichan==10) {
    amp = _asigma /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_msigma ,_wsigma )+
      _asigmap/GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_msigmap,_wsigmap);
  }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newME(0,0,0,0)=amp;
  ME(newME);
  return real(amp*conj(amp));
}

void DtoKPiPiBaBar::dataBaseOutput(ofstream &, bool) const {}

Complex DtoKPiPiBaBar::F1(Energy2 s) const {
  // rho values
  vector<Complex> rho(5);
  if(s>=4.*sqr(_mpi))      rho[0] =            sqrt( 1.-sqr(_mpi  + _mpi  )/s);
  else                     rho[0] = Complex(0.,sqrt(-1.+sqr(_mpi  + _mpi  )/s));
  if(s>=4.*sqr(_mK))       rho[1] = sqrt(1.-sqr(_mK   + _mK   )/s);
  else                     rho[1] = Complex(0.,sqrt(-1.+sqr(_mK   + _mK   )/s));
  if(s<=16.*sqr(_mpi))     rho[2] = 0.;
  else if(s<GeV2)          rho[2] = (*_rho3inter)(s/GeV2);
  else                     rho[2] = 1.-16.*sqr(_mpi)/s;
  if(s>=4.*sqr(_meta))     rho[3] =            sqrt( 1.-sqr(_meta + _meta )/s);
  else                     rho[3] = Complex(0.,sqrt(-1.+sqr(_meta + _meta )/s));
  if(s>=sqr(_meta+_metap)) rho[4] =            sqrt( 1.-sqr(_meta + _metap)/s);
  else                     rho[4] = Complex(0.,sqrt(-1.+sqr(_meta + _metap)/s));
  // K -matrix
  Complex K[5][5];
  for(unsigned int i=0;i<5;++i) {
    for(unsigned int j=0;j<5;++j) {
      K[i][j]=0.;
      for(unsigned int alpha=0;alpha<5;++alpha) {
	K[i][j]+=_gialpha[i][alpha]*_gialpha[j][alpha]/(sqr(_malpha[alpha])-s);
      }
      if(i==0) K[i][j]+= _fscatt[j]*(GeV2-_s0scatt)/(s-_s0scatt);
      else if(j==0) K[i][j]+= _fscatt[i]*(GeV2-_s0scatt)/(s-_s0scatt);
      K[i][j] *= (GeV2-_sA0)/(s-_sA0)*(s-0.5*_sA*sqr(_mpi)/GeV2)/GeV2;
      // multiply in rho
      K[i][j] *= -Complex(0.,1.)*rho[j];
      // add identity
      if(i==j) K[i][j]+=1.;
    }
  }
  vector<Complex> m(25);
  for(unsigned int i=0;i<5;++i) {
    for(unsigned int j=0;j<5;++j) {
      m[i*5+j]=K[i][j];
    }
  }
  // Find all NECESSARY 2x2 dets:  (30 of them)
  Complex Det2_23_01 = m[10]*m[16] - m[11]*m[15];
  Complex Det2_23_02 = m[10]*m[17] - m[12]*m[15];
  Complex Det2_23_03 = m[10]*m[18] - m[13]*m[15];
  Complex Det2_23_04 = m[10]*m[19] - m[14]*m[15];
  Complex Det2_23_12 = m[11]*m[17] - m[12]*m[16];
  Complex Det2_23_13 = m[11]*m[18] - m[13]*m[16];
  Complex Det2_23_14 = m[11]*m[19] - m[14]*m[16];
  Complex Det2_23_23 = m[12]*m[18] - m[13]*m[17];
  Complex Det2_23_24 = m[12]*m[19] - m[14]*m[17];
  Complex Det2_23_34 = m[13]*m[19] - m[14]*m[18];
  Complex Det2_24_01 = m[10]*m[21] - m[11]*m[20];
  Complex Det2_24_02 = m[10]*m[22] - m[12]*m[20];
  Complex Det2_24_03 = m[10]*m[23] - m[13]*m[20];
  Complex Det2_24_04 = m[10]*m[24] - m[14]*m[20];
  Complex Det2_24_12 = m[11]*m[22] - m[12]*m[21];
  Complex Det2_24_13 = m[11]*m[23] - m[13]*m[21];
  Complex Det2_24_14 = m[11]*m[24] - m[14]*m[21];
  Complex Det2_24_23 = m[12]*m[23] - m[13]*m[22];
  Complex Det2_24_24 = m[12]*m[24] - m[14]*m[22];
  Complex Det2_24_34 = m[13]*m[24] - m[14]*m[23];
  Complex Det2_34_01 = m[15]*m[21] - m[16]*m[20];
  Complex Det2_34_02 = m[15]*m[22] - m[17]*m[20];
  Complex Det2_34_03 = m[15]*m[23] - m[18]*m[20];
  Complex Det2_34_04 = m[15]*m[24] - m[19]*m[20];
  Complex Det2_34_12 = m[16]*m[22] - m[17]*m[21];
  Complex Det2_34_13 = m[16]*m[23] - m[18]*m[21];
  Complex Det2_34_14 = m[16]*m[24] - m[19]*m[21];
  Complex Det2_34_23 = m[17]*m[23] - m[18]*m[22];
  Complex Det2_34_24 = m[17]*m[24] - m[19]*m[22];
  Complex Det2_34_34 = m[18]*m[24] - m[19]*m[23];
  // Find all NECESSARY 3x3 dets:   (40 of them)
  Complex Det3_123_012 = m[ 5]*Det2_23_12 - m[ 6]*Det2_23_02 + m[ 7]*Det2_23_01;
  Complex Det3_123_013 = m[ 5]*Det2_23_13 - m[ 6]*Det2_23_03 + m[ 8]*Det2_23_01;
  Complex Det3_123_014 = m[ 5]*Det2_23_14 - m[ 6]*Det2_23_04 + m[ 9]*Det2_23_01;
  Complex Det3_123_023 = m[ 5]*Det2_23_23 - m[ 7]*Det2_23_03 + m[ 8]*Det2_23_02;
  Complex Det3_123_024 = m[ 5]*Det2_23_24 - m[ 7]*Det2_23_04 + m[ 9]*Det2_23_02;
  Complex Det3_123_034 = m[ 5]*Det2_23_34 - m[ 8]*Det2_23_04 + m[ 9]*Det2_23_03;
  Complex Det3_123_123 = m[ 6]*Det2_23_23 - m[ 7]*Det2_23_13 + m[ 8]*Det2_23_12;
  Complex Det3_123_124 = m[ 6]*Det2_23_24 - m[ 7]*Det2_23_14 + m[ 9]*Det2_23_12;
  Complex Det3_123_134 = m[ 6]*Det2_23_34 - m[ 8]*Det2_23_14 + m[ 9]*Det2_23_13;
  Complex Det3_123_234 = m[ 7]*Det2_23_34 - m[ 8]*Det2_23_24 + m[ 9]*Det2_23_23;
  Complex Det3_124_012 = m[ 5]*Det2_24_12 - m[ 6]*Det2_24_02 + m[ 7]*Det2_24_01;
  Complex Det3_124_013 = m[ 5]*Det2_24_13 - m[ 6]*Det2_24_03 + m[ 8]*Det2_24_01;
  Complex Det3_124_014 = m[ 5]*Det2_24_14 - m[ 6]*Det2_24_04 + m[ 9]*Det2_24_01;
  Complex Det3_124_023 = m[ 5]*Det2_24_23 - m[ 7]*Det2_24_03 + m[ 8]*Det2_24_02;
  Complex Det3_124_024 = m[ 5]*Det2_24_24 - m[ 7]*Det2_24_04 + m[ 9]*Det2_24_02;
  Complex Det3_124_034 = m[ 5]*Det2_24_34 - m[ 8]*Det2_24_04 + m[ 9]*Det2_24_03;
  Complex Det3_124_123 = m[ 6]*Det2_24_23 - m[ 7]*Det2_24_13 + m[ 8]*Det2_24_12;
  Complex Det3_124_124 = m[ 6]*Det2_24_24 - m[ 7]*Det2_24_14 + m[ 9]*Det2_24_12;
  Complex Det3_124_134 = m[ 6]*Det2_24_34 - m[ 8]*Det2_24_14 + m[ 9]*Det2_24_13;
  Complex Det3_124_234 = m[ 7]*Det2_24_34 - m[ 8]*Det2_24_24 + m[ 9]*Det2_24_23;
  Complex Det3_134_012 = m[ 5]*Det2_34_12 - m[ 6]*Det2_34_02 + m[ 7]*Det2_34_01;
  Complex Det3_134_013 = m[ 5]*Det2_34_13 - m[ 6]*Det2_34_03 + m[ 8]*Det2_34_01;
  Complex Det3_134_014 = m[ 5]*Det2_34_14 - m[ 6]*Det2_34_04 + m[ 9]*Det2_34_01;
  Complex Det3_134_023 = m[ 5]*Det2_34_23 - m[ 7]*Det2_34_03 + m[ 8]*Det2_34_02;
  Complex Det3_134_024 = m[ 5]*Det2_34_24 - m[ 7]*Det2_34_04 + m[ 9]*Det2_34_02;
  Complex Det3_134_034 = m[ 5]*Det2_34_34 - m[ 8]*Det2_34_04 + m[ 9]*Det2_34_03;
  Complex Det3_134_123 = m[ 6]*Det2_34_23 - m[ 7]*Det2_34_13 + m[ 8]*Det2_34_12;
  Complex Det3_134_124 = m[ 6]*Det2_34_24 - m[ 7]*Det2_34_14 + m[ 9]*Det2_34_12;
  Complex Det3_134_134 = m[ 6]*Det2_34_34 - m[ 8]*Det2_34_14 + m[ 9]*Det2_34_13;
  Complex Det3_134_234 = m[ 7]*Det2_34_34 - m[ 8]*Det2_34_24 + m[ 9]*Det2_34_23;
  Complex Det3_234_012 = m[10]*Det2_34_12 - m[11]*Det2_34_02 + m[12]*Det2_34_01;
  Complex Det3_234_013 = m[10]*Det2_34_13 - m[11]*Det2_34_03 + m[13]*Det2_34_01;
  Complex Det3_234_014 = m[10]*Det2_34_14 - m[11]*Det2_34_04 + m[14]*Det2_34_01;
  Complex Det3_234_023 = m[10]*Det2_34_23 - m[12]*Det2_34_03 + m[13]*Det2_34_02;
  Complex Det3_234_024 = m[10]*Det2_34_24 - m[12]*Det2_34_04 + m[14]*Det2_34_02;
  Complex Det3_234_034 = m[10]*Det2_34_34 - m[13]*Det2_34_04 + m[14]*Det2_34_03;
  Complex Det3_234_123 = m[11]*Det2_34_23 - m[12]*Det2_34_13 + m[13]*Det2_34_12;
  Complex Det3_234_124 = m[11]*Det2_34_24 - m[12]*Det2_34_14 + m[14]*Det2_34_12;
  Complex Det3_234_134 = m[11]*Det2_34_34 - m[13]*Det2_34_14 + m[14]*Det2_34_13;
  Complex Det3_234_234 = m[12]*Det2_34_34 - m[13]*Det2_34_24 + m[14]*Det2_34_23;
  // Find all NECESSARY 4x4 dets:   (25 of them)
  Complex Det4_0123_0123 = m[0]*Det3_123_123 - m[1]*Det3_123_023 
				+ m[2]*Det3_123_013 - m[3]*Det3_123_012;
  Complex Det4_0123_0124 = m[0]*Det3_123_124 - m[1]*Det3_123_024 
				+ m[2]*Det3_123_014 - m[4]*Det3_123_012;
  Complex Det4_0123_0134 = m[0]*Det3_123_134 - m[1]*Det3_123_034 
				+ m[3]*Det3_123_014 - m[4]*Det3_123_013;
  Complex Det4_0123_0234 = m[0]*Det3_123_234 - m[2]*Det3_123_034 
				+ m[3]*Det3_123_024 - m[4]*Det3_123_023;
  Complex Det4_0123_1234 = m[1]*Det3_123_234 - m[2]*Det3_123_134 
				+ m[3]*Det3_123_124 - m[4]*Det3_123_123;
  Complex Det4_0124_0123 = m[0]*Det3_124_123 - m[1]*Det3_124_023 
				+ m[2]*Det3_124_013 - m[3]*Det3_124_012;
  Complex Det4_0124_0124 = m[0]*Det3_124_124 - m[1]*Det3_124_024 
				+ m[2]*Det3_124_014 - m[4]*Det3_124_012;
  Complex Det4_0124_0134 = m[0]*Det3_124_134 - m[1]*Det3_124_034 
				+ m[3]*Det3_124_014 - m[4]*Det3_124_013;
  Complex Det4_0124_0234 = m[0]*Det3_124_234 - m[2]*Det3_124_034 
				+ m[3]*Det3_124_024 - m[4]*Det3_124_023;
  Complex Det4_0124_1234 = m[1]*Det3_124_234 - m[2]*Det3_124_134 
				+ m[3]*Det3_124_124 - m[4]*Det3_124_123;
  Complex Det4_0134_0123 = m[0]*Det3_134_123 - m[1]*Det3_134_023 
				+ m[2]*Det3_134_013 - m[3]*Det3_134_012;
  Complex Det4_0134_0124 = m[0]*Det3_134_124 - m[1]*Det3_134_024 
				+ m[2]*Det3_134_014 - m[4]*Det3_134_012;
  Complex Det4_0134_0134 = m[0]*Det3_134_134 - m[1]*Det3_134_034 
				+ m[3]*Det3_134_014 - m[4]*Det3_134_013;
  Complex Det4_0134_0234 = m[0]*Det3_134_234 - m[2]*Det3_134_034 
				+ m[3]*Det3_134_024 - m[4]*Det3_134_023;
  Complex Det4_0134_1234 = m[1]*Det3_134_234 - m[2]*Det3_134_134 
				+ m[3]*Det3_134_124 - m[4]*Det3_134_123;
  Complex Det4_0234_0123 = m[0]*Det3_234_123 - m[1]*Det3_234_023 
				+ m[2]*Det3_234_013 - m[3]*Det3_234_012;
  Complex Det4_0234_0124 = m[0]*Det3_234_124 - m[1]*Det3_234_024 
				+ m[2]*Det3_234_014 - m[4]*Det3_234_012;
  Complex Det4_0234_0134 = m[0]*Det3_234_134 - m[1]*Det3_234_034 
				+ m[3]*Det3_234_014 - m[4]*Det3_234_013;
  Complex Det4_0234_0234 = m[0]*Det3_234_234 - m[2]*Det3_234_034 
				+ m[3]*Det3_234_024 - m[4]*Det3_234_023;
  Complex Det4_0234_1234 = m[1]*Det3_234_234 - m[2]*Det3_234_134 
				+ m[3]*Det3_234_124 - m[4]*Det3_234_123;
  Complex Det4_1234_0123 = m[5]*Det3_234_123 - m[6]*Det3_234_023 
				+ m[7]*Det3_234_013 - m[8]*Det3_234_012;
  Complex Det4_1234_0124 = m[5]*Det3_234_124 - m[6]*Det3_234_024 
				+ m[7]*Det3_234_014 - m[9]*Det3_234_012;
  Complex Det4_1234_0134 = m[5]*Det3_234_134 - m[6]*Det3_234_034 
				+ m[8]*Det3_234_014 - m[9]*Det3_234_013;
  Complex Det4_1234_0234 = m[5]*Det3_234_234 - m[7]*Det3_234_034 
				+ m[8]*Det3_234_024 - m[9]*Det3_234_023;
  Complex Det4_1234_1234 = m[6]*Det3_234_234 - m[7]*Det3_234_134 
				+ m[8]*Det3_234_124 - m[9]*Det3_234_123;
  // Find the 5x5 det:
  Complex det =   m[0]*Det4_1234_1234 
	 	- m[1]*Det4_1234_0234 
		+ m[2]*Det4_1234_0134 
		- m[3]*Det4_1234_0124 
		+ m[4]*Det4_1234_0123;
  if ( det == Complex(0.,0.) ) {
    throw Exception() << "Matrix inversion fails in DtoKPiPiBaBar::F1()"
		      << Exception::runerror;
  } 
  Complex oneOverDet = 1.0/det;
  Complex mn1OverDet = - oneOverDet;
  Complex inverse[5][5];
  inverse[0][0] =  Det4_1234_1234 * oneOverDet;
  inverse[0][1] =  Det4_0234_1234 * mn1OverDet;
  inverse[0][2] =  Det4_0134_1234 * oneOverDet;
  inverse[0][3] =  Det4_0124_1234 * mn1OverDet;
  inverse[0][4] =  Det4_0123_1234 * oneOverDet;
  inverse[1][0] =  Det4_1234_0234 * mn1OverDet;
  inverse[1][1] =  Det4_0234_0234 * oneOverDet;
  inverse[1][2] =  Det4_0134_0234 * mn1OverDet;
  inverse[1][3] =  Det4_0124_0234 * oneOverDet;
  inverse[1][4] =  Det4_0123_0234 * mn1OverDet;
  inverse[2][0] =  Det4_1234_0134 * oneOverDet;
  inverse[2][1] =  Det4_0234_0134 * mn1OverDet;
  inverse[2][2] =  Det4_0134_0134 * oneOverDet;
  inverse[2][3] =  Det4_0124_0134 * mn1OverDet;
  inverse[2][4] =  Det4_0123_0134 * oneOverDet;
  inverse[3][0] =  Det4_1234_0124 * mn1OverDet;
  inverse[3][1] =  Det4_0234_0124 * oneOverDet;
  inverse[3][2] =  Det4_0134_0124 * mn1OverDet;
  inverse[3][3] =  Det4_0124_0124 * oneOverDet;
  inverse[3][4] =  Det4_0123_0124 * mn1OverDet;
  inverse[4][0] =  Det4_1234_0123 * oneOverDet;
  inverse[4][1] =  Det4_0234_0123 * mn1OverDet;
  inverse[4][2] =  Det4_0134_0123 * oneOverDet;
  inverse[4][3] =  Det4_0124_0123 * mn1OverDet;
  inverse[4][4] =  Det4_0123_0123 * oneOverDet;
  // calculate the production vector
  Complex p[5];
  for(unsigned int i=0;i<5;++i) {
    p[i]=0.;
    for(unsigned int alpha=0;alpha<5;++alpha) {
      p[i]+=_gialpha[i][alpha]*_beta[alpha]/(sqr(_malpha[alpha])-s);
    }
    if(i==0) p[i]+=_fprod*(GeV2-_s0scatt)/(s-_s0scatt);
  }
  // finally compute the answer
  Complex output(0.);
  for(unsigned int j=0;j<5;++j) output += inverse[1][j]*p[j];
  return output;
}
