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
  _gpipi.push_back( 0.229);
  _gpipi.push_back( 0.941);
  _gpipi.push_back( 0.369);
  _gpipi.push_back( 0.337);
  _gpipi.push_back( 0.182);
  // g_KK
  _gKK.push_back(-0.554);
  _gKK.push_back( 0.551);
  _gKK.push_back( 0.239);
  _gKK.push_back( 0.409);
  _gKK.push_back(-0.176);
  // g_4pi
  _g4pi.push_back( 0.000);
  _g4pi.push_back( 0.000);
  _g4pi.push_back( 0.556);
  _g4pi.push_back( 0.857);
  _g4pi.push_back(-0.797);
  // g_etaeta
  _getaeta.push_back(-0.399);
  _getaeta.push_back( 0.391);
  _getaeta.push_back( 0.183);
  _getaeta.push_back( 0.199);
  _getaeta.push_back(-0.004);
  // g_etaetap
  _getaetap.push_back(-0.346);
  _getaetap.push_back( 0.315);
  _getaetap.push_back( 0.187);
  _getaetap.push_back(-0.010);
  _getaetap.push_back( 0.224);
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
  _betare.push_back( -3.78);_betaim.push_back(1.23);
  _betare.push_back(  9.55);_betaim.push_back(3.43);
  _betare.push_back(  0.00);_betaim.push_back(0.00);
  _betare.push_back( 12.97);_betaim.push_back(1.27);
  _betare.push_back(  0.00);_betaim.push_back(0.00);
  // f^prod
  _fprodre.push_back(-10.22);_fprodim.push_back( -6.35);
  _fprodre.push_back(  0.00);_fprodim.push_back(  0.00);
  _fprodre.push_back(  0.00);_fprodim.push_back(  0.00);
  _fprodre.push_back(  0.00);_fprodim.push_back(  0.00);
  _fprodre.push_back(  0.00);_fprodim.push_back(  0.00);
  // The amplitudes and phases for the K-matrix fit
  _k892mre   =-1.159;
  _k892mim   = 1.361;
  _k1430mre0 = 2.482*GeV2;
  _k1430mim0 =-0.653*GeV2;
  _k1430mre2 = 0.852/GeV2;
  _k1430mim2 =-0.729/GeV2;
  _k1410mre  =-0.402;
  _k1410mim  = 0.050;
  _k1680mre  =-1.000;
  _k1680mim  = 1.690;
  _k892pre   = 0.133;
  _k892pim   =-0.132;
  _k1430pre0 = 0.375*GeV2;
  _k1430pim0 =-0.143*GeV2;
  _k1430pre2 = 0.088/GeV2;
  _k1430pim2 =-0.057/GeV2;
  _rho770re  = 1.000;
  _rho770im  = 0.000;
  _omegare   =-0.0182;
  _omegaim   =-0.0367;
  _f2re      = 0.787/GeV2;
  _f2im      =-0.397/GeV2;
  _rho1450re = 0.405;
  _rho1450im =-0.458;
  // amplitudes and phase for the normal fit
  _k892mamp     = 1.781;
  _k892mphase   = 131.0;
  _k1430mamp0   = 2.45*GeV2;
  _k1430mphase0 =-8.3;
  _k1430mamp2   = 1.05/GeV2;
  _k1430mphase2 =-54.3;
  _k1410mamp    = 0.52;
  _k1410mphase  = 154;
  _k1680mamp    = 0.89;
  _k1680mphase  =-139;
  _k892pamp     = 0.180;
  _k892pphase   =-44.1;
  _k1430pamp0   = 0.37*GeV2;
  _k1430pphase0 = 18;
  _k1430pamp2   = 0.075/GeV2;
  _k1430pphase2 =-104;
  _rho770amp    = 1.0;
  _rho770phase  = 0.;
  _omegaamp     = 0.0391;
  _omegaphase   = 115.3;
  _f2amp        = 0.922/GeV2;
  _f2phase      =-21.3;
  _rho1450amp   = 0.52;
  _rho1450phase = 38;
  _f980amp      = 0.482*GeV2;
  _f980phase    =-141.8;
  _f1370amp     = 2.25*GeV2;
  _f1370phase   = 113.2;
  _sigmaamp     = 1.36*GeV2;
  _sigmaphase   =-177.9;
  _sigmapamp    = 0.34*GeV2;
  _sigmapphase  = 153;
  _nonamp       = 3.53;
  _nonphase     = 128;
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
  cerr << "testing _aKm892   " <<_aKm892   << " " << _k892mphase   << "\n";  
  cerr << "testing _aKm14300 " <<_aKm14300 << " " << _k1430mphase0 << "\n";
  cerr << "testing _aKm14302 " <<_aKm14302 << " " << _k1430mphase2 << "\n";
  cerr << "testing _aKm1410  " <<_aKm1410  << " " << _k1410mphase  << "\n";
  cerr << "testing _aKm1680  " <<_aKm1680  << " " << _k1680mphase  << "\n";
  cerr << "testing _aKp892   " <<_aKp892   << " " << _k892pphase   << "\n";  
  cerr << "testing _aKp14300 " <<_aKp14300 << " " << _k1430pphase0 << "\n";
  cerr << "testing _aKp14302 " <<_aKp14302 << " " << _k1430pphase2 << "\n";
  cerr << "testing _arho770  " <<_arho770  << " " << _rho770phase  << "\n";
  cerr << "testing _aomega   " <<_aomega   << " " << _omegaphase   << "\n";  
  cerr << "testing _af2      " <<_af2      << " " << _f2phase      << "\n";  
  cerr << "testing _arho1450 " <<_arho1450 << " " << _rho1450phase << "\n";
  cerr << "testing _af980    " <<_af980    << " " << _f980phase    << "\n";  
  cerr << "testing _af1370   " <<_af1370   << " " << _f1370phase   << "\n";  
  cerr << "testing _asigma   " <<_asigma   << " " << _sigmaphase   << "\n";  
  cerr << "testing _asigmap  " <<_asigmap  << " " << _sigmapphase  << "\n";
  cerr << "testing _aNR      " <<_aNR      << " " << _nonphase     << "\n";  
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
     << _weights;
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
     >> _weights;
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
    amp = 
      _aKm892       *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK892   ,_wK892   )+
      _aKm14300/GeV2*amplitude(0,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14300 ,_wK14300 )+
      _aKm14302*GeV2*amplitude(2,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK14302 ,_wK14302 )+
      _aKm1410      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1410  ,_wK1410  )+
      _aKm1680      *amplitude(1,false,mD,mA,mC,mB,mAC,mAB,mBC,_mK1680  ,_wK1680  )+
      _aKp892       *amplitude(1,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK892   ,_wK892   )+
      _aKp14300/GeV2*amplitude(0,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14300 ,_wK14300 )+
      _aKp14302*GeV2*amplitude(2,false,mD,mA,mB,mC,mAB,mAC,mBC,_mK14302 ,_wK14302 )+
      _arho770      *amplitude(1,true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrho770 ,_wrho770 )+
      _aomega       *amplitude(1,false,mD,mB,mC,mA,mBC,mAB,mAC,_momega  ,_womega  )+
      _af2     *GeV2*amplitude(2,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf2     ,_wf2     )+
      _arho1450     *amplitude(1,true ,mD,mB,mC,mA,mBC,mAB,mAC,_mrho1450,_wrho1450);
    if(_imodel==0) {
      amp+= 
	_af980  /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf980  ,_wf980  )+
	_af1370 /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_mf1370 ,_wf1370 )+
	_asigma /GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_msigma ,_wsigma )+
	_asigmap/GeV2*amplitude(0,false,mD,mB,mC,mA,mBC,mAB,mAC,_msigmap,_wsigmap)+
	_aNR;
    }
    else {
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

void DtoKPiPiBaBar::dataBaseOutput(ofstream & output,
					       bool header) const {
}
