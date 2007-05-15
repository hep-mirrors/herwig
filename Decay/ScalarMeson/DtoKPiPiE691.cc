// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiE691 class.
//

#include "DtoKPiPiE691.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

DtoKPiPiE691::DtoKPiPiE691() {
  // amplitudes and phases for D+ -> K-pi+pi+
  _a1NR    = 1.  ; _phi1NR    =    0.;
  _a1K892  = 0.78; _phi1K892  = - 60.;
  _a1K1430 = 0.53; _phi1K1430 =  132.;
  _a1K1680 = 0.47; _phi1K1680 = - 51.;
  // amplitudes and phases for D0 -> K-pi+pi0
  _a2NR    = 1.00; _phi2NR    =    0.;
  _a2K8920 = 3.19; _phi2K8920 =  167.;
  _a2K892m = 2.96; _phi2K892m = -112.;
  _a2rho   = 8.56; _phi2rho   =   40.;
  // amplitudes and phases for D0 -> Kbar0 pi+pi-
  _a3NR    = 1.00; _phi3NR    =    0.;
  _a3K892  = 2.31; _phi3K892  =  109.;
  _a3rho   = 1.59; _phi3rho   = -123.;
  // masses and widths
  _localparameters=true;
  _mK8920 = 0.8961 *GeV; _wK8920 = 0.0505*GeV;
  _mK892m = 0.89159*GeV; _wK892m = 0.0498*GeV;
  _mK1680 = 1.714  *GeV; _wK1680 = 0.323 *GeV;
  _mK1430 = 1.429  *GeV; _wK1430 = 0.287 *GeV;
  _mrho0  = 0.7681 *GeV; _wrho0  = 0.1515*GeV;
  _mrhop  = 0.7681 *GeV; _wrhop  = 0.1515*GeV;
}

void DtoKPiPiE691::doinit() throw(InitException) {
  cerr << "testing in doinit\n";
  DecayIntegrator::doinit();
  // complex amplitudes calculated from magnitudes and phases
  double fact = pi/180.;
  // D+ -> K-pi+pi+
  _c1NR   = _a1NR   *Complex(cos(_phi1NR   *fact),sin(_phi1NR   *fact)); 
  _c1K892 = _a1K892 *Complex(cos(_phi1K892 *fact),sin(_phi1K892 *fact));
  _c1K1430= _a1K1430*Complex(cos(_phi1K1430*fact),sin(_phi1K1430*fact));
  _c1K1680= _a1K1680*Complex(cos(_phi1K1680*fact),sin(_phi1K1680*fact));
  // D0 -> K-pi+pi0
  _c2NR   = _a2NR   *Complex(cos(_phi2NR   *fact),sin(_phi2NR   *fact));
  _c2K8920= _a2K8920*Complex(cos(_phi2K8920*fact),sin(_phi2K8920*fact));
  _c2K892m= _a2K892m*Complex(cos(_phi2K892m*fact),sin(_phi2K892m*fact));
  _c2rho  = _a2rho  *Complex(cos(_phi2rho  *fact),sin(_phi2rho  *fact));
  // D0 -> Kbar0pi+pi-
  _c3NR   = _a3NR   *Complex(cos(_phi3NR   *fact),sin(_phi3NR   *fact));
  _c3K892 = _a3K892 *Complex(cos(_phi3K892 *fact),sin(_phi3K892 *fact));
  _c3rho  = _a3rho  *Complex(cos(_phi3rho  *fact),sin(_phi3rho  *fact));
  // intermediate resonances
  tPDPtr k8920 = getParticleData(ParticleID::Kstarbar0 );
  tPDPtr k892m = getParticleData(ParticleID::Kstarminus);
  tPDPtr k1430 = getParticleData(ParticleID::Kstar_0bar0);
  tPDPtr k1680 = getParticleData(-30313);
  tPDPtr rho0  = getParticleData(ParticleID::rho0);
  tPDPtr rhop  = getParticleData(ParticleID::rhoplus);
  // D+ -> K-pi+pi+
  PDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::Dplus);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piplus);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int ix=0;
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k8920,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k1430) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1430,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1430,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k1680) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1680,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k1680,0,0., 1,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  // add the mode
  vector<double> wtemp;
  if(ix<=int(_weights.size())) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit,wit+ix);
  }
  else {
    wtemp=vector<double>(ix,1./double(ix));
  }
  if(_maxwgt.empty()) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[0],wtemp);
  // D0 -> K-pi+pi0
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kminus);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::pi0);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int iy=ix;
  if(k8920) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k8920,0,0., 1,2);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(rhop) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rhop,0,0., 2,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  // add the mode
  if(iy<=int(_weights.size())) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<2) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[1],wtemp);
  // D0 -> Kbar0 pi+ pi-
  extpart[0]=getParticleData(ParticleID::D0);
  extpart[1]=getParticleData(ParticleID::Kbar0);
  extpart[2]=getParticleData(ParticleID::piplus);
  extpart[3]=getParticleData(ParticleID::piminus);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  ix=iy;
  if(rho0) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho0,0,0., 2,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k892m) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(k892m,0,0., 1,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  // add the mode
  if(iy<=int(_weights.size())) {
    vector<double>::const_iterator wit=_weights.begin();
    wtemp=vector<double>(wit+ix,wit+iy);
  }
  else {
    wtemp=vector<double>(iy-ix,1./double(iy-ix));
  }
  if(_maxwgt.size()<3) _maxwgt.push_back(1.);
  addMode(mode,_maxwgt[2],wtemp);
  // reset the resonance parameters in the integration if needed
  if(_localparameters) {
    resetIntermediate(k8920,_mK8920,_wK8920);
    resetIntermediate(k892m,_mK892m,_wK892m);
    resetIntermediate(k1430,_mK1680,_wK1680);
    resetIntermediate(k1680,_mK1430,_wK1430);
    resetIntermediate(rho0 ,_mrho0 ,_wrho0 );
    resetIntermediate(rhop ,_mrhop ,_wrhop );
  }
  // get values from the ParticleData objects if needed
  else {
    _mK8920 = k8920->mass();
    _mK892m = k892m->mass();
    _mK1680 = k1430->mass();
    _mK1430 = k1680->mass();
    _mrho0  = rho0 ->mass();
    _mrhop  = rhop ->mass();
    _wK8920 = k8920->width();
    _wK892m = k892m->width();
    _wK1680 = k1430->width();
    _wK1430 = k1680->width();
    _wrho0  = rho0 ->width();
    _wrhop  = rhop ->width();
  }
}

void DtoKPiPiE691::persistentOutput(PersistentOStream & os) const {
  os << _a1NR << _phi1NR << _a1K892 << _phi1K892 << _a1K1430 << _phi1K1430 
     << _a1K1680 << _phi1K1680 << _a2NR << _phi2NR << _a2K8920 << _phi2K8920 
     << _a2K892m << _phi2K892m << _a2rho << _phi2rho << _a3NR << _phi3NR 
     << _a3K892 << _phi3K892 << _a3rho << _phi3rho << _c1NR << _c1K892 << _c1K1430 
     << _c1K1680 << _c2NR << _c2K8920 << _c2K892m << _c2rho << _c3NR << _c3K892 
     << _c3rho << _localparameters << _mK8920 << _wK8920 << _mK892m << _wK892m 
     << _mK1680 << _wK1680 << _mK1430 << _wK1430 << _mrho0 << _wrho0 << _mrhop 
     << _wrhop << _maxwgt << _weights;
}

void DtoKPiPiE691::persistentInput(PersistentIStream & is, int) {
  is >> _a1NR >> _phi1NR >> _a1K892 >> _phi1K892 >> _a1K1430 >> _phi1K1430 
     >> _a1K1680 >> _phi1K1680 >> _a2NR >> _phi2NR >> _a2K8920 >> _phi2K8920 
     >> _a2K892m >> _phi2K892m >> _a2rho >> _phi2rho >> _a3NR >> _phi3NR 
     >> _a3K892 >> _phi3K892 >> _a3rho >> _phi3rho >> _c1NR >> _c1K892 >> _c1K1430 
     >> _c1K1680 >> _c2NR >> _c2K8920 >> _c2K892m >> _c2rho >> _c3NR >> _c3K892 
     >> _c3rho >> _localparameters >> _mK8920 >> _wK8920 >> _mK892m >> _wK892m 
     >> _mK1680 >> _wK1680 >> _mK1430 >> _wK1430 >> _mrho0 >> _wrho0 >> _mrhop 
     >> _wrhop >> _maxwgt >> _weights;
}

ClassDescription<DtoKPiPiE691> DtoKPiPiE691::initDtoKPiPiE691;
// Definition of the static class description member.

void DtoKPiPiE691::Init() {

  static ClassDocumentation<DtoKPiPiE691> documentation
    ("There is no documentation for the DtoKPiPiE691 class");


//   /**
//    *  Amplitudes and phases for the different components
//    */
//   //@{
//   /**
//    *  Amplitude of the non-resonant component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _a1NR;

//   /**
//    *  Phase of the non-resonant component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _phi1NR;

//   /**
//    *  Amplitude of the \f$\bar{K}^*(892)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _a1K892;

//   /**
//    *  Phase of the \f$\bar{K}^*(892)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _phi1K892;

//   /**
//    *  Amplitude of the \f$\bar{K}^*_0(1430)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _a1K1430;

//   /**
//    *  Phase of the \f$\bar{K}^*_0(1430)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _phi1K1430;

//   /**
//    *  Amplitude of the \f$\bar{K}^*_0(1680)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _a1K1680;

//   /**
//    *  Phase of the \f$\bar{K}^*_0(1680)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   double _phi1K1680;

//   /**
//    *  Amplitude of the non-resonant component for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _a2NR;

//   /**
//    *  Phase of the non-resonant component for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _phi2NR;

//   /**
//    *  Amplitude of the \f$\bar{K}^*(892)^0\f$  component for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _a2K8920;

//   /**
//    *  Phase of the \f$\bar{K}^*(892)^0\f$  component for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _phi2K8920;

//   /**
//    *  Amplitude of the \f$K^*(892)^-\f$ component  for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _a2K892m;

//   /**
//    *  Phase of the \f$K^*(892)^-\f$ component  for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _phi2K892m;

//   /**
//    *  Amplitude of the \f$\rho^+\f$ component  for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _a2rho;

//   /**
//    *  Phase of the \f$\rho^+\f$  component for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   double _phi2rho;

//   /**
//    *  Amplitude of the non-resonant component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   double _a3NR;

//   /**
//    *  Phase of the non-resonant component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   double _phi3NR;

//   /**
//    *  Amplitude of the  \f$K^*(892)^-\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   double _a3K892;

//   /**
//    *  Phase of the  \f$K^*(892)^-\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   double _phi3K892;

//   /**
//    *  Amplitude of the  \f$\rho^0\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   double _a3rho;

//   /**
//    *  Phase of the  \f$\rho^0\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   double _phi3rho;
//   //@}

//   /**
//    *  Complex amplitudes for use in the matrix element
//    */
//   //@{
//   /**
//    *  Amplitude of the non-resonant component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   Complex _c1NR;

//   /**
//    *  Amplitude of the \f$\bar{K}^*(892)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   Complex _c1K892;

//   /**
//    *  Amplitude of the \f$\bar{K}^*_0(1430)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   Complex _c1K1430;

//   /**
//    *  Amplitude of the \f$\bar{K}^*_0(1680)^0 component for \f$D^+\to K^-\pi^+\pi^+\f$
//    */
//   Complex _c1K1680;

//   /**
//    *  Amplitude of the non-resonant component for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   Complex _c2NR;

//   /**
//    *  Amplitude of the \f$\bar{K}^*(892)^0\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   Complex _c2K8920;

//   /**
//    *  Amplitude of the \f$K^*(892)^-\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   Complex _c2K892m;

//   /**
//    *  Amplitude of the \f$\rho^+\f$ for \f$D^0\to K^-\pi^+\pi^0\f$
//    */
//   Complex _c2rho;

//   /**
//    *  Amplitude of the non-resonant component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   Complex _c3NR;

//   /**
//    *  Amplitude of the  \f$K^*(892)^-\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   Complex _c3K892;

//   /**
//    *  Amplitude of the  \f$\rho^0\f$ component for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
//    */
//   Complex _c3rho;
//   //@}

//   /**
//    *  Masses and widths of the various resonances
//    */
//   //@{
//   /**
//    *  Use local values for the masses and widths
//    */
//   bool _localparameters;

//   /**
//    * Mass of the \f$K^*(892)^0\f$
//    */
//   Energy _mK8920;

//   /**
//    * Width of the \f$K^*(892)^0\f$
//    */
//   Energy _wK8920;

//   /**
//    * Mass of the \f$K^*(892)^-\f$
//    */
//   Energy _mK892m;

//   /**
//    * Width of the \f$K^*(892)^-\f$
//    */
//   Energy _wK892m;

//   /**
//    * Mass of the \f$K^*(1680)^0\f$
//    */
//   Energy _mK1680;

//   /**
//    * Width of the \f$K^*(1680)^0\f$
//    */
//   Energy _wK1680;

//   /**
//    * Mass of the \f$K^*_0(1430)^0\f$
//    */
//   Energy _mK1430;

//   /**
//    * Width of the \f$K^*_0(1430)^0\f$
//    */
//   Energy _wK1430;

//   /**
//    * Mass of the \f$\rho^0\f$
//    */
//   Energy _mrho0;

//   /**
//    * Width of the \f$\rho^0\f$
//    */
//   Energy _wrho0;

//   /**
//    * Mass of the \f$\rho^+\f$
//    */
//   Energy _mrhop;

//   /**
//    * Width of the \f$\rho^+\f$
//    */
//   Energy _wrhop;
//   //@}

//   /**
//    *  Parameters for the phase-space integration
//    */
//   //@{
//   /**
//    *  Maximum weights for the various modes
//    */
//   vector<double> _maxwgt;

//   /**
//    *  Weights for the different integration channels
//    */
//   vector<double> _weights;
//   //@}
}

int DtoKPiPiE691::modeNumber(bool & cc,const DecayMode & dm) const {
  int id0(dm.parent()->id());
  // incoming particle must be D0 or D+
  if(abs(id0)!=ParticleID::D0&&abs(id0)!=ParticleID::Dplus) return -1;
  cc = id0<0;
  // must be three decay products
  if(dm.products().size()!=3) return -1;
  ParticleMSet::const_iterator pit = dm.products().begin();
  unsigned int npip(0),npim(0),nkm(0),nk0(0),npi0(0);
  int id;
  for( ;pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id          ==ParticleID::piplus)  ++npip;
    else if(id     ==ParticleID::pi0)     ++npi0;
    else if(id     ==ParticleID::piminus) ++npim;
    else if(abs(id)==ParticleID::K0)      ++nk0;
    else if(id     ==ParticleID::K_L0)    ++nk0;
    else if(id     ==ParticleID::K_S0)    ++nk0;
    else if(abs(id)==ParticleID::Kplus)   ++nkm;
  }
  if(abs(id0)==ParticleID::Dplus) {
    if((id0==ParticleID::Dplus &&(nkm==1&&npip==2))||
       (id0==ParticleID::Dminus&&(nkm==1&&npim==2))) return 0;
    else return -1;
  }
  else {
    if(npim==1&&npip==1&&nk0==1) return  2;
    else if(nkm==1&&(npip+npim)==1&&npi0==1) return 1;
    else                        return -1;
  }
}

double DtoKPiPiE691::me2(bool vertex, const int ichan,
			 const Particle & inpart,
			 const ParticleVector & decay) const {
  useMe();
  // wavefunnction for the decaying particle
  tPPtr mytempInpart = const_ptr_cast<tPPtr>(&inpart);
  ScalarWaveFunction(mytempInpart,incoming,true,vertex);
  // wavefunctions for the outgoing particles
  for(unsigned int ix=0;ix<3;++ix) {
    PPtr mytemp = decay[ix]; 
    ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }
  Complex amp;
  // D+ -> K-pi+pi_
  if(imode()==0) {
    double ort(sqrt(0.5));
    Lorentz5Momentum pres1=decay[0]->momentum()+decay[1]->momentum();
    pres1.rescaleMass();
    double ct1 =-decayAngle(inpart.momentum(),pres1,decay[0]->momentum());
    Lorentz5Momentum pres2=decay[0]->momentum()+decay[2]->momentum();
    pres2.rescaleMass();
    double ct2 =-decayAngle(inpart.momentum(),pres2,decay[0]->momentum());
    if(ichan<0) {
      amp = _c1NR
	+ort*_c1K892 *amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920)
	+ort*_c1K892 *amplitude(1,ct2,pres2.mass(),_wK8920,_mK8920)
	+ort*_c1K1430*amplitude(0,ct1,pres1.mass(),_wK1430,_mK1430)
	+ort*_c1K1430*amplitude(0,ct2,pres2.mass(),_wK1430,_mK1430)
	+ort*_c1K1680*amplitude(1,ct1,pres1.mass(),_wK1680,_mK1680)
	+ort*_c1K1680*amplitude(1,ct2,pres2.mass(),_wK1680,_mK1680);
    }
    else if(ichan==0) {
      amp=ort*_c1K892 *amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920);
    }
    else if(ichan==1) {
      amp=ort*_c1K892 *amplitude(1,ct2,pres2.mass(),_wK8920,_mK8920);
    }
    else if(ichan==2) {
      amp=ort*_c1K1430*amplitude(1,ct1,pres1.mass(),_wK1430,_mK1430);
    }
    else if(ichan==3) {
      amp=ort*_c1K1430*amplitude(1,ct2,pres2.mass(),_wK1430,_mK1430);
    }
    else if(ichan==4) {
      amp=ort*_c1K1680*amplitude(1,ct1,pres1.mass(),_wK1680,_mK1680);
    }
    else if(ichan==5) {
      amp=ort*_c1K1680*amplitude(1,ct2,pres2.mass(),_wK1680,_mK1680);
    }
  }
  // D0 -> K-pi+pi0
  else if(imode()==1) {
    Lorentz5Momentum pres1=decay[0]->momentum()+decay[1]->momentum();
    pres1.rescaleMass();
    double ct1 = decayAngle(inpart.momentum(),pres1,decay[0]->momentum());
    Lorentz5Momentum pres2=decay[0]->momentum()+decay[2]->momentum();
    pres2.rescaleMass();
    double ct2 = decayAngle(inpart.momentum(),pres2,decay[2]->momentum());
    Lorentz5Momentum pres3=decay[1]->momentum()+decay[2]->momentum();
    pres3.rescaleMass();
    double ct3 = decayAngle(inpart.momentum(),pres3,decay[2]->momentum());
    if(ichan<0) {
      amp = _c2NR
	+_c2K8920*amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920)
	+_c2K892m*amplitude(1,ct2,pres2.mass(),_wK892m,_mK892m)
	+_c2rho  *amplitude(1,ct3,pres3.mass(),_wrhop ,_mrhop );
    }
    else if(ichan==0) {
      amp = _c2K8920*amplitude(1,ct1,pres1.mass(),_wK8920,_mK8920);
    }
    else if(ichan==1) {
      amp = _c2K892m*amplitude(1,ct2,pres2.mass(),_wK892m,_mK892m);
    }
    else if(ichan==2) {
      amp = _c2rho  *amplitude(1,ct3,pres3.mass(),_wrhop ,_mrhop);
    }
  }
  // D0 -> Kbar0pi+pi-
  else if(imode()==2) {
    Lorentz5Momentum pres1=decay[0]->momentum()+decay[2]->momentum();
    pres1.rescaleMass();
    double ct1 = decayAngle(inpart.momentum(),pres1,decay[0]->momentum());
    Lorentz5Momentum pres2=decay[1]->momentum()+decay[2]->momentum();
    pres2.rescaleMass();
    double ct2 = decayAngle(inpart.momentum(),pres2,decay[1]->momentum());
    if(ichan<0) {
      amp = _c3NR
	+_c3K892*amplitude(1,ct1,pres1.mass(),_wK892m,_mK892m)
	+_c3rho *amplitude(1,ct2,pres2.mass(),_wrho0 ,_mrho0 );
    }
  }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newME(0,0,0,0)=amp;
  ME(newME);
  return real(amp*conj(amp));
}

void DtoKPiPiE691::dataBaseOutput(ofstream & output, bool header) const {
}
