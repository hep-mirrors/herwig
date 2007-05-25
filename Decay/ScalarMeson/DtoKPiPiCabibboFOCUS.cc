// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DtoKPiPiCabibboFOCUS class.
//

#include "DtoKPiPiCabibboFOCUS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using ThePEG::Helicity::RhoDMatrix;
using Herwig::Helicity::ScalarWaveFunction;
using Herwig::Helicity::incoming;
using Herwig::Helicity::outgoing;

DtoKPiPiCabibboFOCUS::DtoKPiPiCabibboFOCUS() {
  // Masses and widths of the resonances
  _localparameters=true;
  _mK892    =  896.1*MeV; _wK892    =  50.7*MeV;
  _mK1410   = 1414. *MeV; _wK1410   = 232. *MeV;
  _mK14300  = 1412. *MeV; _wK14300  = 294. *MeV;
  _mK14302  = 1432.4*MeV; _wK14302  = 109. *MeV;
  _mrho770  =  768.5*MeV; _wrho770  = 150.7*MeV;
  _mf980    =  972. *MeV; _wf980    =  59. *MeV;
  _mrho1450 = 1465. *MeV; _wrho1450 = 310. *MeV;
  // Amplitudes and Phases for the D+ ->  K+pi-pi+
  _aDrho770 = 1.000     ; _phiDrho770 =    0.0;
  _aDK892   = 1.151     ; _phiDK892   = -167.1;
  _aDf980   = 0.476*GeV2; _phiDf980   = -134.5;
  _aDK1430  = 0.451/GeV2; _phiDK1430  =   54.4;
  // Amplitudes and Phases for the D_s+ -> K+pi-pi+
  _aDsNR      = 0.640; _phiDsNR      =   43.1;
  _aDsrho770  = 1.000; _phiDsrho770  =    0.0;
  _aDsK892    = 0.747; _phiDsK892    =  161.7;
  _aDsK1410   = 0.696; _phiDsK1410   = - 34.8;
  _aDsK1430   = 0.444*GeV2; _phiDsK1430   =   59.3;
  _aDsrho1450 = 0.523; _phiDsrho1450 = -151.7;
  // radial sizes
  _rD0  = 5.0/GeV;
  _rres = 1.5/GeV;
}

void DtoKPiPiCabibboFOCUS::persistentOutput(PersistentOStream & os) const {
  os << _localparameters << _mK892 << _wK892 << _mK1410 << _wK1410 << _mK14300 
     << _wK14300 << _mK14302 << _wK14302 << _mrho770 << _wrho770 << _mf980 
     << _wf980 << _mrho1450 << _wrho1450 << _aDrho770 << _phiDrho770 
     << _aDK892 << _phiDK892 << _aDf980 << _phiDf980 << _aDK1430 << _phiDK1430 
     << _cDrho770 << _cDK892 << _cDf980 << _cDK1430 << _aDsNR << _phiDsNR 
     << _aDsrho770 << _phiDsrho770 << _aDsK892 << _phiDsK892 << _aDsK1410 
     << _phiDsK1410 << _aDsK1430 << _phiDsK1430 << _aDsrho1450 << _phiDsrho1450 
     << _cDsNR << _cDsrho770 << _cDsK892 << _cDsK1410 << _cDsK1430 
     << _cDsrho1450 << _maxweight << _weights << _rD0 << _rres;
}

void DtoKPiPiCabibboFOCUS::persistentInput(PersistentIStream & is, int) {
  is >> _localparameters >> _mK892 >> _wK892 >> _mK1410 >> _wK1410 >> _mK14300 
     >> _wK14300 >> _mK14302 >> _wK14302 >> _mrho770 >> _wrho770 >> _mf980 
     >> _wf980 >> _mrho1450 >> _wrho1450 >> _aDrho770 >> _phiDrho770 
     >> _aDK892 >> _phiDK892 >> _aDf980 >> _phiDf980 >> _aDK1430 >> _phiDK1430 
     >> _cDrho770 >> _cDK892 >> _cDf980 >> _cDK1430 >> _aDsNR >> _phiDsNR 
     >> _aDsrho770 >> _phiDsrho770 >> _aDsK892 >> _phiDsK892 >> _aDsK1410 
     >> _phiDsK1410 >> _aDsK1430 >> _phiDsK1430 >> _aDsrho1450 >> _phiDsrho1450 
     >> _cDsNR >> _cDsrho770 >> _cDsK892 >> _cDsK1410 >> _cDsK1430 
     >> _cDsrho1450 >> _maxweight >> _weights >> _rD0 >> _rres;
}

ClassDescription<DtoKPiPiCabibboFOCUS> DtoKPiPiCabibboFOCUS::initDtoKPiPiCabibboFOCUS;
// Definition of the static class description member.

void DtoKPiPiCabibboFOCUS::Init() {

  static ClassDocumentation<DtoKPiPiCabibboFOCUS> documentation
    ("There is no documentation for the DtoKPiPiCabibboFOCUS class");

}

void DtoKPiPiCabibboFOCUS::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // complex amplitudes for D+ -> K+pi-pi+
  double fact = pi/180.;
  _cDrho770   = _aDrho770*Complex(cos(fact*_phiDrho770),sin(fact*_phiDrho770));
  _cDK892     = _aDK892  *Complex(cos(fact*_phiDK892  ),sin(fact*_phiDK892  ));
  _cDf980     = _aDf980  *Complex(cos(fact*_phiDf980  ),sin(fact*_phiDf980  ));
  _cDK1430    = _aDK1430 *Complex(cos(fact*_phiDK1430 ),sin(fact*_phiDK1430 ));
  // complex amplitudes for D_s+ -> K+pi-pi+
  _cDsNR      = _aDsNR     *Complex(cos(fact*_phiDsNR     ),sin(fact*_phiDsNR     ));
  _cDsrho770  = _aDsrho770 *Complex(cos(fact*_phiDsrho770 ),sin(fact*_phiDsrho770 ));
  _cDsK892    = _aDsK892   *Complex(cos(fact*_phiDsK892   ),sin(fact*_phiDsK892   ));
  _cDsK1410   = _aDsK1410  *Complex(cos(fact*_phiDsK1410  ),sin(fact*_phiDsK1410  ));
  _cDsK1430   = _aDsK1430  *Complex(cos(fact*_phiDsK1430  ),sin(fact*_phiDsK1430  ));
  _cDsrho1450 = _aDsrho1450*Complex(cos(fact*_phiDsrho1450),sin(fact*_phiDsrho1450));
  // resonances for the phase-space modes
  tPDPtr k892    = getParticleData(ParticleID::Kstar0);
  tPDPtr k1410   = getParticleData(100313);
  tPDPtr k14300  = getParticleData(ParticleID::Kstar_00);
  tPDPtr k14302  = getParticleData(ParticleID::Kstar_20);
  tPDPtr rho770  = getParticleData(ParticleID::rho0);
  tPDPtr f980    = getParticleData(ParticleID::f_0);
  tPDPtr rho1450 = getParticleData(100113);
  // D+ -> K+ pi- pi+
  PDVector extpart(4);
  extpart[0]=getParticleData(ParticleID::Dplus  );
  extpart[1]=getParticleData(ParticleID::Kplus  );
  extpart[2]=getParticleData(ParticleID::piminus);
  extpart[3]=getParticleData(ParticleID::piplus );
  DecayPhaseSpaceModePtr mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  DecayPhaseSpaceChannelPtr newchannel;
  int ix=0;
  if(rho770) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho770,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k892) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k892,0,0., 1,2);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(f980) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(f980,0,0., 2,3);
    mode->addChannel(newchannel);
    ++ix;
  }
  if(k14302) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k14302,0,0., 1,2);
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
  if(_maxweight.empty()) _maxweight.push_back(1.);
  addMode(mode,_maxweight[0],wtemp);
  // D_s -> K+ pi- pi+ mode
  extpart[0]=getParticleData(ParticleID::D_splus);
  mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
  int iy=ix;
  if(rho770) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho770,0,0., 2,3);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k892) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k892,0,0., 1,2);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k1410) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k1410,0,0., 1,2);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(k14300) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(k14300,0,0., 1,2);
    mode->addChannel(newchannel);
    ++iy;
  }
  if(rho1450) {
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,1);
    newchannel->addIntermediate(rho1450,0,0., 2,3);
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
    for(unsigned int iz=0;iz<wtemp.size();++iz) {
    }
  }
  if(_maxweight.size()<2) _maxweight.push_back(1.);
  addMode(mode,_maxweight[1],wtemp);
  // reset the masses if needed
  if(!_localparameters) {
    _mK892    = k892   ->mass();
    _wK892    = k892   ->width();
    _mK1410   = k1410  ->mass();
    _wK1410   = k1410  ->width();
    _mK14300  = k14300 ->mass();
    _wK14300  = k14300 ->width();
    _mK14302  = k14302 ->mass();
    _wK14302  = k14302 ->width();
    _mrho770  = rho770 ->mass();
    _wrho770  = rho770 ->width();
    _mf980    = f980   ->mass();
    _wf980    = f980   ->width();
    _mrho1450 = rho1450->mass();
    _wrho1450 = rho1450->width();
  }
  else {
    resetIntermediate(k892   ,_mK892   ,_wK892   );
    resetIntermediate(k1410  ,_mK1410  ,_wK1410  );
    resetIntermediate(k14300 ,_mK14300 ,_wK14300 );
    resetIntermediate(k14302 ,_mK14302 ,_wK14302 );
    resetIntermediate(rho770 ,_mrho770 ,_wrho770 );
    resetIntermediate(f980   ,_mf980   ,_wf980   );
    resetIntermediate(rho1450,_mrho1450,_wrho1450);
  } 
}

int DtoKPiPiCabibboFOCUS::modeNumber(bool & cc,const DecayMode & dm) const {
  int id0(dm.parent()->id());
  // incoming particle must be D+/- or D_s+/-
  if(abs(id0)!=ParticleID::Dplus  &&
     abs(id0)!=ParticleID::D_splus) return -1;
  if(dm.products().size()!=3) return -1;
  ParticleMSet::const_iterator pit = dm.products().begin();
  cc = id0<0;
  int isign = cc ? -1 : 1;
  unsigned int npip(0),npim(0),nkp(0);
  for( ;pit!=dm.products().end();++pit) {
    int id=isign*(**pit).id();
    if(     id==ParticleID::piplus ) ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::Kplus  ) ++nkp;
  }
  if(!(npip==1&&npim==1&&nkp==1)) return -1;
  return abs(id0)==ParticleID::Dplus ? 0 : 1;
}

double DtoKPiPiCabibboFOCUS::me2(bool vertex, const int ichan,const Particle & inpart,
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
  // masses of the particles
  Energy mD  = inpart.mass();
  Energy mA  = decay[0]->mass();
  Energy mB  = decay[1]->mass();
  Energy mC  = decay[2]->mass();
  // resonances and angles
  double ct1,ct2;
  Energy2 E1,E2;
  Lorentz5Momentum pres1=decay[1]->momentum()+decay[2]->momentum();
  Lorentz5Momentum pres2=decay[1]->momentum()+decay[0]->momentum();
  pres1.rescaleMass();
  pres2.rescaleMass();
  decayAngle(inpart.momentum(),pres1,decay[2]->momentum(),ct1,E1);
  decayAngle(inpart.momentum(),pres2,decay[0]->momentum(),ct2,E2);
  // matrix element
  Complex amp(0.);
  if(imode()==0) {
    amp = 
      _cDrho770           *amplitude(1,mD,mB,mC,mA,pres1.mass(),_mrho770,_wrho770,E1,ct1)
      +_cDK892 *0.975     *amplitude(1,mD,mA,mB,mC,pres2.mass(),_mK892  ,_wK892  ,E2,ct2)
      +_cDf980 *0.124/GeV2*amplitude(0,mD,mB,mC,mA,pres1.mass(),_mf980  ,_wf980  ,E1,ct1)
      +_cDK1430*32.3 *GeV2*amplitude(2,mD,mA,mB,mC,pres2.mass(),_mK14302,_wK14302,E2,ct2)
      ;
  }
  else {
    amp =
      _cDsNR*0.707
      +_cDsrho770      *amplitude(1,mD,mB,mC,mA,pres1.mass(),_mrho770 ,_wrho770 ,E1,ct1)
      +_cDsK892  *0.981*amplitude(1,mD,mA,mB,mC,pres2.mass(),_mK892   ,_wK892   ,E2,ct2)
      +_cDsK1410 *3.03 *amplitude(1,mD,mA,mB,mC,pres2.mass(),_mK1410  ,_wK1410  ,E2,ct2)
      +_cDsK1430*0.467  /GeV2*amplitude(0,mD,mA,mB,mC,pres2.mass(),_mK14300 ,_wK14300 ,E2,ct2)
      +_cDsrho1450*4.59*amplitude(1,mD,mB,mC,mA,pres1.mass(),_mrho1450,_wrho1450,E1,ct1);
  }
  // now compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0);
  newME(0,0,0,0)=amp;
  ME(newME);
  return real(amp*conj(amp));
}

void DtoKPiPiCabibboFOCUS::dataBaseOutput(ofstream & os,bool header) const {
}
