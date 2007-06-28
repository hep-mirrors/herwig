// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSH2SFSFVertex class.
//

#include "SSH2SFSFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

SSH2SFSFVertex::SSH2SFSFVertex():_sinAlpha(0.),_cosAlpha(0.),_sb(0.),
				 _cb(0.),_mz(0.),_mw(0.),_sw(0.),
				 _cw(0.),_mu(0.),_trilin(9,0.),
				 _cosAB(0.),_glast(0.),_q2last(0.),
				 _hfact(0.),_id1last(0),_id2last(0) {
  vector<int> first,second,third;
  for(unsigned int ix=1000001;ix<1000007;++ix) {
    first.push_back(35);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(unsigned int ix=2000001;ix<2000007;++ix) {
    first.push_back(35);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(unsigned int ix=1000001;ix<1000007;++ix) {
    first.push_back(35);
    second.push_back(ix);
    third.push_back(-(ix+1000000));
  }
  for(unsigned int ix=1000001;ix<1000007;++ix) {
    first.push_back(35);
    second.push_back(-ix);
    third.push_back(ix+1000000);
  }
  for(unsigned int ix=1000011;ix<1000017;++ix) {
    first.push_back(35);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(unsigned int ix=2000011;ix<2000016;ix+=2) {
    first.push_back(35);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(unsigned int ix=1000011;ix<1000016;ix+=2) {
    first.push_back(35);
    second.push_back(ix);
    third.push_back(-(ix+1000000));
  }
  for(unsigned int ix=1000011;ix<1000016;ix+=2) {
    first.push_back(35);
    second.push_back(-ix);
    third.push_back(ix+1000000);
  }
  setList(first,second,third);
}

void SSH2SFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sinAlpha << _cosAlpha << _sb << _cb << _mz 
     << _mw << _sw << _cw << _mu << _trilin <<  _cosAB << _stop << _sbottom << _stau;
}

void SSH2SFSFVertex::persistentInput(PersistentIStream & is, int) {
 is >> _theSS >> _sinAlpha >> _cosAlpha >> _sb >> _cb >> _mz 
    >> _mw >> _sw >> _cw >> _mu >> _trilin >> _cosAB >> _stop >> _sbottom 
    >> _stau;
  _glast=0.;
  _q2last=0.;
  _hfact=0.;
  _id1last=0;
  _id2last=0;
}

void SSH2SFSFVertex::doinit() throw(InitException) {
  SSSVertex::doinit();

  _theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  _sinAlpha = sin(_theSS->higgsMixingAngle());
  _cosAlpha = sqrt(1. - _sinAlpha*_sinAlpha);
  double tb = _theSS->tanBeta();
  _sb = tb/sqrt(1.+tb*tb);
  _cb = sqrt(1.-_cb*_cb);
  _mz = getParticleData(23)->mass();
  _mw = getParticleData(24)->mass();
  _sw = _theSS->sin2ThetaW();
  _cw = sqrt(1. - _sw*_sw);
  _mu = _theSS->muParameter();
  _cosAB = _cosAlpha*_cb - _sinAlpha*_sb;
  
  _trilin[0] = 0.*GeV; 
  _trilin[1] = 0.*GeV;
  _trilin[2] = 0.*GeV; 
  _trilin[3] = 0.*GeV;
  _trilin[4] = _theSS->bottomTrilinear().real();
  _trilin[5] = _theSS->topTrilinear().real();
  _trilin[6] = 0.*GeV;
  _trilin[7] = 0.*GeV;
  _trilin[8] = _theSS->tauTrilinear().real();
 
  _stop = _theSS->stopMix();
  _sbottom = _theSS->sbottomMix();
  _stau = _theSS->stauMix();

  if(!_stop || !_stau || !_sbottom)
    throw InitException() << "SSNFSVertex::doinit -  "
			  << "A null mixing matrix pointer. stop: " << _stop 
			  <<  " sbottom: " << _sbottom << " stau: " 
			  << _stau << Exception::abortnow;

  orderInGem(1);
  orderInGs(0);
}

ClassDescription<SSH2SFSFVertex> SSH2SFSFVertex::initSSH2SFSFVertex;
// Definition of the static class description member.

void SSH2SFSFVertex::Init() {

  static ClassDocumentation<SSH2SFSFVertex> documentation
    ("This is the implementation of the coupling of the H_2 neutral "
     "higgs of the MSSM to the sfermion pairs");

}

void SSH2SFSFVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				 tcPDPtr part2,tcPDPtr part3) {
  long isc1(0),isc2(0);
  if(abs(part1->id()) == 35 ) {
    isc1 = abs(part2->id());
    isc2 = abs(part3->id());
  }
  else if(abs(part2->id()) == 35 ) {
    isc1 = abs(part1->id());
    isc2 = abs(part3->id());
  }
  else if(abs(part3->id()) == 35 ) {
    isc1 = abs(part1->id());
    isc2 = abs(part2->id());
  }
  else {
    throw HelicityConsistencyError() << "SSH2SSVertex::setCoupling "
				     << "No H0 particle found!\n"
				     << Exception::warning;
    setNorm(0.);
  }
  if((isc1>=1000001 && isc1<=1000006)||(isc1>=2000001 && isc1<=2000006)||
     (isc1>=1000011 && isc1<=1000016)||(isc1>=2000011 && isc1<=2000016)) {
    if(q2 != _q2last) {
      double alpha = _theSS->alphaEM(q2);
      _glast = sqrt(4.*Constants::pi*alpha)/_sw;
      _q2last = q2;
    }
    if(isc1 != _id1last && isc2!=_id2last) {
      _id1last=isc1;
      _id2last=isc2;
      unsigned int eig1(isc1/1000000 -1),eig2(isc2/1000000 -1);
      unsigned int imass(0);
      if(eig1==0) {imass=isc1-1000000;}
      else {imass=isc1-2000000;}
      if(isc1==1000001||isc1==1000003||isc1==1000005||
	 isc1==2000001||isc1==2000003||isc1==2000005) {
	Energy dmass = getParticleData(imass)->mass();
	Energy d1 = _mz*_cosAB/_cw;
	Energy d2 = dmass*dmass*_cosAlpha/_mw/_cb;
	double d3 = dmass/2./_mw/_cb;
	double d4 = (_theSS->ed())*_sw*_sw;
	if(isc1==1000005||isc1==2000005) {
	  _hfact = d1*( (*_sbottom)(0,eig1)*(*_sbottom)(0,eig2)*(0.5+d4)
			 - d4*(*_sbottom)(1,eig1)*(*_sbottom)(1,eig2) );
	  _hfact -= d2*( (*_sbottom)(0,eig1)*(*_sbottom)(0,eig2) + 
			 (*_sbottom)(1,eig1)*(*_sbottom)(1,eig2) );
	  _hfact += d3*(_trilin[4]*_cosAlpha + _mu*_sinAlpha)*
	    ( (*_sbottom)(1,eig1)*(*_sbottom)(0,eig2) + 
	      (*_sbottom)(0,eig1)*(*_sbottom)(1,eig2) );
	}
	else {
	  if(eig1==0 && eig2==0) {
	    _hfact = d1*(0.5+d4) - d2;
	  }
	  else if(eig1==1 && eig2==1) {
	    _hfact = d1*d4 - d2;
	  }
	  else {
	    unsigned int tri;
	    if(eig1==1) {tri=isc1-1000001;}
	    else {tri=isc1-2000001;}
	    _hfact = d3*(_mu*_sinAlpha - _trilin[tri]*_cosAlpha);
	  }
	}
      }
      else if(isc1==1000002||isc1==1000004||isc1==1000006||
	      isc1==2000002||isc1==2000004||isc1==2000006) {
	unsigned int imass = isc1-1000000;
	Energy umass = getParticleData(imass)->mass();
	Energy u1 = _mz*_cosAB/_cw;
	Energy u2 = umass*umass*_sinAlpha/_mw/_sb;
	double u3 = umass/2./_mw/_sb;
	double u4 = (_theSS->eu()*_sw*_sw);
	if(isc1==1000006||isc1==2000006) {	
	  _hfact = -u1*( (*_stop)(0,eig1)*(*_stop)(0,eig2)*(0.5 - u4) 
			+ u4*(*_stop)(1,eig1)*(*_stop)(1,eig2));
	  
	  _hfact -= u2*( (*_stop)(0,eig1)*(*_stop)(0,eig2) 
			 + (*_stop)(1,eig1)*(*_stop)(1,eig2) );
	  
	  _hfact -= u3*( (*_stop)(1,eig1)*(*_stop)(0,eig2) 
			 + (*_stop)(0,eig1)*(*_stop)(1,eig2) )*
	  (_trilin[5]*_sinAlpha - _mu*_cosAlpha);
	}
	else {
	  if(eig1==0 && eig2==0) {
	    _hfact = -u1*(0.5 - u4) - u2;
	  }
	  else if(eig1==1 && eig2==1){
	    _hfact = -u1*u4 - u2;
	  }
	  else {
	    unsigned int tri;
	    if(eig1==1) {tri=isc1-1000001;}
	    else {tri=isc1-2000001;}
	    _hfact = -u3*(_trilin[tri]*_sinAlpha - _mu*_cosAlpha);
	  }
	}
      }
      else if(isc1==1000011||isc1==1000013||isc1==1000015||
	      isc1==2000011||isc1==2000013||isc1==2000015){
	unsigned int imass = isc1-1000000;
	Energy l1 = _mz*_cosAB/_cw;
	Energy lmass = getParticleData(imass)->mass();
	Energy l2 = lmass*lmass*_cosAlpha/_mw/_cb;
	double l3 = lmass/2./_mw/_cb;
	if(isc1==1000015||isc1==2000015) {
	  _hfact = l1*( (*_stau)(0,eig1)*(*_stau)(0,eig2)*(0.5-_sw*_sw) 
			 + _sw*_sw*(*_stau)(1,eig1)*(*_stau)(1,eig2)  );
	  _hfact -= l2*( (*_stau)(0,eig1)*(*_stau)(0,eig2) 
			 + (*_stau)(1,eig1)*(*_stau)(1,eig2) );
	  _hfact += l3*(_mu*_sinAlpha - _trilin[8]*_cosAlpha)*
	    ( (*_stau)(0,eig1)*(*_stau)(0,eig2)
	      +(*_stau)(1,eig1)*(*_stau)(1,eig2) );
	}
	else {
	  if(eig1==0 && eig2==0){
	    _hfact = l1*(0.5 - _sw*_sw) - l2;
	  }
	  else if(eig1==1 && eig2==1) {
	    _hfact = l1*_sw*_sw - l2;
	  }
	  else {
	    unsigned int tri;
	    if(eig1==0) {tri=(isc1-999999)/2;}
	    else {tri=(isc1-1999999)/2;}
	    _hfact = -l3*(_mu*_sinAlpha - _trilin[tri]*_cosAlpha );
	  }
	}
      }
      else {
	_hfact = -_mz*_cosAB/2./_cw;
      }
    }
    setNorm(_glast*_hfact);
  }
  else {
    throw HelicityConsistencyError() << "Incorrect particle found "
				     << "in SSH2SFSFVertex."
				     << Exception::warning;
    setNorm(0.);
  }
  
}

