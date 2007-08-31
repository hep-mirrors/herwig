// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNFSVertex class.
//

#include "SSNFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNFSVertex::SSNFSVertex():_tanB(0.), _sw(0.), _cw(0.), _mw(),
			   _sb(0.), _cb(0.), _leftlast(0.), 
			   _rightlast(0), _id1last(0), _id2last(0) {
  vector<int> first,second,third;
  for(unsigned int ineu=0;ineu<5;++ineu) {
    int neu(-1);
    if(ineu==0){neu=1000022;}
    if(ineu==1){neu=1000023;}
    if(ineu==2){neu=1000025;}
    if(ineu==3){neu=1000035;}
    if(ineu==4){neu=1000045;}
    for(unsigned int ix=1;ix<7;++ix){
      first.push_back(neu);
      second.push_back(ix);
      third.push_back(-(1000000+ix));
    }
    for(unsigned int ix=1;ix<7;++ix){
      first.push_back(neu);
      second.push_back(ix);
      third.push_back(-(2000000+ix));
    }
    for(unsigned int ix=1;ix<7;++ix){
      first.push_back(neu);
      second.push_back(-ix);
      third.push_back((1000000+ix));
    }
    for(unsigned int ix=1;ix<7;++ix){
      first.push_back(neu);
      second.push_back(-ix);
      third.push_back((2000000+ix));
    }
    for(unsigned int ix=11;ix<17;++ix){
      first.push_back(neu);
      second.push_back(ix);
      third.push_back(-(1000000+ix));
    }
    for(unsigned int ix=11;ix<17;ix += 2){
      first.push_back(neu);
      second.push_back(ix);
      third.push_back(-(2000000+ix));
    }
    for(unsigned int ix=11;ix<17;++ix){
      first.push_back(neu);
      second.push_back(-ix);
      third.push_back((1000000+ix));
    }
    for(unsigned int ix=11;ix<17;ix +=2){
      first.push_back(neu);
      second.push_back(-ix);
      third.push_back((2000000+ix));
    }
  }
  setList(first,second,third);
}

void SSNFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theStop << _theSbottom << _theStau << _theN << _theSS 
     << _tanB << _sw << _cw << ounit(_mw,GeV) << _sb << _cb;
}

void SSNFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theStop >> _theSbottom >> _theStau >> _theN >> _theSS 
     >> _tanB >> _sw >> _cw >> iunit(_mw,GeV) >> _sb >> _cb;
  _leftlast = 0.;
  _rightlast = 0.;
  _id1last = 0;
  _id2last = 0;
  _id3last = 0;
}

void SSNFSVertex::doinit() throw(InitException) {
  FFSVertex::doinit();
  _theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!_theSS)
    throw InitException() << "SSGSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;

  _theStop = _theSS->stopMix();
  _theSbottom = _theSS->sbottomMix();
  _theStau = _theSS->stauMix();
  _theN = _theSS->neutralinoMix();
  if(!_theStop || !_theStau || !_theSbottom)
    throw InitException() << "SSNFSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _theStop << " sbottom: "
			  << _theSbottom << " stau: " << _theStau 
			  << " N: " << _theN << Exception::abortnow;
  _sw = sqrt(_theSS->sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _tanB = _theSS->tanBeta();
  _cw = sqrt(1. - _sw*_sw);
  _sb = _tanB/sqrt(1 + _tanB*_tanB);
  _cb = sqrt(1 - _sb*_sb);
  
  orderInGem(1);
  orderInGs(0);
}

ClassDescription<SSNFSVertex> SSNFSVertex::initSSNFSVertex;
// Definition of the static class description member.

void SSNFSVertex::Init() {

  static ClassDocumentation<SSNFSVertex> documentation
    ("The SSNFSVertex implements the coupling of a neutralino to "
     "a fermion-sfermion");

}

void SSNFSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3, int iinc) {
  long isc(0), iferm(0), ineut(0);
  if(part1->id() == 1000022 || part1->id() == 1000023 ||
     part1->id() == 1000025 || part1->id() == 1000035 || 
     part1->id() == 1000045) {
    ineut = part1->id();
    if(part2->iSpin() == PDT::Spin1Half) {
      iferm = part2->id();
      isc = part3->id();
    }
    else {
      iferm = part3->id();
      isc = part2->id();
    }
  }
  else if(part2->id() == 1000022 || part2->id() == 1000023 ||
	  part2->id() == 1000025 || part2->id() == 1000035 || 
	  part2->id() == 1000045) {
    ineut = part2->id();
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = part1->id();
      isc = part3->id();
    }
    else {
      iferm = part3->id();
      isc = part1->id();
    } 
  }
  else if(part3->id() == 1000022 || part3->id() == 1000023 ||
	  part3->id() == 1000025 || part3->id() == 1000035 || 
	  part3->id() == 1000045){
    ineut = part3->id();
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = part1->id();
      isc = part2->id();
    }
    else {
      iferm = part2->id();
      isc = part1->id();
    }
  }
  else {
    throw HelicityConsistencyError() 
      << "SSNFSVertex::setCoupling() - There is no neutralino in this vertex!"
      << part1->id() << " " << part2->id() << " " << part3->id()
      << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
  }
  if((abs(iferm) >=1 && abs(iferm) <= 6) || 
     (abs(iferm) >= 11 && abs(iferm) <= 16)) {
    if(abs(ineut) != _id1last || abs(iferm) != _id2last ||
       abs(isc) != _id3last) {
      _id1last = abs(ineut);
      _id2last = abs(iferm);
      _id3last = abs(isc);
      int neu(-1);
      if(ineut == 1000022) neu = 0;
      if(ineut == 1000023) neu = 1;
      if(ineut == 1000025) neu = 2;
      if(ineut == 1000035) neu = 3;
      if(ineut == 1000045) neu = 3;
      double gE = sqrt(4.*Constants::pi*_theSS->alphaEM(q2));
      Complex n1prime = (*_theN)(neu,0)*_cw + (*_theN)(neu,1)*_sw;
      Complex n2prime = (*_theN)(neu,1)*_cw - (*_theN)(neu,0)*_sw;
      Energy fmass = getParticleData(abs(iferm))->mass();
      Complex fact1 = double(gE*fmass/2./_mw/_sw);
      Complex fact2 = (0.5 - _sw*_sw);
      double eu = _theSS->eu();
      Complex fact3 = (0.5 - eu*_sw*_sw);
      double ed = _theSS->ed(); 
      Complex fact4 = (0.5 + ed*_sw*_sw);
      unsigned int eig = (abs(isc)/1000000) - 1;
      if(abs(isc) == 1000012 || abs(isc) == 1000014 || abs(isc) == 1000016) {
	_rightlast = gE*n2prime/2./_cw/_sw;
	_leftlast=0.;
      }
      else if(abs(isc) == 1000011 || abs(isc) == 1000013 || 
	      abs(isc) == 1000015 || abs(isc) == 2000011 ||
	      abs(isc) == 2000013 || abs(isc) == 2000015 ) {
	Complex l1 = fact1*(*_theN)(neu,2)/_cb;
	Complex l2 = gE*(n1prime - (_sw*n2prime/_cw));
	Complex l3 = gE*(n1prime + (n2prime*fact2/_cw/_sw));
	if(abs(iferm) == 15) {
	  _leftlast = (conj(l1)*(*_theStau)(eig,0)) + (*_theStau)(eig,1)*conj(l2);
	  _rightlast = (l1*(*_theStau)(eig,1)) - (*_theStau)(eig,0)*l3;
	}
	else {
	  if(eig == 0) {
	    _leftlast = conj(l1);
	    _rightlast = -l3;
	  }
	  else {
	    _leftlast = conj(l2);
	    _rightlast = l1;
	  }
	}
      }
      else if(abs(isc) == 1000002 || abs(isc) == 1000004 || 
	      abs(isc) == 1000006 || abs(isc) == 2000002 || 
	      abs(isc) == 2000004 || abs(isc) == 2000006 ) {
	Complex u1 = fact1*(*_theN)(neu,3)/_sb;
	Complex u2 = gE*(eu*n1prime + (n2prime*fact3/_sw/_cw)); 
	Complex u3 = gE*eu*(n1prime - (_sw*n2prime/_cw));
	if(abs(iferm) == 6) {
	  _leftlast = (conj(u1)*(*_theStop)(eig,0)) - 
	    ((*_theStop)(eig,1)*conj(u3));
	  _rightlast = (u1*(*_theStop)(eig,1)) + ((*_theStop)(eig,0)*u3);
	}
	else {
	  if(eig == 0) {
	    _leftlast = conj(u1);
	    _rightlast = u2;
	  }
	  else {
	    _leftlast = -conj(u3);
	    _rightlast = u1;
	  }
	}
      }
      else {
	Complex d1 = fact1*(*_theN)(neu,2)/_cb;
	Complex d2 = gE*(ed*n1prime - (n2prime*fact4/_sw/_cw));
	Complex d3 = gE*(ed*n1prime - (ed*_sw*n2prime/_cw) );
	if(abs(iferm) == 5) {
	  _leftlast = (conj(d1)*(*_theSbottom)(eig,0)) - 
	    ((*_theSbottom)(eig,1)*conj(d3));
	  _rightlast = (d1*(*_theSbottom)(eig,1)) + ((*_theSbottom)(eig,0)*d2);
	}
	else {
	  if(eig == 0) {
	    _leftlast = conj(d1);
	    _rightlast = d2;
	  }
	  else {
	    _leftlast = -conj(d3);
	    _rightlast = d1;
	    
	  }
	}
      }
    }
    setNorm(-sqrt(2));
    tcPDPtr incoming, fermion;
    //find incoming
    if(iinc == 1) {
      incoming = part1;
      if(part2->iSpin() == PDT::Spin1Half) fermion = part2;
      else fermion = part3;
    }
    else if(iinc == 2) {
      incoming = part2;
      if(part1->iSpin() == PDT::Spin1Half) fermion = part1;
      else fermion = part3;
    } 
    else {
      incoming = part3;
      if(part1->iSpin() == PDT::Spin1Half) fermion = part1;
      else fermion = part2;
    }
    //determine whether to flip couplings 
    if(incoming->iSpin() == PDT::Spin0) {
      if(incoming->id() > 0) {
	setLeft(_leftlast);
	setRight(_rightlast);
      }
      else {
	setLeft(conj(_rightlast));
	setRight(conj(_leftlast));
      }
    }
    else if(incoming->iSpin() == PDT::Spin1Half && incoming->CC()) {
      if(incoming->id() > 0) {
	setLeft(conj(_rightlast));
	setRight(conj(_leftlast));
      }
      else {
	setLeft(_leftlast);
	setRight(_rightlast);
      }
    }
    else {
      if(fermion->id() < 0) {
	setLeft(conj(_rightlast));
	setRight(conj(_leftlast));
      }
      else {
	setLeft(_leftlast);
	setRight(_rightlast);
      }
    }
  }
  else {
    throw HelicityConsistencyError() 
      << "SSNFSVertex::setCoupling() - Incorrect particle found in vertex. "
      << part1->id() << " " << part2->id() << " " << part3->id()
      << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
  }
}
