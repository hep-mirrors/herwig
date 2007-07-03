// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCFSVertex class.
//

#include "SSCFSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig::Helicity;

SSCFSVertex::SSCFSVertex(): _tanB(0.),_sb(0.),_cb(0.),_mw(0.*MeV),
			    _sw(0.),_q2last(0.*sqr(MeV)),_couplast(0.),
			    _leftlast(0.),_rightlast(0.),
			    _id1last(0), _id2last(0), _id3last(0) {
  vector<int> first,second,third;
  for(unsigned int ic=0;ic<2;++ic) {
    int chargino(-1);
    if(ic==0) {chargino=1000024;}
    else {chargino=1000037;}
    //outgoing chi-
    for(unsigned int ix=2;ix<7;ix+=2) {
      first.push_back(-chargino);
      second.push_back(ix);
      third.push_back(-(999999+ix));
    }
    for(unsigned int ix=2;ix<7;ix+=2) {
      first.push_back(-chargino);
      second.push_back(ix);
      third.push_back(-(1999999+ix));
    }
    for(unsigned int ix=1;ix<6;ix+=2) {
      first.push_back(-chargino);
      second.push_back(-ix);
      third.push_back((1000001+ix));
    }
    for(unsigned int ix=1;ix<6;ix+=2) {
      first.push_back(-chargino);
      second.push_back(-ix);
      third.push_back(2000001+ix);
    }
    for(unsigned int ix=12;ix<17;ix+=2) {
      first.push_back(-chargino);
      second.push_back(ix);
      third.push_back(-(999999+ix));
    }
    for(unsigned int ix=12;ix<17;ix+=2) {
      first.push_back(-chargino);
      second.push_back(ix);
      third.push_back(-(1999999+ix));
    }
    for(unsigned int ix=11;ix<16;ix+=2) {
      first.push_back(-chargino);
      second.push_back(-ix);
      third.push_back(1000001+ix);
    }
    //outgoing chi+
    for(unsigned int ix=2;ix<7;ix+=2) {
      first.push_back(chargino);
      second.push_back(-ix);
      third.push_back((999999+ix));
    }
    for(unsigned int ix=2;ix<7;ix+=2) {
      first.push_back(chargino);
      second.push_back(-ix);
      third.push_back((1999999+ix));
    }
    for(unsigned int ix=1;ix<6;ix+=2) {
      first.push_back(chargino);
      second.push_back(ix);
      third.push_back(-(1000001+ix));
    }
    for(unsigned int ix=1;ix<6;ix+=2) {
      first.push_back(chargino);
      second.push_back(ix);
      third.push_back(-(2000001+ix));
    }
    for(unsigned int ix=12;ix<17;ix+=2) {
      first.push_back(chargino);
      second.push_back(-ix);
      third.push_back((999999+ix));
    }
    for(unsigned int ix=12;ix<17;ix+=2) {
      first.push_back(chargino);
      second.push_back(-ix);
      third.push_back((1999999+ix));
    }
    for(unsigned int ix=11;ix<16;ix+=2) {
      first.push_back(chargino);
      second.push_back(ix);
      third.push_back(-(1000001+ix));
    }
  }
   setList(first,second,third);
}

void SSCFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _tanB << _sb << _cb << ounit(_mw,GeV) << _sw << _stop <<
    _sbottom << _stau << _chargU << _chargV;
}

void SSCFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _tanB >> _sb >> _cb >> iunit(_mw,GeV) >> _sw >> _stop >>
    _sbottom >> _stau >> _chargU >> _chargV;
  _q2last=0.*sqr(MeV);
  _couplast=0.;
  _leftlast=0.;
  _rightlast=0.;
  _id1last=0;
  _id2last=0;
  _id3last = 0;
}


ClassDescription<SSCFSVertex> SSCFSVertex::initSSCFSVertex;
// Definition of the static class description member.

void SSCFSVertex::Init() {

  static ClassDocumentation<SSCFSVertex> documentation
    ("The implementation of the coupling of the charginos to fermion-"
     "sfermions.");
}

void SSCFSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3, int iinc) {
  long iferm(0),icharg(0),isc(0);
  if(abs(part1->id()) == 1000024 || abs(part1->id()) == 1000037) {
    icharg = abs(part1->id());
    if(part2->iSpin() == PDT::Spin1Half) {
      iferm = part2->id();
      isc = part3->id();
    }
    else {
      iferm = part3->id();
      isc = part2->id();
    }
  }
  else if(abs(part2->id()) == 1000024 || abs(part2->id()) == 1000037) {
    icharg = abs(part2->id());
    if(part1->iSpin() == PDT::Spin1Half) {
      iferm = part1->id();
      isc = part3->id();
    }
    else {
      iferm = part3->id();
      isc = part1->id();
    } 
  }
  else if(abs(part3->id()) == 1000024 || abs(part3->id()) == 1000037){
    icharg = abs(part3->id());
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
    throw HelicityConsistencyError() << "SSNFSVertex::setCoupling() - "
				     << "There is no chargino in "
				     << "SSNFSVertex!\n"
				     << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
  }
  if((abs(iferm)>=1 && abs(iferm)<=6)||(abs(iferm)>=11 && abs(iferm) <=16)){
    if(q2 != _q2last){
      double gEw = sqrt(4.*Constants::pi*(_theSS->alphaEM(q2)))/_sw;
      _couplast=-gEw;
      _q2last = q2;
    }
    if(abs(iferm) != _id1last || abs(isc) != _id2last || 
       icharg != _id3last) {
      _id1last = abs(iferm);
      _id2last = abs(isc);
      _id3last = icharg;
      unsigned int chargino;
      if(icharg==1000024)chargino=0;
      else chargino=1;
      unsigned int eig = (abs(isc)/1000000) - 1;
      
      if(abs(iferm)==12||abs(iferm)==14||abs(iferm)==16) {
	Energy lmass = getParticleData(abs(iferm)-1)->mass();
	_leftlast = 0.;
	Complex l1 = lmass*(*_chargU)(chargino,1)/sqrt(2)/_mw/_cb;
	if(abs(iferm)==16) {
	  _rightlast = (*_chargU)(chargino,0)*(*_stau)(0,eig)
	    - l1*(*_stau)(1,eig);
	}
	else {
	  if(eig==0) {
	    _rightlast = (*_chargU)(chargino,0);
	  }
	  else {
	    _rightlast = -l1*(*_chargU)(chargino,1);
	  }
	}
      }
      else if(abs(iferm)==11||abs(iferm)==13||abs(iferm)==15) {
	Energy lmass = getParticleData(abs(iferm))->mass();
	_leftlast = -lmass*conj((*_chargU)(chargino,1))/sqrt(2)/_cb/_mw;
	_rightlast = (*_chargV)(chargino,0);
      }
      else if(abs(iferm)==2||abs(iferm)==4||abs(iferm)==6) {
	Energy massu = getParticleData(abs(iferm))->mass();
	Energy massd = getParticleData(abs(iferm) - 1)->mass();
	Complex u1 = massu*(*_chargV)(chargino,1)/sqrt(2)/_mw/_sb;
	Complex u2 = massd*(*_chargU)(chargino,1)/sqrt(2)/_mw/_cb;
	if(abs(iferm)==6) {
	  _leftlast = -conj(u1)*(*_sbottom)(0,eig);
	  _rightlast = (*_chargU)(chargino,0)*(*_sbottom)(0,eig) 
	    - u2*(*_sbottom)(1,eig);	  
	}
	else {
	  if(eig==0) {
	    _leftlast = -conj(u1);
	    _rightlast = (*_chargU)(chargino,0);
	  }
	  else {
	    _leftlast = 0.;
	    _rightlast = -u2;
	  }
	}
      }
      else {
	Energy massu = getParticleData(abs(iferm)+1)->mass();
	Energy massd = getParticleData(abs(iferm))->mass();
	
	Complex d1 = massd*(*_chargU)(chargino,1)/sqrt(2)/_mw/_cb;
	Complex d2 = massu*(*_chargV)(chargino,1)/sqrt(2)/_mw/_sb;
	if(abs(iferm)==5) {
	  _leftlast = -conj(d1)*(*_stop)(0,eig);
	  _rightlast = (*_chargV)(0,eig) - d2*(*_stop)(1,eig);
	}
	else {
	  if(eig==0) {
	    _leftlast = -conj(d1);
	    _rightlast = (*_chargV)(chargino,0); 
	  }
	  else {
	    _leftlast = 0.;
	    _rightlast = -d2;
	  }
	}
      }
    }
    setNorm(_couplast);
    //work out the correct left and right
    tcPDPtr incoming, fermion;
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
    else if(incoming->iSpin() == PDT::Spin1Half && 
	    abs(incoming->id()) != 1000024 && abs(incoming->id()) != 1000037) {
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
    throw HelicityConsistencyError() << "SSCFSVertex::setcoupling  - "
				     << "Incorrect particle found in "
				     << "SSCFSVertex.\n"
				     << Exception::warning;
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
  }
}
