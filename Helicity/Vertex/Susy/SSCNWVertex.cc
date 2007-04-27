// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCNWVertex class.
//

#include "SSCNWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig::Helicity;

SSCNWVertex::SSCNWVertex() : _sw(0.),  _couplast(0.), _q2last(0.), 
			     _id1last(0), _id2last(0), _leftlast(0.),
			     _rightlast(0.) {
  vector<int> first, second, third;
  //iw == -1 outgoing W-, iw == +1 outgoing W+
  for(int iw = -1; iw < 2; iw += 2) {
    for(unsigned int ine = 0; ine < 4; ++ine) {
      int neu(1000022);
      if(ine == 1) neu = 1000023;
      if(ine == 2) neu = 1000025;
      if(ine == 3) neu = 1000035;
      for(unsigned int ic = 0; ic < 2; ++ic) {
	int cha(1000024);
	if(ic == 1) cha = 1000037;
	first.push_back(-iw*cha);
	second.push_back(neu);
	third.push_back(iw*24);
      }
    }
  }
  setList(first, second, third);
}

void SSCNWVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sw << _theN << _theU << _theV;
}

void SSCNWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _sw >> _theN >> _theU >> _theV;
  _id1last = 0;
  _id2last = 0;
  _q2last = 0.;
  _couplast = 0.;
  _leftlast = 0.;
  _rightlast = 0.;
}

ClassDescription<SSCNWVertex> SSCNWVertex::initSSCNWVertex;
// Definition of the static class description member.

void SSCNWVertex::Init() {

  static ClassDocumentation<SSCNWVertex> documentation
    ("This class implements the coupling of a W boson to a "
     "neutralino and a chargino");

}

void SSCNWVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			      tcPDPtr part3) {
  long neu, cha;
  if(abs(part1->id()) == ParticleID::Wplus) {
    if(part2->charged()) {
      cha = part2->id();
      neu = part3->id();
    }
    else {
      cha = part3->id();
      neu = part2->id();
    }
  }
  else if(abs(part2->id()) == ParticleID::Wplus) {
    if(part1->charged()) {
      cha = part1->id();
      neu = part3->id();
    }
    else {
      cha = part3->id();
      neu = part1->id();
    }
  }
  else if(abs(part3->id()) == ParticleID::Wplus) {
    if(part1->charged()) {
      cha = part1->id();
      neu = part2->id();
    }
    else {
      cha = part2->id();
      neu = part1->id();
    }
  }
  else
    throw HelicityConsistencyError() << "SSCNWVertex::setCoupling() - There is "
				     << "no W-boson in this vertex! "
				     << Exception::warning;
  if((abs(cha) == 1000024 || abs(cha) == 1000037) && 
     (abs(neu) == 1000022 || abs(neu) == 1000023 || 
      abs(neu) == 1000025 || abs(neu) == 1000035) ) {
    if(q2 != _q2last) {
      _q2last = q2;
      _couplast = sqrt(4.*Constants::pi*_theSS->alphaEM(q2))/_sw;
    }
    setNorm(_couplast);
    if(abs(cha) != _id1last || abs(neu) != _id2last) {
      _id1last = abs(cha);
      _id2last = abs(neu);
      unsigned int eigc(0);
      if(abs(cha) == 1000037) eigc = 1;
      unsigned int eign(0);
      if(abs(neu) == 1000023) eign = 1;
      if(abs(neu) == 1000025) eign = 2;
      if(abs(neu) == 1000035) eign = 3;
      _leftlast = (*_theN)(eign, 1)*conj((*_theV)(eigc, 0)) - 
	( (*_theN)(eign, 3)*conj((*_theV)(eigc, 1))/sqrt(2));

      _rightlast = conj((*_theN)(eign, 1))*(*_theU)(eigc, 0) +
	( conj((*_theN)(eign, 2))*(*_theU)(eigc, 1)/sqrt(2));
    }
    setLeft(_leftlast);
    setRight(_rightlast);
  }
  else
    throw HelicityConsistencyError() << "SSCNWVertex::setCoupling() - "
				     << "Something other than a neutralino or a "
				     << "chargino has appeared in this vertex "
				     << neu << "  " << cha
				     << Exception::warning;

}

