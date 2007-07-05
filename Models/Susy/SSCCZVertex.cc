// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCCZVertex class.
//

#include "SSCCZVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCCZVertex::SSCCZVertex() : _sw2(0.), _cw(0.), _couplast(0.),
			     _q2last(), _id1last(0), _id2last(0) {
  vector<int> first, second, third(4, 23);
  for(unsigned int ix = 0; ix < 2; ++ix) {
    int ic1(1000024);
    if(ix == 1) ic1 = 1000037;
    for(unsigned int iy = 0; iy < 2; ++iy) {
      int ic2(1000024);
      if(iy == 1) ic2 = 1000037;
      first.push_back(-ic1);
      second.push_back(ic2);
    }
  }
  setList(first, second, third);
}

void SSCCZVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sw2 << _cw << _theU << _theV;
}

void SSCCZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS >> _sw2 >> _cw >> _theU >> _theV;
  _couplast = 0.;
  _q2last = 0.*GeV2;
  _id1last = 0;
  _id2last = 0;
  _leftlast = 0.;
  _rightlast = 0.;
}

ClassDescription<SSCCZVertex> SSCCZVertex::initSSCCZVertex;
// Definition of the static class description member.

void SSCCZVertex::Init() {

  static ClassDocumentation<SSCCZVertex> documentation
    ("This class implements the coupling of a Z-boson to a pair of"
     " charginos. ");

}

void SSCCZVertex::setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			      tcPDPtr part3) {
  long ichar1, ichar2;
  if(part1->id() == ParticleID::Z0) {
    ichar1 = abs(part2->id());
    ichar2 = abs(part3->id());
  }
  else if(part2->id() == ParticleID::Z0) {
    ichar1 = abs(part1->id());
    ichar2 = abs(part3->id());
  }
  else if(part3->id() == ParticleID::Z0) {
    ichar1 = abs(part1->id());
    ichar2 = abs(part2->id());
  }
  else 
    throw HelicityConsistencyError() << "SSCCZVertex::setCoupling() - There is "
				     << "no Z0 boson this vertex! "
				     << Exception::warning;
  if((ichar1 == 1000024 || ichar1 == 1000037) && 
     (ichar2 == 1000024 || ichar2 == 1000037) ) {
    if(_q2last != q2) {
      _q2last = q2;
      _couplast = sqrt(4.*Constants::pi*_theSS->alphaEM(q2)/_sw2)/_cw;
    }
    setNorm(_couplast);
    if(ichar1 != _id1last || ichar2 != _id2last) {
      _id1last = ichar1;
      _id2last = ichar2;
      unsigned int ic1(0), ic2(0);
      if(ichar1 == 1000037) ic1 = 1;
      if(ichar2 == 1000037) ic2 = 1;
      _leftlast = -(*_theV)(ic1, 0)*conj((*_theV)(ic2, 0)) - 
	0.5*(*_theV)(ic1, 1)*conj((*_theV)(ic2, 1));
      _rightlast = -conj((*_theU)(ic1, 0))*(*_theU)(ic2, 0) - 
	0.5*conj((*_theU)(ic1, 1))*(*_theU)(ic2, 1);
      if(ichar1 == ichar2) {
	_leftlast += _sw2;
	_rightlast += _sw2;
      }
    }
    setLeft(_leftlast);
    setRight(_rightlast);
  }
  else {
    setNorm(0.);
    setLeft(0.);
    setRight(0.);
    throw HelicityConsistencyError() << "SSCCZVertex::setCoupling() - There are "
				     << "no charginos in this vertex! "
				     << Exception::warning;
  }
     
}
