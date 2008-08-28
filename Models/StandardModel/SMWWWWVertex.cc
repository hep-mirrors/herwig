// -*- C++ -*-
//
// SMWWWWVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWWWVertex class.
//

#include "SMWWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMWWWWVertex::SMWWWWVertex() 
  : _couplast(0.0), _q2last(sqr(Constants::MaxEnergy)), 
    _vfact(4,0.0), _sw2(0.), _cw2(0.) {
  // particles
  vector<long> first,second,third,fourth;
  first.push_back(24);
  first.push_back(23);
  first.push_back(22);
  first.push_back(22);
  second.push_back(-24);
  second.push_back(24);
  second.push_back(24);
  second.push_back(24);
  third.push_back(24);
  third.push_back(23);
  third.push_back(22);
  third.push_back(23);
  fourth.push_back(-24);
  fourth.push_back(-24);
  fourth.push_back(-24);
  fourth.push_back(-24);
  setList(first,second,third,fourth);
}

void SMWWWWVertex::doinit() throw(InitException) {
  orderInGem(2);
  orderInGs(0);
  VVVVVertex::doinit();
  // couplings
  _sw2 = generator()->standardModel()->sin2ThetaW();
  _cw2 = 1.-_sw2;
  double sw = sqrt(_sw2);
  double cw = sqrt(_cw2);
  _vfact[0] = -1./_sw2;
  _vfact[1] = _cw2/_sw2;
  _vfact[2] = 1.;
  _vfact[3] = cw/sw;
  // pointer for intermediate particles
  _gamma  = getParticleData(ThePEG::ParticleID::gamma);
  _Z0     = getParticleData(ThePEG::ParticleID::Z0);
  _wplus  = getParticleData(ThePEG::ParticleID::Wplus);
  _wminus = getParticleData(ThePEG::ParticleID::Wminus);
}

void SMWWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _gamma << _Z0 << _wplus << _wminus
     << _vfact  << _sw2 << _cw2;
}

void SMWWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gamma >> _Z0 >> _wplus >> _wminus
     >> _vfact >> _sw2 >> _cw2;
}

ClassDescription<SMWWWWVertex>SMWWWWVertex::initSMWWWWVertex;
// Definition of the static class description member.

void SMWWWWVertex::Init() {
  static ClassDocumentation<SMWWWWVertex> documentation
    ("The SMWWWWVertex class is the implementation of the"
     " Standard Model quartic electroweka gauge boson coupling.");
  
}


// couplings for the WWWW vertex
void SMWWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
    				      tcPDPtr c,tcPDPtr d) {
  // id's of the particles
  int id[4]={a->id(),b->id(),c->id(),d->id()};
  // order the particles
  int ngamma(0),nz(0);
  int iorder[4];
  for(int ix=0;ix<4;++ix) {
    if      (id[ix]==22) ++ngamma;
    else if (id[ix]==23) ++nz;
  }
  // if photons or Z's
  if(ngamma!=0 || nz!=0) {
    int iy=0;
    // put the photons first
    for(int ix=0;iy<ngamma&&ix<4;++ix) {
      if(id[ix]==22) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the Z bosons
    for(int ix=0;iy<ngamma+nz&&ix<4;++ix) {
      if(id[ix]==23) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    if (iy!=3) throw HelicityConsistencyError() 
      << "SMWWWWVertex::setCoupling Error setting order" 
      << Exception::runerror;
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    if(iy!=4) 
      throw HelicityConsistencyError() 
	<< "SMWWWWVertex::setCoupling Error setting order" 
	<< Exception::warning;
  }
  else {
    int iy=0;
    // first the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    if(iy!=2) throw HelicityConsistencyError() 
      << "SMWWWWVertex::setCoupling Error setting order" 
      << Exception::warning;
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    if(iy!=4) throw HelicityConsistencyError() 
      << "SMWWWWVertex::setCoupling Error setting order" 
      << Exception::warning;
    setIntermediate(_gamma,_Z0,_sw2,_cw2);
  }
  setOrder(iorder[0],iorder[1],iorder[2],iorder[3]);
  setType(2);
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = sqr(electroMagneticCoupling(q2));
    _q2last=q2;
  }
  // id's of the first two particles
  int ida(0),idb(0);
  if(iorder[0]==0)      ida = abs(a->id());
  else if(iorder[0]==1) ida = abs(b->id());
  else if(iorder[0]==2) ida = abs(c->id());
  else if(iorder[0]==3) ida = abs(d->id());
  if(iorder[1]==0)      idb = abs(a->id());
  else if(iorder[1]==1) idb = abs(b->id());
  else if(iorder[1]==2) idb = abs(c->id());
  else if(iorder[1]==3) idb = abs(d->id());
  // WWWW coupling
  if(ida==24)               setNorm(_vfact[0]*_couplast);
  // ZZWW coupling
  else if(ida==23&&idb==23) setNorm(_vfact[1]*_couplast);
  // gamma gamma WW coupling
  else if(ida==22&&idb==22) setNorm(_couplast);
  // gamma  Z WW coupling
  else if(ida==22&&idb==23) setNorm(_vfact[3]*_couplast);
  else throw HelicityConsistencyError() 
    << "SMWWWWVertex::setCoupling unknown particles" 
    << Exception::runerror;
}
