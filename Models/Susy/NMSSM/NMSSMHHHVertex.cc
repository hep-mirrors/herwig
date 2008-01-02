// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMHHHVertex class.
//

#include "NMSSMHHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;

NMSSMHHHVertex::NMSSMHHHVertex() : _mw(0.*MeV), _mz(0.*MeV), _sw2(0.),
				   _cw(0.), _lambda(0.), _kappa(0.) ,
				   _lambdaVEV(0.*MeV), _theAl(0.*MeV),
				   _theAk(0.*MeV), _sb(0.), _cb(0.),
				   _s2b(0.), _c2b(0.), _v1(0.*MeV),
				   _v2(0.*MeV), _q2last(0.*MeV2),
				   _glast(0.) {
  vector<int> first, second, third;
  //CP-even Higgs
  first.push_back(25); second.push_back(35); third.push_back(45);
  for( unsigned int i = 25; i < 46; i += 10 ) {
    first.push_back(i); second.push_back(i); third.push_back(25);
    first.push_back(i); second.push_back(i); third.push_back(35);
    first.push_back(i); second.push_back(i); third.push_back(45);
    //Charged Higgs
    first.push_back(i); second.push_back(37); third.push_back(-37);
    //CP-odd Higgs
    first.push_back(i); second.push_back(36); third.push_back(36);
    first.push_back(i); second.push_back(36); third.push_back(46);
    first.push_back(i); second.push_back(46); third.push_back(36);
    first.push_back(i); second.push_back(46); third.push_back(46);
  }
  setList(first, second, third);
}

void NMSSMHHHVertex::doinit() throw(InitException) {
  _theSM = generator()->standardModel();
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if( !nmssm ) 
    throw InitException() << "NMSSMHHHVertex::doinit - The model object is"
			  << "not the NMSSM object."
			  << Exception::runerror;
  
  //SM parameters
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _sw2 = _theSM->sin2ThetaW();
  _cw = sqrt(1. - _sw2);
  
  //NMSSM parameters
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  if( !_mixS || !_mixP ) 
    throw InitException() << "NMSSMHHHVertex::doinit - One of the mixing matrix "
			  << "pointers is null, cannot continue. S: "
			  << _mixS << "  P: " << _mixP << Exception::runerror;
  _lambda = nmssm->lambda();
  _kappa = nmssm->kappa();
  _lambdaVEV = nmssm->lambdaVEV();
  _theAl = nmssm->trilinearLambda();
  _theAk = nmssm->trilinearKappa();
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);
  _s2b = 2.*_sb*_cb;
  _c2b = sqr(_cb) - sqr(_sb);
  _v1 = 2.*_mw*_cb;
  _v2 = 2.*_mw*_sb;
  
  orderInGem(1);
  orderInGs(0);
  SSSVertex::doinit();
}

void NMSSMHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw, GeV) << ounit(_mz,GeV) << _sw2 << _cw <<  _lambda 
     << _kappa <<  ounit(_lambdaVEV,GeV) <<  ounit(_theAl, GeV) 
     << ounit(_theAk,GeV) <<  _sb <<  _cb << _s2b <<  _c2b
     << ounit(_v1,GeV) << ounit(_v2,GeV) << _theSM << _mixS << _mixP;
}

void NMSSMHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw, GeV) >> iunit(_mz,GeV) >> _sw2 >> _cw >>  _lambda 
     >> _kappa >>  iunit(_lambdaVEV,GeV) >>  iunit(_theAl, GeV) 
     >> iunit(_theAk,GeV) >>  _sb >>  _cb >> _s2b >>  _c2b
     >> iunit(_v1,GeV) >> iunit(_v2,GeV) >> _theSM >> _mixS >> _mixP;
}

ClassDescription<NMSSMHHHVertex> NMSSMHHHVertex::initNMSSMHHHVertex;
// Definition of the static class description member.

void NMSSMHHHVertex::Init() {

  static ClassDocumentation<NMSSMHHHVertex> documentation
    ("There is the triple Higgs coupling in the NMSSM.");

}

void NMSSMHHHVertex::setCoupling(Energy2 q2,tcPDPtr p1,tcPDPtr p2, 
				 tcPDPtr p3) {
  long higgs[3] = {p1->id(), p2->id(), p3->id()};
  unsigned int ns(0), np(0), nc(0);
  for( int i = 0; i < 3; ++i ) {
    if( higgs[i] == 25 || higgs[i] == 35 || higgs[i] == 45 )
      ++ns;
    else if( higgs[i] == 36 || higgs[i] == 46 )
      ++np;
    else if( abs(higgs[i]) ==  37 )
      ++nc;
  }
  assert( ns + np + nc == 3 );
  if( q2 != _q2last ) {
    _q2last = q2;
    _glast = sqrt(4.*Constants::pi*_theSM->alphaEM(q2)/_sw2);
  }
  complex<Energy> coupling;
  double rt = sqrt(2);
  if( ns == 3 ) {
    unsigned int a = (higgs[0] - 25)/10;
    unsigned int b = (higgs[1] - 25)/10;
    unsigned int c = (higgs[2] - 25)/10;

    double gf = 0.5*sqr(_glast)*(1. + sqr(_sw2))/rt;

    coupling = 
      -3.*gf*_glast*(_v1*(*_mixS)(a,0)*(*_mixS)(b,0)*(*_mixS)(c,0) + 
		     _v2*(*_mixS)(a,1)*(*_mixS)(b,1)*(*_mixS)(c,1) )
      + (gf - rt*sqr(_lambda))*_glast *
      ( _v1*usMix(a,0,b,1,c,1) + _v2*usMix(a,1,b,0,c,0) )
      + rt*_lambda*_glast*(_kappa*_v2 - _lambda*_v1)*usMix(a,0,b,2,c,2)
      + rt*_lambda*_glast*(_kappa*_v1 - _lambda*_v2)*usMix(a,1,b,2,c,2)
      - rt*_lambda*_lambdaVEV*( usMix(a,2,b,0,c,0) + usMix(a,2,b,1,c,1) )
      + (_lambda*_theAl/rt + rt*_lambdaVEV*_kappa)*( usMix(a,0,b,1,c,2) + 
						 usMix(a,2,b,1,c,0) )
      + (rt*_kappa*_theAk - 6.*rt*sqr(_kappa)*_lambdaVEV/_lambda) * 
      (*_mixS)(a,2)*(*_mixS)(b,2)*(*_mixS)(c,2);
  }
  else if(ns == 1 && np == 2) {
    unsigned int a(0), b(0), c(0);
    if( higgs[0] == 25 || higgs[0] == 35 || higgs[0] == 45 ) {
      a = (higgs[0] - 25)/10;
      b = (higgs[1] - 36)/10;
      c = (higgs[2] - 36)/10;
    }
    else if(higgs[1] == 25 || higgs[1] == 35 || higgs[1] == 45 ) {
      a = (higgs[1] - 25)/10;
      b = (higgs[0] - 36)/10;
      c = (higgs[2] - 36)/10;
    }
    else {
      a = (higgs[2] - 25)/10;
      b = (higgs[0] - 36)/10;
      c = (higgs[1] - 36)/10;
    }
    double rt = sqrt(2);
    double gf = 0.5*sqr(_glast)*(1. + sqr(_sw2))/rt;
    
    coupling =
      -gf*_glast*(_v1*upMix(a,0,b,0,c,0) + _v2*upMix(a,1,b,1,c,1))
      + (gf - rt*sqr(_lambda))*_glast*(_v1*upMix(a,0,b,1,c,1) + 
				       _v2*upMix(a,1,b,0,c,0))
      - rt*_lambda*_glast*(_kappa*_v1 + _lambda*_v2)*upMix(a,1,b,2,c,2)
      - rt*_lambda*_glast*(_kappa*_v2 + _lambda*_v1)*upMix(a,0,b,2,c,2)
      - rt*_lambda*_lambdaVEV*(upMix(a,2,b,0,c,0) + upMix(a,2,b,1,c,1))
      - (2.*rt*sqr(_kappa)*_lambdaVEV/_lambda + rt*_kappa*_theAk) * 
        upMix(a,2,b,2,c,2)
      + rt*_lambda*_kappa*( _v1*(upMix(a,2,b,1,c,2) + upMix(a,2,b,2,c,1) )
			    + _v2*(upMix(a,2,b,0,c,2) + upMix(a,2,b,2,c,0)) )
      + (rt*_lambdaVEV*_kappa - _lambda*_theAl/rt)* 
        (upMix(a,0,b,1,c,2) + upMix(a,0,b,2,c,1) + 
	 upMix(a,1,b,0,c,2) + upMix(a,1,b,2,c,0))
      - (rt*_lambdaVEV*_kappa + _lambda*_theAl/rt) * 
      (upMix(a,2,b,0,c,1) + upMix(a,2,b,1,c,0));
  }
  else {
    unsigned int a(0);
    if( higgs[0] == 25 || higgs[0] == 35 || higgs[0] == 45 )
      a = (higgs[0] - 25)/10;
    else if(higgs[1] == 25 || higgs[1] == 35 || higgs[1] == 45 )
      a = (higgs[1] - 25)/10;
    else
      a = (higgs[2] - 25)/10;
    
    coupling = 
      -_glast*_mw*((*_mixS)(a,0)*_cb + (*_mixS)(a,1)*_sb)
      - 0.5*_glast*_mz*_c2b*((*_mixS)(a,1)*_sb - (*_mixS)(a,0)*_cb)/_cw
      + _glast*sqr(_lambda)*((*_mixS)(a,1)*_v1 - (*_mixS)(a,0)*_v2)/rt
      - (*_mixS)(a,2)*( (2.*_kappa*_lambdaVEV + _theAl)*_s2b 
			+ 2.*_lambdaVEV)/rt;
  }
  setNorm(coupling * UnitRemoval::InvE);
}
