// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMHSFSFVertex class.
//

#include "NMSSMHSFSFVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMHSFSFVertex::NMSSMHSFSFVertex() : 
  _triTp(0.*MeV), _triBt(0.*MeV), _triTa(0.*MeV), _lambda(0.),
  _lambdaVEV(0.*MeV), _v1(0.*MeV), _v2(0.*MeV), _sw(0.), _cw(0.), 
  _mw(0.*MeV), _mz(0.*MeV), _sb(0.), _cb(0.), _tb(0.), _q2last(0.*MeV2), 
  _couplast(0.), _masslast(make_pair(0.*MeV,0.*MeV)), _idlast(make_pair(0,0)) {
  orderInGem(1);
  orderInGs(0);
}

void NMSSMHSFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _mixS << _mixP << _mixTp << _mixBt << _mixTa
     << ounit(_triTp,GeV) << ounit(_triBt,GeV) << ounit(_triTa,GeV) 
     << _lambda << ounit(_lambdaVEV,GeV) << ounit(_v1,GeV) << ounit(_v2,GeV)
     << _sw << _cw << ounit(_mw,GeV) << ounit(_mz,GeV) << _sb << _cb
     << _tb;
}


void NMSSMHSFSFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _mixS >> _mixP >> _mixTp >> _mixBt >> _mixTa
     >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) >> iunit(_triTa,GeV) 
     >> _lambda >> iunit(_lambdaVEV,GeV) >> iunit(_v1,GeV) >> iunit(_v2,GeV) 
     >> _sw >> _cw >> iunit(_mw,GeV) >> iunit(_mz,GeV) >> _sb >> _cb >> _tb;
}

ClassDescription<NMSSMHSFSFVertex> NMSSMHSFSFVertex::initNMSSMHSFSFVertex;
// Definition of the static class description member.

void NMSSMHSFSFVertex::Init() {

  static ClassDocumentation<NMSSMHSFSFVertex> documentation
    ("The coupling of Higgs bosons to sfermions in the MSSM.");

}

void NMSSMHSFSFVertex::doinit() {
  //CP even
  int even[3] = {25, 35, 45};
  for(size_t h = 0; h < 3; ++h ) {
    //squarks
    for(long q = 1; q < 7; ++q) {
      //11
      addToList(even[h], -1000000 - q, 1000000 + q);
      //22
      addToList(even[h], -2000000 - q, 2000000 + q);
      //12
      addToList(even[h], -1000000 - q, 2000000 + q);
      //21
      addToList(even[h], -2000000 - q, 1000000 + q);
    }
    //sleptons
    for(long l = 11; l < 17; ++l) {
      //11
      addToList(even[h], -1000000 - l, 1000000 + l);
      //no right handed sneutrinos
      if( l % 2 != 0 ) {
	//22
	addToList(even[h], -2000000 - l, 2000000 + l);
	//12
	addToList(even[h], -1000000 - l, 2000000 + l);
	//21
	addToList(even[h], -2000000 - l, 1000000 + l);
      }
    }
  }
  //CP odd
  int odd[2] = {36, 46};
  for(size_t h = 0; h < 2; ++h ) {
    //squarks
    for(long q = 1; q < 7; ++q) {
      //12
      addToList(odd[h], -1000000 - q, 2000000 + q);
      //21
      addToList(odd[h], -2000000 - q, 1000000 + q);
    }
    //sleptons
    for(long l = 11; l < 16; l += 2) {
      //12
      addToList(odd[h], -1000000 - l, 2000000 + l);
      //21
      addToList(odd[h], -2000000 - l, 1000000 + l);
    }
  }
  //charged higgs
  //squarks
  for(long q = 1; q < 4; ++q ) {
    //H-
    //LL
    addToList(-37, -2*q - 999999, 2*q + 1000000);
    //RR
    addToList(-37, -2*q - 1999999, 2*q + 2000000);
    //LR
    addToList(-37, -2*q - 999999, 2*q + 2000000);
    //RL
    addToList(-37, -2*q - 1999999, 2*q + 1000000);
    //H+
    //LL
    addToList(37, -2*q - 1000000, 2*q + 999999);
    //RR
    addToList(37, -2*q - 2000000, 2*q + 1999999);
    //LR
    addToList(37, -2*q - 1000000, 2*q + 1999999);
    //RL
    addToList(37, -2*q - 2000000, 2*q + 999999);
  }
  //sleptons
  //easier as there are no right handed sneutrinos
  for(long l = 11; l <= 15; l +=2 ) {
    //H-
    //LL
    addToList(-37, -l - 1000000, l + 1000001);
    //RL
    addToList(-37, -l - 2000000, l + 1000001);
    //H+
    //LL
    addToList(+37, -l - 1000001, l + 1000000);
    //RL
    addToList(+37, -l - 1000001, l + 2000000);
  }
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if( !nmssm )
    throw InitException() << "NMSSMHSFSFVertex::doinit() - The model pointer "
			  << "in this vertex is not an NMSSM one as it "
			  << "should be." << Exception::runerror;
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixTp = nmssm->stopMix();
  _mixBt = nmssm->sbottomMix();
  _mixTa = nmssm->stauMix();
  if( !_mixS || !_mixP || !_mixTp || !_mixBt || !_mixTa ) 
    throw InitException() 
      << "NMSSMHSFSFVertex::doinit() - One of the mixing matrix pointers is "
      << "null, cannot continue. CP-even: " << _mixS << "  CP-odd: " << _mixP
      << "  ~t: " << _mixTp << "  ~b: " << _mixBt << "  ~tau: " << _mixTa
      << Exception::runerror;

  _triTp = nmssm->topTrilinear();
  _triBt = nmssm->bottomTrilinear();
  _triTa = nmssm->tauTrilinear();

  _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();

  _sw = sin2ThetaW();
  _cw = sqrt( 1. - _sw);
  _sw = sqrt(_sw);
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();

  _tb = nmssm->tanBeta();
  double beta = atan(_tb);
  _sb = sin(beta);
  _cb = cos(beta);
  
  _v1 = sqrt(2.)*_mw*_cb;
  _v2 = sqrt(2.)*_mw*_sb;
  SSSVertex::doinit();
}

void NMSSMHSFSFVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				   tcPDPtr part2, tcPDPtr part3) {
  // extract particle ids
  long higgs(part1->id()), isf1(part2->id()), isf2(part3->id());
  // higgs first
  if(abs(isf1)<100) swap(higgs,isf1);
  if(abs(isf2)<100) swap(higgs,isf2);
  // squark second
  if(isf1<0) swap(isf1,isf2);
  // check higgs
  assert( higgs == 25 || higgs == 35 || higgs == 45 || 
	  higgs == 36 || higgs == 46 || abs(higgs) == 37 );
  // abs of antisquark and check
  isf2 *=-1;
  assert(isf1>0&&isf2>0);
  // running coupling
  if( q2 != _q2last ) {
    _q2last = q2;
    _couplast = weakCoupling(q2);
  }
 
  //charged higgs
  if( abs(higgs) == 37 ) {
    norm(_couplast*chargedHiggs(q2, isf1, isf2));
    return;
  }
  // neutral higgs
  // L/R states of the sfermions
  unsigned int alpha = ( isf1 > 2000000 ) ? 1 : 0;
  unsigned int beta  = ( isf2 > 2000000 ) ? 1 : 0;
  // id nad mass of corresponding SM fermion
  long smid = ( alpha == 0 ) ? isf1 - 1000000 : isf1 - 2000000;
  if( q2 != _q2last || smid != _idlast.first) {
    _idlast.first = smid;
    _masslast.first = _theSM->mass(q2, getParticleData(smid));
  }
  double f1 = _masslast.first/_mw;
  complex<Energy> af(ZERO);
  Complex m1a(0.), m1b(0.), m2a(0.), m2b(0.);
  // mixing for down type squarks and charged sleptons
  if( smid % 2 != 0 ) {
    f1 /= _cb;
    // sbottom
    if( smid == 5 ) {
      m1a = (*_mixBt)(alpha, 0);
      m1b = (*_mixBt)(alpha, 1);
      m2a = (*_mixBt)(beta , 0) ;
      m2b = (*_mixBt)(beta , 1);
      af = _triBt;
    }
    // stau
    else if( smid == 15 ) {
      m1a = (*_mixTa)(alpha, 0);
      m1b = (*_mixTa)(alpha, 1);
      m2a = (*_mixTa)(beta , 0) ;
      m2b = (*_mixTa)(beta , 1);
      af = _triTa;
    }
    // 1st 2 generations
    else {
      m1a = (alpha == 0) ? 1. : 0.;
      m1b = (alpha == 0) ? 0. : 1.;
      m2a = (beta  == 0) ? 1. : 0.;
      m2b = (beta  == 0) ? 0. : 1.;
      af = ZERO;
    }
  }
  // mixing for up type squarks and sneutrions
  else {
    f1 /= _sb;
    // stop
    if( smid == 6 ) {	  
      m1a = (*_mixTp)(alpha, 0);
      m1b = (*_mixTp)(alpha, 1);
      m2a = (*_mixTp)(beta , 0);
      m2b = (*_mixTp)(beta , 1);
      af = _triTp;
    }
    // everything else
    else {
      m1a = (alpha == 0) ? 1. : 0.;
      m1b = (alpha == 0) ? 0. : 1.;
      m2a = (beta == 0) ? 1. : 0.;
      m2b = (beta == 0) ? 0. : 1.;
      af = 0.*MeV;
    }
  }
  // scalar higgs bosons
  complex<Energy> fact(ZERO);
  if( higgs == 25 || higgs == 35 || higgs == 45 ) {
    int iloc = (higgs - 25)/10;
    complex<Energy> f2 = 0.5*_mz*( - _cb*(*_mixS)(iloc,0) 
				   + _sb*(*_mixS)(iloc,1))/_cw;

    // down type squarks and charged sleptons
    if( smid % 2 != 0 ) {
      double ef = (smid < 7) ? -1./3. : -1.;
      fact = - f2*( (1. + 2.*ef*sqr(_sw))*m1a*m2a - 2.*ef*sqr(_sw)*m1b*m2b)
	- f1*_masslast.first*(*_mixS)(iloc,0)*(m1a*m2a + m1b*m2b) 
	- 0.5*f1*(( - _lambdaVEV*(*_mixS)(iloc,1) 
		    - _lambda*_v2*(*_mixS)(iloc,2)/_couplast 
		    +  af*(*_mixS)(iloc,0)) * 
		  (m2a*m1b + m1a*m2b) );
    }
    // up type squarks and sneutrinos
    else {
      double ef = (smid < 7) ? 2./3. : 0.;
      fact =  +f2*( (1. - 2.*ef*sqr(_sw))*m1a*m2a + 2.*ef*sqr(_sw)*m1b*m2b )
	- f1*_masslast.first*(*_mixS)(iloc,1)*(m1a*m2a + m1b*m2b)  
	-  0.5*f1*(( - _lambdaVEV*(*_mixS)(iloc,0) 
		     - _lambda*_v1*(*_mixS)(iloc,2)/_couplast
		     +  af*(*_mixS)(iloc,1) ) *
		   (m2a*m1b + m1a*m2b));
    }
  }
  // pseudo scalar
  else if( higgs == 36 || higgs == 46 ) {
    int iloc = (higgs - 36)/10;
    // down type squarks and charged sleptons
    if( smid % 2 != 0 ) {
      fact = -0.5*f1*Complex(0.0,1.0)*
	( _lambdaVEV*(*_mixP)(iloc,1) +
	  _lambda*_v2*(*_mixP)(iloc,2)/_couplast +
	  af*(*_mixP)(iloc,0) );
    }
    // up-type squarks and sneutrinos
    else {
      fact =-0.5*f1*Complex(0.0,1.0)*
	( _lambdaVEV *(*_mixP)(iloc,0) +
	  _lambda*_v1*(*_mixP)(iloc,2)/_couplast +
	  af*(*_mixP)(iloc,1));
    }
    if(alpha<beta) fact *= -1.;
  }
  norm(_couplast*fact*UnitRemoval::InvE);
}

Complex NMSSMHSFSFVertex::chargedHiggs(Energy2 q2, long id1, long id2) {
  //have id1 as up-type
  if( id1 % 2 != 0) swap(id1, id2);
  // sfermion L/R states
  unsigned int alpha = ( id1/1000000 == 2 ) ? 1 : 0;
  unsigned int beta  = ( id2/1000000 == 2 ) ? 1 : 0;
  // type of quarks
  long utype = (alpha == 0) ? id1 - 1000000 : id1 - 2000000;
  long dtype = ( beta == 0) ? id2 - 1000000 : id2 - 2000000;
  // compute the running masses
  if( q2 != _q2last || id1 != _idlast.first || id2 != _idlast.second) {
    _idlast.first  = id1;
    _idlast.second = id2;
    _masslast.first  = _theSM->mass(q2, getParticleData(utype) );
    _masslast.second = _theSM->mass(q2, getParticleData(dtype) );
  }
  Energy2 facta = 2.*sqr(_mw)*_sb*_cb;
  complex<Energy2> coupling(ZERO);
  // sleptons
  if( dtype == 11 || dtype == 13 || dtype == 15) {
    Complex l1b = 0., l2b = 0.;
    complex<Energy> tri(ZERO);
    // 1st 2 generations
    if (dtype == 11 || dtype == 13) {
      l1b = (beta == 0) ? 1.0 : 0.0;
      l2b = (beta == 0) ? 0.0 : 1.0;
    }
    // stau
    else {
      l1b = (*_mixTa)(beta, 0) ;
      l2b = (*_mixTa)(beta, 1);
      tri = _triTa;
    }
    coupling = ( l1b*(sqr(_masslast.second)*_tb - facta) +
		 l2b*_masslast.second*(tri*_tb + _lambdaVEV) );
  }
  // squarks
  else {
    Complex q1a(0.0), q1b(0.0), q2a(0.0), q2b(0.0);
    complex<Energy> triD(ZERO), triU(ZERO);
    // up-type bit
    // stop
    if(utype == 6){
      q1a =  (*_mixTp)(alpha, 0) ;
      q2a =  (*_mixTp)(alpha, 1);
      triU = _triTp;
    }
    // light
    else{
      q1a = (alpha == 0) ? 1.0 : 0.0;
      q2a = (alpha == 0) ? 0.0 : 1.0;
    }
    // down-type bit
    // sbottom
    if(dtype == 5){
      q1b =  (*_mixBt)(beta, 0) ;
      q2b =  (*_mixBt)(beta, 1);
      triD = _triBt;
    }
    // light
    else{
      q1b = (beta == 0) ? 1.0 : 0.0;
      q2b = (beta == 0) ? 0.0 : 1.0;
    }
    Energy mfu = _masslast.first;
    Energy mfd = _masslast.second;
    coupling = ( q1a*q1b*((sqr(mfd)*_tb + sqr(mfu)/_tb) - facta)
		 + q2a*q2b*mfu*mfd*(_tb + (1./_tb))
		 + q1a*q2b*mfd*(triD*_tb + _lambdaVEV)
		 + q2a*q1b*mfu*(triU/_tb + _lambdaVEV));
  }
  return coupling * UnitRemoval::InvE/_mw/sqrt(2.);
}
