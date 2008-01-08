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
  vector<int> first, second, third;
  //CP even
  int even[3] = {25, 35, 45};
  for( unsigned int h = 0; h < 3; ++h ) {
    //squarks
    for( unsigned int q = 1; q < 7; ++q) {
      //11
      first.push_back(even[h]);
      second.push_back(-1000000 - q);
      third.push_back(1000000 + q);
      //22
      first.push_back(even[h]);
      second.push_back(-2000000 - q);
      third.push_back(2000000 + q);
      //12
      first.push_back(even[h]);
      second.push_back(-1000000 - q);
      third.push_back(2000000 + q);
      //21
      first.push_back(even[h]);
      second.push_back(-2000000 - q);
      third.push_back(1000000 + q);
    }
    //sleptons
    for( unsigned int l = 11; l < 17; ++l) {
      //11
      first.push_back(even[h]);
      second.push_back(-1000000 - l);
      third.push_back(1000000 + l);
      //no right handed sneutrinos
      if( l % 2 != 0 ) {
	//22
	first.push_back(even[h]);
	second.push_back(-2000000 - l);
	third.push_back(2000000 + l);
	//12
	first.push_back(even[h]);
	second.push_back(-1000000 - l);
	third.push_back(2000000 + l);
	//21
	first.push_back(even[h]);
	second.push_back(-2000000 - l);
	third.push_back(1000000 + l);
      }
    }
  }
  //CP odd
  int odd[2] = {36, 46};
  for( unsigned int h = 0; h < 2; ++h ) {
    //squarks
    for( unsigned int q = 1; q < 7; ++q) {
      //12
      first.push_back(odd[h]);
      second.push_back(-1000000 - q);
      third.push_back(2000000 + q);
      //21
      first.push_back(odd[h]);
      second.push_back(-2000000 - q);
      third.push_back(1000000 + q);
    }
    //sleptons
    for( unsigned int l = 11; l < 17; l += 2) {
      //12
      first.push_back(odd[h]);
      second.push_back(-1000000 - l);
      third.push_back(2000000 + l);
      //21
      first.push_back(odd[h]);
      second.push_back(-2000000 - l);
      third.push_back(1000000 + l);
    }
  }
  //charged higgs
  //squarks
  for( unsigned int q = 1; q < 4; ++q ) {
    //H-
    //LL
    first.push_back(-37);
    second.push_back(-2*q - 999999);
    third.push_back(2*q + 1000000);
    //RR
    first.push_back(-37);
    second.push_back(-2*q - 1999999);
    third.push_back(2*q + 2000000);
    //LR
    first.push_back(-37);
    second.push_back(-2*q - 999999);
    third.push_back(2*q + 2000000);
    //RL
    first.push_back(-37);
    second.push_back(-2*q - 1999999);
    third.push_back(2*q + 1000000);

    //H+
    //LL
    first.push_back(37);
    second.push_back(-2*q - 1000000);
    third.push_back(2*q + 999999);
    //RR
    first.push_back(37);
    second.push_back(-2*q - 2000000);
    third.push_back(2*q + 1999999);
    //LR
    first.push_back(37);
    second.push_back(-2*q - 1000000);
    third.push_back(2*q + 1999999);
    //RL
    first.push_back(37);
    second.push_back(-2*q - 2000000);
    third.push_back(2*q + 999999);
  }
  //sleptons
  //easier as there are no right handed sneutrinos
  for( unsigned int l = 11; l < 16; l +=2 ) {
    //H-
    //LL
    first.push_back(-37);
    second.push_back(-l - 1000000);
    third.push_back(l + 1000001);
    //H+
    //LL
    first.push_back(+37);
    second.push_back(-l - 1000001);
    third.push_back(l + 1000000);
  }
  setList(first, second, third);
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

void NMSSMHSFSFVertex::doinit() throw(InitException) {
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


  _sw = sqrt(_theSM->sin2ThetaW());
  _cw = sqrt( 1. - _theSM->sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();

  _v1 = 2.*_mw*_cb;
  _v2 = 2.*_mw*_sb;

  _tb = nmssm->tanBeta();
  double beta = atan(_tb);
  _sb = sin(beta);
  _cb = cos(beta);

  orderInGs(0);
  orderInGs(1);
  SSSVertex::doinit();
}

void NMSSMHSFSFVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				   tcPDPtr part2, tcPDPtr part3) {
  long higgs(0), isf1(0), isf2(0);
  if( part1->id() == 25 || part1->id() == 35 || part1->id() == 45 ||
      part1->id() == 36 || part1->id() == 46 || 
      abs(part1->id()) == 37 ) {
    higgs = part1->id();
    isf1 = part2->id();
    isf2 = part3->id();
  }
  else if( part2->id() == 25 || part2->id() == 35 || part2->id() == 45 ||
	   part2->id() == 36 || part2->id() == 46 || 
	   abs(part2->id()) == 37 ) { 
    higgs = part2->id();
    isf1 = part1->id();
    isf2 = part3->id();
  }
  else if( part3->id() == 25 || part3->id() == 35 || part3->id() == 45 ||
	   part3->id() == 36 || part3->id() == 46 || 
	   abs(part3->id()) == 37 ) { 
    higgs = part3->id();
    isf1 = part1->id();
    isf2 = part2->id();
  }
  else {
    throw HelicityConsistencyError()
	  << "NMSSMHSFSFVertex::setCoupling - There is no higgs particle "
	  << "in this vertex. " << part1->id() << " " << part2->id() << " "
	  << part3->id() << Exception::warning;
    setNorm(Complex(0.));
    return;
  }
  if( isf2 < 0 ) swap(isf1, isf2);
  int alpha = (abs(isf1) < 2000000) ? 0 : 1;
  int beta = (isf2 < 2000000) ? 0 : 1;
  Complex fact(0.);
    //charged higgs
  if( abs(higgs) == 37 ) {
    fact = chargedHiggs(q2, abs(isf1), abs(isf2));
  }
  else { 
    long smid = (alpha == 0) ? isf2 - 1000000 : isf2 - 2000000;
    if( q2 != _q2last || isf2 != _idlast.first) {
      _idlast.first = isf2;
      tcPDPtr p = getParticleData(smid);
      _masslast.first = _theSM->mass(q2, p);
    }
    double f1 = _masslast.first/_mw;
    complex<Energy> af(0.*MeV);
    Complex m1a(0.), m1b(0.), m2a(0.), m2b(0.);
    if( smid % 2 != 0 ) {
      f1 /= _cb;
      if( smid == 5 ) {
	m1a = (*_mixBt)(alpha, 0);
	m2a = (*_mixBt)(alpha, 1);
	m1b = (*_mixBt)(beta, 0);
	m2b = (*_mixBt)(beta, 1);
	af = _triBt;
      }
      else if( smid == 15 ) {
	m1a = (*_mixTa)(alpha, 0);
	m2a = (*_mixTa)(alpha, 1);
	m1b = (*_mixTa)(beta, 0);
	m2b = (*_mixTa)(beta, 1);
	af = _triTa;
      }
      else {
	m1a = (alpha == 0) ? 1. : 0.;
      	m2a = (alpha == 0) ? 0. : 1.;
	m1b = (beta == 0) ? 1. : 0.;
	m2b = (beta == 0) ? 0. : 1.;
	af = 0.*MeV;
      }
    }
    else {
      f1 /= _sb;
      if( smid == 6 ) {
	m1a = (*_mixTp)(alpha, 0);
	m2a = (*_mixTp)(alpha, 1);
	m1b = (*_mixTp)(beta, 0);
	m2b = (*_mixTp)(beta, 1);
	af = _triTp;
      }
      else {
	m1a = (alpha == 0) ? 1. : 0.;
      	m2a = (alpha == 0) ? 0. : 1.;
	m1b = (beta == 0) ? 1. : 0.;
	m2b = (beta == 0) ? 0. : 1.;
	af = 0.*MeV;
      }
    }
    if( higgs == 25 || higgs == 35 || higgs == 45 ) {
      int iloc = (higgs - 25)/10;
      Complex f2 = 0.5*_mz*(_cb*(*_mixS)(iloc,0) - _sb*(*_mixS)(iloc,1))/_cw
	* UnitRemoval::InvE;
      if( smid % 2 != 0 ) {
	double ef = (smid < 7) ? 1./3. : 1.;
	fact = f2*( (1. - 2.*ef*sqr(_sw))*m1a*m1b - 2.*ef*sqr(_sw)*m2a*m2b)
	  - f1*_masslast.first*(*_mixS)(iloc,0)*(m1a*m1b + m2a*m2b) * 
	  UnitRemoval::InvE;
	if( alpha != beta ) 
	  fact -= 0.5*f1*( _lambdaVEV*(*_mixS)(iloc,1) 
			   + _lambda*_couplast*_v2*(*_mixS)(iloc,2) 
			   +  af*(*_mixS)(iloc,1)) * 
	    (m2a*m1b + m1a*m2b) * UnitRemoval::InvE;
      }
      else {
	f1 /= _sb;
	double ef = (smid < 7) ? 2./3. : 0.;
	fact = -f2*( (1. - 2.*ef*sqr(_sw))*m1a*m1b + 2*ef*sqr(_sw)*m2a*m2b )
	  - f1*_masslast.first*(*_mixS)(iloc,1)*(m1a*m1b + m2a*m2b) * 
	  UnitRemoval::InvE;
	if( alpha != beta )
	  fact -= 0.5*f1*( _lambdaVEV*(*_mixS)(iloc,0) 
			   + _lambda*_couplast*_v1*(*_mixS)(iloc,2)
			   + af*(*_mixS)(iloc,1) ) *
	    (m2a*m1b + m1a*m2b) * UnitRemoval::InvE;
      }
    }
    //CP-odd
    else if( higgs == 36 || higgs == 46 ) {
      int iloc = (higgs - 36)/10;
      if( smid % 2 != 0 ) {
	fact = Complex(0.,0.5)*f1*(_lambdaVEV*(*_mixP)(iloc,1) 
				   + _lambda*_couplast*_v2*(*_mixP)(iloc,2)
				   - af*(*_mixP)(iloc,0)) *
	  (m2a*m1b - m1a*m2b) * UnitRemoval::InvE;
      }
      else {
	fact = Complex(0., 0.5)*f1*( _lambdaVEV*(*_mixP)(iloc,0)
				     + _lambda*_couplast*_v1*(*_mixP)(iloc,2)
				     - af*(*_mixP)(iloc,1) ) * 
	  (m2a*m1b - m1a*m2b) * UnitRemoval::InvE;
      }
    }
  }

  if( q2 != _q2last ) {
    _q2last = q2;
    _couplast = weakCoupling(q2);
  }
  setNorm(_couplast*fact);
}

Complex NMSSMHSFSFVertex::chargedHiggs(Energy2 q2, long id1, long id2) {
  //have id1 as up-type
  if( id1 % 2 != 0)
    swap(id1, id2);
  unsigned int alpha = ( id1/1000000 == 2 ) ? 1 : 0;
  unsigned int beta = ( id2/1000000 == 2 ) ? 1 : 0;
  
  long utype = (beta == 0) ? id1 - 1000000 : id1 - 2000000;
  long dtype = (alpha == 0) ? id2 - 1000000 : id2 - 2000000;
  if( q2 != _q2last || id1 != _idlast.first || id2 != _idlast.second) {
    _idlast.first = id1;
    _idlast.second = id2;
    tcPDPtr p = getParticleData(utype);
    _masslast.first = _theSM->mass(q2, p);
    p = getParticleData(dtype);
    _masslast.second = _theSM->mass(q2, p);
  }
  Energy2 facta = 2.*sqr(_mw)*_sb*_cb;
  complex<Energy> coupling(0.*MeV);
  if( dtype == 11 || dtype == 13 || dtype == 15) {
    Complex l1b = (beta == 0) ? 1.0 : 0.0;
    Complex l2b = (beta == 0) ? 0.0 : 1.0;
    complex<Energy> tri(0.*MeV);
    if( dtype == 15 ) {
      l1b = (*_mixTa)(beta, 0);
      l2b = (*_mixTa)(beta, 1);
      tri = _triTa;
    }
    coupling = ( l1b*(sqr(_masslast.second)*_tb - facta) 
		 + l2b*_masslast.second*(tri*_tb + _lambdaVEV) )/_mw/sqrt(2.);
  }
  else {
    Complex q1a(0.0), q1b(0.0), q2a(0.0), q2b(0.0);
    complex<Energy> triD(0.*MeV), triU(0.*MeV);
    if( utype == 2 || utype == 4 ) {
      q1a = (alpha == 0) ? 1.0 : 0.0;
      q2a = (alpha == 0) ? 0.0 : 1.0;
    }
    else {
      q1a = (*_mixTp)(alpha, 0);
      q2a = (*_mixTp)(alpha, 1);
      triU = _triTp;
    }
    if( dtype == 1 || dtype == 3 ) {
      q1b = (beta == 0) ? 1.0 : 0.0;
      q2b = (beta == 0) ? 0.0 : 1.0;
    }
    else {
      q1b = (*_mixBt)(beta, 0);
      q2b = (*_mixBt)(beta, 1);
      triD = _triBt;
    }
    Energy mfu = _masslast.first;
    Energy mfd = _masslast.second;
    coupling = ( q1a*q1b*(sqr(mfd)*_tb +sqr(mfu)/_tb - facta)
		 + q2a*q2b*mfu*mfd*(_tb + (1./_tb))
		 + q1a*q1b*mfd*(triD*_tb + _lambdaVEV)
		 + q2a*q1b*mfu*(triU/_tb + _lambdaVEV))/_mw/sqrt(2.);
  }
  return coupling * UnitRemoval::InvE;
}
