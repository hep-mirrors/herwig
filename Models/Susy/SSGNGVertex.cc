// -*- C++ -*-
//
// SSGNGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGNGVertex class.
//

#include "SSGNGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/Susy/MixingMatrix.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

namespace {

  unsigned int neutralinoIndex(long id) {
    if(id> 1000000)
      return id<1000025 ? id-1000022 : (id-1000005)/10;
    else if(abs(id)<=16) 
      return (abs(id)-4)/2;
    else
      return id-13;
  }
}

SSGNGVertex::SSGNGVertex() : _includeOnShell(false), _realIntegral(false), 
			     _omitLightQuarkYukawas(false),
			     _sw(0.), _cw(0.), _idlast(0), 
			     _q2last(ZERO), _couplast(0.),
			     _leftlast(ZERO), _rightlast(ZERO),
			     _initLoops(false) {
  orderInGem(1);
  orderInGs(2);
}

void SSGNGVertex::doinit() {
  if(!_initLoops) {
    Looptools::ltini();
    _initLoops = true;
  }
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if(!model)
    throw InitException() << "SSGNGVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _theN  = model->neutralinoMix();
  if(!_theN )
    throw InitException() << "SSGNGVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  vector<long> ineu(4);
  ineu[0] = 1000022; ineu[1] = 1000023; 
  ineu[2] = 1000025; ineu[3] = 1000035;
  if(_theN->size().first==5)
    ineu.push_back(1000045);
  else if(_theN->size().first==7) {
    if(model->majoranaNeutrinos()) {
      ineu.push_back(17);
      ineu.push_back(18);
      ineu.push_back(19);
    }
    else {
      ineu.push_back(12);
      ineu.push_back(14);
      ineu.push_back(16);
    }
  }
  for(unsigned int i = 0; i < ineu.size(); ++i) {
    addToList(1000021, ineu[i], 21);
  }
  GeneralFFVVertex::doinit();
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1 - _sw*_sw);
  _mw = getParticleData(ParticleID::Wplus)->mass();
  double tb = model->tanBeta();
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
  _stop = model->stopMix();
  _sbot = model->sbottomMix();
  Looptools::ltexi();
}

void SSGNGVertex::dofinish() {
  Looptools::ltexi();
  GeneralFFVVertex::dofinish();
}

void SSGNGVertex::doinitrun() {
  if(!_initLoops) {
    Looptools::ltini();
    _initLoops = true;
  }
  GeneralFFVVertex::doinitrun();
}

void SSGNGVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw << _cw << _theN << ounit(_mw,GeV) << _sb << _cb
     << _stop << _sbot << _includeOnShell << _omitLightQuarkYukawas
     << _realIntegral;
}

void SSGNGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _theN >> iunit(_mw,GeV) >> _sb >> _cb
     >> _stop >> _sbot >> _includeOnShell >> _omitLightQuarkYukawas
     >> _realIntegral;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGNGVertex,GeneralFFVVertex>
describeHerwigSSGNGVertex("Herwig::SSGNGVertex", "HwSusy.so");

void SSGNGVertex::Init() {

  static ClassDocumentation<SSGNGVertex> documentation
    ("The loop-mediated coupling of the gluino to a gluon and a neutralino.");

  static Switch<SSGNGVertex,bool> interfaceIncludeOnShellIntermediates
    ("IncludeOnShellIntermediates",
     "Whether or not to include on-shell intermediate states",
     &SSGNGVertex::_includeOnShell, false, false, false);
  static SwitchOption interfaceIncludeOnShellIntermediatesYes
    (interfaceIncludeOnShellIntermediates,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeOnShellIntermediatesNo
    (interfaceIncludeOnShellIntermediates,
     "No",
     "Don't incldue them",
     false);

  static Switch<SSGNGVertex,bool> interfaceOmitLightQuarkYukawas
    ("OmitLightQuarkYukawas",
     "Omit the yukawa type couplings for down, up, strange"
     " and charm quarks, mainly for testing vs ISAJET",
     &SSGNGVertex::_omitLightQuarkYukawas, false, false, false);
  static SwitchOption interfaceOmitLightQuarkYukawasNo
    (interfaceOmitLightQuarkYukawas,
     "No",
     "Include the Yukawas",
     false);
  static SwitchOption interfaceOmitLightQuarkYukawasYes
    (interfaceOmitLightQuarkYukawas,
     "Yes",
     "Omit Yukawas",
     true);

  static Switch<SSGNGVertex,bool> interfaceRealIntegral
    ("RealIntegral",
     "Only include the real parts of the integrals",
     &SSGNGVertex::_realIntegral, false, false, false);
  static SwitchOption interfaceRealIntegralYes
    (interfaceRealIntegral,
     "Yes",
     "Only include the real part",
     true);
  static SwitchOption interfaceRealIntegralNo
    (interfaceRealIntegral,
     "No",
     "Don't include the real part",
     false);
}

void SSGNGVertex::setCoupling(Energy2 q2, tcPDPtr part1,
#ifndef NDEBUG
			      tcPDPtr part2,tcPDPtr part3) {
#else
                              tcPDPtr part2,tcPDPtr) {
#endif
  int o[2]={1,0};
  long in1 = part1->id();
  long in2 = part2->id();
  Energy Mj = part1->mass();
  Energy Mi = part2->mass();
  if(in1!=ParticleID::SUSY_g) {
    swap(in1,in2);
    swap(Mj,Mi);
  }
  // checks of the particle ids
  assert(part3->id()==ParticleID::g);
  assert(in1 == ParticleID::SUSY_g);
  assert(in2 == ParticleID::SUSY_chi_10 || in2 == ParticleID::SUSY_chi_20 ||
	 in2 == ParticleID::SUSY_chi_30 || in2 == ParticleID::SUSY_chi_40 ||
	 in2 == 1000045 || in2==12 || in2==14 || in2==16 || in2==17 || in2==18 || in2==19 );
  // normal couplings are zero
  setLeft (0.);
  setRight(0.);
  if(in2 != _idlast || q2 !=_q2last) {
    if(!_initLoops) {
      Looptools::ltini();
      _initLoops = true;
    }
    Looptools::clearcache();
    _leftlast  = ZERO;
    _rightlast = ZERO;
    _idlast = in2;
    unsigned int neu = neutralinoIndex(in2);
    Complex n1prime = (*_theN)(neu,0)*_cw + (*_theN)(neu,1)*_sw;
    Complex n2prime = (*_theN)(neu,1)*_cw - (*_theN)(neu,0)*_sw;
    // squark/quark loops
    for(long iferm=1;iferm<7;++iferm) {
      tcPDPtr smf = getParticleData(iferm);
      Energy mf = smf->mass();
      double qf = smf->charge()/eplus;
      double y = (!(iferm<=4&&_omitLightQuarkYukawas)) ?
	0.5*double(dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel())->mass(q2,smf)/_mw) : 0.;
      Complex bracketl = qf*_sw*( conj(n1prime) - _sw*conj(n2prime)/_cw );
      double lambda(0.);
      // neutralino mixing element
      Complex nlf(0.);
      if( iferm % 2 == 0 ) {
  	y /= _sb;
  	lambda = -0.5 + qf*sqr(_sw);
   	nlf = (*_theN)(neu,3);
      }
      else { 
  	y /= _cb;
  	lambda = 0.5 + qf*sqr(_sw);
  	nlf = (*_theN)(neu,2);
      }
      Complex bracketr = _sw*qf*n1prime - n2prime*lambda/_cw;
      for(long iy=0;iy<2;++iy) {
	long isf = 1000000*(1+iy)+iferm;
	Energy msf = getParticleData(isf)->mass();
	if(!_includeOnShell&&(mf+msf<Mj||mf+msf<Mi)) continue;
	Complex g[2][2];
	Complex ma1(0.), ma2(0.);
	// heavy fermions
	if( iferm == 5 || iferm == 6 ) {
	  if( iferm == 5 ) {
	    ma1 = (*_sbot)(iy,0);
	    ma2 = (*_sbot)(iy,1); 
	  } 
	  else if( iferm == 6 ) {
	    ma1 = (*_stop)(iy,0);
	    ma2 = (*_stop)(iy,1);
	  }
	}
	else if(iy==0) {
	  ma1 = 1.;
	}
	else {
	  ma2 = 1.;
	}
	g[0][0] = y*conj(nlf)*ma1 - ma2*bracketl;
	g[0][1] = y*nlf      *ma2 + ma1*bracketr;
	g[1][0] =  - ma2;
	g[1][1] =  + ma1;
 	swap(g[0][0],g[0][1]);
	complex<InvEnergy2> I,J,K,I2;
	loopIntegrals(Mi,Mj,msf,mf,I,J,K,I2);
	complex<InvEnergy> coup[2];
	for(unsigned int ix=0;ix<2;++ix) {
	  coup[ix] = Mj*(I2-K)*(g[0][ix]*g[1][o[ix]]-conj(g[0][o[ix]]*g[1][ix]))
	    +Mi*K*(g[0][o[ix]]*g[1][ix]-conj(g[0][ix]*g[1][o[ix]]))
	    +mf*I*(g[0][ix]*g[1][ix]-conj(g[0][o[ix]]*g[1][o[ix]]));
	}
	_leftlast  += 2.*coup[0];
	_rightlast += 2.*coup[1];
      }
    }
  }
  if(q2 != _q2last || _couplast==0.) {
    _q2last = q2;
    _couplast = weakCoupling(q2)*sqr(strongCoupling(q2))/
      32./sqr(Constants::pi);
  }
  norm(_couplast);
  setLeftSigma ( _leftlast);
  setRightSigma(_rightlast);
}

void SSGNGVertex::loopIntegrals(Energy Mi, Energy Mj, Energy M, Energy m,
				complex<InvEnergy2> & I, complex<InvEnergy2> & J,
				complex<InvEnergy2> & K, complex<InvEnergy2> & I2) {
  Energy2 m2(sqr(m)),M2(sqr(M)),Mi2(sqr(Mi)),Mj2(sqr(Mj));
  double min2  = Mj2*UnitRemoval::InvE2;
  double mout2 = Mi2*UnitRemoval::InvE2;
  double mf2   = m2 *UnitRemoval::InvE2;
  double ms2   = M2 *UnitRemoval::InvE2;
  I  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  J  = Looptools::C0i(Looptools::cc0,min2,mout2,0.,ms2,mf2,ms2)*UnitRemoval::InvE2;
  I2 =-Looptools::C0i(Looptools::cc1,min2,mout2,0.,mf2,ms2,mf2)*UnitRemoval::InvE2;
  K  = (1.+Complex(m2*I+M2*J-Mj2*I2))/(Mi2-Mj2);
  if(_realIntegral) {
    I  = I .real();
    J  = J .real();
    I2 = I2.real();
    K  = K .real();
  }
}
