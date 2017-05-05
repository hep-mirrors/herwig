// -*- C++ -*-
//
// SSNFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSNFSVertex class.
//

#include "SSNFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSNFSVertex::SSNFSVertex() :  _sw(0.), _cw(0.), _mw(), 
			     _sb(0.), _cb(0.), _q2last(), _couplast(0.),
			     _leftlast(0.), _rightlast(0.), _id1last(0), 
			     _id2last(0),
			      yukawa_(1) {
  orderInGem(1);
  orderInGs(0);
}

void SSNFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _stop << _sbot << _stau << _nmix << _theSS  << _sw << _cw 
     << ounit(_mw,GeV) << _sb << _cb << yukawa_;
}

void SSNFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _stop >> _sbot >> _stau >> _nmix >> _theSS >> _sw >> _cw 
     >> iunit(_mw,GeV) >> _sb >> _cb >> yukawa_;
}

void SSNFSVertex::doinit() {
  long neut[5] = {1000022, 1000023, 1000025, 1000035, 1000045};
  for(unsigned int nl = 0; nl < 5; ++nl) {
    //quarks
    for(long ix=1;ix<7;++ix){
      addToList( neut[nl],  ix, -(1000000+ix) );
      addToList( neut[nl],  ix, -(2000000+ix) );
      addToList( neut[nl], -ix,  (1000000+ix) );
      addToList( neut[nl], -ix,  (2000000+ix) );
    }
    //leptons
    for(long ix=11;ix<17;++ix) {
      addToList( neut[nl],  ix, -(1000000+ix) );
      addToList( neut[nl], -ix,  (1000000+ix) );
     
      if( ix % 2 != 0 ) {
	addToList( neut[nl],  ix, -(2000000+ix) );
	addToList( neut[nl], -ix,  (2000000+ix) );
      }
    }
    
  }
  FFSVertex::doinit();
  _theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!_theSS)
    throw InitException() << "SSGSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;

  _stop = _theSS->stopMix();
  _sbot = _theSS->sbottomMix();
  _stau = _theSS->stauMix();
  _nmix = _theSS->neutralinoMix();
  if(!_stop || !_stau || !_sbot || !_nmix)
    throw InitException() << "SSNFSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: "
			  << _sbot << " stau: " << _stau 
			  << " N: " << _nmix << Exception::abortnow;

  _sw = sqrt(sin2ThetaW());
  _mw = getParticleData(24)->mass();
  double tb = _theSS->tanBeta();
  _cw = sqrt(1. - sqr(_sw));
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1 - sqr(_sb));
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSNFSVertex,FFSVertex>
describeHerwigSSNFSVertex("Herwig::SSNFSVertex", "HwSusy.so");

void SSNFSVertex::Init() {

  static ClassDocumentation<SSNFSVertex> documentation
    ("The SSNFSVertex implements the coupling of a neutralino to "
     "a fermion-sfermion");

  static Switch<SSNFSVertex,unsigned int> interfaceYukawa
    ("Yukawa",
     "Whether or not to include the Yukawa type couplings",
     &SSNFSVertex::yukawa_, 1, false, false);
  static SwitchOption interfaceYukawaYes
    (interfaceYukawa,
     "Yes",
     "Include the terms",
     1);
  static SwitchOption interfaceYukawaNo
    (interfaceYukawa,
     "No",
     "Don't include them",
     0);
  static SwitchOption interfaceYukawa3rdGen
    (interfaceYukawa,
     "ThirdGeneration",
     "Only include for the third generation",
     2);
}

void SSNFSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long isc(abs(part3->id())), ism(abs(part1->id())),
    ineut(abs(part2->id()));
  tcPDPtr smfermion = part1;
  if( ism / 1000000 == 1 )  {
    swap( ism, ineut);
    smfermion = part2;
  }
  
  if(q2!=_q2last || _couplast==0.) {
    _couplast = -sqrt(2)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);

  if( ineut != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ineut;
    _id2last = ism;
    _id3last = isc;
    // determine neutralino and squark eigenstates
    unsigned int alpha(isc/1000000 - 1), nl(0);
    switch( ineut ) {
    case 1000022 : nl = 0;
      break;
    case 1000023 : nl = 1;
      break;
    case 1000025 : nl = 2;
      break;
    case 1000035 : nl = 3;
      break;
    case 1000045 : nl = 4;
      break;
    default : assert(false);
    }
    // common primed neutralino matrices
    Complex n2prime = (*_nmix)(nl,1)*_cw - (*_nmix)(nl,0)*_sw;
    //handle neutrinos first
    if( ism == 12 || ism == 14 || ism == 16 ) {
      _leftlast = Complex(0., 0.);
      _rightlast = n2prime/2./_cw;
    }
    else {
      Complex n1prime = (*_nmix)(nl,0)*_cw + (*_nmix)(nl,1)*_sw;
      tcPDPtr smf = getParticleData(ism);
      double qf = smf->charge()/eplus;
      Complex bracketl = qf*_sw*( conj(n1prime) - _sw*conj(n2prime)/_cw );
      double y = 0.;
      if(yukawa_==1 || ((ism==5 || ism==6 || ism==15) && yukawa_==2))
	y = double(_theSS->mass(q2, smf)/2./_mw);
      double lambda(0.);
      //neutralino mixing element
      Complex nlf(0.);
      if( ism % 2 == 0 ) {
	y /= _sb;
	lambda = -0.5 + qf*sqr(_sw);
	nlf = (*_nmix)(nl,3);
      }
      else { 
	y /= _cb;
	lambda = 0.5 + qf*sqr(_sw);
	nlf = (*_nmix)(nl,2);
      }
      Complex bracketr = _sw*qf*n1prime - n2prime*lambda/_cw;
      
      //heavy quarks/sleptons
      if( ism == 5 || ism == 6 || ism == 15 ) {
	Complex ma1(0.), ma2(0.);
	if( ism == 5 ) {
	  ma1 = (*_sbot)(alpha,0);
	  ma2 = (*_sbot)(alpha,1); 
	} 
	else if( ism == 6 ) {
	  ma1 = (*_stop)(alpha,0);
	  ma2 = (*_stop)(alpha,1);
	} 
	else {
	  ma1 = (*_stau)(alpha,0);
	  ma2 = (*_stau)(alpha,1);
	}
	_leftlast = y*conj(nlf)*ma1 - ma2*bracketl;
	_rightlast = y*nlf*ma2 + ma1*bracketr;
      }
      else {
	if( alpha == 0 ) {
	  _leftlast = y*conj(nlf);
	  _rightlast = bracketr;
	} 
	else {
	  _leftlast = -bracketl;
	  _rightlast = y*nlf;
	}
      }
    }
  }
  //determine the helicity order of the vertex
  if( smfermion->id() < 0 ) {
    left(conj(_rightlast));
    right(conj(_leftlast));
  }
  else {
    left(_leftlast);
    right(_rightlast);
  }
}
