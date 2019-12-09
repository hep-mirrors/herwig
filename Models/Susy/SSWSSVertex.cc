// -*- C++ -*-
//
// SSWSSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWSSVertex class.
//

#include "SSWSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWSSVertex::SSWSSVertex():_sw(0.), _cw(0.), _q2last(),_couplast(0.), 
			   _ulast(0), _dlast(0), _gblast(0),
			   _factlast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SSWSSVertex::doinit() {
  //W-
  //LL-squarks
  for(long ix=1000001;ix<1000006;ix+=2) {
    addToList(-24,ix+1,-ix);
  }
  //1-2 stop sbottom
  addToList(-24,1000006,-2000005);
  //2-1 stop sbottom
  addToList(-24,2000006,-1000005);
  //2-2 stop sbottom
  addToList(-24,2000006,-2000005);
 
  //LL-sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    addToList(-24,-ix,ix+1);
  }
  //2-L stau
  addToList(-24,-2000015,1000016);
  //W+
  for(long ix=1000001;ix<1000006;ix+=2) {
    addToList(24,-(ix+1),ix);
  }

//1-2 stop sbottom
  addToList(24,-1000006,2000005);
  //2-1 stop sbottom
  addToList(24,-2000006,1000005);
  //2-2 stop sbottom
  addToList(24,-2000006,2000005);

  //LL-sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    addToList(24,ix,-ix-1);
  }
  //2-L stau
  addToList(24,2000015,-1000016);
  
  //---Z0----
//LL-sleptons
  for(long ix=1000011;ix<1000017;++ix) {
    addToList(23,ix,-ix);
  }
  //RR-sleptons
  for(long ix=2000011;ix<2000016;ix+=2) {
    addToList(23,ix,-ix);
  }
  //L-Rbar stau
  addToList(23,1000015,-2000015);
  //Lbar-R stau
  addToList(23,-1000015,2000015);
   
  //LL squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(23,ix,-ix);
  }
  //RR squarks
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(23,ix,-ix);
  }
  //L-Rbar stop
  addToList(23,1000006,-2000006);
  //Lbar-R stop
  addToList(23,-1000006,2000006);

  //L-Rbar sbottom
  addToList(23,1000005,-2000005);
  //Lbar-R sbottom
  addToList(23,-1000005,2000005);
  
  //----gamma----
  //sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    addToList(22,ix,-ix);
  }
  for(long ix=2000011;ix<2000016;ix+=2) {
    addToList(22,ix,-ix);
  }
  //squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(22,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(22,ix,-ix);
  }
  VSSVertex::doinit();
  tMSSMPtr theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSWSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt( 1. - sqr(_sw) );
  _stop = theSS->stopMix();
  _sbottom = theSS->sbottomMix();
  _stau = theSS->stauMix();
  if(!_stop || !_stau || !_sbottom)
    throw InitException() << "SSWSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbottom
			  << " stau: " << _stau << Exception::abortnow;
}

void SSWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw  << _cw << _stau << _stop << _sbottom;
}

void SSWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _stau >> _stop >> _sbottom;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWSSVertex,VSSVertex>
describeHerwigSSWSSVertex("Herwig::SSWSSVertex", "HwSusy.so");

void SSWSSVertex::Init() {

  static ClassDocumentation<SSWSSVertex> documentation
    ("This is the implementation of the coupling of an SM boson "
     "a pair of sfermions");
  
}

void SSWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long boson(abs(part1->id()));
  assert( boson == ParticleID::Wplus || boson == ParticleID::Z0 ||
	  boson == ParticleID::gamma );
  long sf1(abs(part2->id())),sf2(abs(part3->id()));

  assert( (sf1 >= 1000001 && sf1 <= 1000006) 
	  || (sf1 >= 1000011 && sf1 <= 1000016)
	  || (sf1 >= 2000001 && sf1 <= 2000006)
	  || (sf1 >= 2000011 && sf1 <= 2000016) );
  
  assert( (sf2 >= 1000001 && sf2 <= 1000006) 
	  || (sf2 >= 1000011 && sf2 <= 1000016)
	  || (sf2 >= 2000001 && sf2 <= 2000006)
	  || (sf2 >= 2000011 && sf2 <= 2000016) );

  if( sf1 % 2 != 0 ) swap(sf1, sf2);
  if( sf1 != _ulast || sf2 != _dlast || boson != _gblast) {
    _gblast = boson;
    _ulast = sf1;
    _dlast = sf2;
    //photon is simplest
    if( boson == ParticleID::gamma )
      _factlast = getParticleData(sf1)->charge()/eplus;
    else {
      //determine which helicity state
      unsigned int alpha(sf1/1000000 - 1), beta(sf2/1000000 - 1);
      //mixing factors
      Complex m1a(0.), m1b(0.);
      if( sf1 == ParticleID::SUSY_t_1 || sf1 == ParticleID::SUSY_t_2 )
	m1a = (*_stop)(alpha, 0);
      else if( sf1 == ParticleID::SUSY_b_1 || sf1 == ParticleID::SUSY_b_2 )
	m1a = (*_sbottom)(alpha, 0);
      else if( sf1 == ParticleID::SUSY_tau_1minus || 
	       sf1 == ParticleID::SUSY_tau_2minus )
	m1a = (*_stau)(alpha, 0);
      else
	m1a = (alpha == 0) ? Complex(1.) : Complex(0.);

      if( sf2 == ParticleID::SUSY_t_1 || sf2 == ParticleID::SUSY_t_2 )
	m1b = (*_stop)(beta, 0);
      else if( sf2 == ParticleID::SUSY_b_1 || sf2 == ParticleID::SUSY_b_2 )
	m1b = (*_sbottom)(beta, 0);
      else if( sf2 == ParticleID::SUSY_tau_1minus || 
	       sf2 == ParticleID::SUSY_tau_2minus )
	m1b = (*_stau)(beta, 0);
      else
	m1b = (beta == 0) ? Complex(1.) : Complex(0.);
    
      //W boson
      if( boson == ParticleID::Wplus ) {
	_factlast = m1a*m1b/sqrt(2)/_sw;
      }
      //Z boson
      else {
	if( sf1 == ParticleID::SUSY_nu_eL || sf1 == ParticleID::SUSY_nu_muL ||
	    sf1 == ParticleID::SUSY_nu_tauL ) {
	  _factlast = 1./_cw/2./_sw;
	}
	else {
	  double lmda(1.);
	  if( sf2 % 2 == 0 ) lmda = -1.;
	  _factlast = lmda*m1a*m1b;
	  if( alpha == beta) {
	    double ef = getParticleData(sf1)->charge()/eplus;
	    _factlast += 2.*ef*sqr(_sw);
	  }
	  _factlast *= -1./2./_cw/_sw; 
	}
      }
    }
  }
  if( q2 != _q2last || _couplast==0. ) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  if(part2->id()>0)
    norm(-_couplast*_factlast);
  else
    norm(+_couplast*_factlast);
}

