// -*- C++ -*-
//
// SSWSSVertex.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWSSVertex class.
//

#include "SSWSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWSSVertex::SSWSSVertex():_sw(0.), _cw(0.), _q2last(),_couplast(0.), 
				  _ulast(0), _dlast(0), _gblast(0),
				  _factlast(0.) {
  vector<long> first,second,third;
  //W-
  //LL-squarks
  for(long ix=1000001;ix<1000006;ix+=2) {
    first.push_back(-24);
    second.push_back(ix+1);
    third.push_back(-ix);
  }
  //1-2 stop sbottom
  first.push_back(-24);
  second.push_back(1000006);
  third.push_back(-2000005);
  //2-1 stop sbottom
  first.push_back(-24);
  second.push_back(2000006);
  third.push_back(-1000005);
  //2-2 stop sbottom
  first.push_back(-24);
  second.push_back(2000006);
  third.push_back(-2000005);
 
  //LL-sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    first.push_back(-24);
    second.push_back(-ix);
    third.push_back(ix+1);
  }
  //2-L stau
  first.push_back(-24);
  second.push_back(-2000015);
  third.push_back(1000016);
  //W+
  for(long ix=1000001;ix<1000006;ix+=2) {
    first.push_back(24);
    second.push_back(-(ix+1));
    third.push_back(ix);
  }

//1-2 stop sbottom
  first.push_back(24);
  second.push_back(-1000006);
  third.push_back(2000005);
  //2-1 stop sbottom
  first.push_back(24);
  second.push_back(-2000006);
  third.push_back(1000005);
  //2-2 stop sbottom
  first.push_back(24);
  second.push_back(-2000006);
  third.push_back(2000005);

  //LL-sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    first.push_back(24);
    second.push_back(ix);
    third.push_back(-ix-1);
  }
  //2-L stau
  first.push_back(24);
  second.push_back(2000015);
  third.push_back(-1000016);
  
  //---Z0----
//LL-sleptons
  for(long ix=1000011;ix<1000017;++ix) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //RR-sleptons
  for(long ix=2000011;ix<2000016;ix+=2) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //L-Rbar stau
  first.push_back(23);
  second.push_back(1000015);
  third.push_back(-2000015);
  //Lbar-R stau
  first.push_back(23);
  second.push_back(-1000015);
  third.push_back(2000015);
   
  //LL squarks
  for(long ix=1000001;ix<1000007;++ix) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //RR squarks
  for(long ix=2000001;ix<2000007;++ix) {
    first.push_back(23);
    second.push_back(ix);
    third.push_back(-ix);
  }
 //L-Rbar stop
  first.push_back(23);
  second.push_back(1000006);
  third.push_back(-2000006);
  //Lbar-R stop
  first.push_back(23);
  second.push_back(-1000006);
  third.push_back(2000006);

  //L-Rbar sbottom
  first.push_back(23);
  second.push_back(1000005);
  third.push_back(-2000005);
  //Lbar-R sbottom
  first.push_back(23);
  second.push_back(-1000005);
  third.push_back(2000005);
  
  //----gamma----
  //sleptons
  for(long ix=1000011;ix<1000016;ix+=2) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(long ix=2000011;ix<2000016;ix+=2) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  //squarks
  for(long ix=1000001;ix<1000007;++ix) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    first.push_back(22);
    second.push_back(ix);
    third.push_back(-ix);
  }
  setList(first,second,third);
}

void SSWSSVertex::doinit() {
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
  orderInGem(1);
  orderInGs(0);
}

void SSWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw  << _cw << _stau << _stop << _sbottom;
}

void SSWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _stau >> _stop >> _sbottom;
}

ClassDescription<SSWSSVertex> SSWSSVertex::initSSWSSVertex;
// Definition of the static class description member.

void SSWSSVertex::Init() {

  static ClassDocumentation<SSWSSVertex> documentation
    ("This is the implementation of the coupling of an SM boson "
     "a pair of sfermions");
  
}

void SSWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3){
  long boson(abs(part1->id()));
  if( boson != ParticleID::Wplus && boson != ParticleID::Z0 && 
      boson != ParticleID::gamma ) {
    throw HelicityConsistencyError()
      << "SSWSSVertex::setCoupling() - The vector particle in this "
      << "vertex is not a W/Z. " << boson << Exception::warning;
    setNorm(0.);
  }
  long sf1(abs(part2->id())),sf2(abs(part3->id()));
  if( (sf1 > 1000006 && sf1 < 1000011 && sf1 > 1000016) ||
      (sf1 > 2000006 && sf1 < 2000011 && sf1 > 2000016) ||
      (sf2 > 1000006 && sf2 < 1000011 && sf2 > 1000016) ||
      (sf2 > 2000006 && sf2 < 2000011 && sf2 > 2000016) )
    throw HelicityConsistencyError()
      << "SSWSSVertex::setCoupling() - There are no sfermions in "
      << "this vertex! " << part2->id() << " " << part3->id() 
      << Exception::warning;

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
  setNorm(_couplast*_factlast);
}

