// -*- C++ -*-
//
// SSWGSSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSWGSSVertex class.
//

#include "SSWGSSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSWGSSVertex::SSWGSSVertex() : _sw(0.), _cw(0.), _q2last(),_emcouplast(0.),
			       _scouplast(0.), _ulast(0), _dlast(0),
			       _gblast(0), _factlast(0.)  {
  colourStructure(ColourStructure::SU3TFUND);
}

void SSWGSSVertex::doinit() {
  //W-
  //LL-squarks
  for(long ix=1000001;ix<1000006;ix+=2) {
    addToList(-24,21,ix+1,-ix);
  }
  //1-2 stop sbottom
  addToList(-24,21,1000006,-2000005);
  //2-1 stop sbottom
  addToList(-24,21,2000006,-1000005);
  //2-2 stop sbottom
  addToList(-24,21,2000006,-2000005);

  //W+
  for(long ix=1000001;ix<1000006;ix+=2) {
    addToList(24,21,-(ix+1),ix);
  }
  //1-2 stop sbottom
  addToList(24,21,-1000006,2000005);
  //2-1 stop sbottom
  addToList(24,21,-2000006,1000005);
  //2-2 stop sbottom
  addToList(24,21,-2000006,2000005);
  
  //---Z0----   
  //LL squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(23,21,ix,-ix);
  }
  //RR squarks
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(23,21,ix,-ix);
  }
  //L-Rbar stop
  addToList(23,21,1000006,-2000006);
  //Lbar-R stop
  addToList(23,21,-1000006,2000006);

  //L-Rbar sbottom
  addToList(23,21,1000005,-2000005);
  //Lbar-R sbottom
  addToList(23,21,-1000005,2000005);
  
  //----gamma----
  //squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(22,21,ix,-ix);
  }
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(22,21,ix,-ix);
  }
  orderInGem(1);
  orderInGs(1);

  VVSSVertex::doinit();
  tMSSMPtr theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  if(!theSS)
    throw InitException() << "SSWGSSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt( 1. - sqr(_sw) );
  _stop = theSS->stopMix();
  _sbottom = theSS->sbottomMix();
  if(!_stop || !_sbottom)
    throw InitException() << "SSWGSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbottom
			  << Exception::abortnow;
}


void SSWGSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _sw  << _cw << _stop << _sbottom;
}

void SSWGSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sw >> _cw >> _stop >> _sbottom;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSWGSSVertex,VVSSVertex>
describeHerwigSSWGSSVertex("Herwig::SSWGSSVertex", "HwSusy.so");


void SSWGSSVertex::Init() {

  static ClassDocumentation<SSWGSSVertex> documentation
    ("This implements the gluon-gluon-squark-squark vertex.");
}

void SSWGSSVertex::setCoupling(Energy2 q2,   tcPDPtr part1,
			       tcPDPtr part2,tcPDPtr part3,tcPDPtr part4) { 

  long boson(abs(part1->id()));
  long gluon(abs(part2->id()));
  if (gluon > boson) swap(gluon, boson);

  if( boson != ParticleID::Wplus && boson != ParticleID::Z0 && 
      boson != ParticleID::gamma ) {
    throw HelicityConsistencyError()
      << "SSWGSSVertex::setCoupling() - Vector particle in this "
      << "vertex is not a W/Z/gamma. " << boson << Exception::warning;
    norm(0.);
  }
  
  if( gluon != ParticleID::g ) {
    throw HelicityConsistencyError()
      << "SSWGSSVertex::setCoupling() - Vector particle in this "
      << "vertex is not a gluon. " << gluon << Exception::warning;
    norm(0.);
  }
  long sq1(abs(part3->id())),sq2(abs(part4->id()));
  if( (sq1 < 1000001 && sq1 > 1000006 && sq1 < 2000001 && sq1 > 2000006) ||
      (sq2 < 1000001 && sq2 > 1000006 && sq2 < 2000001 && sq2 > 2000006))
    throw HelicityConsistencyError()
      << "SSWGSSVertex::setCoupling() - There are no squarks in "
      << "this vertex! " << part3->id() << " " << part4->id() 
      << Exception::warning;

  if( sq1 % 2 != 0 ) swap(sq1, sq2);
  if( sq1 != _ulast || sq2 != _dlast || boson != _gblast) {
    _gblast = boson;
    _ulast = sq1;
    _dlast = sq2;
    //photon is simplest
    if( boson == ParticleID::gamma )
      _factlast = -2.*getParticleData(sq1)->charge()/eplus;

    else {
     //determine which helicity state
      unsigned int alpha(sq1/1000000 - 1), beta(sq2/1000000 - 1);
      //mixing factors
      Complex m1a(0.), m1b(0.);
      if( sq1 == ParticleID::SUSY_t_1 || sq1 == ParticleID::SUSY_t_2 )
	m1a = (*_stop)(alpha, 0);
      else if( sq1 == ParticleID::SUSY_b_1 || sq1 == ParticleID::SUSY_b_2 ) 
	m1a = (*_sbottom)(alpha, 0);
      else
	m1a = (alpha == 0) ? Complex(1.) : Complex(0.);

      if( sq2 == ParticleID::SUSY_t_1 || sq2 == ParticleID::SUSY_t_2 )
	m1b = (*_stop)(beta, 0);
      else if( sq2 == ParticleID::SUSY_b_1 || sq2 == ParticleID::SUSY_b_2 )
	m1b = (*_sbottom)(beta, 0);
      else
	m1b = (beta == 0) ? Complex(1.) : Complex(0.);
    
      //W boson
      if( boson == ParticleID::Wplus ) {
	_factlast = -1.*m1a*m1b*sqrt(2)/_sw;
      }
      //Z boson
      else {
	double lmda(1.);
	if( sq2 % 2 == 0 ) lmda = -1.;
	_factlast = lmda*m1a*m1b;
	if( alpha == beta) {
	    double ef = getParticleData(sq1)->charge()/eplus;
	    _factlast += 2.*ef*sqr(_sw);
	}
	_factlast *= 1./_cw/_sw; 	
      }
    }
  }
  if( q2 != _q2last || _emcouplast==0. || _scouplast==0. ) {
    _q2last = q2;
    _emcouplast = electroMagneticCoupling(q2);
    _scouplast = strongCoupling(q2);
  }
  norm(-_emcouplast*_scouplast*_factlast);
}
