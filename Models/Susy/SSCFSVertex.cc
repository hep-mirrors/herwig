// -*- C++ -*-
//
// SSCFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSCFSVertex class.
//

#include "SSCFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSCFSVertex::SSCFSVertex(): _sb(0.),_cb(0.),_mw(ZERO),
			    _q2last(0.*GeV2), _couplast(0.),
			    _leftlast(0.),_rightlast(0.),
			    _id1last(0), _id2last(0), _id3last(0),
			    yukawa_(1) {
  orderInGem(1);
  orderInGs(0);
}

void SSCFSVertex::doinit() {
  long chargino[2] = {1000024, 1000037};
  for(unsigned int ic = 0; ic < 2; ++ic) {
    //quarks 
    for(long ix = 1; ix < 7; ++ix) {
      if( ix % 2 == 0 ) {
	addToList(-chargino[ic],ix,-(999999+ix));
	
	addToList(-chargino[ic],ix,-(1999999+ix));
	
	addToList(-ix,chargino[ic],(999999+ix));
	
	addToList(-ix,chargino[ic],(1999999+ix));
      }
      else {
	addToList(-chargino[ic],-ix,(1000001+ix));

	addToList(-chargino[ic],-ix,2000001+ix);

	addToList(chargino[ic],ix,-(1000001+ix));

	addToList(chargino[ic],ix,-(2000001+ix));
      }
    }
    //leptons
    for(long ix = 11; ix < 17; ++ix) {
      if( ix % 2 == 0 ) {
	addToList(-chargino[ic],ix,-(999999+ix));
      
	addToList(-chargino[ic],ix,-(1999999+ix));

	addToList(-ix,chargino[ic],(999999+ix));

	addToList(-ix,chargino[ic],(1999999+ix));	
      }
      else {
	addToList(-chargino[ic],-ix,1000001+ix);

	addToList(chargino[ic],ix,-(1000001+ix));
      }
    }
  } 
  FFSVertex::doinit();
  _theSS = dynamic_ptr_cast<MSSMPtr>(generator()->standardModel());
  //mixing matrices
  _stop = _theSS->stopMix();
  _sbot = _theSS->sbottomMix();
  _stau = _theSS->stauMix();
  _umix = _theSS->charginoUMix();
  _vmix = _theSS->charginoVMix();

  if(!_stop || !_stau || !_sbot || !_umix || !_vmix)
    throw InitException() << "SSCFSVertex:: doinit  - " 
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbot
			  << " stau: " << _stau << " U: " << _umix
			  << " V:" << _vmix
			  << Exception::abortnow;

  _mw = getParticleData(24)->mass();
  double tb = _theSS->tanBeta();
  _sb = tb/sqrt(1 + sqr(tb));
  _cb = sqrt(1.- sqr(_sb));
}


void SSCFSVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSS << _sb << _cb << ounit(_mw,GeV) << _stop 
     << _sbot << _stau << _umix << _vmix << yukawa_;
}

void SSCFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSS  >> _sb >> _cb >> iunit(_mw,GeV) >> _stop
     >> _sbot >> _stau >> _umix >> _vmix >> yukawa_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSCFSVertex,FFSVertex>
describeHerwigSSCFSVertex("Herwig::SSCFSVertex", "HwSusy.so");

void SSCFSVertex::Init() {

  static ClassDocumentation<SSCFSVertex> documentation
    ("The implementation of the coupling of the charginos to fermion-"
     "sfermions.");

  static Switch<SSCFSVertex,unsigned int> interfaceYukawa
    ("Yukawa",
     "Whether or not to include the Yukawa type couplings",
     &SSCFSVertex::yukawa_, true, false, false);
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
     "Only include them for the third generation",
     2);

}

void SSCFSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long isc(abs(part3->id())), ism(abs(part1->id())), 
    ichg(abs(part2->id()));
  tcPDPtr smfermion = part1;
  if( ism / 1000000 == 1 )  {
    swap( ism, ichg);
    smfermion = part2;
  }
  //overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _q2last=q2;
    _couplast = -weakCoupling(q2);
  }
  norm(_couplast);


  if( ichg != _id1last || ism != _id2last || isc != _id3last ) {
    _id1last = ichg;
    _id2last = ism;
    _id3last = isc;
    // determine chargino and sfermion eigenstates
    unsigned int alpha(isc/1000000 - 1);
    unsigned int ch = (ichg == 1000024 ) ? 0 : 1;

    Complex ul1 = (*_umix)(ch,0);
    Complex ul2 = (*_umix)(ch,1);
    Complex vl1 = (*_vmix)(ch,0);
    Complex vl2 = (*_vmix)(ch,1);

    if( ism >= 11 && ism <= 16 ) {
      long lept = ( ism % 2 == 0 ) ? ism - 1 : ism;
      double y = 0.;
      if(yukawa_==1 || (lept==15 && yukawa_==2))
	y = double(_theSS->mass(q2, getParticleData(lept))/_mw/sqrt(2)/_cb);


      if( ism == 12 || ism == 14 ) {
	_leftlast = Complex(0., 0.);
	if( alpha == 0 )
	  _rightlast = ul1;
	else
	  _rightlast = -y*ul2;
      }
      else if( ism == 16 ) {
	_leftlast = Complex(0., 0.);
	_rightlast = ul1*(*_stau)(alpha, 0) - y*(*_stau)(alpha,1)*ul2;
      }
      else if( ism == 11 || ism == 13 || ism == 15 ) {
	_leftlast = -y*conj(ul2);
	_rightlast = vl1;
      }
    }
    else {
      double yd(0.), yu(0.);
      if(yukawa_==1 || ((ism==5 || ism==6 ) && yukawa_==2)) {
	if( ism % 2 == 0) {
	  yu = _theSS->mass(q2, getParticleData(ism))/_mw/sqrt(2)/_sb;
	  yd = _theSS->mass(q2, getParticleData(ism - 1))/_mw/sqrt(2)/_cb;
	}
	else {
	  yu = _theSS->mass(q2, getParticleData(ism + 1))/_mw/sqrt(2)/_sb;
	  yd = _theSS->mass(q2, getParticleData(ism))/_mw/sqrt(2)/_cb;
	}
      }
      //heavy quarks
      if( ism == 5 ) {
	_leftlast = -yd*conj(ul2)*(*_stop)(alpha,0);
	_rightlast = vl1*(*_stop)(alpha, 0) - yu*vl2*(*_stop)(alpha,1);
      }
      else if( ism == 6 ) {
	_leftlast = -yu*conj(vl2)*(*_sbot)(alpha,0);
	_rightlast = ul1*(*_sbot)(alpha, 0) - yd*ul2*(*_sbot)(alpha,1);
      }
      else {
	if( alpha == 0 ) {
	  _leftlast = (ism % 2 == 0) ? -yu*conj(vl2) : -yd*conj(ul2);
	  _rightlast = (ism % 2 == 0) ? ul1 : vl1;
	}
	else {
	  _leftlast = Complex(0.);
	  _rightlast = (ism % 2 == 0) ? -yd*ul2 : -yu*vl2;
	}
      }
    }
  }//end of coupling calculation

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
