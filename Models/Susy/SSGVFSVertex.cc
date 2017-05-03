// -*- C++ -*-
//
// SSGVFSVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGVFSVertex class.
//

#include "SSGVFSVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGVFSVertex::SSGVFSVertex() : MPlanck_(2.4e18*GeV) {
  orderInGem(1);
  orderInGs(0);
}

void SSGVFSVertex::persistentOutput(PersistentOStream & os) const {
  os << stop_ << sbot_ << stau_ << ounit(MPlanck_,GeV);
}

void SSGVFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> stop_ >> sbot_ >> stau_ >> iunit(MPlanck_,GeV);
}

void SSGVFSVertex::doinit() {
  //quarks
  for(long ix=1;ix<7;++ix){
    addToList( ParticleID::SUSY_Gravitino,  ix, -(1000000+ix) );
    addToList( ParticleID::SUSY_Gravitino,  ix, -(2000000+ix) );
    addToList( ParticleID::SUSY_Gravitino, -ix,  (1000000+ix) );
    addToList( ParticleID::SUSY_Gravitino, -ix,  (2000000+ix) );
  }
  //leptons
  for(long ix=11;ix<17;++ix) {
    addToList( ParticleID::SUSY_Gravitino,  ix, -(1000000+ix) );
    addToList( ParticleID::SUSY_Gravitino, -ix,  (1000000+ix) );
    
    if( ix % 2 != 0 ) {
      addToList( ParticleID::SUSY_Gravitino,  ix, -(2000000+ix) );
      addToList( ParticleID::SUSY_Gravitino, -ix,  (2000000+ix) );
    }
  }
  RFSVertex::doinit();
  tMSSMPtr model = dynamic_ptr_cast<tMSSMPtr>(generator()->standardModel());
  if( !model )
    throw InitException() << "SSGVFSVertex::doinit() - "
			  << "The model pointer is null."
			  << Exception::abortnow;
  MPlanck_ = model->MPlanck();
  stop_ = model->stopMix();
  sbot_ = model->sbottomMix();
  stau_ = model->stauMix();
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SSGVFSVertex,RFSVertex>
describeHerwigSSGVFSVertex("Herwig::SSGVFSVertex", "libHwSusy.so");

void SSGVFSVertex::Init() {

  static ClassDocumentation<SSGVFSVertex> documentation
    ("The SSGVFSVertex implements the coupling of  the gravitino to "
     "a fermion-sfermion");
}

void SSGVFSVertex::setCoupling(Energy2 ,
#ifndef NDEBUG
			       tcPDPtr part1,
#else
			       tcPDPtr,
#endif
			       tcPDPtr part2,tcPDPtr part3) {
  assert(part1->id()==ParticleID::SUSY_Gravitino);
  assert(part3->iSpin()==PDT::Spin0);
  norm(double(sqrt(2.)/MPlanck_*UnitRemoval::E));
  // sfermion mass eigenstate
  unsigned int alpha(abs(part3->id())/1000000 - 1);
  unsigned int ism(abs(part2->id()));
  Complex lc,rc;
  //heavy quarks/sleptons
  if( ism == 5 || ism == 6 || ism == 15 ) {
    Complex ma1(0.), ma2(0.);
    if( ism == 5 ) {
      ma1 = (*sbot_)(alpha,0);
      ma2 = (*sbot_)(alpha,1); 
    } 
    else if( ism == 6 ) {
      ma1 = (*stop_)(alpha,0);
      ma2 = (*stop_)(alpha,1);
    } 
    else {
      ma1 = (*stau_)(alpha,0);
      ma2 = (*stau_)(alpha,1);
    }
    lc = - ma2;
    rc = + ma1;
  }
  else {
    if( alpha == 0 ) {
      lc =  0.;
      rc =  1.;
    } 
    else {
      lc = -1.;
      rc =  0.;
    }
  }
  // determine the helicity order of the vertex
  if( part2->id() < 0 ) {
    left (conj(rc));
    right(conj(lc));
  }
  else {
    left (lc);
    right(rc);
  }
}
