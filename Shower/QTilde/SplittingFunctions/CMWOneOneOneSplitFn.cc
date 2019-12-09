  // -*- C++ -*-
  //
  // CMWOneOneOneSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the CMWOneOneOneSplitFn class.
  //

#include "CMWOneOneOneSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"


using namespace Herwig;

DescribeClass<CMWOneOneOneSplitFn,Herwig::OneOneOneSplitFn>
describeCMWOneOneOneSplitFn ("Herwig::CMWOneOneOneSplitFn","HwShower.so");

void CMWOneOneOneSplitFn::Init() {
  
  static ClassDocumentation<CMWOneOneOneSplitFn> documentation
  ("The CMWOneOneOneSplitFn class implements the g -> gg splitting function");
  
  static Reference<CMWOneOneOneSplitFn,ShowerAlpha>
  interfaceAlpha("Alpha",
                 "A reference to the Alpha object",
                 &Herwig::CMWOneOneOneSplitFn::alpha_,
                 false, false, true, false);
  
  static Switch<CMWOneOneOneSplitFn, bool> interfaceIsIS
  ("isInititalState",
   "Switch on if this kernel is used for initial state emission.",
   &CMWOneOneOneSplitFn::isIS_, 0, false, false);
  static SwitchOption interfaceIsISYes
  (interfaceIsIS,"No","The kernel is used for final state emissions.", 0);
  static SwitchOption interfaceIsISNo
  (interfaceIsIS,"Yes","The kernel is used for final state emissions.", 1);
  
}

void CMWOneOneOneSplitFn::persistentOutput(PersistentOStream & os) const {
  os << alpha_ << isIS_ ;
}

void CMWOneOneOneSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> alpha_ >> isIS_ ;
}



double CMWOneOneOneSplitFn::P(const double z, const Energy2 t,
                              const IdList & ids, const bool, const RhoDMatrix &)const {
  
  auto scale2=t;
  if (!isIS_){
    scale2*=pTScale() ? z*(1.-z):1.;
  }else{
    scale2*=pTScale() ? z*(1.-z):z;
  }
  
  return colourFactor(ids) * alpha_->value(scale2) * Kg(scale2)/2./Constants::pi/(z*(1.-z));
  
}

double CMWOneOneOneSplitFn::ratioP(const double z, const Energy2 t,
                                   const IdList & , const bool, const RhoDMatrix &) const {
  
  auto scale2=t;
  if (!isIS_){
    scale2*=pTScale() ? z*(1.-z):1.;
  }else{
    scale2*=pTScale() ? z*(1.-z):z;
  }
  
  return alpha_->value(scale2)  * Kg(scale2)/2./Constants::pi;
}

