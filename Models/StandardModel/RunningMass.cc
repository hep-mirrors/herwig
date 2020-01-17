// -*- C++ -*-
//
// RunningMass.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RunningMass class.
//

#include "RunningMass.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"

namespace Herwig {
using namespace ThePEG;

void RunningMass::persistentOutput(PersistentOStream & os) const {
  os << _theQCDOrder << _thePower << _theCoefficient << _theMaxFlav
     << _theStandardModel << _lightOption << _heavyOption;
}

void RunningMass::persistentInput(PersistentIStream & is, int) {
  is >> _theQCDOrder >> _thePower >> _theCoefficient >> _theMaxFlav 
     >> _theStandardModel >> _lightOption >> _heavyOption;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RunningMass,RunningMassBase>
describeHerwigRunningMass("Herwig::RunningMass", "Herwig.so");

void RunningMass::Init() {

  static Parameter<RunningMass,unsigned int> interfaceQCDOrder
    ("QCDOrder",
     "The order in alpha_S",
     &RunningMass::_theQCDOrder, 1, 1, 2,
     false, false, true);
    
  static Parameter<RunningMass,unsigned int> interfaceMaxFlav
    ("MaxFlav",
     "maximum number of flavours",
     &RunningMass::_theMaxFlav, 6, 3, 6,
     false, false, true);
  
  static ClassDocumentation<RunningMass> documentation
    ("The RunningMass class is the implementation of the"
     " QCD running mass to one or two loop in alpha_S");

  static Switch<RunningMass,unsigned int> interfaceLightQuarkMass
    ("LightQuarkMass",
     "Option for the treatment of light mass masses",
     &RunningMass::_lightOption, 1, false, false);
  static SwitchOption interfaceLightQuarkMassRunning
    (interfaceLightQuarkMass,
     "Running",
     "Use a running, probably zero mass",
     0);
  static SwitchOption interfaceLightQuarkMassPole
    (interfaceLightQuarkMass,
     "Pole",
     "Use the pole mass",
     1);

  static Switch<RunningMass,unsigned int> interfaceBottomCharmMass
    ("TopBottomCharmMass",
     "Option for using a running or pole mass for the top, bottom and charm quarks",
     &RunningMass::_heavyOption, 0, false, false);
  static SwitchOption interfaceBottomCharmMassRunning
    (interfaceBottomCharmMass,
     "Running",
     "Use the running mass",
     0);
  static SwitchOption interfaceBottomCharmMassPole
    (interfaceBottomCharmMass,
     "Pole",
     "Use the pole mass",
     1);

}
// Return the masses used.
vector<Energy> RunningMass::mass() const {
  using Constants::pi;
  vector<Energy> masses;
  for ( long f = 1; f <= long(_theMaxFlav); ++f ) {
    PDPtr p = getParticleData(f);
    Energy massf = p ? p->mass() : ZERO;
    if ( !_theStandardModel->alphaSPtr()->quarkMasses().empty() ) {
      if ( f < 3 ) {
	massf = 
	  f == 1 ?
	  _theStandardModel->alphaSPtr()->quarkMasses()[1] :
	  _theStandardModel->alphaSPtr()->quarkMasses()[0];
      } else {
	massf = _theStandardModel->alphaSPtr()->quarkMasses()[f-1];
      }
    }
    if((f<=3&&_lightOption==0) ||
       (f>3 &&_heavyOption==0)) {
      double coeff = _theQCDOrder==2 ? _theCoefficient[f-1]+4./3./pi : 0.;
      double as = _theStandardModel->alphaS(massf*massf);
      massf = as>0 ? massf/(1.+coeff*as)/pow(as,_thePower[f-1]) : ZERO;
    }
    masses.push_back(massf);
  }
  return masses;
}

// Return the running mass for a given scale and particle type.
Energy RunningMass::value(Energy2 scale, tcPDPtr part) const {
  Energy output;
  unsigned int id=abs(part->id());
  // calculate the running mass
  if(id<=_theMaxFlav) {
    if( (id <= 3 && _lightOption == 1 ) ||
	(id >= 4 && _heavyOption == 1 ) ) {
      output= part->mass();
    }
    else {
      // calculate the value of alphaS and number of flavours
      //unsigned int nf=  _theStandardModel->Nf(scale);
      unsigned int nf=id;
      double as = _theStandardModel->alphaS(scale);
      id=id-1;
      output = massElement(id)*pow(as,_thePower[nf-1]);
      if(_theQCDOrder==2){output*=(1.+as*_theCoefficient[nf-1]);}
    }
  }
  // 
  else {
    output= part->mass();
  }
  return output;
}

void RunningMass::doinit() {
  using Constants::pi;
  _theStandardModel = generator()->standardModel();
  _theStandardModel->alphaSPtr()->init();
  // coefficients for the calculation
  double c = 1./pi,cprime,b,bprime,power,coeff;
  for(unsigned int f=1;f<=_theMaxFlav;++f) {
    // the basic parameters for the running mass
    cprime =     c*(303.-10.*f)/72.;
    b      =     c*(33. -2. *f)/12.;
    bprime = 0.5*c*(153.-19.*f)/(33.-2.*f);
    power = c/b;
    coeff = c*(cprime-bprime)/b;
    _thePower.push_back(power);
    _theCoefficient.push_back(coeff);
  }
  RunningMassBase::doinit();
}
}
