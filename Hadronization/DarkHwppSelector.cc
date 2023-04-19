// -*- C++ -*-
//
// DarkHwppSelector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DarkHwppSelector class.
//

#include "DarkHwppSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Repository/UseRandom.h"
#include <cassert>
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<DarkHwppSelector,DarkHadronSpectrum>
describeDarkHwppSelector("Herwig::DarkHwppSelector","Herwig.so");

IBPtr DarkHwppSelector::clone() const {
  return new_ptr(*this);
}

IBPtr DarkHwppSelector::fullclone() const {
  return new_ptr(*this);
}

void DarkHwppSelector::doinit() {
  DarkHadronSpectrum::doinit();
}

void DarkHwppSelector::persistentOutput(PersistentOStream & os) const {
  os << _mode << ounit(_m0Decay,GeV) << _massMeasure;
}

void DarkHwppSelector::persistentInput(PersistentIStream & is, int) {
  is >> _mode >> iunit(_m0Decay,GeV) >> _massMeasure;
}

void DarkHwppSelector::Init() {

  static ClassDocumentation<DarkHwppSelector> documentation
    ("The DarkHwppSelector class implements the Herwig algorithm for selecting"
     " the hadrons",
     "The hadronization used the selection algorithm described in \\cite{Kupco:1998fx}.",
     "%\\cite{Kupco:1998fx}\n"
     "\\bibitem{Kupco:1998fx}\n"
     "  A.~Kupco,\n"
     "  ``Cluster hadronization in HERWIG 5.9,''\n"
     "  arXiv:hep-ph/9906412.\n"
     "  %%CITATION = HEP-PH/9906412;%%\n"
     );
    // put useMe() only in correct place!

  static Switch<DarkHwppSelector,unsigned int> interfaceMode
    ("Mode",
     "Which algorithm to use",
     &DarkHwppSelector::_mode, 1, false, false);
  static SwitchOption interfaceModeKupco
    (interfaceMode,
     "Kupco",
     "Use the Kupco approach",
     0);
  static SwitchOption interfaceModeHwpp
    (interfaceMode,
     "Hwpp",
     "Use the Herwig approach",
     1);

   static Switch<DarkHwppSelector,int> interfaceMassMeasure
     ("MassMeasure",
      "Option to use different mass measures",
      &DarkHwppSelector::_massMeasure,0,false,false);
   static SwitchOption interfaceMassMeasureMass
     (interfaceMassMeasure,
      "Mass",
      "Mass Measure",
      0);
   static SwitchOption interfaceMassMeasureLambda
     (interfaceMassMeasure,
      "Lambda",
      "Lambda Measure",
      1);

  static Parameter<DarkHwppSelector,Energy> interfaceDecayMassScale
    ("DecayMassScale",
     "Cluster decay mass scale",
     &DarkHwppSelector::_m0Decay, GeV, 1.0*GeV, 0.1*GeV, 50.*GeV,
     false, false, Interface::limited);

}

double DarkHwppSelector::baryonWeight(long id) const {
  const int pspin = id % 10;
  if(pspin == 2) {
    // Singlet (Lambda-like) baryon
    if( (id/100)%10 < (id/10 )%10 ) return sqr(_sngWt);
  }
  // Decuplet baryon
  else if (pspin == 4)              return sqr(_decWt);
  return 1.;
}

std::tuple<bool,bool,bool> DarkHwppSelector::selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const {
  useMe();
  std::tuple<bool,bool,bool> output(true,true,true);
  if(_mode ==1) {
    if(UseRandom::rnd() > 1./(1.+_pwtDIquark) && cluMass > massLightestBaryonPair(par1,par2)) {
      std::get<0>(output)  = false;
    }
    else {
      std::get<1>(output)  = false;
      std::get<2>(output)  = false;
    }
  }
  return output;
}

tcPDPair DarkHwppSelector::lightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const {
  // Make sure that we don't have any diquarks as input, return arbitrarily
  // large value if we do
  Energy currentSum = Constants::MaxEnergy;
  tcPDPair output;
  for(unsigned int ix=0; ix<partons().size(); ++ix) {
    if(!DiquarkMatcher::Check(partons()[ix]->id())) continue;
    HadronTable::const_iterator
      tit1=table().find(make_pair(abs(ptr1->id()),partons()[ix]->id())),
      tit2=table().find(make_pair(partons()[ix]->id(),abs(ptr2->id())));
    if( tit1==table().end() || tit2==table().end()) continue;
    if(tit1->second.empty()||tit2->second.empty()) continue;
    Energy s = tit1->second.begin()->mass + tit2->second.begin()->mass;
    if(currentSum > s) {
      currentSum = s;
      output.first  = tit1->second.begin()->ptrData;
      output.second = tit2->second.begin()->ptrData;
    }
  }
  return output;
}
