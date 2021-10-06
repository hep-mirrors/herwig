// -*- C++ -*-
//
// HwppSelector.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwppSelector class.
//

#include "HwppSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Repository/UseRandom.h"
#include "CheckId.h"
#include <cassert>
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<HwppSelector,HadronSelector>
describeHwppSelector("Herwig::HwppSelector","Herwig.so");

IBPtr HwppSelector::clone() const {
  return new_ptr(*this);
}

IBPtr HwppSelector::fullclone() const {
  return new_ptr(*this);
}

void HwppSelector::doinit() {
  // the default partons allowed
  // the quarks
  for ( int ix=1; ix<=5; ++ix ) {
    partons().push_back(getParticleData(ix));
  }
  // the diquarks
  for(unsigned int ix=1;ix<=5;++ix) {
    for(unsigned int iy=1; iy<=ix;++iy) {
      if(ix==iy)
	partons().push_back(getParticleData(CheckId::makeDiquarkID(ix,iy,long(3))));
      else
        partons().push_back(getParticleData(CheckId::makeDiquarkID(ix,iy,long(1))));
    }
  }
  // weights for the different quarks etc
  for(unsigned int ix=0; ix<partons().size(); ++ix) {
    pwt()[partons()[ix]->id()]=0.;
  }
  pwt()[1]  = _pwtDquark;
  pwt()[2]  = _pwtUquark;
  pwt()[3]  = _pwtSquark;
  pwt()[4]  = _pwtCquark;
  pwt()[5]  = _pwtBquark;
  pwt()[1103] =       _pwtDIquark * _pwtDquark * _pwtDquark;
  pwt()[2101] = 0.5 * _pwtDIquark * _pwtUquark * _pwtDquark;
  pwt()[2203] =       _pwtDIquark * _pwtUquark * _pwtUquark;
  pwt()[3101] = 0.5 * _pwtDIquark * _pwtSquark * _pwtDquark;
  pwt()[3201] = 0.5 * _pwtDIquark * _pwtSquark * _pwtUquark;
  pwt()[3303] =       _pwtDIquark * _pwtSquark * _pwtSquark;
  HadronSelector::doinit();
  // lightest members (baryons)
  for(const PDPtr & p1 : partons()) {
    if(DiquarkMatcher::Check(p1->id())) continue;
    for(const PDPtr & p2 : partons()) {
      if(DiquarkMatcher::Check(p2->id())) continue;
      lightestBaryons_[make_pair(p1->id(),p2->id())] = lightestBaryonPair(p1,p2);
    }
  }
}

void HwppSelector::persistentOutput(PersistentOStream & os) const {
  os << _pwtDquark  << _pwtUquark << _pwtSquark
     << _pwtCquark << _pwtBquark << _pwtDIquark
     << _sngWt << _decWt 
     << _mode << _enhanceSProb << ounit(_m0Decay,GeV) << _massMeasure
     << _scHadronWtFactor << _sbHadronWtFactor << lightestBaryons_;
}

void HwppSelector::persistentInput(PersistentIStream & is, int) {
  is >> _pwtDquark  >> _pwtUquark >> _pwtSquark
     >> _pwtCquark >> _pwtBquark >> _pwtDIquark
     >> _sngWt >> _decWt 
     >> _mode >> _enhanceSProb >> iunit(_m0Decay,GeV) >> _massMeasure
     >> _scHadronWtFactor >> _sbHadronWtFactor >> lightestBaryons_;
}

void HwppSelector::Init() {

  static ClassDocumentation<HwppSelector> documentation
    ("The HwppSelector class implements the Herwig algorithm for selecting"
     " the hadrons",
     "The hadronization used the selection algorithm described in \\cite{Kupco:1998fx}.",
     "%\\cite{Kupco:1998fx}\n"
     "\\bibitem{Kupco:1998fx}\n"
     "  A.~Kupco,\n"
     "  ``Cluster hadronization in HERWIG 5.9,''\n"
     "  arXiv:hep-ph/9906412.\n"
     "  %%CITATION = HEP-PH/9906412;%%\n"
     );
  
  static Parameter<HwppSelector,double>
    interfacePwtDquark("PwtDquark","Weight for choosing a quark D",
		       &HwppSelector::_pwtDquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<HwppSelector,double>
    interfacePwtUquark("PwtUquark","Weight for choosing a quark U",
		       &HwppSelector::_pwtUquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<HwppSelector,double>
    interfacePwtSquark("PwtSquark","Weight for choosing a quark S",
		       &HwppSelector::_pwtSquark, 0, 1.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<HwppSelector,double>
    interfacePwtCquark("PwtCquark","Weight for choosing a quark C",
		       &HwppSelector::_pwtCquark, 0, 0.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<HwppSelector,double>
    interfacePwtBquark("PwtBquark","Weight for choosing a quark B",
		       &HwppSelector::_pwtBquark, 0, 0.0, 0.0, 10.0,
		       false,false,false);

  static Parameter<HwppSelector,double>
    interfacePwtDIquark("PwtDIquark","Weight for choosing a DIquark",
			&HwppSelector::_pwtDIquark, 0, 1.0, 0.0, 100.0,
			false,false,false);

  static Parameter<HwppSelector,double>
    interfaceSngWt("SngWt","Weight for singlet baryons",
                  &HwppSelector::_sngWt, 0, 1.0, 0.0, 10.0,
		   false,false,false);

  static Parameter<HwppSelector,double>
    interfaceDecWt("DecWt","Weight for decuplet baryons",
                  &HwppSelector::_decWt, 0, 1.0, 0.0, 10.0,
		   false,false,false);
  
  static Switch<HwppSelector,unsigned int> interfaceMode
    ("Mode",
     "Which algorithm to use",
     &HwppSelector::_mode, 1, false, false);

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

  static Switch<HwppSelector,int> interfaceEnhanceSProb
    ("EnhanceSProb",
     "Option for enhancing strangeness",
     &HwppSelector::_enhanceSProb, 0, false, false);

  static SwitchOption interfaceEnhanceSProbNo
    (interfaceEnhanceSProb,
     "No",
     "No strangeness enhancement.",
     0);

  static SwitchOption interfaceEnhanceSProbScaled
    (interfaceEnhanceSProb,
     "Scaled",
     "Scaled strangeness enhancement",
     1);

  static SwitchOption interfaceEnhanceSProbExponential
    (interfaceEnhanceSProb,
     "Exponential",
     "Exponential strangeness enhancement",
     2);

   static Switch<HwppSelector,int> interfaceMassMeasure
     ("MassMeasure",
      "Option to use different mass measures",
      &HwppSelector::_massMeasure,0,false,false);

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

   static Parameter<HwppSelector,double> interfacescHadronWtFactor
     ("scHadronWtFactor",
     "Wight factor for strenge-charm heavy hadrns",
     &HwppSelector::_scHadronWtFactor, 1., 0., 10.,
     false, false, Interface::limited);

   static Parameter<HwppSelector,double> interfacesbHadronWtFactor
     ("sbHadronWtFactor",
     "Wight factor for strenge-bottom heavy hadrns",
     &HwppSelector::_sbHadronWtFactor, 1., 0., 10.,
     false, false, Interface::limited);

  static Parameter<HwppSelector,Energy> interfaceDecayMassScale
    ("DecayMassScale",
     "Cluster decay mass scale",
     &HwppSelector::_m0Decay, GeV, 1.0*GeV, 0.1*GeV, 50.*GeV,
     false, false, Interface::limited);

}

double HwppSelector::baryonWeight(long id) const {
  const int pspin = id % 10;
  if(pspin == 2) {
    // Singlet (Lambda-like) baryon
    if( (id/100)%10 < (id/10 )%10 ) return sqr(_sngWt);
  }
  // Decuplet baryon
  else if (pspin == 4)              return sqr(_decWt);
  return 1.;
}

std::tuple<bool,bool,bool> HwppSelector::selectBaryon(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const {
  useMe();
  std::tuple<bool,bool,bool> output(true,true,true);
  useMe();
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

double HwppSelector::strangeWeight(const Energy cluMass, tcPDPtr par1, tcPDPtr par2) const {
  // Decoupling the weight of heavy strenge hadrons
  if(_enhanceSProb == 0 && abs(par1->id()) == 4) {
    return pwt(3)*_scHadronWtFactor;
  }
  else if(_enhanceSProb == 0 && abs(par1->id()) == 5) {
    return pwt(3)*_sbHadronWtFactor;
  }
  // Scaling strangeness enhancement
  else if(_enhanceSProb == 1) {
    double scale = double(sqr(_m0Decay/cluMass));
    return (_maxScale < scale) ? 0. : pow(pwt(3),scale);
  }
  // Exponential strangeness enhancement
  else if(_enhanceSProb == 2) {
    Energy2 mass2;
    Energy endpointmass = par1->mass() + par2->mass();
    // Choose to use either the cluster mass
    // or to use the lambda measure
    mass2 = (_massMeasure == 0) ? sqr(cluMass) :
      sqr(cluMass) - sqr(endpointmass);
    double scale = double(sqr(_m0Decay)/mass2);
    return (_maxScale < scale) ? 0. : exp(-scale);
  }
  return pwt(3);
}

tcPDPair HwppSelector::lightestBaryonPair(tcPDPtr ptr1, tcPDPtr ptr2) const {
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
