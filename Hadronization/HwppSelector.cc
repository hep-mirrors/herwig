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
describeHwppSelector("Herwig::HwppSelector","");

IBPtr HwppSelector::clone() const {
  return new_ptr(*this);
}

IBPtr HwppSelector::fullclone() const {
  return new_ptr(*this);
}

void HwppSelector::doinit() {
  HadronSelector::doinit();
}

void HwppSelector::persistentOutput(PersistentOStream & os) const {
  os << _mode << _enhanceSProb << ounit(_m0Decay,GeV) << _massMeasure;
}

void HwppSelector::persistentInput(PersistentIStream & is, int) {
  is >> _mode >> _enhanceSProb >> iunit(_m0Decay,GeV) >> _massMeasure;
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
    // put useMe() only in correct place!


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

  static Parameter<HwppSelector,Energy> interfaceDecayMassScale
    ("DecayMassScale",
     "Cluster decay mass scale",
     &HwppSelector::_m0Decay, GeV, 1.0*GeV, 0.1*GeV, 50.*GeV,
     false, false, Interface::limited);

}

pair<tcPDPtr,tcPDPtr> HwppSelector::chooseHadronPair(const Energy cluMass,tcPDPtr par1,
						     tcPDPtr par2,tcPDPtr ) const {
  // if either of the input partons is a diquark don't allow diquarks to be
  // produced
  bool diquark = !(DiquarkMatcher::Check(par1->id()) || DiquarkMatcher::Check(par2->id()));
  bool quark = true;
  // if the Herwig algorithm
  if(_mode ==1) {
    if(UseRandom::rnd() > 1./(1.+pwtDIquark())
       &&cluMass > massLightestBaryonPair(par1,par2)) {
      diquark = true;
      quark = false;
    }
    else {
      useMe();
      diquark = false;
      quark = true;
    }
  }
  // weights for the different possibilities
  Energy weight, wgtsum(ZERO);
  // loop over all hadron pairs with the allowed flavours
  static vector<Kupco> hadrons;
  hadrons.clear();
  for(unsigned int ix=0;ix<partons().size();++ix) {
    tcPDPtr quarktopick  = partons()[ix];
    if(!quark  &&  abs(int(quarktopick->iColour())) == 3
       && !DiquarkMatcher::Check(quarktopick->id())) continue;
    if(!diquark && abs(int(quarktopick->iColour())) == 3
       && DiquarkMatcher::Check(quarktopick->id())) continue;
    HadronTable::const_iterator
      tit1 = table().find(make_pair(abs(par1->id()),quarktopick->id()));
    HadronTable::const_iterator
      tit2 = table().find(make_pair(quarktopick->id(),abs(par2->id())));
    // If not in table skip
    if(tit1 == table().end()||tit2==table().end()) continue;
    // tables empty skip
    const KupcoData & T1 = tit1->second;
    const KupcoData & T2 = tit2->second;
    if(T1.empty()||T2.empty()) continue;
    // if too massive skip
    if(cluMass <= T1.begin()->mass +
                  T2.begin()->mass) continue;
    // quark weight
    double quarkWeight =  pwt(quarktopick->id());
    if (quarktopick->id() == 3) {
      // Scaling strangeness enhancement
      if (_enhanceSProb == 1){
	double scale = double(sqr(_m0Decay/cluMass));
	quarkWeight = (_maxScale < scale) ? 0. : pow(quarkWeight,scale);
      }
      // Exponential strangeness enhancement
      else if (_enhanceSProb == 2) {
	Energy2 mass2;
	Energy endpointmass = par1->mass() + par2->mass();
	// Choose to use either the cluster mass
	// or to use the lambda measure
	mass2 = (_massMeasure == 0) ? sqr(cluMass) :
	  sqr(cluMass) - sqr(endpointmass);
	double scale = double(sqr(_m0Decay)/mass2);
	quarkWeight = (_maxScale < scale) ? 0. : exp(-scale);
      }
    }
    // loop over the hadrons
    KupcoData::const_iterator H1,H2;
    for(H1 = T1.begin();H1 != T1.end(); ++H1) {
      for(H2 = T2.begin();H2 != T2.end(); ++H2) {
 	// break if cluster too light
 	if(cluMass < H1->mass + H2->mass) break;
	weight = quarkWeight * H1->overallWeight * H2->overallWeight *
	  Kinematics::pstarTwoBodyDecay(cluMass, H1->mass, H2->mass );

	int signQ = 0;
	assert (par1 && quarktopick);
	assert (par2);

	assert(quarktopick->CC());

	if(CheckId::canBeHadron(par1, quarktopick->CC())
	   && CheckId::canBeHadron(quarktopick, par2))
	   signQ = +1;
	else if(CheckId::canBeHadron(par1, quarktopick)
		&& CheckId::canBeHadron(quarktopick->CC(), par2))
	   signQ = -1;
	else {
	  cerr << "Could not make sign for" << par1->id()<< " " << quarktopick->id()
	       << " " << par2->id() << "\n";
	  assert(false);
	}

	if (signQ  == -1)
	  quarktopick = quarktopick->CC();
	// construct the object with the info
	Kupco a(quarktopick, H1->ptrData, H2->ptrData, weight);
	hadrons.push_back(a);
	wgtsum += weight;
      }
    }
  }
  if (hadrons.empty())
    return make_pair(tcPDPtr(),tcPDPtr());
  // select the hadron
  wgtsum *= UseRandom::rnd();
  unsigned int ix=0;
  do {
    wgtsum-= hadrons[ix].weight;
    ++ix;
  }
  while(wgtsum > ZERO && ix < hadrons.size());
  if(ix == hadrons.size() && wgtsum > ZERO)
      return make_pair(tcPDPtr(),tcPDPtr());
  --ix;
  assert(hadrons[ix].idQ);
  int signHad1 = signHadron(par1, hadrons[ix].idQ->CC(), hadrons[ix].hadron1);
  int signHad2 = signHadron(par2, hadrons[ix].idQ, hadrons[ix].hadron2);
  assert( signHad1 != 0 && signHad2 != 0 );
  return make_pair
    ( signHad1 > 0 ? hadrons[ix].hadron1 : tcPDPtr(hadrons[ix].hadron1->CC()),
      signHad2 > 0 ? hadrons[ix].hadron2 : tcPDPtr(hadrons[ix].hadron2->CC()));
}
