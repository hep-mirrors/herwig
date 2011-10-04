// -*- C++ -*-
//
// HwppSelector.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwppSelector class.
//

#include "HwppSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Kinematics.h"
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



namespace {
  int abs(PDT::Colour c) {
    return c > 0 ? c : -c;
  }
}

void HwppSelector::doinit() {
  HadronSelector::doinit();
}

void HwppSelector::persistentOutput(PersistentOStream & os) const {
  os << _mode;
}

void HwppSelector::persistentInput(PersistentIStream & is, int) {
  is >> _mode;
}

void HwppSelector::Init() {

  static ClassDocumentation<HwppSelector> documentation
    ("The HwppSelector class implements the Herwig++ algorithm for selecting"
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
     "Use the Herwig++ approach",
     1);

}

pair<tcPDPtr,tcPDPtr> HwppSelector::chooseHadronPair(const Energy cluMass,tcPDPtr par1, 
						     tcPDPtr par2,tcPDPtr )
  {
  // if either of the input partons is a diquark don't allow diquarks to be 
  // produced
  bool diquark = !(DiquarkMatcher::Check(par1->id()) || DiquarkMatcher::Check(par2->id()));
  bool quark = true;
  // if the Herwig++ algorithm 
  if(_mode ==1) {
    if(cluMass > massLightestBaryonPair(par1,par2) && 
       UseRandom::rnd() > 1./(1.+pwtDIquark())) {
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
  vector<Kupco> hadrons;
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
    if(tit1->second.empty()||tit2->second.empty()) continue;
    // if too massive skip
    if(cluMass <= tit1->second.begin()->mass + 
                  tit2->second.begin()->mass) continue; 
    // loop over the hadrons
    KupcoData::iterator H1,H2;
    for(H1 = tit1->second.begin();H1 != tit1->second.end(); ++H1) {
      for(H2 = tit2->second.begin();H2 != tit2->second.end(); ++H2) {
 	// break if cluster too light
 	if(cluMass < H1->mass + H2->mass) break;
 	// calculate the weight
 	weight = pwt(quarktopick) * H1->overallWeight * H2->overallWeight *
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
  assert(!( signHad1 == 0 || signHad2 == 0));
//     throw Exception() << "HwppSelector::selectPair " 
// 		      << "***Inconsistent Hadron " 
// 		      << hadrons[ix].idQ->id() << " " 
// 		      << hadrons[ix].hadron1->id() << " " 
// 		      << hadrons[ix].hadron2->id() << " " 
// 		      << signHad1 << " " << signHad2
// 		      << Exception::runerror;
  return make_pair
    ( signHad1 > 0 ? hadrons[ix].hadron1 : tcPDPtr(hadrons[ix].hadron1->CC()),
      signHad2 > 0 ? hadrons[ix].hadron2 : tcPDPtr(hadrons[ix].hadron2->CC()));
}
